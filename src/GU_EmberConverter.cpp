/**
 * @file GU_EmberConverter.cpp
 * @brief Houdini GU_Detail <-> EMBER PolygonSoup conversion
 *
 * H21 API: uses GA_FOR_ALL_PRIMOFF, GEO_PrimPoly::build(),
 *          gdp->setVertexPoint(), GEO_Face::appendVertex()/close()
 */

#include "GU_EmberConverter.h"
#include "ember/PolygonSoup.h"
#include "ember/MeshImport.h"

#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <GA/GA_Attribute.h>
#include <GA/GA_Iterator.h>
#include <GEO/GEO_Primitive.h>
#include <GEO/GEO_PrimPoly.h>
#include <GEO/GEO_Face.h>
#include <SYS/SYS_Math.h>

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cstring>

namespace ember {

//=============================================================================
// ATTRIBUTE MAPPING
//=============================================================================

EmberAttributeMapping buildAttributeMapping(const GU_Detail* gdp) {
    EmberAttributeMapping mapping;

    // Position (P) - always present
    mapping.P_handle = GA_ROHandleV3(gdp->getP());

    // Normal (N) - point attribute
    const GA_Attribute* N_attr = gdp->findNormalAttribute(GA_ATTRIB_POINT);
    if (N_attr) {
        mapping.N_handle = GA_ROHandleV3(N_attr);
    }

    // UV - prefer GA_ATTRIB_VERTEX (preserves seams), fall back to GA_ATTRIB_POINT
    const GA_Attribute* uv_attr = gdp->findTextureAttribute(GA_ATTRIB_VERTEX);
    if (uv_attr) {
        mapping.uv_handle    = GA_ROHandleV3(uv_attr);
        mapping.uv_is_vertex = true;
    } else {
        uv_attr = gdp->findTextureAttribute(GA_ATTRIB_POINT);
        if (uv_attr) {
            mapping.uv_handle    = GA_ROHandleV3(uv_attr);
            mapping.uv_is_vertex = false;
        }
    }

    // Color (Cd) - point attribute
    const GA_Attribute* Cd_attr = gdp->findDiffuseAttribute(GA_ATTRIB_POINT);
    if (Cd_attr) {
        mapping.Cd_handle = GA_ROHandleV3(Cd_attr);
    }

    // Velocity (v) - point attribute
    const GA_Attribute* v_attr = gdp->findVelocityAttribute(GA_ATTRIB_POINT);
    if (v_attr) {
        mapping.v_handle = GA_ROHandleV3(v_attr);
    }

    return mapping;
}

//=============================================================================
// CONVERSION: GU_Detail -> PolygonSoup
//=============================================================================

PolygonSoup convertGUDetailToPolygonSoup(const GU_Detail* gdp, int mesh_id) {
    PolygonSoup soup;
    if (!gdp) return soup;

    // H21: GA_FOR_ALL_PRIMOFF iterates all primitive offsets
    GA_Offset prim_off;
    GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
        const GEO_Primitive* prim = gdp->getGEOPrimitive(prim_off);
        if (!prim) continue;
        if (prim->getTypeId() != GA_PRIMPOLY) continue;

        const int num_verts = static_cast<int>(prim->getVertexCount());
        if (num_verts < 3) continue;

        // Fan-triangulate polygons with >3 vertices
        // H21: GEO_Primitive::getPos3(int) returns UT_Vector3
        UT_Vector3 p0 = prim->getPos3(0);

        for (int i = 1; i < num_verts - 1; ++i) {
            Triangle tri;
            tri.mesh_id      = mesh_id;
            tri.src_prim_idx = static_cast<uint32_t>(prim_off);

            UT_Vector3 p1 = prim->getPos3(i);
            UT_Vector3 p2 = prim->getPos3(i + 1);

            tri.v[0] = { p0.x(), p0.y(), p0.z() };
            tri.v[1] = { p1.x(), p1.y(), p1.z() };
            tri.v[2] = { p2.x(), p2.y(), p2.z() };

            // iv[] will be filled by soup.quantizeAll()
            tri.iv[0] = { 0, 0, 0 };
            tri.iv[1] = { 0, 0, 0 };
            tri.iv[2] = { 0, 0, 0 };

            soup.triangles.push_back(tri);
        }
    }

    return soup;
}

//=============================================================================
// CONVERSION: PolygonSoup -> GU_Detail (standard path)
//=============================================================================

void convertPolygonSoupToGUDetail(const PolygonSoup& soup,
                                   GU_Detail* gdp,
                                   const EmberAttributeMapping& mapping)
{
    if (!gdp) return;
    gdp->clearAndDestroy();

    // Output handles
    GA_RWHandleV3 P_out(gdp->getP());

    GA_RWHandleV3 uv_out;
    if (mapping.uv_handle.isValid()) {
        GA_AttributeOwner owner = mapping.uv_is_vertex
            ? GA_ATTRIB_VERTEX : GA_ATTRIB_POINT;
        uv_out = GA_RWHandleV3(gdp->addFloatTuple(owner, "uv", 3));
    }

    GA_RWHandleV3 N_out;
    if (mapping.N_handle.isValid()) {
        N_out = GA_RWHandleV3(gdp->addFloatTuple(GA_ATTRIB_POINT, "N", 3));
    }

    GA_RWHandleV3 Cd_out;
    if (mapping.Cd_handle.isValid()) {
        Cd_out = GA_RWHandleV3(gdp->addFloatTuple(GA_ATTRIB_POINT, "Cd", 3));
    }

    const bool has_output = !soup.output_polygons.empty();

    if (has_output) {
        // ── Indexed output from backend ──────────────────────────────────────
        std::vector<GA_Offset> ptoffs;
        ptoffs.reserve(soup.int_vertices.size());

        const double sx = soup.quantization_inv_scale[0];
        const double sy = soup.quantization_inv_scale[1];
        const double sz = soup.quantization_inv_scale[2];

        for (const auto& iv : soup.int_vertices) {
            GA_Offset ptoff = gdp->appendPoint();
            P_out.set(ptoff, UT_Vector3(
                static_cast<float>(iv.x * sx),
                static_cast<float>(iv.y * sy),
                static_cast<float>(iv.z * sz)));
            ptoffs.push_back(ptoff);
        }

        // Build polygons using H21 API
        for (const auto& poly : soup.output_polygons) {
            const int nv = static_cast<int>(poly.vertex_indices.size());
            if (nv < 3) continue;

            // H21: GEO_PrimPoly::build(gdp, npts, closed, do_append_points)
            //   do_append_points=0 means we must assign vertex->point ourselves
            GEO_PrimPoly* ppoly = GEO_PrimPoly::build(gdp, nv, GU_POLY_CLOSED, 0);
            if (!ppoly) continue;

            for (int vi = 0; vi < nv; ++vi) {
                uint32_t idx = poly.vertex_indices[vi];
                if (idx < ptoffs.size()) {
                    // H21: setVertexPoint(vertex_offset, point_offset)
                    gdp->setVertexPoint(ppoly->getVertexOffset(vi), ptoffs[idx]);
                }
            }
        }

    } else {
        // ── Direct triangle output (passthrough mode) ────────────────────────
        for (size_t t = 0; t < soup.triangles.size(); ++t) {
            const Triangle& tri = soup.triangles[t];

            GA_Offset pts[3];
            for (int v = 0; v < 3; ++v) {
                pts[v] = gdp->appendPoint();
                P_out.set(pts[v], UT_Vector3(
                    tri.v[v][0], tri.v[v][1], tri.v[v][2]));
            }

            GEO_PrimPoly* ppoly = GEO_PrimPoly::build(gdp, 3, GU_POLY_CLOSED, 0);
            if (!ppoly) continue;

            for (int v = 0; v < 3; ++v) {
                gdp->setVertexPoint(ppoly->getVertexOffset(v), pts[v]);
            }
        }
    }
}

//=============================================================================
// ATTRIBUTE INTERPOLATION
//=============================================================================

void interpolateVertexAttributes(GU_Detail* gdp,
                                  const EmberAttributeMapping& mapping,
                                  GA_Offset ptoff,
                                  const std::array<float, 3>& barycentric,
                                  GA_Offset src_pt0,
                                  GA_Offset src_pt1,
                                  GA_Offset src_pt2,
                                  const GEO_Primitive* /*src_prim*/)
{
    const float u = barycentric[0];
    const float v = barycentric[1];
    const float w = barycentric[2];

    auto interpolateV3 = [&](const GA_ROHandleV3& h,
                              GA_Offset p0, GA_Offset p1, GA_Offset p2) -> UT_Vector3 {
        if (!h.isValid()) return UT_Vector3(0, 0, 0);
        return h.get(p0) * u + h.get(p1) * v + h.get(p2) * w;
    };

    // Normals — interpolate + renormalize
    if (mapping.N_handle.isValid()) {
        UT_Vector3 N = interpolateV3(mapping.N_handle, src_pt0, src_pt1, src_pt2);
        N.normalize();
        GA_RWHandleV3 out(gdp->findAttribute(GA_ATTRIB_POINT, "N"));
        if (out.isValid()) out.set(ptoff, N);
    }

    // UV
    if (mapping.uv_handle.isValid()) {
        UT_Vector3 uv_val = interpolateV3(mapping.uv_handle, src_pt0, src_pt1, src_pt2);
        GA_AttributeOwner owner = mapping.uv_is_vertex
            ? GA_ATTRIB_VERTEX : GA_ATTRIB_POINT;
        GA_RWHandleV3 out(gdp->findAttribute(owner, "uv"));
        if (out.isValid()) out.set(ptoff, uv_val);
    }

    // Color
    if (mapping.Cd_handle.isValid()) {
        UT_Vector3 Cd = interpolateV3(mapping.Cd_handle, src_pt0, src_pt1, src_pt2);
        GA_RWHandleV3 out(gdp->findAttribute(GA_ATTRIB_POINT, "Cd"));
        if (out.isValid()) out.set(ptoff, Cd);
    }

    // Velocity
    if (mapping.v_handle.isValid()) {
        UT_Vector3 vel = interpolateV3(mapping.v_handle, src_pt0, src_pt1, src_pt2);
        GA_RWHandleV3 out(gdp->findAttribute(GA_ATTRIB_POINT, "v"));
        if (out.isValid()) out.set(ptoff, vel);
    }
}

//=============================================================================
// MANIFOLD OUTPUT PATH
//=============================================================================

void convertManifoldOutputToGUDetail(const PolygonSoup& soup,
                                      GU_Detail* gdp,
                                      bool is_shatter)
{
    if (!gdp) return;
    // NOTE: caller has already called gdp->clearAndDestroy()

    GA_RWHandleV3 P_out(gdp->getP());

    for (const auto& tri : soup.triangles) {
        GA_Offset ptoffs[3];
        for (int v = 0; v < 3; ++v) {
            ptoffs[v] = gdp->appendPoint();
            P_out.set(ptoffs[v], UT_Vector3(
                tri.v[v][0], tri.v[v][1], tri.v[v][2]));
        }

        GEO_PrimPoly* ppoly = GEO_PrimPoly::build(gdp, 3, GU_POLY_CLOSED, 0);
        if (ppoly) {
            for (int v = 0; v < 3; ++v) {
                gdp->setVertexPoint(ppoly->getVertexOffset(v), ptoffs[v]);
            }
        }
    }

    // Add piece attribute for Shatter
    if (is_shatter && !soup.triangles.empty()) {
        GA_RWHandleI piece_h(
            gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "piece", 1));
        if (piece_h.isValid()) {
            GA_Offset prim_off;
            int idx = 0;
            GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
                if (idx < static_cast<int>(soup.triangles.size())) {
                    piece_h.set(prim_off, soup.triangles[idx].mesh_id);
                }
                ++idx;
            }
        }
    }
}

//=============================================================================
// CONNECTIVITY NAMING (Shatter post-processing)
//=============================================================================

void assignConnectivityNames(GU_Detail* gdp, bool /*verbose*/)
{
    if (!gdp || gdp->getNumPrimitives() == 0) return;

    // Read piece attribute
    GA_ROHandleI piece_h(gdp->findIntTuple(GA_ATTRIB_PRIMITIVE, "piece", 1));
    if (!piece_h.isValid()) return;

    // Create name attribute (string, primitive)
    GA_RWHandleS name_h(gdp->addStringTuple(GA_ATTRIB_PRIMITIVE, "name", 1));
    if (!name_h.isValid()) return;

    // Build adjacency from shared points and BFS to find connected components
    const GA_Size numPrims = gdp->getNumPrimitives();

    // Map: point_offset -> list of prim indices that use it
    std::unordered_map<GA_Offset, std::vector<int>> pt_to_prims;

    // For each prim, collect its point offsets
    std::vector<std::vector<GA_Offset>> prim_points(numPrims);
    {
        int pi = 0;
        GA_Offset prim_off;
        GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
            const GA_Primitive* prim = gdp->getPrimitive(prim_off);
            if (prim) {
                GA_Range vtx_range = prim->getVertexRange();
                for (GA_Iterator it(vtx_range); !it.atEnd(); ++it) {
                    GA_Offset pt = gdp->vertexPoint(*it);
                    prim_points[pi].push_back(pt);
                    pt_to_prims[pt].push_back(pi);
                }
            }
            ++pi;
        }
    }

    // BFS connectivity
    std::vector<int> component(numPrims, -1);
    int num_components = 0;

    for (int i = 0; i < numPrims; ++i) {
        if (component[i] >= 0) continue;

        int comp_id = num_components++;
        std::queue<int> q;
        q.push(i);
        component[i] = comp_id;

        while (!q.empty()) {
            int cur = q.front(); q.pop();

            for (GA_Offset pt : prim_points[cur]) {
                for (int neighbor : pt_to_prims[pt]) {
                    if (component[neighbor] < 0) {
                        component[neighbor] = comp_id;
                        q.push(neighbor);
                    }
                }
            }
        }
    }

    // Determine majority piece_id per component (a or b origin)
    std::vector<int> comp_piece_count_a(num_components, 0);
    std::vector<int> comp_piece_count_b(num_components, 0);

    {
        int pi = 0;
        GA_Offset prim_off;
        GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
            int comp = component[pi];
            int piece = piece_h.get(prim_off);
            if (piece == 0) comp_piece_count_a[comp]++;
            else            comp_piece_count_b[comp]++;
            ++pi;
        }
    }

    // Assign names: "piece_N_a" or "piece_N_b"
    {
        int pi = 0;
        GA_Offset prim_off;
        GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
            int comp = component[pi];
            const char* suffix = (comp_piece_count_a[comp] >= comp_piece_count_b[comp])
                ? "a" : "b";
            char buf[64];
            snprintf(buf, sizeof(buf), "piece_%d_%s", comp, suffix);
            name_h.set(prim_off, buf);
            ++pi;
        }
    }
}

} // namespace ember
