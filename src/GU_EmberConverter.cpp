/**
 * @file GU_EmberConverter.cpp
 * @brief Houdini GU_Detail <-> EMBER PolygonSoup conversion
 */

#include "GU_EmberConverter.h"
#include "ember/PolygonSoup.h"
#include "ember/MeshImport.h"

#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <GA/GA_Attribute.h>
#include <GEO/GEO_Primitive.h>
#include <GEO/GEO_Face.h>        // FIX: appendVertex / close live on GEO_Face
#include <SYS/SYS_Math.h>

namespace ember {

//=============================================================================
// ATTRIBUTE MAPPING
//=============================================================================

EmberAttributeMapping buildAttributeMapping(const GU_Detail* gdp) {
    EmberAttributeMapping mapping;

    // Position (P) - always present
    mapping.P_handle = GA_ROHandleV3(gdp->getP());

    // Normal (N) - point attribute
    // FIX: findNormalAttribute returns const GA_Attribute* on a const gdp
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

    GA_Offset prim_off;
    GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
        const GEO_Primitive* prim = gdp->getGEOPrimitive(prim_off);
        if (!prim) continue;
        if (prim->getTypeId() != GA_PRIMPOLY) continue;

        const int num_verts = static_cast<int>(prim->getVertexCount());
        if (num_verts < 3) continue;

        // FIX: GEO_Primitive::getPos3(int) replaces the non-existent getPos(int)
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

            // iv[] will be filled by ember::quantizeAll() in the pipeline
            tri.iv[0] = { 0, 0, 0 };
            tri.iv[1] = { 0, 0, 0 };
            tri.iv[2] = { 0, 0, 0 };

            soup.triangles.push_back(tri);
        }
    }

    return soup;
}

//=============================================================================
// CONVERSION: PolygonSoup -> GU_Detail
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

    // Use output_polygons if populated by the backend, else fall back to
    // raw triangle list (passthrough / unprocessed soup).
    const bool has_output = !soup.output_polygons.empty();

    if (has_output) {
        // ── Indexed output from backend ──────────────────────────────────────
        // Build point array from int_vertices
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

        // Build polygons
        for (const auto& poly : soup.output_polygons) {
            if (poly.vertex_indices.size() < 3) continue;

            // FIX: cast to GEO_Face* — appendVertex / close live there,
            // not on the GEO_Primitive base class.
            GA_Primitive* raw = gdp->appendPrimitive(GA_PRIMPOLY);
            GEO_Face* face    = static_cast<GEO_Face*>(raw);
            if (!face) continue;

            for (uint32_t idx : poly.vertex_indices) {
                if (idx < ptoffs.size()) {
                    face->appendVertex(ptoffs[idx]);
                }
            }
            face->close();
        }

    } else {
        // ── Direct triangle output (passthrough mode) ────────────────────────
        for (size_t t = 0; t < soup.triangles.size(); ++t) {
            const Triangle& tri = soup.triangles[t];

            GA_Offset ptoffs[3];
            for (int v = 0; v < 3; ++v) {
                ptoffs[v] = gdp->appendPoint();
                P_out.set(ptoffs[v], UT_Vector3(
                    tri.v[v][0], tri.v[v][1], tri.v[v][2]));
            }

            GA_Primitive* raw = gdp->appendPrimitive(GA_PRIMPOLY);
            GEO_Face* face    = static_cast<GEO_Face*>(raw);
            if (!face) continue;

            for (int v = 0; v < 3; ++v) {
                face->appendVertex(ptoffs[v]);
            }
            face->close();
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
                                  const GEO_Primitive* src_prim)
{
    const float u = barycentric[0];
    const float v = barycentric[1];
    const float w = barycentric[2];

    auto interpolateV3 = [&](const GA_ROHandleV3& h,
                              GA_Offset p0, GA_Offset p1, GA_Offset p2) -> UT_Vector3 {
        if (!h.isValid()) return UT_Vector3(0, 0, 0);
        return h.get(p0) * u + h.get(p1) * v + h.get(p2) * w;
    };

    auto interpolateF = [&](const GA_ROHandleF& h,
                             GA_Offset p0, GA_Offset p1, GA_Offset p2) -> fpreal {
        if (!h.isValid()) return 0.0f;
        return h.get(p0) * u + h.get(p1) * v + h.get(p2) * w;
    };

    // Normals — interpolate + renormalize
    if (mapping.N_handle.isValid()) {
        UT_Vector3 N = interpolateV3(mapping.N_handle, src_pt0, src_pt1, src_pt2);
        N.normalize();
        GA_RWHandleV3 out(gdp->findAttribute(GA_ATTRIB_POINT, "N"));
        if (out.isValid()) out.set(ptoff, N);
    }

    // Colors
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

    // UVs
    if (mapping.uv_handle.isValid()) {
        UT_Vector3 uv_val(0, 0, 0);

        if (mapping.uv_is_vertex && src_prim != nullptr) {
            // Per-vertex UVs: read from the source primitive's vertex offsets
            GA_Offset vtx0 = src_prim->getVertexOffset(0);
            GA_Offset vtx1 = src_prim->getVertexOffset(1);
            GA_Offset vtx2 = src_prim->getVertexOffset(2);
            uv_val = mapping.uv_handle.get(vtx0) * u
                   + mapping.uv_handle.get(vtx1) * v
                   + mapping.uv_handle.get(vtx2) * w;
        } else {
            // Per-point UVs
            uv_val = interpolateV3(mapping.uv_handle, src_pt0, src_pt1, src_pt2);
        }

        GA_AttributeOwner owner = mapping.uv_is_vertex
            ? GA_ATTRIB_VERTEX : GA_ATTRIB_POINT;
        GA_RWHandleV3 out(gdp->findAttribute(owner, "uv"));
        if (out.isValid()) out.set(ptoff, uv_val);
    }

    // FIX: Removed loops over mapping.customPositionHandles,
    // mapping.customVectorHandles, mapping.customFloatHandles —
    // those fields live in SOP_EmberBoolean::AttributeState, not
    // in EmberAttributeMapping. If custom attribute propagation is
    // needed, it should be done in SOP_EmberBoolean::propagateAttributes().
}

} // namespace ember
