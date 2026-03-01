/**
 * @file GU_EmberConverter.h
 * @brief Houdini GU_Detail <-> EMBER PolygonSoup conversion
 * 
 * @author EMBER Boolean Engine
 * @version 1.2.0 (Manifold integration)
 */

#pragma once

#include <GA/GA_Handle.h>
#include <UT/UT_Vector3.h>

// Forward declarations
class GU_Detail;
class GEO_Primitive;

namespace ember {

// Forward declarations
struct PolygonSoup;
struct Triangle;

/**
 * @brief Attribute mapping between Houdini and EMBER
 */
struct EmberAttributeMapping {
    GA_ROHandleV3 P_handle;      // Position (always present)
    GA_ROHandleV3 N_handle;      // Normal
    GA_ROHandleV3 uv_handle;     // UV/Texture coordinates
    GA_ROHandleV3 Cd_handle;     // Diffuse color
    GA_ROHandleV3 v_handle;      // Velocity
    bool uv_is_vertex = false;   // true if UV is GA_ATTRIB_VERTEX
};

/**
 * @brief Build attribute mapping from GU_Detail
 */
EmberAttributeMapping buildAttributeMapping(const GU_Detail* gdp);

/**
 * @brief Convert GU_Detail to EMBER PolygonSoup (fan-triangulates polygons)
 * @param gdp     Input Houdini geometry
 * @param mesh_id Mesh identifier (0 = target, 1+ = cutters)
 */
PolygonSoup convertGUDetailToPolygonSoup(const GU_Detail* gdp, int mesh_id);

/**
 * @brief Convert EMBER PolygonSoup to GU_Detail (standard path)
 * Uses output_polygons[] + int_vertices[] if available, else raw triangles[].
 */
void convertPolygonSoupToGUDetail(const PolygonSoup& soup, 
                                   GU_Detail* gdp,
                                   const EmberAttributeMapping& mapping);

/**
 * @brief Convert Manifold output to GU_Detail
 * Uses triangles[] with float positions written directly by ManifoldBackend.
 * NOTE: caller must have already called gdp->clearAndDestroy().
 */
void convertManifoldOutputToGUDetail(const PolygonSoup& soup,
                                      GU_Detail* gdp,
                                      bool is_shatter);

/**
 * @brief Interpolate vertex attributes for a new point
 */
void interpolateVertexAttributes(GU_Detail* gdp,
                                  const EmberAttributeMapping& mapping,
                                  GA_Offset ptoff,
                                  const std::array<float, 3>& barycentric,
                                  GA_Offset src_pt0,
                                  GA_Offset src_pt1,
                                  GA_Offset src_pt2,
                                  const GEO_Primitive* src_prim);

/**
 * @brief Assign connectivity-based names for Shatter output
 * BFS flood-fill connected components, name as "piece_N_a" or "piece_N_b"
 */
void assignConnectivityNames(GU_Detail* gdp, bool verbose = false);

} // namespace ember
