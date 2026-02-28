/**
 * @file GU_EmberConverter.h
 * @brief Houdini GU_Detail <-> EMBER PolygonSoup conversion
 * 
 * @author EMBER Boolean Engine
 * @version 1.1.0 (Integrated)
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
 * 
 * Tracks which Houdini attributes are present and their owner class
 * (point vs vertex) for proper interpolation.
 */
struct EmberAttributeMapping {
    GA_ROHandleV3 P_handle;      // Position (always present)
    GA_ROHandleV3 N_handle;      // Normal
    GA_ROHandleV3 uv_handle;     // UV/Texture coordinates
    GA_ROHandleV3 Cd_handle;     // Diffuse color
    GA_ROHandleV3 v_handle;      // Velocity
    
    // ═══════════════════════════════════════════════════════════════════════════════
    // FIX 3.1: Track UV owner class for proper interpolation
    // ═══════════════════════════════════════════════════════════════════════════════
    bool uv_is_vertex = false;   // true if UV is GA_ATTRIB_VERTEX, false if GA_ATTRIB_POINT
};

/**
 * @brief Build attribute mapping from GU_Detail
 * 
 * Scans the input geometry for relevant attributes and creates handles.
 * FIX 3.1: Prioritizes GA_ATTRIB_VERTEX for UVs to preserve seams.
 * 
 * @param gdp  Input Houdini geometry
 * @return Attribute mapping structure
 */
EmberAttributeMapping buildAttributeMapping(const GU_Detail* gdp);

/**
 * @brief Convert GU_Detail to EMBER PolygonSoup
 * 
 * @param gdp  Input Houdini geometry
 * @param mesh_id  Mesh identifier (0 = target, 1+ = cutters)
 * @return EMBER PolygonSoup representation
 */
PolygonSoup convertGUDetailToPolygonSoup(const GU_Detail* gdp, int mesh_id);

/**
 * @brief Convert EMBER PolygonSoup to GU_Detail
 * 
 * @param soup  EMBER PolygonSoup
 * @param gdp  Output Houdini geometry (will be cleared)
 * @param mapping  Attribute mapping for interpolation
 */
void convertPolygonSoupToGUDetail(const PolygonSoup& soup, 
                                   GU_Detail* gdp,
                                   const EmberAttributeMapping& mapping);

/**
 * @brief Interpolate vertex attributes for a new point (P0 FIX: Full implementation)
 * 
 * Performs barycentric interpolation of all attributes (P, N, uv, Cd, v, rest, etc.)
 * from a source triangle to a new point created by Boolean intersection.
 * 
 * @param gdp  Output geometry
 * @param mapping  Attribute mapping from source geometry
 * @param ptoff  Point offset in output geometry
 * @param barycentric  Barycentric coordinates (u, v, w) where u+v+w=1
 * @param src_pt0  Source triangle vertex 0 (point offset)
 * @param src_pt1  Source triangle vertex 1 (point offset)
 * @param src_pt2  Source triangle vertex 2 (point offset)
 * @param src_prim  Source primitive (for vertex attributes like UVs)
 */
void interpolateVertexAttributes(GU_Detail* gdp,
                                  const EmberAttributeMapping& mapping,
                                  GA_Offset ptoff,
                                  const std::array<float, 3>& barycentric,
                                  GA_Offset src_pt0,
                                  GA_Offset src_pt1,
                                  GA_Offset src_pt2,
                                  const GEO_Primitive* src_prim = nullptr);

} // namespace ember
