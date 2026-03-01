// ═══════════════════════════════════════════════════════════════════════════════
// ADD THESE DECLARATIONS TO THE END OF src/GU_EmberConverter.h
// (inside the ember namespace, before the closing brace)
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * @brief Convert Manifold backend output to GU_Detail
 *
 * Manifold outputs raw triangles with float positions in
 * soup.triangles[].v[][] and piece IDs in mesh_id.
 * Cap faces are encoded with CAP_FLAG (bit 31) in src_prim_idx.
 *
 * Welds vertices by exact float equality (bitwise memcmp).
 * Creates piece/name/side primitive attributes for Shatter.
 *
 * @param soup       Output soup from ManifoldBackend::execute()
 * @param gdp        Target GU_Detail (NOT cleared — caller does clearAndDestroy)
 * @param is_shatter If true, creates piece/side attributes
 */
void convertManifoldOutputToGUDetail(const PolygonSoup& soup,
                                      GU_Detail* gdp,
                                      bool is_shatter);

/**
 * @brief Assign connectivity-based names to Shatter pieces
 *
 * BFS flood-fill within each piece_id to find connected components.
 * Creates/updates "name" string prim attribute:
 *   - Single component: "piece_0", "piece_1"
 *   - Multiple components: "piece_0_a", "piece_0_b", ...
 *
 * @param gdp        Geometry with "piece" prim attribute
 * @param verbose    If true, also adds "side" attribute (0=cap, 1=shell)
 */
void assignConnectivityNames(GU_Detail* gdp, bool verbose);
