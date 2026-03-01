// ═══════════════════════════════════════════════════════════════════════════════
// APPEND THIS CODE TO THE END OF src/GU_EmberConverter.cpp
// (inside the ember namespace, before the closing brace)
// ═══════════════════════════════════════════════════════════════════════════════

// Required additional includes (add at the top of the file if not present):
// #include <unordered_map>
// #include <unordered_set>
// #include <queue>
// #include <cstring>    // memcmp

//=============================================================================
// §1  MANIFOLD OUTPUT → GU_DETAIL
//=============================================================================
//
// Manifold backend writes output triangles with FLOAT positions into
// soup.triangles[].v[][].  We weld by exact float-bit equality (memcmp)
// so shared seam vertices get ONE GA_Offset → manifold output.
//
// CAP_FLAG = 0x80000000u in src_prim_idx marks inside/cap faces.
// mesh_id stores piece_id for Shatter (0 = A-B, 1 = A∩B).

static constexpr uint32_t CAP_FLAG = 0x80000000u;

// FNV-1a hash on raw float bits for vertex dedup
struct ManifoldVtxHash {
    size_t operator()(const std::array<float,3>& p) const {
        size_t h = 14695981039346656037ULL;
        const unsigned char* raw = reinterpret_cast<const unsigned char*>(p.data());
        for (size_t i = 0; i < sizeof(float)*3; ++i) {
            h ^= raw[i];
            h *= 1099511628211ULL;
        }
        return h;
    }
};

struct ManifoldVtxEq {
    bool operator()(const std::array<float,3>& a,
                    const std::array<float,3>& b) const {
        return std::memcmp(a.data(), b.data(), sizeof(float)*3) == 0;
    }
};

void convertManifoldOutputToGUDetail(const PolygonSoup& soup,
                                      GU_Detail* gdp,
                                      bool is_shatter)
{
    if (!gdp) return;
    // NOTE: caller already did clearAndDestroy()

    if (soup.triangles.empty()) return;

    GA_RWHandleV3 P_out(gdp->getP());

    // ── Optional Shatter attributes ──────────────────────────────────────
    GA_RWHandleI piece_h;
    GA_RWHandleI side_h;

    if (is_shatter) {
        piece_h = GA_RWHandleI(
            gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "piece", 1));
        side_h  = GA_RWHandleI(
            gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "side", 1));
    }

    // ── Vertex welding by exact float equality ───────────────────────────
    std::unordered_map<std::array<float,3>, GA_Offset,
                       ManifoldVtxHash, ManifoldVtxEq> vtx_map;
    vtx_map.reserve(soup.triangles.size() * 2);  // rough estimate

    for (const auto& tri : soup.triangles) {
        GA_Offset ptoffs[3];

        for (int v = 0; v < 3; ++v) {
            std::array<float,3> key = { tri.v[v][0], tri.v[v][1], tri.v[v][2] };

            auto it = vtx_map.find(key);
            if (it != vtx_map.end()) {
                ptoffs[v] = it->second;
            } else {
                GA_Offset ptoff = gdp->appendPoint();
                P_out.set(ptoff, UT_Vector3(key[0], key[1], key[2]));
                vtx_map[key] = ptoff;
                ptoffs[v] = ptoff;
            }
        }

        // ── Build closed triangle primitive ──────────────────────────────
        GA_Primitive* raw = gdp->appendPrimitive(GA_PRIMPOLY);
        GEO_Face* face = static_cast<GEO_Face*>(raw);
        if (!face) continue;

        for (int v = 0; v < 3; ++v)
            face->appendVertex(ptoffs[v]);
        face->close();

        // ── Shatter attributes ───────────────────────────────────────────
        if (is_shatter) {
            GA_Offset prim_off = face->getMapOffset();
            piece_h.set(prim_off, tri.mesh_id);

            bool is_cap = (tri.src_prim_idx & CAP_FLAG) != 0;
            side_h.set(prim_off, is_cap ? 0 : 1);  // 0 = cap, 1 = shell
        }
    }
}

//=============================================================================
// §2  BFS CONNECTIVITY NAMES
//=============================================================================
//
// For Shatter output: flood-fill connected components within each piece_id.
// Assigns string "name" attribute:
//   piece_0, piece_1          (single component per piece)
//   piece_0_a, piece_0_b      (multiple components per piece)

void assignConnectivityNames(GU_Detail* gdp, bool /*verbose*/)
{
    if (!gdp) return;

    GA_ROHandleI piece_h(gdp->findIntTuple(GA_ATTRIB_PRIMITIVE, "piece", 1));
    if (!piece_h.isValid()) return;

    // ── Create "name" string attribute ───────────────────────────────────
    GA_RWHandleS name_h(gdp->addStringTuple(GA_ATTRIB_PRIMITIVE, "name", 1));
    if (!name_h.isValid()) return;

    const GA_IndexMap& prim_map = gdp->getPrimitiveMap();
    const GA_Size num_prims = gdp->getNumPrimitives();
    if (num_prims == 0) return;

    // ── Build point → primitive adjacency ────────────────────────────────
    // For each point, which primitives use it?
    std::unordered_map<GA_Offset, std::vector<GA_Index>> pt_to_prims;

    GA_Offset prim_off;
    GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
        const GEO_Primitive* prim = gdp->getGEOPrimitive(prim_off);
        if (!prim) continue;

        GA_Index prim_idx = prim_map.indexFromOffset(prim_off);
        for (GA_Size v = 0; v < prim->getVertexCount(); ++v) {
            GA_Offset vtx = prim->getVertexOffset(v);
            GA_Offset pt  = gdp->vertexPoint(vtx);
            pt_to_prims[pt].push_back(prim_idx);
        }
    }

    // ── BFS flood per piece ──────────────────────────────────────────────
    std::vector<int> visited(num_prims, -1);  // -1 = unvisited

    // Group primitives by piece_id
    std::unordered_map<int, std::vector<GA_Index>> piece_prims;
    GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
        GA_Index idx = prim_map.indexFromOffset(prim_off);
        int pid = piece_h.get(prim_off);
        piece_prims[pid].push_back(idx);
    }

    // For each piece, find connected components
    std::unordered_map<int, int> piece_component_count;

    for (auto& [pid, indices] : piece_prims) {
        int component = 0;

        for (GA_Index seed : indices) {
            if (visited[seed] >= 0) continue;

            // BFS from this seed, staying within same piece_id
            std::queue<GA_Index> Q;
            Q.push(seed);
            visited[seed] = component;

            while (!Q.empty()) {
                GA_Index cur = Q.front(); Q.pop();
                GA_Offset cur_off = prim_map.offsetFromIndex(cur);
                const GEO_Primitive* prim = gdp->getGEOPrimitive(cur_off);
                if (!prim) continue;

                // Find neighbors via shared points
                for (GA_Size v = 0; v < prim->getVertexCount(); ++v) {
                    GA_Offset vtx = prim->getVertexOffset(v);
                    GA_Offset pt  = gdp->vertexPoint(vtx);

                    for (GA_Index nbr : pt_to_prims[pt]) {
                        if (visited[nbr] >= 0) continue;

                        // Stay within same piece
                        GA_Offset nbr_off = prim_map.offsetFromIndex(nbr);
                        if (piece_h.get(nbr_off) != pid) continue;

                        visited[nbr] = component;
                        Q.push(nbr);
                    }
                }
            }
            component++;
        }
        piece_component_count[pid] = component;
    }

    // ── Assign names ─────────────────────────────────────────────────────
    GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
        GA_Index idx = prim_map.indexFromOffset(prim_off);
        int pid = piece_h.get(prim_off);
        int comp = visited[idx];

        char buf[64];
        if (piece_component_count[pid] <= 1) {
            // Single component: "piece_0"
            snprintf(buf, sizeof(buf), "piece_%d", pid);
        } else {
            // Multiple components: "piece_0_a", "piece_0_b"
            char suffix = 'a' + static_cast<char>(comp % 26);
            snprintf(buf, sizeof(buf), "piece_%d_%c", pid, suffix);
        }

        name_h.set(prim_off, buf);
    }
}
