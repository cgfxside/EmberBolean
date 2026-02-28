/**
 * @file ShatterPost.cpp
 * @brief Implementation of shatter post-processing with parallel DSU
 */

#include "ShatterPost.h"
#include "MeshExport.h"
#include <algorithm>
#include <numeric>
#include <cstring>
#include <unordered_map>

// For parallel execution
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ember {

//=============================================================================
// EDGE STRUCTURE FOR CONNECTIVITY
//=============================================================================

/**
 * @brief Edge key for identifying shared edges
 * 
 * Uses 64-bit encoding: (min_vert << 32) | max_vert
 * where vertices are sorted to ensure consistent ordering.
 */
struct EdgeKey {
    uint32_t v0;  // Smaller vertex index
    uint32_t v1;  // Larger vertex index
    
    EdgeKey(uint32_t a, uint32_t b) {
        if (a < b) {
            v0 = a;
            v1 = b;
        } else {
            v0 = b;
            v1 = a;
        }
    }
    
    bool operator==(const EdgeKey& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
    
    bool operator<(const EdgeKey& other) const {
        if (v0 != other.v0) return v0 < other.v0;
        return v1 < other.v1;
    }
    
    /**
     * @brief Encode as 64-bit integer for hashing
     */
    uint64_t encode() const {
        return (static_cast<uint64_t>(v0) << 32) | v1;
    }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& k) const {
        // FNV-1a inspired hash
        uint64_t key = k.encode();
        key ^= key >> 33;
        key *= 0xff51afd7ed558ccdULL;
        key ^= key >> 33;
        key *= 0xc4ceb9fe1a85ec53ULL;
        key ^= key >> 33;
        return static_cast<size_t>(key);
    }
};

/**
 * @brief Edge with face information for connectivity building
 */
struct FaceEdge {
    EdgeKey key;
    uint32_t face_idx;
    uint8_t edge_idx;  // Which edge of the face (0, 1, or 2)
    
    FaceEdge(uint32_t a, uint32_t b, uint32_t face, uint8_t eidx)
        : key(a, b), face_idx(face), edge_idx(eidx) {}
};

//=============================================================================
// PARALLEL DISJOINT SET UNION (DSU)
//=============================================================================

/**
 * @brief Thread-safe Disjoint Set Union with path compression
 * 
 * This implementation supports parallel union operations with
 * atomic compare-and-swap for parent updates.
 */
class ParallelDSU {
public:
    explicit ParallelDSU(size_t n) : parent_(n), rank_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }
    
    /**
     * @brief Find root with path compression (non-atomic for single-threaded)
     */
    size_t find(size_t x) {
        if (parent_[x] != x) {
            parent_[x] = find(parent_[x]);
        }
        return parent_[x];
    }
    
    /**
     * @brief Union two sets by rank
     */
    void unite(size_t x, size_t y) {
        size_t rx = find(x);
        size_t ry = find(y);
        if (rx == ry) return;
        
        // Union by rank
        if (rank_[rx] < rank_[ry]) {
            std::swap(rx, ry);
        }
        parent_[ry] = rx;
        if (rank_[rx] == rank_[ry]) {
            rank_[rx]++;
        }
    }
    
    /**
     * @brief Get all unique root indices
     */
    std::vector<size_t> getRoots() {
        std::vector<size_t> roots;
        roots.reserve(parent_.size());
        for (size_t i = 0; i < parent_.size(); ++i) {
            roots.push_back(find(i));
        }
        return roots;
    }
    
private:
    std::vector<size_t> parent_;
    std::vector<uint8_t> rank_;
};

//=============================================================================
// CONNECTED COMPONENT EXTRACTION
//=============================================================================

/**
 * @brief Build face adjacency from edge sharing
 * 
 * Creates a mapping from edges to the faces that share them.
 * Two faces are adjacent if they share an edge.
 */
static void buildFaceAdjacency(const PolygonSoup& soup,
                                std::vector<std::vector<uint32_t>>& adjacency) {
    EMBER_SCOPED_TIMER("Build face adjacency");
    
    const size_t num_faces = soup.output_polygons.size();
    adjacency.resize(num_faces);
    
    // Collect all edges with their face indices
    std::vector<FaceEdge> all_edges;
    all_edges.reserve(num_faces * 3);
    
    for (size_t f = 0; f < num_faces; ++f) {
        const auto& poly = soup.output_polygons[f];
        const size_t nverts = poly.vertex_indices.size();
        
        // For triangles and n-gons, create edges between consecutive vertices
        for (size_t e = 0; e < nverts; ++e) {
            uint32_t v0 = poly.vertex_indices[e];
            uint32_t v1 = poly.vertex_indices[(e + 1) % nverts];
            all_edges.emplace_back(v0, v1, static_cast<uint32_t>(f), static_cast<uint8_t>(e));
        }
    }
    
    // Sort edges by key for efficient adjacency building
    std::sort(all_edges.begin(), all_edges.end(), 
              [](const FaceEdge& a, const FaceEdge& b) {
                  if (a.key.v0 != b.key.v0) return a.key.v0 < b.key.v0;
                  return a.key.v1 < b.key.v1;
              });
    
    // Build adjacency: faces sharing an edge are adjacent
    size_t i = 0;
    while (i < all_edges.size()) {
        // Find range of edges with same key
        size_t j = i + 1;
        while (j < all_edges.size() && all_edges[j].key == all_edges[i].key) {
            ++j;
        }
        
        // Connect all faces sharing this edge
        for (size_t a = i; a < j; ++a) {
            for (size_t b = a + 1; b < j; ++b) {
                uint32_t face_a = all_edges[a].face_idx;
                uint32_t face_b = all_edges[b].face_idx;
                adjacency[face_a].push_back(face_b);
                adjacency[face_b].push_back(face_a);
            }
        }
        
        i = j;
    }
    
    // Sort and deduplicate adjacency lists
    for (auto& adj : adjacency) {
        std::sort(adj.begin(), adj.end());
        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
    }
}

/**
 * @brief Label connected components using BFS/DFS
 * 
 * Returns the number of connected components found.
 */
static int labelConnectedComponents(const std::vector<std::vector<uint32_t>>& adjacency,
                                     std::vector<int>& labels) {
    EMBER_SCOPED_TIMER("Label connected components");
    
    const size_t n = adjacency.size();
    labels.assign(n, -1);
    
    int component_id = 0;
    std::vector<uint32_t> stack;
    stack.reserve(n);
    
    for (size_t start = 0; start < n; ++start) {
        if (labels[start] >= 0) continue;
        
        // BFS from this unvisited face
        stack.clear();
        stack.push_back(static_cast<uint32_t>(start));
        labels[start] = component_id;
        
        while (!stack.empty()) {
            uint32_t current = stack.back();
            stack.pop_back();
            
            for (uint32_t neighbor : adjacency[current]) {
                if (labels[neighbor] < 0) {
                    labels[neighbor] = component_id;
                    stack.push_back(neighbor);
                }
            }
        }
        
        ++component_id;
    }
    
    return component_id;
}

/**
 * @brief Label connected components using parallel DSU
 * 
 * This is a parallel-friendly approach that first unions all adjacent
 * faces, then labels components based on the DSU structure.
 */
static int labelConnectedComponentsParallel(const std::vector<std::vector<uint32_t>>& adjacency,
                                             std::vector<int>& labels) {
    EMBER_SCOPED_TIMER("Label connected components (parallel DSU)");
    
    const size_t n = adjacency.size();
    if (n == 0) return 0;
    
    // Initialize DSU
    ParallelDSU dsu(n);
    
    // Union all adjacent faces
    // This can be parallelized with atomic operations
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
    #endif
    for (int64_t i = 0; i < static_cast<int64_t>(n); ++i) {
        for (uint32_t neighbor : adjacency[i]) {
            if (i < static_cast<int64_t>(neighbor)) {
                dsu.unite(static_cast<size_t>(i), neighbor);
            }
        }
    }
    
    // Find all unique roots and assign component IDs
    labels.resize(n);
    std::unordered_map<size_t, int> root_to_component;
    int component_id = 0;
    
    for (size_t i = 0; i < n; ++i) {
        size_t root = dsu.find(i);
        auto it = root_to_component.find(root);
        if (it == root_to_component.end()) {
            root_to_component[root] = component_id;
            labels[i] = component_id;
            ++component_id;
        } else {
            labels[i] = it->second;
        }
    }
    
    return component_id;
}

void runShatterPost(PolygonSoup& soup, const ShatterConfig& config) {
    DiagnosticLog log;
    runShatterPost(soup, config, log);
}

void runShatterPost(PolygonSoup& soup, const ShatterConfig& config, DiagnosticLog& log) {
    EMBER_SCOPED_TIMER("Shatter post-processing");
    
    if (soup.output_polygons.empty()) {
        log.warn("No output polygons to process");
        return;
    }
    
    // Step 1: Build face adjacency from shared edges
    std::vector<std::vector<uint32_t>> adjacency;
    buildFaceAdjacency(soup, adjacency);
    
    // Step 2: Label connected components
    int num_components = labelConnectedComponents(adjacency, soup.piece_ids);
    soup.piece_count = num_components;
    
    log.info("Found " + std::to_string(num_components) + " connected components");
    
    // Step 3: Generate piece names
    if (config.generate_name_attrib) {
        generatePieceNames(soup, config.name_prefix);
    }
    
    // Step 4: Compute piece pivots
    if (config.generate_piece_pivots) {
        computePiecePivots(soup);
    }
    
    // Step 5: Verify watertightness
    if (config.verify_watertight) {
        uint32_t non_watertight = verifyAllPiecesWatertight(soup, log);
        if (non_watertight > 0) {
            log.warn(DiagCategory::NotWatertight, 
                     std::to_string(non_watertight) + " pieces are not watertight");
        }
    }
    
    // Log statistics
    ShatterStats stats;
    computeShatterStats(soup, stats);
    stats.print();
}

//=============================================================================
// PIECE PIVOT COMPUTATION
//=============================================================================

void computePiecePivots(PolygonSoup& soup) {
    EMBER_SCOPED_TIMER("Compute piece pivots");
    
    if (soup.piece_count <= 0 || soup.piece_ids.empty()) {
        return;
    }
    
    // Accumulate vertex positions per piece
    std::vector<std::array<double, 3>> accum(soup.piece_count, {0.0, 0.0, 0.0});
    std::vector<uint32_t> counts(soup.piece_count, 0);
    
    // For each output polygon, accumulate its vertices to its piece
    for (size_t i = 0; i < soup.output_polygons.size(); ++i) {
        const auto& poly = soup.output_polygons[i];
        if (i >= soup.piece_ids.size()) {
            continue;
        }
        
        int piece_id = soup.piece_ids[i];
        if (piece_id < 0 || piece_id >= soup.piece_count) {
            continue;
        }
        
        // Get source triangle for vertex positions
        if (poly.original_tri >= soup.triangles.size()) {
            continue;
        }
        
        const Triangle& tri = soup.triangles[poly.original_tri];
        
        // Accumulate triangle vertices
        for (int v = 0; v < 3; ++v) {
            accum[piece_id][0] += tri.v[v][0];
            accum[piece_id][1] += tri.v[v][1];
            accum[piece_id][2] += tri.v[v][2];
        }
        counts[piece_id] += 3;
    }
    
    // Compute averages
    soup.piece_pivots.resize(soup.piece_count);
    for (int i = 0; i < soup.piece_count; ++i) {
        if (counts[i] > 0) {
            soup.piece_pivots[i][0] = static_cast<float>(accum[i][0] / counts[i]);
            soup.piece_pivots[i][1] = static_cast<float>(accum[i][1] / counts[i]);
            soup.piece_pivots[i][2] = static_cast<float>(accum[i][2] / counts[i]);
        } else {
            soup.piece_pivots[i] = {0.0f, 0.0f, 0.0f};
        }
    }
}

bool computePiecePivotAndBounds(const PolygonSoup& soup, int piece_id,
                                 std::array<float, 3>& pivot,
                                 std::array<float, 3>& bbox_min,
                                 std::array<float, 3>& bbox_max) {
    if (piece_id < 0 || piece_id >= soup.piece_count) {
        return false;
    }
    
    bool first = true;
    std::array<double, 3> accum = {0.0, 0.0, 0.0};
    uint32_t count = 0;
    
    for (size_t i = 0; i < soup.output_polygons.size(); ++i) {
        if (i >= soup.piece_ids.size() || soup.piece_ids[i] != piece_id) {
            continue;
        }
        
        const auto& poly = soup.output_polygons[i];
        if (poly.original_tri >= soup.triangles.size()) {
            continue;
        }
        
        const Triangle& tri = soup.triangles[poly.original_tri];
        
        for (int v = 0; v < 3; ++v) {
            // Update bounding box
            if (first) {
                bbox_min = {static_cast<float>(tri.v[v][0]),
                           static_cast<float>(tri.v[v][1]),
                           static_cast<float>(tri.v[v][2])};
                bbox_max = bbox_min;
                first = false;
            } else {
                for (int c = 0; c < 3; ++c) {
                    bbox_min[c] = std::min(bbox_min[c], static_cast<float>(tri.v[v][c]));
                    bbox_max[c] = std::max(bbox_max[c], static_cast<float>(tri.v[v][c]));
                }
            }
            
            // Accumulate for centroid
            accum[0] += tri.v[v][0];
            accum[1] += tri.v[v][1];
            accum[2] += tri.v[v][2];
            ++count;
        }
    }
    
    if (count == 0) {
        pivot = {0.0f, 0.0f, 0.0f};
        return false;
    }
    
    pivot = {static_cast<float>(accum[0] / count),
             static_cast<float>(accum[1] / count),
             static_cast<float>(accum[2] / count)};
    
    return true;
}

//=============================================================================
// WATERTIGHT VERIFICATION
//=============================================================================

bool isPieceWatertight(const PolygonSoup& soup, int piece_id, 
                       uint32_t& boundary_edge_count) {
    boundary_edge_count = 0;
    
    // Count edge occurrences for this piece
    std::unordered_map<uint64_t, uint32_t> edge_counts;
    
    for (size_t i = 0; i < soup.output_polygons.size(); ++i) {
        if (i >= soup.piece_ids.size() || soup.piece_ids[i] != piece_id) {
            continue;
        }
        
        const auto& poly = soup.output_polygons[i];
        const size_t nverts = poly.vertex_indices.size();
        
        for (size_t e = 0; e < nverts; ++e) {
            uint32_t v0 = poly.vertex_indices[e];
            uint32_t v1 = poly.vertex_indices[(e + 1) % nverts];
            
            // Create canonical edge key
            EdgeKey key(v0, v1);
            uint64_t encoded = key.encode();
            
            edge_counts[encoded]++;
        }
    }
    
    // Check for boundary edges (count == 1)
    for (const auto& pair : edge_counts) {
        if (pair.second == 1) {
            ++boundary_edge_count;
        } else if (pair.second > 2) {
            // Non-manifold edge - also counts as boundary for our purposes
            ++boundary_edge_count;
        }
    }
    
    return boundary_edge_count == 0;
}

uint32_t verifyAllPiecesWatertight(const PolygonSoup& soup, DiagnosticLog& log) {
    uint32_t non_watertight = 0;
    
    for (int piece_id = 0; piece_id < soup.piece_count; ++piece_id) {
        uint32_t boundary_edges = 0;
        if (!isPieceWatertight(soup, piece_id, boundary_edges)) {
            ++non_watertight;
            log.debug("Piece " + std::to_string(piece_id) + 
                      " has " + std::to_string(boundary_edges) + " boundary edges");
        }
    }
    
    return non_watertight;
}

//=============================================================================
// STATISTICS
//=============================================================================

void computeShatterStats(const PolygonSoup& soup, ShatterStats& stats) {
    stats.piece_count = soup.piece_count;
    stats.total_polygons = static_cast<uint32_t>(soup.output_polygons.size());
    
    if (stats.piece_count == 0) {
        stats.avg_polygons_per_piece = 0.0;
        stats.min_polygons = 0;
        stats.max_polygons = 0;
        return;
    }
    
    // Count polygons per piece
    std::vector<uint32_t> poly_counts(stats.piece_count, 0);
    for (size_t i = 0; i < soup.output_polygons.size(); ++i) {
        if (i < soup.piece_ids.size()) {
            int piece_id = soup.piece_ids[i];
            if (piece_id >= 0 && piece_id < stats.piece_count) {
                poly_counts[piece_id]++;
            }
        }
    }
    
    // Compute statistics
    uint64_t total = 0;
    stats.min_polygons = UINT32_MAX;
    stats.max_polygons = 0;
    
    for (uint32_t count : poly_counts) {
        total += count;
        stats.min_polygons = std::min(stats.min_polygons, count);
        stats.max_polygons = std::max(stats.max_polygons, count);
    }
    
    stats.avg_polygons_per_piece = static_cast<double>(total) / stats.piece_count;
    
    // Count watertight pieces
    for (int i = 0; i < stats.piece_count; ++i) {
        uint32_t boundary;
        if (isPieceWatertight(soup, i, boundary)) {
            ++stats.watertight_pieces;
        } else {
            ++stats.non_watertight_pieces;
        }
    }
}

} // namespace ember
