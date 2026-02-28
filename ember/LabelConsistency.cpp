/**
 * @file LabelConsistency.cpp
 * @brief Implementation of Sheng-style label propagation
 */

#include "LabelConsistency.h"
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace ember {

//=============================================================================
// WINDING NUMBER COMPUTATION
//=============================================================================

/**
 * @brief Compute solid angle of a triangle from a point
 * 
 * The solid angle is used for winding number computation.
 * Uses the formula from Jacobson et al. "Robust Inside-Outside Segmentation"
 */
static double computeSolidAngle(const double p[3],
                                 const double a[3],
                                 const double b[3],
                                 const double c[3]) {
    // Vectors from p to triangle vertices
    double pa[3] = {a[0] - p[0], a[1] - p[1], a[2] - p[2]};
    double pb[3] = {b[0] - p[0], b[1] - p[1], b[2] - p[2]};
    double pc[3] = {c[0] - p[0], c[1] - p[1], c[2] - p[2]};
    
    // Normalize vectors
    auto normalize = [](double v[3]) {
        double len = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (len > 1e-10) {
            v[0] /= len;
            v[1] /= len;
            v[2] /= len;
        }
    };
    
    normalize(pa);
    normalize(pb);
    normalize(pc);
    
    // Compute scalar triple product
    double det = pa[0] * (pb[1] * pc[2] - pb[2] * pc[1]) -
                 pa[1] * (pb[0] * pc[2] - pb[2] * pc[0]) +
                 pa[2] * (pb[0] * pc[1] - pb[1] * pc[0]);
    
    // Compute dot products for denominator
    double dot_ab = pa[0]*pb[0] + pa[1]*pb[1] + pa[2]*pb[2];
    double dot_bc = pb[0]*pc[0] + pb[1]*pc[1] + pb[2]*pc[2];
    double dot_ca = pc[0]*pa[0] + pc[1]*pa[1] + pc[2]*pa[2];
    
    double denom = 1.0 + dot_ab + dot_bc + dot_ca;
    
    if (std::abs(denom) < 1e-10) {
        return 0.0;
    }
    
    // Solid angle = 2 * atan2(det, denom)
    return 2.0 * std::atan2(det, denom);
}

WindingClassification computeWindingNumber(const PolygonSoup& soup,
                                            int mesh_id,
                                            const double point[3]) {
    double total_angle = 0.0;
    
    // Get triangle range for this mesh
    auto range = soup.getMeshTriangleRange(mesh_id);
    
    for (size_t i = range.first; i < range.second; ++i) {
        const Triangle& tri = soup.triangles[i];
        
        double a[3] = {tri.v[0][0], tri.v[0][1], tri.v[0][2]};
        double b[3] = {tri.v[1][0], tri.v[1][1], tri.v[1][2]};
        double c[3] = {tri.v[2][0], tri.v[2][1], tri.v[2][2]};
        
        total_angle += computeSolidAngle(point, a, b, c);
    }
    
    // Winding number = total_angle / (4*pi)
    int wn = static_cast<int>(std::round(total_angle / (4.0 * M_PI)));
    
    return WindingClassification::fromWindingNumber(wn);
}

WindingClassification computeWindingNumberForTriangle(const PolygonSoup& soup,
                                                       int mesh_id,
                                                       uint32_t tri_idx) {
    if (tri_idx >= soup.triangles.size()) {
        return WindingClassification::fromWindingNumber(0);
    }
    
    const Triangle& tri = soup.triangles[tri_idx];
    
    // Use centroid as query point
    double centroid[3] = {
        (tri.v[0][0] + tri.v[1][0] + tri.v[2][0]) / 3.0,
        (tri.v[0][1] + tri.v[1][1] + tri.v[2][1]) / 3.0,
        (tri.v[0][2] + tri.v[1][2] + tri.v[2][2]) / 3.0
    };
    
    return computeWindingNumber(soup, mesh_id, centroid);
}

//=============================================================================
// FACE ADJACENCY COMPUTATION
//=============================================================================

/**
 * @brief Edge key for adjacency building
 */
struct EdgeKeyLC {
    uint32_t v0, v1;
    
    EdgeKeyLC(uint32_t a, uint32_t b) {
        if (a < b) {
            v0 = a;
            v1 = b;
        } else {
            v0 = b;
            v1 = a;
        }
    }
    
    bool operator==(const EdgeKeyLC& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
    
    bool operator<(const EdgeKeyLC& other) const {
        if (v0 != other.v0) return v0 < other.v0;
        return v1 < other.v1;
    }
};

struct EdgeKeyLCHash {
    size_t operator()(const EdgeKeyLC& k) const {
        return (static_cast<uint64_t>(k.v0) << 32) | k.v1;
    }
};

void computeFaceAdjacencyForPropagation(const PolygonSoup& soup,
                                        std::vector<std::vector<uint32_t>>& adjacency) {
    EMBER_SCOPED_TIMER("Compute face adjacency");
    
    const size_t num_faces = soup.output_polygons.size();
    adjacency.resize(num_faces);
    
    // Map edges to faces
    std::unordered_map<EdgeKeyLC, std::vector<uint32_t>, EdgeKeyLCHash> edge_faces;
    
    for (size_t f = 0; f < num_faces; ++f) {
        const auto& poly = soup.output_polygons[f];
        const size_t nverts = poly.vertex_indices.size();
        
        for (size_t e = 0; e < nverts; ++e) {
            uint32_t v0 = poly.vertex_indices[e];
            uint32_t v1 = poly.vertex_indices[(e + 1) % nverts];
            
            EdgeKeyLC key(v0, v1);
            edge_faces[key].push_back(static_cast<uint32_t>(f));
        }
    }
    
    // Build adjacency: faces sharing an edge are adjacent
    for (const auto& pair : edge_faces) {
        const auto& faces = pair.second;
        
        // Skip seam edges (boundary or intersection edges)
        // For now, we assume edges shared by exactly 2 faces are interior
        if (faces.size() != 2) {
            continue;
        }
        
        uint32_t f0 = faces[0];
        uint32_t f1 = faces[1];
        
        adjacency[f0].push_back(f1);
        adjacency[f1].push_back(f0);
    }
    
    // Sort and deduplicate
    for (auto& adj : adjacency) {
        std::sort(adj.begin(), adj.end());
        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
    }
}

//=============================================================================
// LABEL PROPAGATION
//=============================================================================

LabelResult propagateLabels(const PolygonSoup& soup) {
    DiagnosticLog log;
    return propagateLabels(soup, log);
}

LabelResult propagateLabels(const PolygonSoup& soup, DiagnosticLog& log) {
    EMBER_SCOPED_TIMER("Propagate labels");
    
    LabelResult result;
    const size_t num_faces = soup.output_polygons.size();
    
    if (num_faces == 0) {
        return result;
    }
    
    result.face_labels.assign(num_faces, -1);
    
    // Build face adjacency
    std::vector<std::vector<uint32_t>> adjacency;
    computeFaceAdjacencyForPropagation(soup, adjacency);
    
    // Select seed faces
    std::vector<uint32_t> seeds;
    selectSeedFaces(soup, SeedStrategy::SpatialGrid, 
                    static_cast<int>(std::min(size_t(100), num_faces / 10 + 1)), 
                    seeds);
    
    result.num_seeds = static_cast<int>(seeds.size());
    
    // Label seeds and propagate
    int next_label = 0;
    for (uint32_t seed : seeds) {
        if (result.face_labels[seed] >= 0) {
            continue;  // Already labeled
        }
        
        // Compute winding number for seed
        const auto& poly = soup.output_polygons[seed];
        int mesh_to_test = (poly.mesh_id == 0) ? 1 : 0;  // Test against other mesh
        
        WindingClassification wc = computeWindingNumberForTriangle(soup, mesh_to_test, 
                                                                    poly.original_tri);
        
        int label = wc.is_inside ? next_label++ : next_label++;
        
        // BFS propagate
        int propagated = propagateLabelBFS(result.face_labels, static_cast<int>(seed),
                                           label, adjacency);
        result.num_propagated += propagated;
        
        log.debug("Label " + std::to_string(label) + " propagated to " + 
                  std::to_string(propagated) + " faces");
    }
    
    result.num_labels = next_label;
    
    // Handle unlabeled faces (compute directly)
    int unlabeled = 0;
    for (size_t i = 0; i < num_faces; ++i) {
        if (result.face_labels[i] < 0) {
            ++unlabeled;
            
            const auto& poly = soup.output_polygons[i];
            int mesh_to_test = (poly.mesh_id == 0) ? 1 : 0;
            
            WindingClassification wc = computeWindingNumberForTriangle(soup, mesh_to_test,
                                                                        poly.original_tri);
            
            // Assign to existing label or create new
            result.face_labels[i] = wc.is_inside ? 0 : 1;
        }
    }
    
    log.info("Label propagation: " + std::to_string(result.num_propagated) + 
             " faces propagated, " + std::to_string(unlabeled) + 
             " faces computed directly, " + std::to_string(result.num_labels) + 
             " labels total");
    
    return result;
}

int propagateLabelBFS(std::vector<int>& labels, int start_face, 
                      int label, const PolygonSoup& soup) {
    std::vector<std::vector<uint32_t>> adjacency;
    computeFaceAdjacencyForPropagation(soup, adjacency);
    return propagateLabelBFS(labels, start_face, label, adjacency);
}

int propagateLabelBFS(std::vector<int>& labels, int start_face, 
                      int label, const std::vector<std::vector<uint32_t>>& adjacency) {
    if (start_face < 0 || start_face >= static_cast<int>(labels.size())) {
        return 0;
    }
    
    if (labels[start_face] >= 0) {
        return 0;  // Already labeled
    }
    
    int count = 0;
    std::queue<int> queue;
    
    labels[start_face] = label;
    queue.push(start_face);
    ++count;
    
    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();
        
        if (current >= static_cast<int>(adjacency.size())) {
            continue;
        }
        
        for (uint32_t neighbor : adjacency[current]) {
            if (labels[neighbor] < 0) {
                labels[neighbor] = label;
                queue.push(static_cast<int>(neighbor));
                ++count;
            }
        }
    }
    
    return count;
}

//=============================================================================
// SEED SELECTION
//=============================================================================

void selectSeedFaces(const PolygonSoup& soup,
                     SeedStrategy strategy,
                     int num_seeds,
                     std::vector<uint32_t>& seeds) {
    seeds.clear();
    
    const size_t num_faces = soup.output_polygons.size();
    if (num_faces == 0 || num_seeds <= 0) {
        return;
    }
    
    num_seeds = static_cast<int>(std::min(static_cast<size_t>(num_seeds), num_faces));
    
    switch (strategy) {
        case SeedStrategy::Random: {
            // Simple deterministic "random" selection
            size_t step = num_faces / num_seeds;
            for (int i = 0; i < num_seeds; ++i) {
                seeds.push_back(static_cast<uint32_t>(i * step));
            }
            break;
        }
        
        case SeedStrategy::SpatialGrid: {
            // Compute bounding box
            double min[3] = {1e300, 1e300, 1e300};
            double max[3] = {-1e300, -1e300, -1e300};
            
            for (const auto& poly : soup.output_polygons) {
                if (poly.original_tri >= soup.triangles.size()) continue;
                
                const Triangle& tri = soup.triangles[poly.original_tri];
                for (int v = 0; v < 3; ++v) {
                    for (int c = 0; c < 3; ++c) {
                        min[c] = std::min(min[c], static_cast<double>(tri.v[v][c]));
                        max[c] = std::max(max[c], static_cast<double>(tri.v[v][c]));
                    }
                }
            }
            
            // Grid-based selection
            int grid_res = static_cast<int>(std::ceil(std::pow(num_seeds, 1.0/3.0)));
            
            for (int i = 0; i < grid_res && static_cast<int>(seeds.size()) < num_seeds; ++i) {
                for (int j = 0; j < grid_res && static_cast<int>(seeds.size()) < num_seeds; ++j) {
                    for (int k = 0; k < grid_res && static_cast<int>(seeds.size()) < num_seeds; ++k) {
                        // Find face closest to grid point
                        double gx = min[0] + (max[0] - min[0]) * (i + 0.5) / grid_res;
                        double gy = min[1] + (max[1] - min[1]) * (j + 0.5) / grid_res;
                        double gz = min[2] + (max[2] - min[2]) * (k + 0.5) / grid_res;
                        
                        uint32_t closest = 0;
                        double closest_dist = 1e300;
                        
                        for (size_t f = 0; f < num_faces; ++f) {
                            const auto& poly = soup.output_polygons[f];
                            if (poly.original_tri >= soup.triangles.size()) continue;
                            
                            const Triangle& tri = soup.triangles[poly.original_tri];
                            double cx = (tri.v[0][0] + tri.v[1][0] + tri.v[2][0]) / 3.0;
                            double cy = (tri.v[0][1] + tri.v[1][1] + tri.v[2][1]) / 3.0;
                            double cz = (tri.v[0][2] + tri.v[1][2] + tri.v[2][2]) / 3.0;
                            
                            double dist = (cx-gx)*(cx-gx) + (cy-gy)*(cy-gy) + (cz-gz)*(cz-gz);
                            if (dist < closest_dist) {
                                closest_dist = dist;
                                closest = static_cast<uint32_t>(f);
                            }
                        }
                        
                        seeds.push_back(closest);
                    }
                }
            }
            break;
        }
        
        case SeedStrategy::Extrema: {
            // Select faces at bounding box extrema
            double min[3] = {1e300, 1e300, 1e300};
            double max[3] = {-1e300, -1e300, -1e300};
            uint32_t min_face[3] = {0, 0, 0};
            uint32_t max_face[3] = {0, 0, 0};
            
            for (size_t f = 0; f < num_faces; ++f) {
                const auto& poly = soup.output_polygons[f];
                if (poly.original_tri >= soup.triangles.size()) continue;
                
                const Triangle& tri = soup.triangles[poly.original_tri];
                double cx = (tri.v[0][0] + tri.v[1][0] + tri.v[2][0]) / 3.0;
                double cy = (tri.v[0][1] + tri.v[1][1] + tri.v[2][1]) / 3.0;
                double cz = (tri.v[0][2] + tri.v[1][2] + tri.v[2][2]) / 3.0;
                
                for (int c = 0; c < 3; ++c) {
                    double coord = (c == 0) ? cx : (c == 1) ? cy : cz;
                    if (coord < min[c]) {
                        min[c] = coord;
                        min_face[c] = static_cast<uint32_t>(f);
                    }
                    if (coord > max[c]) {
                        max[c] = coord;
                        max_face[c] = static_cast<uint32_t>(f);
                    }
                }
            }
            
            for (int c = 0; c < 3; ++c) {
                seeds.push_back(min_face[c]);
                seeds.push_back(max_face[c]);
            }
            break;
        }
        
        case SeedStrategy::LargestComponent: {
            // For now, fall back to random
            size_t step = num_faces / num_seeds;
            for (int i = 0; i < num_seeds; ++i) {
                seeds.push_back(static_cast<uint32_t>(i * step));
            }
            break;
        }
    }
    
    // Remove duplicates
    std::sort(seeds.begin(), seeds.end());
    seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
}

//=============================================================================
// BOOLEAN CLASSIFICATION
//=============================================================================

int classifyForBoolean(const PolygonSoup& soup,
                       BooleanOp op,
                       std::vector<BooleanClass>& classifications) {
    EMBER_SCOPED_TIMER("Classify for Boolean");
    
    classifications.clear();
    
    // Get labels
    LabelResult labels = propagateLabels(soup);
    
    const size_t num_faces = soup.output_polygons.size();
    classifications.resize(num_faces, BooleanClass::Unknown);
    
    int classified = 0;
    
    for (size_t i = 0; i < num_faces; ++i) {
        const auto& poly = soup.output_polygons[i];
        
        // Determine inside/outside based on label and mesh
        bool is_inside = (labels.getLabel(i) % 2 == 0);  // Simplified
        
        BooleanClass cls = BooleanClass::Unknown;
        
        if (poly.mesh_id == 0) {
            // Mesh A
            cls = is_inside ? BooleanClass::A_in_B : BooleanClass::A_out_B;
        } else {
            // Mesh B
            cls = is_inside ? BooleanClass::B_in_A : BooleanClass::B_out_A;
        }
        
        classifications[i] = cls;
        ++classified;
    }
    
    return classified;
}

void filterByClassification(const PolygonSoup& soup,
                            const std::vector<BooleanClass>& classifications,
                            const std::vector<BooleanClass>& keep_class,
                            std::vector<uint32_t>& filtered_indices) {
    filtered_indices.clear();
    
    const size_t num_faces = soup.output_polygons.size();
    if (classifications.size() != num_faces) {
        return;
    }
    
    for (size_t i = 0; i < num_faces; ++i) {
        for (BooleanClass kc : keep_class) {
            if (classifications[i] == kc) {
                filtered_indices.push_back(static_cast<uint32_t>(i));
                break;
            }
        }
    }
}

//=============================================================================
// STATISTICS
//=============================================================================

void computeLabelStats(const LabelResult& result, LabelStats& stats) {
    stats.total_faces = static_cast<int>(result.face_labels.size());
    stats.num_labels = result.num_labels;
    stats.num_seeds = result.num_seeds;
    
    // Count labeled faces
    stats.labeled_faces = 0;
    for (int label : result.face_labels) {
        if (label >= 0) {
            ++stats.labeled_faces;
        }
    }
    
    stats.unlabeled_faces = stats.total_faces - stats.labeled_faces;
    
    if (stats.total_faces > 0) {
        stats.propagation_ratio = static_cast<double>(stats.labeled_faces) / stats.total_faces;
    }
    
    if (result.num_labels > 0) {
        stats.avg_component_size = static_cast<double>(stats.labeled_faces) / result.num_labels;
    }
}

} // namespace ember
