#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// EMBER MCUTBackend — Boolean Operations via MCUT Library
// ═══════════════════════════════════════════════════════════════════════════════
//
// Backend implementation using the MCUT (Mesh Cutting) library.
// MCUT performs robust mesh-mesh intersection and Boolean operations
// using a fragment-based approach.
//
// Reference: https://github.com/cutdigital/mcut

#include "IBooleanBackend.h"

namespace ember {

class MCUTBackend : public IBooleanBackend {
public:
    // Backend identification
    std::string name() const override { return "mcut"; }
    
    // Capabilities
    bool supportsOpenMesh() const override { return true; }
    
    // Main execution method
    BackendResult execute(const PolygonSoup& soup, BooleanOp op) override;

private:
    // Internal conversion helpers
    void convertToMCUTFormat(const PolygonSoup& soup,
                             std::vector<double>& meshA_vertices,
                             std::vector<uint32_t>& meshA_faces,
                             std::vector<double>& meshB_vertices,
                             std::vector<uint32_t>& meshB_faces);
    
    BackendResult convertFromMCUTResult(const std::vector<double>& vertices,
                                        const std::vector<uint32_t>& faces,
                                        BooleanOp op);
};

} // namespace ember
