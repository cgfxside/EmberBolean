#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// EMBER ManifoldBackend — Boolean Operations via Manifold Library
// ═══════════════════════════════════════════════════════════════════════════════
//
// Internal conversion helpers (ManifoldMesh struct, convertToManifoldFormat,
// convertFromManifoldResult) are defined as file-local functions in
// ManifoldBackend.cpp to avoid exposing <manifold/manifold.h> and
// <glm/glm.hpp> to translation units that don't link those libraries.

#include "IBooleanBackend.h"

namespace ember {

class ManifoldBackend : public IBooleanBackend {
public:
    std::string name() const override { return "manifold"; }
    bool supportsOpenMesh() const override { return false; }
    BackendResult execute(const PolygonSoup& soup, BooleanOp op) override;
};

} // namespace ember
