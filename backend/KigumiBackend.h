#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// EMBER KigumiBackend — Boolean Operations via kigumi Library
// ═══════════════════════════════════════════════════════════════════════════════
//
// Internal conversion helpers (convertToKigumiFormat, convertFromKigumiResult)
// are defined as file-local functions in KigumiBackend.cpp to avoid exposing
// <kigumi/mesh.h> to translation units that don't link the kigumi library.

#include "IBooleanBackend.h"

namespace ember {

class KigumiBackend : public IBooleanBackend {
public:
    std::string name() const override { return "kigumi"; }
    bool supportsOpenMesh() const override { return true; }
    BackendResult execute(const PolygonSoup& soup, BooleanOp op) override;
};

} // namespace ember
