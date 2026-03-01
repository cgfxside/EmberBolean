/**
 * @file ManifoldBackend.h
 * @brief EMBER Boolean backend using the Manifold library (v3.x)
 *
 * Does NOT expose <manifold/manifold.h> to other translation units.
 * No GLM dependency.
 */

#pragma once

#include "IBooleanBackend.h"

namespace ember {

class ManifoldBackend : public IBooleanBackend {
public:
    std::string name() const override;
    std::string description() const override;
    bool supportsOpenMesh() const override;
    bool providesExactTopology() const override;

    BackendResult execute(const PolygonSoup& soup, BooleanOp op) override;
};

} // namespace ember
