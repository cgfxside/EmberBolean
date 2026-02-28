/**
 * @file IBooleanBackend.h
 * @brief Abstract interface for Boolean operation backends
 */

#pragma once

#include <string>
#include <memory>
#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"

namespace ember {

// Forward declarations
struct PolygonSoup;
class DiagnosticLog;

/**
 * @brief Result from a Boolean backend execution
 */
struct BackendResult {
    PolygonSoup output_soup;
    bool success = false;
    std::string error_message;
    DiagnosticLog diagnostics;

    // Execution statistics
    uint32_t num_vertices = 0;
    uint32_t num_triangles = 0;
    double execution_time_ms = 0.0;

    BackendResult() = default;
};

/**
 * @brief Abstract interface for Boolean operation backends
 *
 * Implementations:
 * - CherchiBackend: Exact CDT-based (recommended)
 * - MCUTBackend: MCUT library wrapper
 * - ManifoldBackend: manifold library wrapper
 * - KigumiBackend: kigumi library wrapper
 */
class IBooleanBackend {
public:
    virtual ~IBooleanBackend() = default;

    /// Get backend name
    virtual std::string name() const = 0;

    /// Get backend description (optional, defaults to name)
    virtual std::string description() const { return name(); }

    /// Execute Boolean operation
    virtual BackendResult execute(const PolygonSoup& soup, BooleanOp op) = 0;

    /// Does this backend support open (non-watertight) meshes?
    virtual bool supportsOpenMesh() const { return false; }

    /// Does this backend provide exact topology?
    virtual bool providesExactTopology() const { return false; }
};

/**
 * @brief Factory function to create backends by name
 *
 * Supported names: "mcut", "manifold", "kigumi", "cherchi", "auto"
 */
std::unique_ptr<IBooleanBackend> createBackend(const std::string& name);

} // namespace ember
