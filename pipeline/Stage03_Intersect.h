/**
 * @file Stage03_Intersect.h
 * @brief Stage 03: Intersection Resolution for EMBER pipeline
 * 
 * This stage dispatches to the selected Boolean backend for exact
 * mesh-mesh intersection computation. It:
 * 
 * - Selects the appropriate backend based on configuration
 * - Prepares input data for the backend
 * - Executes the Boolean operation
 * - Validates and processes backend results
 * 
 * Supported backends:
 * - MCUT: Fast approximate Booleans
 * - Manifold: Robust manifold-preserving Booleans
 * - Kigumi: Exact arithmetic Booleans
 * - Cherchi: State-of-the-art exact Booleans (default)
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#pragma once

#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"
#include "backend/IBooleanBackend.h"

#include <string>
#include <memory>
#include <vector>

namespace ember {
namespace stages {

//=============================================================================
// BACKEND SELECTION
//=============================================================================

/**
 * @brief Backend selection result
 */
struct BackendSelection {
    std::string name;                    ///< Selected backend name
    std::string description;             ///< Backend description
    bool available;                      ///< Whether backend is available
    bool exact;                          ///< Whether backend provides exact results
    std::string fallback_reason;         ///< Reason for fallback (if applicable)
};

/**
 * @brief Select the best backend for the given configuration
 * 
 * Considers:
 * - User preference (config.backend)
 * - Backend availability
 * - Input characteristics
 * - Operation requirements
 * 
 * @param config_backend Requested backend name ("auto" for automatic)
 * @param soup Input polygon soup
 * @param op Boolean operation
 * @return Backend selection result
 */
BackendSelection selectBackend(const std::string& config_backend,
                                const PolygonSoup& soup,
                                BooleanOp op);

/**
 * @brief Get list of available backends
 * @return Vector of available backend names
 */
std::vector<std::string> getAvailableBackends();

/**
 * @brief Check if a backend is available
 * @param name Backend name
 * @return true if backend is available
 */
bool isBackendAvailable(const std::string& name);

//=============================================================================
// BACKEND FACTORY
//=============================================================================

/**
 * @brief Create backend instance by name
 * 
 * Factory function that creates the appropriate backend instance.
 * Returns nullptr if backend is not available.
 * 
 * @param name Backend name
 * @return Unique pointer to backend, or nullptr if unavailable
 */
std::unique_ptr<IBooleanBackend> createBackend(const std::string& name);

//=============================================================================
// INPUT PREPARATION
//=============================================================================

/**
 * @brief Prepare input data for backend processing
 * 
 * Converts the polygon soup to the format expected by the backend.
 * This may include:
 * - Vertex deduplication
 * - Index remapping
 * - Attribute extraction
 * 
 * @param soup Input polygon soup
 * @return Prepared input data (backend-specific)
 */
struct PreparedInput {
    std::vector<double> vertices;        ///< Flat vertex array [x0,y0,z0, x1,y1,z1, ...]
    std::vector<uint32_t> indices;       ///< Triangle indices [t0v0,t0v1,t0v2, t1v0,...]
    std::vector<uint32_t> mesh_offsets;  ///< Start index for each mesh in indices
    size_t num_vertices;                 ///< Total vertex count
    size_t num_triangles;                ///< Total triangle count
};

PreparedInput prepareBackendInput(const PolygonSoup& soup);

/**
 * @brief Validate prepared input
 * 
 * Checks for common issues that would cause backend failure:
 * - Empty input
 * - Degenerate triangles
 * - Invalid indices
 * 
 * @param input Prepared input data
 * @param log Diagnostic log for issues
 * @return true if input is valid
 */
bool validatePreparedInput(const PreparedInput& input, DiagnosticLog& log);

//=============================================================================
// RESULT PROCESSING
//=============================================================================

/**
 * @brief Process backend result into polygon soup format
 * 
 * Converts backend output to EMBER's internal format.
 * Handles vertex deduplication and polygon reconstruction.
 * 
 * @param backend_result Result from backend
 * @param soup Output polygon soup
 * @param log Diagnostic log
 * @return true if processing succeeded
 */
bool processBackendResult(const BackendResult& backend_result,
                          PolygonSoup& soup,
                          DiagnosticLog& log);

/**
 * @brief Validate backend result
 * 
 * Checks for issues in the backend output:
 * - Empty result when input was non-empty
 * - Invalid polygons
 * - Inconsistent winding
 * 
 * @param result Backend result
 * @param log Diagnostic log
 * @return true if result is valid
 */
bool validateBackendResult(const BackendResult& result, DiagnosticLog& log);

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

/**
 * @brief Stage 03: Intersection Resolution
 * 
 * Main entry point for exact Boolean computation.
 * 
 * Processing flow:
 *   1. Select appropriate backend
 *   2. Prepare input data
 *   3. Execute backend
 *   4. Process and validate results
 * 
 * @param soup Input polygon soup with candidate pairs
 * @param op Boolean operation to perform
 * @param backend_name Backend identifier ("auto" for automatic)
 * @return BackendResult containing output polygons
 */
BackendResult stage03_intersect(const PolygonSoup& soup, 
                                BooleanOp op, 
                                const std::string& backend_name);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Intersection stage statistics
 */
struct IntersectionStats {
    std::string backend_used;            ///< Name of backend that was used
    double backend_time_ms = 0.0;        ///< Backend execution time
    double preparation_time_ms = 0.0;    ///< Input preparation time
    double processing_time_ms = 0.0;     ///< Result processing time
    size_t input_triangles = 0;          ///< Input triangle count
    size_t output_triangles = 0;         ///< Output triangle count
    bool fallback_occurred = false;      ///< Whether fallback backend was used
    
    void print() const {
        std::printf("[Intersect] Backend: %s\n", backend_used.c_str());
        std::printf("[Intersect] Backend time: %.2f ms\n", backend_time_ms);
        std::printf("[Intersect] Input: %zu tris, Output: %zu tris\n",
                    input_triangles, output_triangles);
        if (fallback_occurred) {
            std::printf("[Intersect] Fallback backend was used\n");
        }
    }
};

/**
 * @brief Get statistics from last intersection execution
 */
const IntersectionStats& getLastIntersectionStats();

} // namespace stages
} // namespace ember
