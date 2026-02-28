/**
 * @file EmberPipeline.h
 * @brief Main pipeline orchestrator for the EMBER Boolean engine
 * 
 * This file defines the main pipeline entry point and configuration structures
 * for the EMBER Boolean engine. The pipeline consists of 6 stages:
 * 
 *   Stage 01: Diagnostic Gatekeeper - Input validation and triangulation
 *   Stage 02: Embree Broad-Phase - BVH-based candidate pair detection
 *   Stage 02b: Plane Fast Path - Analytical plane classification for Tier 1 cutters
 *   Stage 03: Intersection Resolution - Backend Boolean computation
 *   Stage 04: Classification - Winding number evaluation and label propagation
 *   Stage 05: Output Reconstruction - Build result mesh with attributes
 *   Stage 06: Shatter Post-Processing - Connected component extraction
 * 
 * The pipeline supports both exact mesh-mesh Booleans (Tier 2) and fast
 * analytical plane cutting (Tier 1) through automatic cutter dispatch.
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#pragma once

#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"
#include "backend/IBooleanBackend.h"

#include <string>
#include <chrono>

namespace ember {

// Forward declarations
struct DispatchResult;
struct CutterDescriptor;

//=============================================================================
// PIPELINE CONFIGURATION
//=============================================================================

/**
 * @brief Configuration parameters for the EMBER pipeline
 * 
 * These parameters control the behavior of the Boolean operation,
 * including noise generation, grout width, and backend selection.
 */
struct PipelineConfig {
    // Boolean operation type
    BooleanOp operation = BooleanOp::Union;
    
    // Backend selection ("auto", "mcut", "manifold", "kigumi", "cherchi")
    std::string backend = "auto";
    
    // Noise generation parameters (for Tier 2 cutters)
    float noise_amplitude = 0.0f;    ///< Maximum displacement amplitude
    float noise_frequency = 1.0f;    ///< Spatial frequency of noise
    int noise_octaves = 3;           ///< Number of octaves for fBm
    uint32_t noise_seed = 0;         ///< Random seed for reproducibility
    
    // Grout parameters (for Tier 1G planar cutters)
    float grout_width = 0.0f;        ///< Width of grout region around plane
    
    // Shell noise parameters
    float shell_noise = 0.0f;        ///< Amplitude for shell-based noise
    
    // Tessellation parameters (for noised cutters)
    float segment_size = 0.1f;       ///< Target segment length for tessellation
    int max_segments = 64;           ///< Maximum segments per cutter edge
    
    // Pipeline control flags
    bool cull_external = true;       ///< Remove triangles outside result
    bool force_exact = false;        ///< Force Tier 2 even for planar cutters
    
    // ═══════════════════════════════════════════════════════════════════════════════
    // P2 FIX: Configurable aspect ratio threshold for quantization
    // ═══════════════════════════════════════════════════════════════════════════════
    double aspect_ratio_threshold = 1000.0;  ///< Threshold for per-axis quantization
    
    /**
     * @brief Check if noise is enabled
     */
    bool hasNoise() const { return noise_amplitude > 0.0f; }
    
    /**
     * @brief Check if grout is enabled
     */
    bool hasGrout() const { return grout_width > 0.0f; }
    
    /**
     * @brief Check if this is a shatter operation
     */
    bool isShatter() const { return operation == BooleanOp::Shatter; }
};

//=============================================================================
// PIPELINE RESULT
//=============================================================================

/**
 * @brief Result structure from pipeline execution
 * 
 * Contains the output mesh, success status, error information,
 * and diagnostic log from the entire pipeline execution.
 */
struct PipelineResult {
    PolygonSoup output;              ///< Result polygon soup
    bool success = false;            ///< True if operation succeeded
    std::string error_message;       ///< Error description if failed
    DiagnosticLog diagnostics;       ///< Full diagnostic log
    double execution_time_ms = 0.0;  ///< Total execution time
    
    // Per-stage timing breakdown
    double stage01_time_ms = 0.0;    ///< Diagnostic stage time
    double stage02_time_ms = 0.0;    ///< Broad-phase time
    double stage03_time_ms = 0.0;    ///< Intersection time
    double stage04_time_ms = 0.0;    ///< Classification time
    double stage05_time_ms = 0.0;    ///< Reconstruction time
    double stage06_time_ms = 0.0;    ///< Shatter post-processing time
    
    /**
     * @brief Check if result contains valid output
     */
    bool hasOutput() const { return success && !output.output_polygons.empty(); }
    
    /**
     * @brief Get total triangle count in output
     */
    size_t outputTriangleCount() const { return output.output_polygons.size(); }
    
    /**
     * @brief Print timing breakdown
     */
    void printTiming() const {
        std::printf("[EMBER Pipeline] Total: %.2f ms\n", execution_time_ms);
        std::printf("  Stage 01 (Diagnostic):   %.2f ms\n", stage01_time_ms);
        std::printf("  Stage 02 (Broad-Phase):  %.2f ms\n", stage02_time_ms);
        std::printf("  Stage 03 (Intersect):    %.2f ms\n", stage03_time_ms);
        std::printf("  Stage 04 (Classify):     %.2f ms\n", stage04_time_ms);
        std::printf("  Stage 05 (Reconstruct):  %.2f ms\n", stage05_time_ms);
        if (stage06_time_ms > 0.0) {
            std::printf("  Stage 06 (Shatter):      %.2f ms\n", stage06_time_ms);
        }
    }
};

//=============================================================================
// MAIN PIPELINE ENTRY POINT
//=============================================================================

/**
 * @brief Execute the complete EMBER Boolean pipeline
 * 
 * This is the main entry point for Boolean operations. It orchestrates
 * all pipeline stages based on the input configuration and automatically
 * selects between Tier 1 (plane fast-path) and Tier 2 (exact mesh-mesh)
 * processing based on cutter analysis.
 * 
 * @param input Input polygon soup (target mesh + cutters)
 * @param config Pipeline configuration parameters
 * @return PipelineResult containing output mesh and diagnostics
 * 
 * Example usage:
 * @code
 *   PipelineConfig config;
 *   config.operation = BooleanOp::Difference;
 *   config.backend = "cherchi";
 *   
 *   PipelineResult result = runPipeline(input_soup, config);
 *   if (result.success) {
 *       // Use result.output
 *   }
 * @endcode
 */
PipelineResult runPipeline(const PolygonSoup& input, const PipelineConfig& config);

//=============================================================================
// INDIVIDUAL STAGES (Advanced Usage)
//=============================================================================

/**
 * @brief Individual pipeline stages for advanced usage
 * 
 * These functions allow direct access to individual pipeline stages
 * for debugging, testing, or custom pipeline construction.
 */
namespace stages {

/**
 * @brief Stage 01: Diagnostic Gatekeeper
 * 
 * Validates and prepares input mesh:
 * - Triangulates all polygons
 * - Removes degenerate triangles
 * - Normalizes winding order
 * - Checks watertightness
 * - Performs cutter type detection (Tier 1 vs Tier 2)
 * 
 * @param soup Input/output polygon soup (modified in place)
 * @return Diagnostic log with any issues found
 */
DiagnosticLog stage01_diagnostic(PolygonSoup& soup, const PipelineConfig& config, 
                                  DispatchResult& dispatch);

/**
 * @brief Stage 02: Embree Broad-Phase
 * 
 * Builds BVH using Embree and finds candidate triangle pairs:
 * - Constructs BVH over all triangles
 * - Performs self-collision queries for target-cutter pairs
 * - Stores candidate pairs in soup.candidate_pairs
 * 
 * @param soup Input/output polygon soup
 */
void stage02_broadPhase(PolygonSoup& soup);

/**
 * @brief Stage 02b: Plane Fast Path
 * 
 * Analytical plane classification for Tier 1 planar cutters:
 * - Skips BVH construction
 * - Classifies each triangle against plane equations
 * - Splits straddling triangles
 * - Produces output directly for simple cases
 * 
 * @param soup Input/output polygon soup
 * @param config Pipeline configuration
 * @param dispatch Cutter dispatch results with plane equations
 */
void stage02b_planeFastPath(PolygonSoup& soup, const PipelineConfig& config,
                            const DispatchResult& dispatch);

/**
 * @brief Stage 03: Intersection Resolution
 * 
 * Dispatches to selected backend for exact Boolean computation:
 * - Selects backend based on config.backend
 * - Executes Boolean operation
 * - Returns backend result with output polygons
 * 
 * @param soup Input polygon soup with candidate pairs
 * @param op Boolean operation to perform
 * @param backend_name Backend identifier
 * @return BackendResult containing output polygons
 */
BackendResult stage03_intersect(const PolygonSoup& soup, BooleanOp op, 
                                const std::string& backend_name);

/**
 * @brief Stage 04: Classification
 * 
 * Performs winding number evaluation and label propagation:
 * - Computes winding numbers for output polygons
 * - Propagates labels for consistency
 * - Filters polygons based on Boolean operation
 * 
 * @param soup Input/output polygon soup
 * @param op Boolean operation for filtering
 */
void stage04_classify(PolygonSoup& soup, BooleanOp op);

/**
 * @brief Stage 05: Output Reconstruction
 * 
 * Builds final output mesh with attributes:
 * - Constructs output polygons
 * - Interpolates vertex attributes
 * - Materializes floating-point coordinates
 * 
 * @param soup Input/output polygon soup
 */
void stage05_reconstruct(PolygonSoup& soup);

/**
 * @brief Stage 06: Shatter Post-Processing
 * 
 * Connected component extraction for shatter operations:
 * - Extracts connected components from result
 * - Generates piece metadata
 * - Assigns piece IDs and pivot points
 * 
 * @param soup Input/output polygon soup
 */
void stage06_shatterPost(PolygonSoup& soup);

} // namespace stages

//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

/**
 * @brief Get EMBER pipeline version string
 */
const char* getPipelineVersion();

/**
 * @brief Get default pipeline configuration
 */
PipelineConfig getDefaultConfig();

/**
 * @brief Validate pipeline configuration
 * @return True if configuration is valid
 */
bool validateConfig(const PipelineConfig& config, DiagnosticLog& log);

} // namespace ember
