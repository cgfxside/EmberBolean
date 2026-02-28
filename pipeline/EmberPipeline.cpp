/**
 * @file EmberPipeline.cpp
 * @brief Main pipeline orchestrator implementation for EMBER Boolean engine
 * 
 * This file implements the complete EMBER pipeline with automatic dispatch
 * between Tier 1 (plane fast-path) and Tier 2 (exact mesh-mesh) processing.
 * 
 * Pipeline Flow:
 *   1. Stage 01: Validate input, detect cutter types
 *   2. Decision: All planar? → Stage 02b (fast path)
 *              : Mixed/General? → Stage 02 + 03 (exact path)
 *   3. Stage 04-06: Classification, reconstruction, post-processing
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#include "EmberPipeline.h"
#include "Stage01_Diagnostic.h"
#include "Stage02_EmbreeBVH.h"
#include "Stage02b_PlaneFastPath.h"
#include "Stage03_Intersect.h"
#include "Stage04_Classify.h"
#include "Stage05_Reconstruct.h"
#include "Stage06_ShatterPost.h"

#include "backend/IBooleanBackend.h"
#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"

#include <chrono>
#include <memory>
#include <string>

// Cutter dispatch from upload directory
#include "../upload/CutterDispatch.h"

namespace ember {

//=============================================================================
// VERSION AND DEFAULTS
//=============================================================================

static constexpr const char* PIPELINE_VERSION = "EMBER 1.0.0";

const char* getPipelineVersion() {
    return PIPELINE_VERSION;
}

PipelineConfig getDefaultConfig() {
    PipelineConfig config;
    config.operation = BooleanOp::Union;
    config.backend = "auto";
    config.noise_amplitude = 0.0f;
    config.noise_frequency = 1.0f;
    config.noise_octaves = 3;
    config.noise_seed = 0;
    config.grout_width = 0.0f;
    config.shell_noise = 0.0f;
    config.segment_size = 0.1f;
    config.max_segments = 64;
    config.cull_external = true;
    config.force_exact = false;
    return config;
}

bool validateConfig(const PipelineConfig& config, DiagnosticLog& log) {
    bool valid = true;
    
    // Validate noise parameters
    if (config.noise_amplitude < 0.0f) {
        log.error(DiagCategory::Error, "Noise amplitude cannot be negative");
        valid = false;
    }
    
    if (config.noise_frequency <= 0.0f) {
        log.error(DiagCategory::Error, "Noise frequency must be positive");
        valid = false;
    }
    
    if (config.noise_octaves < 1 || config.noise_octaves > 10) {
        log.warn(DiagCategory::Warning, "Noise octaves should be between 1 and 10");
    }
    
    // Validate grout parameters
    if (config.grout_width < 0.0f) {
        log.error(DiagCategory::Error, "Grout width cannot be negative");
        valid = false;
    }
    
    // Validate segment parameters
    if (config.segment_size <= 0.0f) {
        log.error(DiagCategory::Error, "Segment size must be positive");
        valid = false;
    }
    
    if (config.max_segments < 1) {
        log.error(DiagCategory::Error, "Max segments must be at least 1");
        valid = false;
    }
    
    // Validate backend name
    const std::string& backend = config.backend;
    if (backend != "auto" && backend != "mcut" && backend != "manifold" &&
        backend != "kigumi" && backend != "cherchi") {
        log.error(DiagCategory::Error, "Unknown backend: " + backend);
        valid = false;
    }
    
    return valid;
}

//=============================================================================
// MAIN PIPELINE ENTRY POINT
//=============================================================================

PipelineResult runPipeline(const PolygonSoup& input, const PipelineConfig& config) {
    PipelineResult result;
    auto pipeline_start = std::chrono::steady_clock::now();
    
    EMBER_LOG_INFO("Starting EMBER pipeline v%s", getPipelineVersion());
    
    // Validate configuration
    if (!validateConfig(config, result.diagnostics)) {
        result.success = false;
        result.error_message = "Invalid pipeline configuration";
        return result;
    }
    
    // Create mutable copy of input for processing
    PolygonSoup soup = input;
    
    //=========================================================================
    // STAGE 01: Diagnostic Gatekeeper
    //=========================================================================
    EMBER_LOG_INFO("Stage 01: Diagnostic Gatekeeper");
    auto stage01_start = std::chrono::steady_clock::now();
    
    DispatchResult dispatch;
    DiagnosticLog stage01_log = stages::stage01_diagnostic(soup, config, dispatch);
    result.diagnostics.merge(stage01_log);
    
    auto stage01_end = std::chrono::steady_clock::now();
    result.stage01_time_ms = std::chrono::duration<double, std::milli>(
        stage01_end - stage01_start).count();
    
    // Check for fatal errors in stage 01
    if (result.diagnostics.has_fatal) {
        result.success = false;
        result.error_message = "Fatal error in diagnostic stage";
        return result;
    }
    
    EMBER_LOG_INFO("Stage 01 complete: %zu triangles, %zu Tier 1, %zu Tier 2",
                   soup.triangles.size(), 
                   static_cast<size_t>(dispatch.tier1_count),
                   static_cast<size_t>(dispatch.tier2_count));
    
    //=========================================================================
    // PIPELINE BRANCH: Fast Path vs Exact Path
    //=========================================================================
    
    if (dispatch.isFastPathOnly() && !config.force_exact) {
        //=====================================================================
        // TIER 1 FAST PATH: Plane-based cutting
        //=====================================================================
        EMBER_LOG_INFO("Taking Tier 1 fast path (planar cutters)");
        
        auto stage02b_start = std::chrono::steady_clock::now();
        stages::stage02b_planeFastPath(soup, config, dispatch);
        auto stage02b_end = std::chrono::steady_clock::now();
        result.stage02_time_ms = std::chrono::duration<double, std::milli>(
            stage02b_end - stage02b_start).count();
        
        // For fast path, skip Stage 03 (no backend needed)
        result.stage03_time_ms = 0.0;
        
    } else {
        //=====================================================================
        // TIER 2 EXACT PATH: Mesh-mesh Boolean
        //=====================================================================
        EMBER_LOG_INFO("Taking Tier 2 exact path (mesh-mesh Booleans)");
        
        //=====================================================================
        // STAGE 02: Embree Broad-Phase
        //=====================================================================
        EMBER_LOG_INFO("Stage 02: Embree Broad-Phase");
        auto stage02_start = std::chrono::steady_clock::now();
        
        stages::stage02_broadPhase(soup);
        
        auto stage02_end = std::chrono::steady_clock::now();
        result.stage02_time_ms = std::chrono::duration<double, std::milli>(
            stage02_end - stage02_start).count();
        
        EMBER_LOG_INFO("Stage 02 complete: %zu candidate pairs",
                       soup.candidate_pairs.size());
        
        // Check for excessive candidate pairs
        const size_t max_candidates = 100000000;  // 100M pairs
        if (soup.candidate_pairs.size() > max_candidates) {
            result.diagnostics.warn(DiagCategory::ExcessiveCandidates,
                "Excessive candidate pairs: " + 
                std::to_string(soup.candidate_pairs.size()) +
                ", may cause memory issues");
        }
        
        //=====================================================================
        // STAGE 03: Intersection Resolution
        //=====================================================================
        EMBER_LOG_INFO("Stage 03: Intersection Resolution (backend: %s)", 
                       config.backend.c_str());
        auto stage03_start = std::chrono::steady_clock::now();
        
        BackendResult backend_result = stages::stage03_intersect(
            soup, config.operation, config.backend);
        
        auto stage03_end = std::chrono::steady_clock::now();
        result.stage03_time_ms = std::chrono::duration<double, std::milli>(
            stage03_end - stage03_start).count();
        
        // Check backend result
        if (!backend_result.success) {
            result.success = false;
            result.error_message = backend_result.error_message;
            result.diagnostics.error(DiagCategory::Error, 
                "Backend failed: " + backend_result.error_message);
            return result;
        }
        
        // Transfer backend output to soup
        soup.output_polygons = std::move(backend_result.polygons);
        
        EMBER_LOG_INFO("Stage 03 complete: %zu output polygons",
                       soup.output_polygons.size());
    }
    
    //=========================================================================
    // STAGE 04: Classification
    //=========================================================================
    EMBER_LOG_INFO("Stage 04: Classification");
    auto stage04_start = std::chrono::steady_clock::now();
    
    stages::stage04_classify(soup, config.operation);
    
    auto stage04_end = std::chrono::steady_clock::now();
    result.stage04_time_ms = std::chrono::duration<double, std::milli>(
        stage04_end - stage04_start).count();
    
    EMBER_LOG_INFO("Stage 04 complete");
    
    //=========================================================================
    // STAGE 05: Output Reconstruction
    //=========================================================================
    EMBER_LOG_INFO("Stage 05: Output Reconstruction");
    auto stage05_start = std::chrono::steady_clock::now();
    
    stages::stage05_reconstruct(soup);
    
    auto stage05_end = std::chrono::steady_clock::now();
    result.stage05_time_ms = std::chrono::duration<double, std::milli>(
        stage05_end - stage05_start).count();
    
    EMBER_LOG_INFO("Stage 05 complete: %zu output polygons",
                   soup.output_polygons.size());
    
    //=========================================================================
    // STAGE 06: Shatter Post-Processing (if needed)
    //=========================================================================
    if (config.isShatter()) {
        EMBER_LOG_INFO("Stage 06: Shatter Post-Processing");
        auto stage06_start = std::chrono::steady_clock::now();
        
        stages::stage06_shatterPost(soup);
        
        auto stage06_end = std::chrono::steady_clock::now();
        result.stage06_time_ms = std::chrono::duration<double, std::milli>(
            stage06_end - stage06_start).count();
        
        EMBER_LOG_INFO("Stage 06 complete: %d pieces", soup.piece_count);
    }
    
    //=========================================================================
    // FINALIZE RESULT
    //=========================================================================
    auto pipeline_end = std::chrono::steady_clock::now();
    result.execution_time_ms = std::chrono::duration<double, std::milli>(
        pipeline_end - pipeline_start).count();
    
    // Transfer output
    result.output = std::move(soup);
    result.success = true;
    
    EMBER_LOG_INFO("Pipeline complete: %.2f ms total", result.execution_time_ms);
    
    return result;
}

} // namespace ember
