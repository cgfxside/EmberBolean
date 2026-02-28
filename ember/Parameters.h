/**
 * @file ember/Parameters.h
 * @brief Core parameter structure for EMBER Boolean operations
 * 
 * This is a Houdini-agnostic parameter struct used by the core library.
 * The SOP layer (src/SOPParameters.h) provides conversion from Houdini parms.
 */

#pragma once

#include <cstdint>
#include <cmath>

namespace ember {

/**
 * @brief Core parameters for EMBER Boolean operations.
 * 
 * This struct is Houdini-agnostic and can be used by any frontend.
 * All values are in normalized units (not Houdini parameter units).
 */
struct Parameters {
    // Operation type
    enum class Operation {
        Union = 0,
        Intersection = 1,
        DifferenceAB = 2,  // A - B
        DifferenceBA = 3,  // B - A
        Shatter = 4,
        Seam = 5
    };
    
    // Backend selection
    enum class Backend {
        Cherchi = 0,   // Exact (default)
        MCUT = 1,
        Manifold = 2,
        Kigumi = 3     // Exact
    };
    
    Operation operation = Operation::Union;
    Backend backend = Backend::Cherchi;
    
    // Noise parameters (for noised cutters)
    float noise_amplitude = 0.0f;
    float noise_frequency = 1.0f;
    int noise_octaves = 3;
    int noise_seed = 0;
    
    // Grout parameters
    float grout_width = 0.0f;
    float shell_noise = 0.0f;
    
    // Segment/tessellation parameters
    float segment_size = 0.1f;
    int max_segments = 64;
    
    // Quantization parameters
    float aspect_ratio_threshold = 1000.0f;
    
    // Flags
    bool cull_external = true;
    bool force_exact = false;
    
    // Convenience constructors
    Parameters() = default;
    
    explicit Parameters(Operation op) : operation(op) {}
    
    // Check if this is a difference operation
    bool isDifference() const {
        return operation == Operation::DifferenceAB || 
               operation == Operation::DifferenceBA;
    }
    
    // Check if we need Winding Number Volume computation
    bool needsWNV() const {
        return operation == Operation::Union ||
               operation == Operation::Intersection ||
               operation == Operation::DifferenceAB ||
               operation == Operation::DifferenceBA;
    }
};

} // namespace ember
