/**
 * @file SOPParameters.cpp
 * @brief Parameter definitions for EMBER Boolean SOP node
 */

#include "SOPParameters.h"

namespace SOPParameters {

    // Operation type
    PRM_Name operationName("operation", "Operation");
    PRM_Name backendName("backend", "Backend");
    
    // Noise parameters
    PRM_Name noiseAmpName("noise_amp", "Noise Amplitude");
    PRM_Name noiseFreqName("noise_freq", "Noise Frequency");
    PRM_Name noiseOctavesName("noise_octaves", "Noise Octaves");
    PRM_Name noiseSeedName("noise_seed", "Noise Seed");
    
    // Grout parameters
    PRM_Name groutWidthName("grout_width", "Grout Width");
    PRM_Name shellNoiseName("shell_noise", "Shell Noise");
    
    // Segment parameters
    PRM_Name segmentSizeName("segment_size", "Segment Size");
    PRM_Name maxSegmentsName("max_segments", "Max Segments");
    
    // Quantization parameters
    PRM_Name aspectRatioThresholdName("aspect_ratio_threshold", "Aspect Ratio Threshold");
    
    // Flags
    PRM_Name cullExternalName("cull_external", "Cull External Faces");
    PRM_Name forceExactName("force_exact", "Force Exact Mode");

    // Choice list names
    static PRM_Name operationNames[] = {
        PRM_Name("union", "Union"),
        PRM_Name("intersection", "Intersection"),
        PRM_Name("diffab", "Difference (A - B)"),
        PRM_Name("diffba", "Difference (B - A)"),
        PRM_Name("shatter", "Shatter"),
        PRM_Name("seam", "Seam"),
        PRM_Name(0)
    };

    static PRM_Name backendNames[] = {
        PRM_Name("cherchi", "Cherchi (Exact)"),
        PRM_Name("mcut", "MCUT"),
        PRM_Name("manifold", "Manifold"),
        PRM_Name("kigumi", "kigumi (Exact)"),
        PRM_Name(0)
    };

    // Choice lists
    PRM_ChoiceList operationMenu(PRM_CHOICELIST_SINGLE, operationNames);
    PRM_ChoiceList backendMenu(PRM_CHOICELIST_SINGLE, backendNames);

    // Default values
    PRM_Default noiseAmpDefault(0.0f);
    PRM_Default noiseFreqDefault(1.0f);
    PRM_Default noiseOctavesDefault(3);
    PRM_Default noiseSeedDefault(0);
    PRM_Default groutWidthDefault(0.0f);
    PRM_Default shellNoiseDefault(0.0f);
    PRM_Default segmentSizeDefault(0.1f);
    PRM_Default maxSegmentsDefault(64);
    PRM_Default aspectRatioThresholdDefault(1000.0f);

    // Default values for integer parameters (0 = first choice in menu)
    PRM_Default operationDefault(0);
    PRM_Default backendDefault(0);
    
    // Default values for toggle parameters
    PRM_Default cullExternalDefault(1);   // Default ON (1)
    PRM_Default forceExactDefault(0);     // Default OFF (0)

    // Parameter template list
    // PRM_Template(type, tupleSize, namePtr, defaultPtr, menuPtr, rangePtr, callback, spareData, exportFlag)
    static PRM_Template templateList[] = {
        PRM_Template(PRM_INT_E, 1, &operationName, &operationDefault, &operationMenu),
        PRM_Template(PRM_INT_E, 1, &backendName, &backendDefault, &backendMenu),
        
        // Noise parameters
        PRM_Template(PRM_FLT_E, 1, &noiseAmpName, &noiseAmpDefault),
        PRM_Template(PRM_FLT_E, 1, &noiseFreqName, &noiseFreqDefault),
        PRM_Template(PRM_INT_E, 1, &noiseOctavesName, &noiseOctavesDefault),
        PRM_Template(PRM_INT_E, 1, &noiseSeedName, &noiseSeedDefault),
        
        // Grout parameters
        PRM_Template(PRM_FLT_E, 1, &groutWidthName, &groutWidthDefault),
        PRM_Template(PRM_FLT_E, 1, &shellNoiseName, &shellNoiseDefault),
        
        // Segment parameters
        PRM_Template(PRM_FLT_E, 1, &segmentSizeName, &segmentSizeDefault),
        PRM_Template(PRM_INT_E, 1, &maxSegmentsName, &maxSegmentsDefault),
        
        // P2 FIX: Configurable aspect ratio threshold
        PRM_Template(PRM_FLT_E, 1, &aspectRatioThresholdName, &aspectRatioThresholdDefault),
        
        // Flags
        PRM_Template(PRM_TOGGLE_E, 1, &cullExternalName, &cullExternalDefault),
        PRM_Template(PRM_TOGGLE_E, 1, &forceExactName, &forceExactDefault),
        
        PRM_Template()
    };

    PRM_Template* getTemplateList() {
        return templateList;
    }
    
    // Note: toCoreParameters implementation is in SOP_EmberBoolean.cpp
    // because it needs access to the SOP node's eval methods

} // namespace SOPParameters
