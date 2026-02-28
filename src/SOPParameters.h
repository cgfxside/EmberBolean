/**
 * @file SOPParameters.h
 * @brief Houdini SOP parameter definitions for EMBER Boolean
 * 
 * This file provides Houdini-specific parameter handling.
 * The core parameter struct is defined in ember/Parameters.h
 * 
 * EMBER - Exact Mesh Boolean Engine for Robust operations
 * Copyright (c) 2024-2025
 */

#pragma once

#include <PRM/PRM_Template.h>
#include <PRM/PRM_Name.h>
#include <PRM/PRM_Default.h>
#include <PRM/PRM_ChoiceList.h>
#include "ember/Parameters.h"

/**
 * @brief Parameter names for EMBER Boolean SOP
 */
namespace SOPParameters {

    // Operation type
    extern PRM_Name operationName;
    extern PRM_Name backendName;
    
    // Noise parameters
    extern PRM_Name noiseAmpName;
    extern PRM_Name noiseFreqName;
    extern PRM_Name noiseOctavesName;
    extern PRM_Name noiseSeedName;
    
    // Grout parameters
    extern PRM_Name groutWidthName;
    extern PRM_Name shellNoiseName;
    
    // Segment parameters
    extern PRM_Name segmentSizeName;
    extern PRM_Name maxSegmentsName;
    
    // Quantization parameters
    extern PRM_Name aspectRatioThresholdName;
    
    // Flags
    extern PRM_Name cullExternalName;
    extern PRM_Name forceExactName;

    // Choice lists
    extern PRM_ChoiceList operationMenu;
    extern PRM_ChoiceList backendMenu;

    // Default values
    extern PRM_Default noiseAmpDefault;
    extern PRM_Default noiseFreqDefault;
    extern PRM_Default noiseOctavesDefault;
    extern PRM_Default noiseSeedDefault;
    extern PRM_Default groutWidthDefault;
    extern PRM_Default shellNoiseDefault;
    extern PRM_Default segmentSizeDefault;
    extern PRM_Default maxSegmentsDefault;
    extern PRM_Default aspectRatioThresholdDefault;

    // Parameter template list
    PRM_Template* getTemplateList();
    
    // Convert Houdini parameters to core ember::Parameters
    // This should be called from the SOP's cook method
    ember::Parameters toCoreParameters(const class OP_Context& context);

} // namespace SOPParameters
