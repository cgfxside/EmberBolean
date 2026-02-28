/**
 * @file SOP_EmberBoolean.h
 * @brief Houdini SOP node for EMBER Boolean operations
 * 
 * EMBER - Exact Mesh Boolean Engine for Robust operations
 * Copyright (c) 2024-2025
 */

#pragma once

#include <SOP/SOP_Node.h>
#include <GA/GA_Handle.h>
#include <UT/UT_Array.h>

// Parameter definitions
#include "SOPParameters.h"

namespace ember {
    struct PipelineConfig;
    struct PipelineResult;
    struct PolygonSoup;
}

/**
 * @brief EMBER Boolean SOP Node
 * 
 * Performs exact Boolean operations on polygon meshes using the EMBER engine.
 * Supports Union, Intersection, Difference, and Shatter operations with
 * multiple backend options.
 */
class SOP_EmberBoolean : public SOP_Node {
public:
    // Standard Houdini SOP interface
    static OP_Node* myConstructor(OP_Network* net, const char* name, OP_Operator* op);
    static PRM_Template* buildTemplates();
    
    SOP_EmberBoolean(OP_Network* net, const char* name, OP_Operator* op);
    virtual ~SOP_EmberBoolean() = default;
    
protected:
    // Main cook method
    virtual OP_ERROR cookMySop(OP_Context& context) override;
    
    // Input changed notification
    virtual int isRefInput(OP_InputIdx i) const override { return i > 0; }
    
private:
    // Attribute state for each input
    struct AttributeState {
        GA_ROHandleV3 P_handle;
        GA_ROHandleV3 N_handle;
        GA_ROHandleV3 uv_handle;
        GA_ROHandleV3 Cd_handle;
        GA_ROHandleV3 v_handle;
        UT_Array<GA_ROHandleV3> customPositionHandles;
        UT_Array<GA_ROHandleV3> customVectorHandles;
        UT_Array<GA_ROHandleF> customFloatHandles;
        UT_Array<GA_Attribute*> customAttributes;
        bool uv_is_vertex = false;
    };
    
    AttributeState m_attrStateA;
    AttributeState m_attrStateB;
    
    // Parameter evaluation
    void evalParameters(ember::PipelineConfig& config, OP_Context& context);
    
    // Input capture
    void captureInputAttributes(const GU_Detail* gdp,
                                 UT_Array<GA_ROHandleV3>& posHandles,
                                 UT_Array<GA_ROHandleV3>& vecHandles,
                                 UT_Array<GA_ROHandleF>& floatHandles);
    
    // Output generation
    void buildOutputGeometry(const ember::PipelineResult& result,
                              OP_Context& context);
    
    // Attribute propagation
    void propagateAttributes(GU_Detail* output_gdp,
                              const ember::PolygonSoup& soup);
    
    // Parameter conversion (for ember core library)
    ember::Parameters getCoreParameters(OP_Context& context);
    
    // Parameter template list (defined in SOPParameters.cpp)
    static PRM_Template myTemplateList[];
};
