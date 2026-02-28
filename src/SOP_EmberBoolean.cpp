/**
 * @file SOP_EmberBoolean.cpp
 * @brief Implementation of EMBER Boolean SOP node
 *
 * FIXES IN THIS VERSION:
 *
 *   BUG 1 (CRASH — FIXED PREV): evalString() on PRM_INT_E "backend" param.
 *     Use evalInt() + static name table instead.
 *
 *   BUG 2 (CRASH — FIXED PREV): ODR violation, CherchiBackend compiled twice.
 *     CMakeLists.txt must only compile CherchiBackend_stub.cpp.
 *
 *   BUG 3 (CRASH — FIXED PREV): Uncaught C++ exceptions from backend.
 *     try/catch wraps the entire cook body.
 *
 *   BUG 4 (CRASH RISK — FIXED PREV): Manual lock/unlock without RAII.
 *     OP_AutoLockInputs used instead.
 *
 *   BUG 5 (NO OUTPUT): quantizeAll() was never called.
 *     tri.iv[][] stayed all-zeros for every triangle. Any backend that
 *     uses quantized integer coordinates (intersection detection, plane
 *     computation, MCUT vertex dedup) received garbage. Fixed: call
 *     computeQuantization() + quantizeAll() on the combined soup before
 *     passing it to the backend.
 *
 *   BUG 6 (MISLEADING ERROR): Unavailable backend (kigumi, mcut, manifold)
 *     returned success=false with a terse error string, making it look like
 *     a cook failure. Fixed: detect UnavailableBackend result, emit a
 *     Houdini WARNING (not error), fall back to Cherchi so the cook still
 *     produces passthrough geometry.
 */

#include "GU_EmberConverter.h"
#include "SOP_EmberBoolean.h"

#include <GU/GU_Detail.h>
#include <GU/GU_DetailHandle.h>
#include <GA/GA_Primitive.h>
#include <GA/GA_Types.h>
#include <GA/GA_Attribute.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <PRM/PRM_Name.h>
#include <UT/UT_String.h>
#include <UT/UT_WorkBuffer.h>
#include <UT/UT_ErrorManager.h>

#include "ember/PolygonSoup.h"
#include "ember/MeshImport.h"
#include "ember/Diagnostics.h"
#include "backend/IBooleanBackend.h"

// ── Backend name table (index must match backendNames[] in SOPParameters.cpp) ─
static const char* kBackendNames[] = {
    "cherchi",   // 0
    "mcut",      // 1
    "manifold",  // 2
    "kigumi"     // 3
};
static constexpr int kNumBackends = 4;

// Prefix used by UnavailableBackend::error_message so we can detect it.
static constexpr const char* kUnavailablePrefix = "Backend '";

namespace ember {
    struct PipelineConfig {
        BooleanOp   operation              = BooleanOp::Union;
        std::string backend                = "cherchi";
        float       noise_amplitude        = 0.0f;
        float       noise_frequency        = 1.0f;
        int         noise_octaves          = 3;
        uint32_t    noise_seed             = 0;
        float       grout_width            = 0.0f;
        float       shell_noise            = 0.0f;
        float       segment_size           = 0.1f;
        int         max_segments           = 64;
        double      aspect_ratio_threshold = 1000.0;
        bool        cull_external          = true;
        bool        force_exact            = false;

        bool hasNoise() const { return noise_amplitude > 0.0f; }
    };

    struct PipelineResult {
        PolygonSoup   output;
        bool          success = false;
        std::string   error_message;
        DiagnosticLog diagnostics;
        double        execution_time_ms = 0.0;
    };
}

// ═══════════════════════════════════════════════════════════════════════════════
// PARAMETER TEMPLATES
// ═══════════════════════════════════════════════════════════════════════════════

static PRM_Default operationDefault(0);
static PRM_Default backendDefault(0);
static PRM_Default cullExternalDefault(1);
static PRM_Default forceExactDefault(0);

PRM_Template SOP_EmberBoolean::myTemplateList[] = {
    PRM_Template(PRM_INT_E, 1, &SOPParameters::operationName, &operationDefault, &SOPParameters::operationMenu),
    PRM_Template(PRM_INT_E, 1, &SOPParameters::backendName,   &backendDefault,   &SOPParameters::backendMenu),
    PRM_Template(PRM_FLT_E, 1, &SOPParameters::noiseAmpName,             &SOPParameters::noiseAmpDefault),
    PRM_Template(PRM_FLT_E, 1, &SOPParameters::noiseFreqName,            &SOPParameters::noiseFreqDefault),
    PRM_Template(PRM_INT_E, 1, &SOPParameters::noiseOctavesName,         &SOPParameters::noiseOctavesDefault),
    PRM_Template(PRM_INT_E, 1, &SOPParameters::noiseSeedName,            &SOPParameters::noiseSeedDefault),
    PRM_Template(PRM_FLT_E, 1, &SOPParameters::groutWidthName,           &SOPParameters::groutWidthDefault),
    PRM_Template(PRM_FLT_E, 1, &SOPParameters::shellNoiseName,           &SOPParameters::shellNoiseDefault),
    PRM_Template(PRM_FLT_E, 1, &SOPParameters::segmentSizeName,          &SOPParameters::segmentSizeDefault),
    PRM_Template(PRM_INT_E, 1, &SOPParameters::maxSegmentsName,          &SOPParameters::maxSegmentsDefault),
    PRM_Template(PRM_FLT_E, 1, &SOPParameters::aspectRatioThresholdName, &SOPParameters::aspectRatioThresholdDefault),
    PRM_Template(PRM_TOGGLE_E, 1, &SOPParameters::cullExternalName, &cullExternalDefault),
    PRM_Template(PRM_TOGGLE_E, 1, &SOPParameters::forceExactName,   &forceExactDefault),
    PRM_Template()
};

// ═══════════════════════════════════════════════════════════════════════════════
// CONSTRUCTOR / FACTORY
// ═══════════════════════════════════════════════════════════════════════════════

OP_Node* SOP_EmberBoolean::myConstructor(OP_Network* net, const char* name, OP_Operator* op) {
    return new SOP_EmberBoolean(net, name, op);
}

PRM_Template* SOP_EmberBoolean::buildTemplates() {
    return myTemplateList;
}

SOP_EmberBoolean::SOP_EmberBoolean(OP_Network* net, const char* name, OP_Operator* op)
    : SOP_Node(net, name, op)
{}

// ═══════════════════════════════════════════════════════════════════════════════
// PARAMETER EVALUATION
// ═══════════════════════════════════════════════════════════════════════════════

void SOP_EmberBoolean::evalParameters(ember::PipelineConfig& config, OP_Context& context) {
    fpreal t = context.getTime();

    switch (evalInt("operation", 0, t)) {
        case 0: config.operation = ember::BooleanOp::Union;        break;
        case 1: config.operation = ember::BooleanOp::Intersection; break;
        case 2: config.operation = ember::BooleanOp::DiffAB;       break;
        case 3: config.operation = ember::BooleanOp::DiffBA;       break;
        case 4: config.operation = ember::BooleanOp::Shatter;      break;
        case 5: config.operation = ember::BooleanOp::Seam;         break;
        default: config.operation = ember::BooleanOp::Union;       break;
    }

    // BUG 1 FIX: evalInt, not evalString, for PRM_INT_E parameters.
    int backend_idx = static_cast<int>(evalInt("backend", 0, t));
    config.backend = (backend_idx >= 0 && backend_idx < kNumBackends)
        ? kBackendNames[backend_idx]
        : "cherchi";

    config.noise_amplitude = static_cast<float>(evalFloat("noise_amp",     0, t));
    config.noise_frequency = static_cast<float>(evalFloat("noise_freq",    0, t));
    config.noise_octaves   = static_cast<int>  (evalInt  ("noise_octaves", 0, t));
    config.noise_seed      = static_cast<uint32_t>(evalInt("noise_seed",   0, t));
    config.grout_width     = static_cast<float>(evalFloat("grout_width",   0, t));
    config.shell_noise     = static_cast<float>(evalFloat("shell_noise",   0, t));
    config.segment_size    = static_cast<float>(evalFloat("segment_size",  0, t));
    config.max_segments    = static_cast<int>  (evalInt  ("max_segments",  0, t));
    config.aspect_ratio_threshold = evalFloat("aspect_ratio_threshold", 0, t);
    config.cull_external   = evalInt("cull_external", 0, t) != 0;
    config.force_exact     = evalInt("force_exact",   0, t) != 0;
}

// ═══════════════════════════════════════════════════════════════════════════════
// PARAMETER CONVERSION (for ember core library)
// ═══════════════════════════════════════════════════════════════════════════════

ember::Parameters SOP_EmberBoolean::getCoreParameters(OP_Context& context) {
    ember::Parameters params;
    fpreal t = context.getTime();

    switch (evalInt("operation", 0, t)) {
        case 0: params.operation = ember::Parameters::Operation::Union;         break;
        case 1: params.operation = ember::Parameters::Operation::Intersection;  break;
        case 2: params.operation = ember::Parameters::Operation::DifferenceAB;  break;
        case 3: params.operation = ember::Parameters::Operation::DifferenceBA;  break;
        case 4: params.operation = ember::Parameters::Operation::Shatter;       break;
        case 5: params.operation = ember::Parameters::Operation::Seam;          break;
        default: params.operation = ember::Parameters::Operation::Union;        break;
    }
    switch (evalInt("backend", 0, t)) {
        case 0: params.backend = ember::Parameters::Backend::Cherchi;  break;
        case 1: params.backend = ember::Parameters::Backend::MCUT;     break;
        case 2: params.backend = ember::Parameters::Backend::Manifold; break;
        case 3: params.backend = ember::Parameters::Backend::Kigumi;   break;
        default: params.backend = ember::Parameters::Backend::Cherchi; break;
    }

    params.noise_amplitude        = static_cast<float>(evalFloat("noise_amp",              0, t));
    params.noise_frequency        = static_cast<float>(evalFloat("noise_freq",             0, t));
    params.noise_octaves          = static_cast<int>  (evalInt  ("noise_octaves",          0, t));
    params.noise_seed             = static_cast<int>  (evalInt  ("noise_seed",             0, t));
    params.grout_width            = static_cast<float>(evalFloat("grout_width",            0, t));
    params.shell_noise            = static_cast<float>(evalFloat("shell_noise",            0, t));
    params.segment_size           = static_cast<float>(evalFloat("segment_size",           0, t));
    params.max_segments           = static_cast<int>  (evalInt  ("max_segments",           0, t));
    params.aspect_ratio_threshold = static_cast<float>(evalFloat("aspect_ratio_threshold", 0, t));
    params.cull_external          = evalInt("cull_external", 0, t) != 0;
    params.force_exact            = evalInt("force_exact",   0, t) != 0;

    return params;
}

// ═══════════════════════════════════════════════════════════════════════════════
// ATTRIBUTE CAPTURE
// ═══════════════════════════════════════════════════════════════════════════════

void SOP_EmberBoolean::captureInputAttributes(const GU_Detail* gdp,
                                               UT_Array<GA_ROHandleV3>& posHandles,
                                               UT_Array<GA_ROHandleV3>& vecHandles,
                                               UT_Array<GA_ROHandleF>&  floatHandles) {
    const GA_AttributeDict& attrDict = gdp->getAttributeDict(GA_ATTRIB_POINT);
    for (GA_AttributeDict::iterator it = attrDict.begin(); !it.atEnd(); ++it) {
        GA_Attribute* attr = it.attrib();
        if (!attr) continue;
        const char* name = attr->getName();

        GA_ROHandleV3 hvec(attr);
        if (hvec.isValid()) {
            if (strcmp(name, "rest") == 0 || strcmp(name, "pref") == 0)
                posHandles.append(hvec);
            else if (strcmp(name, "v") == 0 || strcmp(name, "vel") == 0)
                vecHandles.append(hvec);
            continue;
        }
        GA_ROHandleF hfloat(attr);
        if (hfloat.isValid())
            floatHandles.append(hfloat);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// MAIN COOK
// ═══════════════════════════════════════════════════════════════════════════════

OP_ERROR SOP_EmberBoolean::cookMySop(OP_Context& context) {

    // BUG 4 FIX: RAII lock — unlockInputs() always called, even on exception.
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    // BUG 3 FIX: Never let a C++ exception escape into Houdini's call stack.
    try {

        const GU_Detail* gdpA = inputGeo(0, context);
        const GU_Detail* gdpB = inputGeo(1, context);

        if (!gdpA) { addWarning(SOP_MESSAGE, "Input A required"); return error(); }
        if (!gdpB) { addWarning(SOP_MESSAGE, "Input B required"); return error(); }

        ember::PipelineConfig config;
        evalParameters(config, context);

        // ── Stage 1: GU_Detail → PolygonSoup ─────────────────────────────────
        ember::PolygonSoup soupA = ember::convertGUDetailToPolygonSoup(gdpA, 0);
        ember::PolygonSoup soupB = ember::convertGUDetailToPolygonSoup(gdpB, 1);

        if (soupA.triangles.empty()) {
            addWarning(SOP_MESSAGE, "Input A has no triangles after conversion");
            return error();
        }
        if (soupB.triangles.empty()) {
            addWarning(SOP_MESSAGE, "Input B has no triangles after conversion");
            return error();
        }

        // ── Stage 2: Merge soups ──────────────────────────────────────────────
        ember::PolygonSoup combined;
        combined.triangles.reserve(soupA.triangles.size() + soupB.triangles.size());
        combined.triangles.insert(combined.triangles.end(),
            soupA.triangles.begin(), soupA.triangles.end());
        combined.triangles.insert(combined.triangles.end(),
            soupB.triangles.begin(), soupB.triangles.end());
        combined.mesh_tri_count = {
            static_cast<uint32_t>(soupA.triangles.size()),
            static_cast<uint32_t>(soupB.triangles.size())
        };

        // ── Stage 3: Quantize all vertices ───────────────────────────────────
        //
        // BUG 5 FIX: quantizeAll() was never called in the original code.
        //
        // Without this, tri.iv[v][0..2] stays {0,0,0} for every triangle.
        // Consequences:
        //   • MCUT vertex dedup sees every vertex as the same point → empty mesh
        //   • computePlaneFromTriangle() gets zero-length edges → zero normals
        //   • Any intersection predicate operates on a degenerate point cloud
        //   • The Cherchi stub works around this because it uses tri.v (float)
        //     directly, but once real intersection code is added it will break too
        //
        // computeQuantization() analyzes the merged bounding box to choose
        // uniform or per-axis scaling. quantizeAll() fills tri.iv[v][0..2]
        // and stores the scale in combined.quantization_scale[]/inv_scale[].
        {
            ember::QuantizationContext qctx =
                ember::computeQuantization(combined, config.aspect_ratio_threshold);
            ember::quantizeAll(combined, qctx);
        }

        // ── Stage 4: Run backend ──────────────────────────────────────────────
        std::unique_ptr<ember::IBooleanBackend> backend =
            ember::createBackend(config.backend);

        if (!backend) {
            addError(SOP_MESSAGE,
                     ("EMBER: could not create backend: " + config.backend).c_str());
            return error();
        }

        ember::BackendResult result = backend->execute(combined, config.operation);

        // BUG 6 FIX: Graceful handling of unavailable backends.
        //
        // kigumi, mcut, manifold all compile as UnavailableBackend stubs when
        // their libraries are not linked. execute() on an UnavailableBackend
        // immediately returns success=false with a message starting with
        // "Backend '...".
        //
        // Previously this caused cookMySop to call addError() and return early
        // with no geometry — the node turned red and the user had to remember
        // which backends are compiled in.
        //
        // Now we detect the "not available" message prefix, emit a yellow
        // warning instead of a red error, and automatically retry with the
        // Cherchi stub so the node still cooks and produces passthrough geometry.
        if (!result.success) {
            const bool is_unavailable =
                result.error_message.rfind(kUnavailablePrefix, 0) == 0;

            if (is_unavailable) {
                UT_String warn;
                warn.sprintf("EMBER: %s Falling back to Cherchi (passthrough).",
                             result.error_message.c_str());
                addWarning(SOP_MESSAGE, warn.c_str());

                auto fallback = ember::createBackend("cherchi");
                if (fallback)
                    result = fallback->execute(combined, config.operation);
            }

            if (!result.success) {
                addError(SOP_MESSAGE, result.error_message.c_str());
                return error();
            }
        }

        // ── Stage 5: PolygonSoup → GU_Detail ─────────────────────────────────
        // convertPolygonSoupToGUDetail() calls gdp->clearAndDestroy() itself.
        ember::EmberAttributeMapping mappingA = ember::buildAttributeMapping(gdpA);
        ember::convertPolygonSoupToGUDetail(result.output_soup, gdp, mappingA);

    } catch (const std::exception& e) {
        UT_String msg;
        msg.sprintf("EMBER internal error: %s", e.what());
        addError(SOP_MESSAGE, msg.c_str());
    } catch (...) {
        addError(SOP_MESSAGE, "EMBER: unknown internal error during cook");
    }

    return error();
}

// ═══════════════════════════════════════════════════════════════════════════════
// OUTPUT / ATTRIBUTE STUBS
// ═══════════════════════════════════════════════════════════════════════════════

void SOP_EmberBoolean::buildOutputGeometry(const ember::PipelineResult& result,
                                            OP_Context& /*context*/) {
    if (!result.success)
        addError(SOP_MESSAGE, result.error_message.c_str());
}

void SOP_EmberBoolean::propagateAttributes(GU_Detail* /*output_gdp*/,
                                            const ember::PolygonSoup& /*soup*/) {
    // Full attribute propagation lives in GU_EmberConverter.cpp.
}
