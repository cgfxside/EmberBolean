/**
 * @file SOP_EmberBoolean.cpp
 * @brief Houdini SOP node — EMBER Boolean operations
 *
 * Phase 1.6 — Manifold backend integration
 */

#include "SOP_EmberBoolean.h"
#include "SOPParameters.h"
#include "GU_EmberConverter.h"

#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPoly.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_String.h>

#include "ember/PolygonSoup.h"
#include "ember/MeshImport.h"
#include "ember/Diagnostics.h"
#include "backend/IBooleanBackend.h"
#include "backend/BackendFactory.h"

#include <chrono>
#include <string>

// ═══════════════════════════════════════════════════════════════════════════════
// Backend name table — index must match SOPParameters::backendMenu
// ═══════════════════════════════════════════════════════════════════════════════

static const char* kBackendNames[] = {
    "cherchi",   // 0
    "mcut",      // 1
    "manifold",  // 2
    "kigumi"     // 3
};
static constexpr int kNumBackends = 4;

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

// NOTE: ~SOP_EmberBoolean() is defined as = default in the header.

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

    int backend_idx = static_cast<int>(evalInt("backend", 0, t));
    config.backend = (backend_idx >= 0 && backend_idx < kNumBackends)
        ? kBackendNames[backend_idx]
        : "cherchi";

    config.noise_amplitude        = static_cast<float>(evalFloat("noise_amp", 0, t));
    config.noise_frequency        = static_cast<float>(evalFloat("noise_freq", 0, t));
    config.noise_octaves          = static_cast<int>(evalInt("noise_octaves", 0, t));
    config.noise_seed             = static_cast<uint32_t>(evalInt("noise_seed", 0, t));
    config.grout_width            = static_cast<float>(evalFloat("grout_width", 0, t));
    config.shell_noise            = static_cast<float>(evalFloat("shell_noise", 0, t));
    config.segment_size           = static_cast<float>(evalFloat("segment_size", 0, t));
    config.max_segments           = static_cast<int>(evalInt("max_segments", 0, t));
    config.aspect_ratio_threshold = evalFloat("aspect_ratio_threshold", 0, t);
    config.cull_external          = evalInt("cull_external", 0, t) != 0;
    config.force_exact            = evalInt("force_exact", 0, t) != 0;
}

ember::Parameters SOP_EmberBoolean::getCoreParameters(OP_Context& context) {
    ember::PipelineConfig config;
    evalParameters(config, context);
    ember::Parameters params;
    // Minimal conversion — extend as needed
    return params;
}

// ═══════════════════════════════════════════════════════════════════════════════
// COOK
// ═══════════════════════════════════════════════════════════════════════════════

OP_ERROR SOP_EmberBoolean::cookMySop(OP_Context& ctx)
{
    // ── RAII input lock ──────────────────────────────────────────────────────
    OP_AutoLockInputs lock(this);
    if (lock.lock(ctx) >= UT_ERROR_ABORT)
        return error();

    gdp->clearAndDestroy();

    // ── Validate inputs ──────────────────────────────────────────────────────
    const GU_Detail* gdpA = inputGeo(0, ctx);
    const GU_Detail* gdpB = inputGeo(1, ctx);

    if (!gdpA || !gdpB) {
        addError(SOP_MESSAGE, "Both polygon inputs required.");
        return error();
    }

    // ── Read parameters ──────────────────────────────────────────────────────
    ember::PipelineConfig config;
    evalParameters(config, ctx);

    std::string backendName = config.backend;
    bool is_shatter = (config.operation == ember::BooleanOp::Shatter);

    try {
        // ── Convert inputs to PolygonSoup ────────────────────────────────────
        // Uses existing GU_EmberConverter API: convertGUDetailToPolygonSoup(gdp, mesh_id)
        ember::PolygonSoup soupA = ember::convertGUDetailToPolygonSoup(gdpA, 0);
        ember::PolygonSoup soupB = ember::convertGUDetailToPolygonSoup(gdpB, 1);

        // Capture triangle counts BEFORE moving
        const uint32_t nA = static_cast<uint32_t>(soupA.triangles.size());
        const uint32_t nB = static_cast<uint32_t>(soupB.triangles.size());

        // Merge into one soup
        ember::PolygonSoup soup = std::move(soupA);
        soup.triangles.insert(soup.triangles.end(),
                              soupB.triangles.begin(),
                              soupB.triangles.end());

        // CRITICAL: backends use mesh_tri_count to split A vs B triangles.
        // Without this, Cherchi stub returns "need two non-empty meshes".
        soup.mesh_tri_count = { nA, nB };
        soup.mesh_count = 2;

        if (soup.triangles.empty()) {
            addWarning(SOP_MESSAGE, "Empty input geometry after triangulation.");
            return error();
        }

        // ── Quantize ─────────────────────────────────────────────────────────
        // Free functions from ember/MeshImport.h
        ember::QuantizationContext qctx =
            ember::computeQuantization(soup, config.aspect_ratio_threshold);
        ember::quantizeAll(soup, qctx);

        // ── Create backend ───────────────────────────────────────────────────
        auto backend = ember::createBackend(backendName);

        // ── Execute Boolean operation ────────────────────────────────────────
        auto t0 = std::chrono::high_resolution_clock::now();
        ember::BackendResult result = backend->execute(soup, config.operation);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // ── Handle unavailable backend — fall back to Cherchi ────────────────
        if (!result.success &&
            result.error_message.find("not available") != std::string::npos)
        {
            addWarning(SOP_MESSAGE,
                (std::string("Backend '") + backendName
                 + "' not available, falling back to Cherchi.").c_str());

            backend = ember::createBackend("cherchi");
            result = backend->execute(soup, config.operation);
            backendName = "cherchi";
        }

        if (!result.success) {
            addError(SOP_MESSAGE, result.error_message.c_str());
            return error();
        }

        // ── Convert output to GU_Detail ──────────────────────────────────────
        //
        // Manifold backend → output in result.output_soup.triangles[]
        //                     with float positions, piece_id in mesh_id
        // Other backends   → output in result.output_soup.output_polygons[]
        //                     with indexed int_vertices
        //
        bool use_manifold_path =
            (backendName == "manifold") ||
            (!result.output_soup.triangles.empty() &&
             result.output_soup.output_polygons.empty());

        if (use_manifold_path) {
            // ── Manifold path: triangles[] has float positions ────────────────
            // Uses H21 API: GEO_PrimPoly::build() + gdp->setVertexPoint()
            GA_RWHandleV3 P_out(gdp->getP());

            for (const auto& tri : result.output_soup.triangles) {
                GA_Offset ptoffs[3];
                for (int v = 0; v < 3; ++v) {
                    ptoffs[v] = gdp->appendPoint();
                    P_out.set(ptoffs[v], UT_Vector3(
                        tri.v[v][0], tri.v[v][1], tri.v[v][2]));
                }

                // H21: GEO_PrimPoly::build(gdp, npts, closed, do_append_points=0)
                GEO_PrimPoly* ppoly = GEO_PrimPoly::build(gdp, 3, GU_POLY_CLOSED, 0);
                if (ppoly) {
                    for (int v = 0; v < 3; ++v) {
                        gdp->setVertexPoint(ppoly->getVertexOffset(v), ptoffs[v]);
                    }
                }
            }

            // For Shatter, add piece + name attributes from mesh_id
            if (is_shatter && !result.output_soup.triangles.empty()) {
                GA_RWHandleI piece_h(
                    gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "piece", 1));
                GA_RWHandleS name_h(
                    gdp->addStringTuple(GA_ATTRIB_PRIMITIVE, "name", 1));

                GA_Offset prim_off;
                int idx = 0;
                GA_FOR_ALL_PRIMOFF(gdp, prim_off) {
                    if (idx < static_cast<int>(result.output_soup.triangles.size())) {
                        int pid = result.output_soup.triangles[idx].mesh_id;
                        if (piece_h.isValid())
                            piece_h.set(prim_off, pid);
                        if (name_h.isValid()) {
                            char buf[32];
                            snprintf(buf, sizeof(buf), "piece%d", pid);
                            name_h.set(prim_off, buf);
                        }
                    }
                    ++idx;
                }
            }
        } else {
            // ── Standard path: use existing converter ────────────────────────
            ember::EmberAttributeMapping mappingA = ember::buildAttributeMapping(gdpA);
            ember::convertPolygonSoupToGUDetail(result.output_soup, gdp, mappingA);
        }

        // ── Performance info ─────────────────────────────────────────────────
        {
            char info[256];
            snprintf(info, sizeof(info),
                "EMBER: %s | %s | %u tris in -> %u tris out | %.1f ms",
                backendName.c_str(),
                is_shatter ? "Shatter" : "Boolean",
                static_cast<unsigned>(soup.triangles.size()),
                result.num_triangles,
                ms);
            addMessage(SOP_MESSAGE, info);
        }

    } catch (const std::exception& e) {
        UT_String msg;
        msg.sprintf("EMBER internal error: %s", e.what());
        addError(SOP_MESSAGE, msg.c_str());
    } catch (...) {
        addError(SOP_MESSAGE, "EMBER: unknown internal error during cook");
    }

    gdp->getP()->bumpDataId();
    gdp->bumpAllDataIds();

    return error();
}

// ═══════════════════════════════════════════════════════════════════════════════
// OUTPUT / ATTRIBUTE STUBS (declared in header)
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

void SOP_EmberBoolean::captureInputAttributes(
    const GU_Detail* /*gdp*/,
    UT_Array<GA_ROHandleV3>& /*posHandles*/,
    UT_Array<GA_ROHandleV3>& /*vecHandles*/,
    UT_Array<GA_ROHandleF>& /*floatHandles*/)
{
    // Stub — attribute capture for interpolation (future)
}
