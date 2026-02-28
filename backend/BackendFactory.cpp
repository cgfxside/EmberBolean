/**
 * @file BackendFactory.cpp
 * @brief Factory implementation for creating Boolean backend instances
 */

#include "IBooleanBackend.h"
#include "CherchiBackend.h"
#include "MCUTBackend.h"
#include "ManifoldBackend.h"
#include "KigumiBackend.h"

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

//=============================================================================
// AVAILABILITY MACROS
//=============================================================================

#define EMBER_HAS_MCUT_BACKEND 0      // Enable when MCUT is linked
#define EMBER_HAS_MANIFOLD_BACKEND 0  // Enable when Manifold is linked
#define EMBER_HAS_KIGUMI_BACKEND 0    // Enable when kigumi is linked
#define EMBER_HAS_CHERCHI_BACKEND 1   // Always available (stub included)

namespace ember {

namespace {

//=============================================================================
// UNAVAILABLE BACKEND STUB
//=============================================================================

class UnavailableBackend : public IBooleanBackend {
    std::string backend_name_;
    std::string dependency_name_;

public:
    UnavailableBackend(const std::string& name, const std::string& dep)
        : backend_name_(name), dependency_name_(dep) {}

    std::string name() const override {
        return backend_name_ + " (unavailable)";
    }

    std::string description() const override {
        return backend_name_ + " backend is not available. "
               "Rebuild with " + dependency_name_ + " support enabled.";
    }

    BackendResult execute(const PolygonSoup& soup, BooleanOp op) override {
        (void)soup;
        (void)op;
        BackendResult result;
        result.success = false;
        result.error_message =
            "Backend '" + backend_name_ + "' is not available. "
            "Please rebuild EMBER with " + dependency_name_ + " support enabled.";
        return result;
    }
};

//=============================================================================
// NAME NORMALIZATION
//=============================================================================

std::string normalizeBackendName(const std::string& name) {
    std::string result;
    result.reserve(name.size());
    for (char c : name) {
        if (!std::isspace(static_cast<unsigned char>(c))) {
            result.push_back(std::tolower(static_cast<unsigned char>(c)));
        }
    }
    return result;
}

//=============================================================================
// PER-BACKEND CREATION HELPERS
//=============================================================================

std::unique_ptr<IBooleanBackend> createMCUTBackendImpl() {
#if EMBER_HAS_MCUT_BACKEND
    return std::make_unique<MCUTBackend>();
#else
    return std::make_unique<UnavailableBackend>("mcut", "MCUT library");
#endif
}

std::unique_ptr<IBooleanBackend> createManifoldBackendImpl() {
#if EMBER_HAS_MANIFOLD_BACKEND
    return std::make_unique<ManifoldBackend>();
#else
    return std::make_unique<UnavailableBackend>("manifold", "Manifold library");
#endif
}

std::unique_ptr<IBooleanBackend> createKigumiBackendImpl() {
#if EMBER_HAS_KIGUMI_BACKEND
    return std::make_unique<KigumiBackend>();
#else
    return std::make_unique<UnavailableBackend>("kigumi", "kigumi library");
#endif
}

std::unique_ptr<IBooleanBackend> createCherchiBackendImpl() {
#if EMBER_HAS_CHERCHI_BACKEND
    return std::make_unique<CherchiBackend>();
#else
    return std::make_unique<UnavailableBackend>("cherchi", "Cherchi library");
#endif
}

std::vector<std::string> getBackendPreferenceOrder() {
    return {"cherchi", "kigumi", "manifold", "mcut"};
}

bool isBackendAvailable(const std::string& n) {
    if (n == "mcut")     return EMBER_HAS_MCUT_BACKEND     != 0;
    if (n == "manifold") return EMBER_HAS_MANIFOLD_BACKEND != 0;
    if (n == "kigumi")   return EMBER_HAS_KIGUMI_BACKEND   != 0;
    if (n == "cherchi")  return EMBER_HAS_CHERCHI_BACKEND  != 0;
    return false;
}

} // anonymous namespace

//=============================================================================
// PUBLIC FACTORY — makeBackend (primary implementation)
//=============================================================================

std::unique_ptr<IBooleanBackend> makeBackend(const std::string& name) {
    const std::string n = normalizeBackendName(name);

    if (n == "mcut")     return createMCUTBackendImpl();
    if (n == "manifold") return createManifoldBackendImpl();
    if (n == "kigumi")   return createKigumiBackendImpl();
    if (n == "cherchi")  return createCherchiBackendImpl();

    if (n == "auto") {
        for (const auto& preferred : getBackendPreferenceOrder()) {
            if (isBackendAvailable(preferred)) {
                return makeBackend(preferred);
            }
        }
        // No backends fully available — fall through to cherchi stub
        return createCherchiBackendImpl();
    }

    throw std::invalid_argument(
        "Unknown backend: '" + name +
        "'. Supported: mcut, manifold, kigumi, cherchi, auto");
}

//=============================================================================
// PUBLIC FACTORY — createBackend (alias, declared in IBooleanBackend.h)
//=============================================================================

std::unique_ptr<IBooleanBackend> createBackend(const std::string& name) {
    return makeBackend(name);
}

//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

std::vector<std::string> getAvailableBackends() {
    std::vector<std::string> available;
    for (const auto& name : getBackendPreferenceOrder()) {
        if (isBackendAvailable(name)) {
            available.push_back(name);
        }
    }
    return available;
}

std::string getDefaultBackendName() {
    for (const auto& name : getBackendPreferenceOrder()) {
        if (isBackendAvailable(name)) {
            return name;
        }
    }
    return "cherchi"; // always has a stub
}

} // namespace ember
