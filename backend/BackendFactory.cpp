// ═══════════════════════════════════════════════════════════════════════════════
// EMBER BackendFactory — Creates IBooleanBackend instances by name
// ═══════════════════════════════════════════════════════════════════════════════
//
// Supported backends:
//   "cherchi"  — CherchiBackend (always available via stub)
//   "manifold" — ManifoldBackend (requires EMBER_HAS_MANIFOLD)
//   "mcut"     — MCUTBackend (requires EMBER_HAS_MCUT)
//   "kigumi"   — KigumiBackend (requires EMBER_HAS_KIGUMI)
//
// Unknown or unavailable backends return UnavailableBackend which
// reports a clear error message and sets success=false.
// ═══════════════════════════════════════════════════════════════════════════════

#include "IBooleanBackend.h"
#include "CherchiBackend.h"

#ifdef EMBER_HAS_MANIFOLD
#include "ManifoldBackend.h"
#endif

#ifdef EMBER_HAS_MCUT
#include "MCUTBackend.h"
#endif

#ifdef EMBER_HAS_KIGUMI
#include "KigumiBackend.h"
#endif

#include <memory>
#include <string>

namespace ember {

// ─────────────────────────────────────────────────────────────────────────────
// UnavailableBackend — graceful fallback when a backend isn't compiled in
// ─────────────────────────────────────────────────────────────────────────────

static const char* kUnavailablePrefix = "[UNAVAILABLE]";

class UnavailableBackend : public IBooleanBackend {
    std::string m_name;
public:
    explicit UnavailableBackend(const std::string& requested)
        : m_name(requested) {}

    std::string name() const override { return m_name; }

    std::string description() const override {
        return m_name + " (not compiled — rebuild with EMBER_HAS_"
             + m_name + ")";
    }

    BackendResult execute(const PolygonSoup& /*soup*/, BooleanOp /*op*/) override {
        BackendResult r;
        r.success       = false;
        r.error_message = std::string(kUnavailablePrefix)
                        + " Backend '" + m_name + "' is not available. "
                        + "Rebuild with -DEMBER_HAS_"
                        + m_name + "=ON.";
        return r;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Factory
// ─────────────────────────────────────────────────────────────────────────────

std::unique_ptr<IBooleanBackend> createBackend(const std::string& name)
{
    if (name == "cherchi")
        return std::make_unique<CherchiBackend>();

#ifdef EMBER_HAS_MANIFOLD
    if (name == "manifold")
        return std::make_unique<ManifoldBackend>();
#endif

#ifdef EMBER_HAS_MCUT
    if (name == "mcut")
        return std::make_unique<MCUTBackend>();
#endif

#ifdef EMBER_HAS_KIGUMI
    if (name == "kigumi")
        return std::make_unique<KigumiBackend>();
#endif

    // Not compiled in or unknown name → graceful stub
    return std::make_unique<UnavailableBackend>(name);
}

/// Detect UnavailableBackend errors so the SOP can fall back
bool isUnavailableBackendError(const std::string& msg)
{
    return msg.rfind(kUnavailablePrefix, 0) == 0;
}

} // namespace ember
