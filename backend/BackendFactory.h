// ═══════════════════════════════════════════════════════════════════════════════
// EMBER BackendFactory.h — Creates IBooleanBackend instances by name
// ═══════════════════════════════════════════════════════════════════════════════
//
// Supported backends: "cherchi", "manifold", "mcut", "kigumi"
// Unknown/unavailable backends return UnavailableBackend stub.

#pragma once

#include "IBooleanBackend.h"
#include <memory>
#include <string>

namespace ember {

/// Create a backend by name.
/// Returns UnavailableBackend if the requested backend isn't compiled in.
/// NOTE: This is also declared in IBooleanBackend.h for backward compat.
///       The definition lives in BackendFactory.cpp.
std::unique_ptr<IBooleanBackend> createBackend(const std::string& name);

/// Returns true if the error was produced by UnavailableBackend.
bool isUnavailableBackendError(const std::string& msg);

} // namespace ember
