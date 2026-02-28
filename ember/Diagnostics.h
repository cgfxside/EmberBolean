/**
 * @file Diagnostics.h
 * @brief Diagnostic logging system for EMBER Boolean engine
 */

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <cstdarg>

namespace ember {

//=============================================================================
// DIAGNOSTIC CATEGORIES
//=============================================================================

enum class DiagCategory {
    Degenerate,         // Degenerate triangles
    NonManifold,        // Non-manifold edges or vertices
    BowTie,             // Self-intersecting (bow-tie) polygons
    InconsistentWinding,// Inconsistent face winding
    NotWatertight,      // Mesh has boundary edges
    SelfIntersecting,   // Mesh intersects itself
    QuantizationWarning,// Quantization may lose precision
    Info,               // Informational message
    Warning,            // General warning
    Error               // Fatal error
};

//=============================================================================
// DIAGNOSTIC ENTRY
//=============================================================================

struct DiagEntry {
    DiagCategory category;
    std::string message;
    int count;  // For aggregated messages (e.g., "N degenerate triangles")
    
    DiagEntry(DiagCategory cat, const std::string& msg, int cnt = 0)
        : category(cat), message(msg), count(cnt) {}
};

//=============================================================================
// DIAGNOSTIC LOG
//=============================================================================

class DiagnosticLog {
public:
    std::vector<DiagEntry> entries;
    bool has_fatal = false;
    
    // Add a warning entry
    void warn(DiagCategory cat, const std::string& msg, int count = 0) {
        entries.emplace_back(cat, msg, count);
    }
    
    // Add a fatal error entry
    void fatal(DiagCategory cat, const std::string& msg) {
        entries.emplace_back(cat, "FATAL: " + msg, 0);
        has_fatal = true;
    }
    
    // Add an info entry
    void info(const std::string& msg) {
        entries.emplace_back(DiagCategory::Info, msg, 0);
    }
    
    // Check if any errors occurred
    bool hasErrors() const {
        return has_fatal;
    }
    
    // Print all entries to stderr
    void print() const {
        for (const auto& entry : entries) {
            const char* prefix = "";
            switch (entry.category) {
                case DiagCategory::Info: prefix = "[EMBER INFO] "; break;
                case DiagCategory::Warning: prefix = "[EMBER WARN] "; break;
                case DiagCategory::Error: prefix = "[EMBER ERROR] "; break;
                case DiagCategory::Degenerate: prefix = "[EMBER DEGENERATE] "; break;
                case DiagCategory::NonManifold: prefix = "[EMBER NON-MANIFOLD] "; break;
                case DiagCategory::BowTie: prefix = "[EMBER BOWTIE] "; break;
                case DiagCategory::InconsistentWinding: prefix = "[EMBER WINDING] "; break;
                case DiagCategory::NotWatertight: prefix = "[EMBER WATERTIGHT] "; break;
                case DiagCategory::SelfIntersecting: prefix = "[EMBER SELF-INTERSECT] "; break;
                case DiagCategory::QuantizationWarning: prefix = "[EMBER QUANTIZE] "; break;
            }
            std::cerr << prefix << entry.message;
            if (entry.count > 0) {
                std::cerr << " (count: " << entry.count << ")";
            }
            std::cerr << std::endl;
        }
    }
};

//=============================================================================
// LOGGING MACROS
//=============================================================================

#define EMBER_LOG_INFO(fmt, ...)  do { printf("[EMBER] " fmt "\n", ##__VA_ARGS__); } while(0)
#define EMBER_LOG_WARN(fmt, ...)  do { printf("[EMBER WARN] " fmt "\n", ##__VA_ARGS__); } while(0)
#define EMBER_LOG_ERROR(fmt, ...) do { printf("[EMBER ERROR] " fmt "\n", ##__VA_ARGS__); } while(0)

} // namespace ember
