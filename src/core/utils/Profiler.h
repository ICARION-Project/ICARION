// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file Profiler.h
 * @brief Performance profiling system for ICARION
 * 
 * Provides low-overhead timing instrumentation to identify bottlenecks
 * in simulation performance. Particularly useful for diagnosing why
 * OpenMP parallelization may not be providing expected speedups.
 */

#pragma once

#include <string>
#include <chrono>
#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <iostream>
#include <mutex>
#include <limits>

namespace ICARION {
namespace profiling {

/**
 * @brief Timing result for a single profiling section
 */
struct TimingResult {
    std::string name;           ///< Section name
    double duration_ms;         ///< Total accumulated time (ms)
    size_t call_count;          ///< Number of times section was called
    double min_ms;              ///< Minimum duration (ms)
    double max_ms;              ///< Maximum duration (ms)
    double avg_ms;              ///< Average duration per call (ms)
};

/**
 * @brief RAII profiling timer
 * 
 * Automatically records timing when going out of scope.
 * Thread-safe for use in parallel regions.
 */
class ProfileSection {
public:
    explicit ProfileSection(const std::string& name);
    ~ProfileSection();
    
    ProfileSection(const ProfileSection&) = delete;
    ProfileSection& operator=(const ProfileSection&) = delete;
    
private:
    std::string name_;
    std::chrono::steady_clock::time_point start_;
};

/**
 * @brief Global profiler for performance measurement
 * 
 * Singleton pattern for easy access throughout the codebase.
 * Thread-safe for concurrent use in OpenMP regions.
 */
class Profiler {
public:
    static Profiler& getInstance();
    
    void enable(bool enable = true);
    bool isEnabled() const { return enabled_; }
    
    void startSection(const std::string& name);
    void endSection(const std::string& name);
    
    void reset();
    std::vector<TimingResult> getResults() const;
    
    void printSummary(std::ostream& os = std::cout) const;
    void exportJSON(const std::string& filename) const;
    void exportCSV(const std::string& filename) const;
    
private:
    Profiler() = default;
    ~Profiler() = default;
    Profiler(const Profiler&) = delete;
    Profiler& operator=(const Profiler&) = delete;
    
    struct SectionData {
        double total_time_ms = 0.0;
        size_t call_count = 0;
        double min_ms = std::numeric_limits<double>::max();
        double max_ms = 0.0;
        std::chrono::steady_clock::time_point start;
    };
    
    bool enabled_ = false;
    std::map<std::string, SectionData> sections_;
    mutable std::mutex mutex_;  ///< Thread-safety for OpenMP regions
};

/**
 * @brief Convenience macro for profiling a scope
 * 
 * Usage:
 * @code
 * void expensive_function() {
 *     PROFILE_SCOPE("expensive_function");
 *     // ... code ...
 * }
 * @endcode
 */
#define PROFILE_SCOPE_CONCAT_(a, b) a##b
#define PROFILE_SCOPE_CONCAT(a, b) PROFILE_SCOPE_CONCAT_(a, b)
#define PROFILE_SCOPE(name) \
    ICARION::profiling::ProfileSection PROFILE_SCOPE_CONCAT(_profile_section_, __LINE__)(name)

/**
 * @brief Conditional profiling (only if profiler is enabled)
 * 
 * Lower overhead than PROFILE_SCOPE when profiler is disabled.
 */
#define PROFILE_SCOPE_IF_ENABLED(name) \
    std::unique_ptr<ICARION::profiling::ProfileSection> PROFILE_SCOPE_CONCAT(_profile_section_, __LINE__); \
    if (ICARION::profiling::Profiler::getInstance().isEnabled()) { \
        PROFILE_SCOPE_CONCAT(_profile_section_, __LINE__) = std::make_unique<ICARION::profiling::ProfileSection>(name); \
    }

} // namespace profiling
} // namespace ICARION
