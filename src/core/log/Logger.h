// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <memory>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <string>

namespace ICARION {
namespace log {

/**
 * @brief Modern structured logging for ICARION
 * 
 * Features:
 * - Automatic category-based loggers (config, integrator, hdf5, etc.)
 * - Console + File output (configurable)
 * - JSON format option for automated analysis
 * - Low overhead when disabled (compile-time checks on LOG_TIMER/LOG_DEBUG_IF)
 * - Thread-safe
 * - Rotating log files (prevents disk filling)
 * 
 * Usage:
 * @code
 * // Initialize once at program start
 * Logger::init("INFO", "simulation.log");
 * 
 * // Use category-specific loggers
 * Logger::config()->info("Loaded {} species", count);
 * Logger::integrator()->debug("Step {}, ion {}", step, ion_id);
 * Logger::hdf5()->error("Write failed: {}", e.what());
 * 
 * // Auto-log performance
 * LOG_TIMER("expensive_function");
 * 
 * // Shutdown before exit
 * Logger::shutdown();
 * @endcode
 */
class Logger {
public:
    /**
     * @brief Initialize logging system
     * 
     * @param level Log level: "DEBUG", "INFO", "WARN", "ERROR"
     * @param log_file Path to log file (empty = console only)
     * @param format Log format: "text" or "json" (default: "text")
     * @param max_file_size_mb Rotate log files after this size (default: 50 MB)
     * 
     * Call this once at program startup, before any logging.
     */
    static void init(
        const std::string& level = "INFO",
        const std::string& log_file = "",
        const std::string& format = "text",
        size_t max_file_size_mb = 50
    );
    
    /**
     * @brief Get or create category-specific logger
     * 
     * @param category Logger category (e.g., "config", "integrator")
     * @return Shared logger instance
     * 
     * Categories are automatically created on first access.
     * All loggers share the same sinks (console + file).
     */
    static std::shared_ptr<spdlog::logger> get(const std::string& category);
    
    /**
     * @brief Convenience shortcuts for common categories
     * 
     * Pre-defined categories:
     * - main: Program lifecycle, startup/shutdown
     * - config: Configuration loading and parsing
     * - integrator: Trajectory integration (use DEBUG sparingly!)
     * - hdf5: HDF5 I/O operations
     * - physics: Physics computations (collisions, reactions, fields)
     * - perf: Performance timing measurements
     * - gpu: GPU-related operations
     * - reactions: Chemical reaction events
     */
    static std::shared_ptr<spdlog::logger> main() { return get("main"); }
    static std::shared_ptr<spdlog::logger> config() { return get("config"); }
    static std::shared_ptr<spdlog::logger> integrator() { return get("integrator"); }
    static std::shared_ptr<spdlog::logger> hdf5() { return get("hdf5"); }
    static std::shared_ptr<spdlog::logger> physics() { return get("physics"); }
    static std::shared_ptr<spdlog::logger> perf() { return get("perf"); }
    static std::shared_ptr<spdlog::logger> gpu() { return get("gpu"); }
    static std::shared_ptr<spdlog::logger> reactions() { return get("reactions"); }
    
    /**
     * @brief Shutdown logging (flush buffers)
     * 
     * Call before program exit to ensure all logs are written.
     * After shutdown, no more logging should occur.
     */
    static void shutdown();
    
    /**
     * @brief Check if debug logging is enabled
     * 
     * Use for expensive debug operations:
     * @code
     * if (Logger::is_debug_enabled()) {
     *     std::string expensive_info = compute_diagnostics();
     *     Logger::integrator()->debug("{}", expensive_info);
     * }
     * @endcode
     */
    static bool is_debug_enabled();

private:
    static spdlog::level::level_enum parse_level(const std::string& level);
    
    static std::unordered_map<std::string, std::shared_ptr<spdlog::logger>> loggers_;
    static spdlog::level::level_enum current_level_;
    static std::vector<spdlog::sink_ptr> sinks_;
};

/**
 * @brief RAII timer for performance measurement
 * 
 * Automatically logs duration on destruction.
 * Logs to "perf" category at DEBUG level.
 * 
 * Usage:
 * @code
 * void expensive_function() {
 *     LOG_TIMER("expensive_function");
 *     // ... code ...
 *     // Duration logged automatically when function exits
 * }
 * @endcode
 */
class Timer {
public:
    /**
     * @brief Create timer and start measurement
     * @param name Name for this timing measurement
     */
    explicit Timer(const std::string& name);
    
    /**
     * @brief Destructor logs elapsed time
     */
    ~Timer();
    
private:
    std::string name_;
    std::chrono::steady_clock::time_point start_;
};

/**
 * @brief Macro for zero-overhead performance timing
 * 
 * In Release builds (NDEBUG defined), this is a no-op.
 * In Debug builds, creates a Timer that auto-logs on scope exit.
 * 
 * @code
 * void my_function() {
 *     LOG_TIMER("my_function");  // Only active in Debug builds
 *     // ... code ...
 * }
 * @endcode
 */
#ifdef NDEBUG
    #define LOG_TIMER(name) ((void)0)
#else
    #define LOG_TIMER(name) ICARION::log::Timer _log_timer_##__LINE__(name)
#endif

/**
 * @brief Conditional debug logging (zero overhead in Release)
 * 
 * @code
 * LOG_DEBUG_IF(step % 1000 == 0, "Checkpoint at step {}", step);
 * @endcode
 */
#ifdef NDEBUG
    #define LOG_DEBUG_IF(cond, msg, ...) ((void)0)
#else
    #define LOG_DEBUG_IF(cond, msg, ...) \
        if (cond) ICARION::log::Logger::integrator()->debug(msg, ##__VA_ARGS__)
#endif

/**
 * @brief Simple debug logger for legacy code compatibility
 * @param msg Message to log
 * 
 * Used by core computation modules for debugging output.
 * Only outputs when log level is DEBUG (--verbose or --log-level DEBUG).
 * 
 * Delegates to structured Logger system.
 * Falls back to stderr if Logger not yet initialized.
 */
void debug_log(const std::string& msg);

}  // namespace log
}  // namespace ICARION
