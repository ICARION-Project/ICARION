// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_NUMERICAL_SAFETY_LOGGER_H
#define ICARION_NUMERICAL_SAFETY_LOGGER_H


#include <fstream>
#include <string>
#include <chrono>
#include <mutex>
#include <memory>
#include <vector>
#include <unordered_map>
#include "core/utils/mathUtils.h"

namespace ICARION {
namespace safety {

// Use the Vec3 from ICARION::core
using ICARION::core::Vec3;

/**
 * @brief Violation type classification for detailed logging
 */
enum class ViolationType {
    NAN_POSITION,
    NAN_VELOCITY, 
    INF_POSITION,
    INF_VELOCITY,
    BOUNDS_POSITION,
    BOUNDS_VELOCITY,
    BOUNDS_ACCELERATION,
    REJECTED_STEP,
    TIMESTEP_TOO_SMALL
};

/**
 * @brief Detailed information about a numerical violation event
 */
struct ViolationEvent {
    ViolationType type;
    std::chrono::steady_clock::time_point timestamp;
    int ion_index;
    int step_number;
    double simulation_time;
    double timestep;
    
    // Ion state at violation
    Vec3 position;
    Vec3 velocity;
    Vec3 acceleration;
    
    // Violation-specific data
    std::string violation_context;
    double violation_magnitude;
    bool recovery_attempted;
    bool recovery_successful;
    
    // Performance impact
    double safety_check_time_ns;
};

/**
 * @brief Statistical summary of numerical safety performance
 */
struct SafetyStatistics {
    size_t total_checks = 0;
    size_t total_violations = 0;
    size_t nan_violations = 0;
    size_t inf_violations = 0;
    size_t bounds_violations = 0;
    size_t rejected_steps = 0;
    size_t recovery_attempts = 0;
    size_t successful_recoveries = 0;
    
    double total_safety_time_ns = 0.0;
    double max_safety_time_ns = 0.0;
    double avg_safety_time_ns = 0.0;
    
    // Per-violation type statistics
    std::unordered_map<ViolationType, size_t> violation_counts;
    std::unordered_map<ViolationType, double> violation_times;
};

/**
 * @brief Numerical safety logger with buffered writes (debug use)
 * 
 * Buffered logging and basic statistics for debugging safety issues. Not intended
 * for performance-critical production runs.
 */
class NumericalSafetyLogger {
private:
    std::string log_file_path_;
    std::ofstream log_file_;
    mutable std::mutex log_mutex_;
    bool enabled_;
    bool verbose_mode_;
    
    // Performance tracking
    SafetyStatistics statistics_;
    std::vector<ViolationEvent> violation_history_;
    size_t max_history_size_;
    
    // Buffered logging for performance
    std::vector<std::string> log_buffer_;
    size_t buffer_size_;
    size_t buffer_count_;
    
    // Singleton pattern for global access
    static std::unique_ptr<NumericalSafetyLogger> instance_;
    static std::mutex instance_mutex_;
    
    NumericalSafetyLogger(const std::string& log_file, bool verbose = false);

public:
    /**
     * @brief Get singleton logger instance
     */
    static NumericalSafetyLogger& getInstance(const std::string& log_file = "numerical_safety.log", 
                                            bool verbose = false);
    
    /**
     * @brief Configure logger settings
     */
    void configure(bool enable_logging, bool verbose_mode, size_t buffer_size = 1000, 
                   size_t max_history = 10000);
    
    /**
     * @brief Log a numerical violation event with full context
     */
    void logViolation(const ViolationEvent& event);
    
    /**
     * @brief Log a rejected integration step
     */
    void logRejectedStep(int ion_index, int step_number, double simulation_time,
                        double attempted_timestep, double new_timestep, 
                        const std::string& reason);
    
    /**
     * @brief Log NaN detection with ion state
     */
    void logNaNDetection(int ion_index, int step_number, double simulation_time,
                        const Vec3& position, const Vec3& velocity,
                        const std::string& context);
    
    /**
     * @brief Log bounds violation with detailed analysis
     */
    void logBoundsViolation(int ion_index, int step_number, double simulation_time,
                           const Vec3& position, const Vec3& velocity,
                           ViolationType violation_type, double magnitude);
    
    /**
     * @brief Update performance statistics
     */
    void updateStatistics(double safety_check_time_ns, bool violation_found);
    
    /**
     * @brief Generate comprehensive safety report
     */
    void generateSafetyReport(const std::string& report_file = "safety_report.txt");
    
    /**
     * @brief Flush buffered log entries
     */
    void flush();
    
    /**
     * @brief Get current safety statistics
     */
    SafetyStatistics getStatistics() const;
    
    /**
     * @brief Reset all statistics and history
     */
    void reset();
    
    /**
     * @brief Enable/disable logging at runtime
     */
    void setEnabled(bool enabled) { enabled_ = enabled; }
    bool isEnabled() const { return enabled_; }
    
    /**
     * @brief Destructor ensures proper cleanup
     */
    ~NumericalSafetyLogger();

private:
    /**
     * @brief Write buffered entries to file
     */
    void writeBuffer();
    
    /**
     * @brief Format violation event for logging
     */
    std::string formatViolationEvent(const ViolationEvent& event);
    
    /**
     * @brief Convert violation type to string
     */
    std::string violationTypeToString(ViolationType type);
    
    /**
     * @brief Add entry to buffer with automatic flushing
     */
    void addToBuffer(const std::string& entry);
};

/**
 * @brief Convenience macros for logging with automatic timing
 */
#define ICARION_SAFETY_LOG_VIOLATION(event) \
    do { \
        if (ICARION::safety::NumericalSafetyLogger::getInstance().isEnabled()) { \
            auto start_time = std::chrono::steady_clock::now(); \
            ICARION::safety::NumericalSafetyLogger::getInstance().logViolation(event); \
            auto end_time = std::chrono::steady_clock::now(); \
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count(); \
            ICARION::safety::NumericalSafetyLogger::getInstance().updateStatistics(duration, true); \
        } \
    } while(0)

#define ICARION_SAFETY_LOG_NAN(ion_idx, step_num, sim_time, pos, vel, context) \
    do { \
        if (ICARION::safety::NumericalSafetyLogger::getInstance().isEnabled()) { \
            ICARION::safety::NumericalSafetyLogger::getInstance().logNaNDetection( \
                ion_idx, step_num, sim_time, pos, vel, context); \
        } \
    } while(0)

#define ICARION_SAFETY_LOG_BOUNDS(ion_idx, step_num, sim_time, pos, vel, type, mag) \
    do { \
        if (ICARION::safety::NumericalSafetyLogger::getInstance().isEnabled()) { \
            ICARION::safety::NumericalSafetyLogger::getInstance().logBoundsViolation( \
                ion_idx, step_num, sim_time, pos, vel, type, mag); \
        } \
    } while(0)

#define ICARION_SAFETY_LOG_REJECTED_STEP(ion_idx, step_num, sim_time, old_dt, new_dt, reason) \
    do { \
        if (ICARION::safety::NumericalSafetyLogger::getInstance().isEnabled()) { \
            ICARION::safety::NumericalSafetyLogger::getInstance().logRejectedStep( \
                ion_idx, step_num, sim_time, old_dt, new_dt, reason); \
        } \
    } while(0)

/**
 * @brief RAII-style performance timer for safety checks
 */
class SafetyTimer {
private:
    std::chrono::steady_clock::time_point start_time_;
    bool violation_found_;
    
public:
    SafetyTimer() : start_time_(std::chrono::steady_clock::now()), violation_found_(false) {}
    
    void markViolation() { violation_found_ = true; }
    
    ~SafetyTimer() {
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time_).count();
        NumericalSafetyLogger::getInstance().updateStatistics(duration, violation_found_);
    }
};

#define ICARION_SAFETY_TIMER() ICARION::safety::SafetyTimer timer_##__LINE__
#define ICARION_SAFETY_MARK_VIOLATION() timer_##__LINE__.markViolation()

} // namespace safety
} // namespace ICARION

#endif // ICARION_NUMERICAL_SAFETY_LOGGER_H
