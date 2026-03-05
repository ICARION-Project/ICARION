// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "numericalSafetyLogger.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace ICARION {
namespace safety {

// Static member initialization
std::unique_ptr<NumericalSafetyLogger> NumericalSafetyLogger::instance_ = nullptr;
std::mutex NumericalSafetyLogger::instance_mutex_;

NumericalSafetyLogger::NumericalSafetyLogger(const std::string& log_file, bool verbose)
    : log_file_path_(log_file), enabled_(true), verbose_mode_(verbose),
      max_history_size_(10000), buffer_size_(1000), buffer_count_(0) {
    
    // Initialize log file (best-effort; disable if it fails)
    log_file_.open(log_file_path_, std::ios::out | std::ios::app);
    if (!log_file_.is_open()) {
        std::cerr << "Warning: Could not open numerical safety log file: " << log_file_path_ << std::endl;
        enabled_ = false;
        return;
    }
    
    // Write header
    log_file_ << "\n=== ICARION Numerical Safety Log Session Started ===\n";
    log_file_ << "Timestamp: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count() << " ms\n";
    log_file_ << "Verbose Mode: " << (verbose_mode_ ? "Enabled" : "Disabled") << "\n";
    log_file_ << "Buffer Size: " << buffer_size_ << "\n";
    log_file_ << "Note: Logging is intended for debugging and may impact performance.\n";
    log_file_ << "=======================================================\n\n";
    log_file_.flush();
    
    // Reserve buffer space
    log_buffer_.reserve(buffer_size_);
    violation_history_.reserve(max_history_size_);
}

NumericalSafetyLogger& NumericalSafetyLogger::getInstance(const std::string& log_file, bool verbose) {
    std::lock_guard<std::mutex> lock(instance_mutex_);
    if (!instance_) {
        instance_ = std::unique_ptr<NumericalSafetyLogger>(new NumericalSafetyLogger(log_file, verbose));
    }
    return *instance_;
}

void NumericalSafetyLogger::configure(bool enable_logging, bool verbose_mode, 
                                    size_t buffer_size, size_t max_history) {
    std::lock_guard<std::mutex> lock(log_mutex_);
    enabled_ = enable_logging;
    verbose_mode_ = verbose_mode;
    buffer_size_ = buffer_size;
    max_history_size_ = max_history;
    
    // Adjust containers
    log_buffer_.reserve(buffer_size_);
    violation_history_.reserve(max_history_size_);
    
    if (enabled_ && log_file_.is_open()) {
        log_file_ << "[CONFIG] Logger reconfigured - Enabled: " << enabled_ 
                 << ", Verbose: " << verbose_mode_ << ", Buffer: " << buffer_size_ << "\n";
        log_file_.flush();
    }
}

void NumericalSafetyLogger::logViolation(const ViolationEvent& event) {
    if (!enabled_) return;
    
    std::lock_guard<std::mutex> lock(log_mutex_);
    
    // Update statistics
    statistics_.total_violations++;
    statistics_.violation_counts[event.type]++;
    statistics_.violation_times[event.type] += event.safety_check_time_ns;
    
    switch (event.type) {
        case ViolationType::NAN_POSITION:
        case ViolationType::NAN_VELOCITY:
        case ViolationType::INF_POSITION:
        case ViolationType::INF_VELOCITY:
            statistics_.nan_violations++;
            break;
        case ViolationType::BOUNDS_POSITION:
        case ViolationType::BOUNDS_VELOCITY:
        case ViolationType::BOUNDS_ACCELERATION:
            statistics_.bounds_violations++;
            break;
        case ViolationType::REJECTED_STEP:
        case ViolationType::TIMESTEP_TOO_SMALL:
            statistics_.rejected_steps++;
            break;
    }
    
    if (event.recovery_attempted) statistics_.recovery_attempts++;
    if (event.recovery_successful) statistics_.successful_recoveries++;
    
    // Store in history (with size management)
    if (violation_history_.size() >= max_history_size_) {
        violation_history_.erase(violation_history_.begin());
    }
    violation_history_.push_back(event);
    
    // Format and buffer log entry
    std::string log_entry = formatViolationEvent(event);
    addToBuffer(log_entry);
}

void NumericalSafetyLogger::logRejectedStep(int ion_index, int step_number, double simulation_time,
                                           double attempted_timestep, double new_timestep, 
                                           const std::string& reason) {
    if (!enabled_) return;
    
    ViolationEvent event;
    event.type = ViolationType::REJECTED_STEP;
    event.timestamp = std::chrono::steady_clock::now();
    event.ion_index = ion_index;
    event.step_number = step_number;
    event.simulation_time = simulation_time;
    event.timestep = attempted_timestep;
    event.violation_context = reason;
    event.violation_magnitude = attempted_timestep / new_timestep; // rejection factor
    event.recovery_attempted = true;
    event.recovery_successful = true; // step size reduction is always successful
    event.safety_check_time_ns = 0.0; // Not applicable for step rejection
    
    logViolation(event);
    
    // Additional specific logging for rejected steps
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6);
    ss << "[REJECTED_STEP] Ion " << ion_index << ", Step " << step_number 
       << ", Time " << simulation_time << "s, dt " << attempted_timestep 
       << " -> " << new_timestep << " (factor: " << (new_timestep/attempted_timestep)
       << "), Reason: " << reason << "\n";
    
    addToBuffer(ss.str());
}

void NumericalSafetyLogger::logNaNDetection(int ion_index, int step_number, double simulation_time,
                                           const Vec3& position, const Vec3& velocity,
                                           const std::string& context) {
    if (!enabled_) return;
    
    ViolationEvent event;
    event.type = std::isnan(position.x) || std::isnan(position.y) || std::isnan(position.z) ? 
                 ViolationType::NAN_POSITION : ViolationType::NAN_VELOCITY;
    event.timestamp = std::chrono::steady_clock::now();
    event.ion_index = ion_index;
    event.step_number = step_number;
    event.simulation_time = simulation_time;
    event.position = position;
    event.velocity = velocity;
    event.violation_context = context;
    event.recovery_attempted = false; // NaN recovery depends on configuration
    event.recovery_successful = false;
    event.safety_check_time_ns = 0.0; // Will be updated by caller
    
    logViolation(event);
}

void NumericalSafetyLogger::logBoundsViolation(int ion_index, int step_number, double simulation_time,
                                              const Vec3& position, const Vec3& velocity,
                                              ViolationType violation_type, double magnitude) {
    if (!enabled_) return;
    
    ViolationEvent event;
    event.type = violation_type;
    event.timestamp = std::chrono::steady_clock::now();
    event.ion_index = ion_index;
    event.step_number = step_number;
    event.simulation_time = simulation_time;
    event.position = position;
    event.velocity = velocity;
    event.violation_magnitude = magnitude;
    event.violation_context = "Bounds exceeded";
    event.recovery_attempted = false; // Depends on configuration
    event.recovery_successful = false;
    event.safety_check_time_ns = 0.0; // Will be updated by caller
    
    logViolation(event);
}

void NumericalSafetyLogger::updateStatistics(double safety_check_time_ns, bool violation_found) {
    (void)violation_found;
    if (!enabled_) return;
    
    std::lock_guard<std::mutex> lock(log_mutex_);
    statistics_.total_checks++;
    statistics_.total_safety_time_ns += safety_check_time_ns;
    statistics_.max_safety_time_ns = std::max(statistics_.max_safety_time_ns, safety_check_time_ns);
    statistics_.avg_safety_time_ns = statistics_.total_safety_time_ns / statistics_.total_checks;
}

void NumericalSafetyLogger::generateSafetyReport(const std::string& report_file) {
    if (!enabled_) return;
    
    std::lock_guard<std::mutex> lock(log_mutex_);
    
    std::ofstream report(report_file);
    if (!report.is_open()) {
        std::cerr << "Error: Could not create safety report file: " << report_file << std::endl;
        return;
    }
    
    report << "=== ICARION Numerical Safety Report ===\n\n";
    
    // Overall statistics
    report << "OVERALL STATISTICS:\n";
    report << "  Total Safety Checks: " << statistics_.total_checks << "\n";
    report << "  Total Violations: " << statistics_.total_violations << "\n";
    report << "  Violation Rate: " << (statistics_.total_checks > 0 ? 
        100.0 * statistics_.total_violations / statistics_.total_checks : 0.0) << "%\n\n";
    
    // Violation breakdown
    report << "VIOLATION BREAKDOWN:\n";
    report << "  NaN/Inf Violations: " << statistics_.nan_violations << "\n";
    report << "  Bounds Violations: " << statistics_.bounds_violations << "\n";
    report << "  Rejected Steps: " << statistics_.rejected_steps << "\n\n";
    
    // Recovery statistics
    report << "RECOVERY STATISTICS:\n";
    report << "  Recovery Attempts: " << statistics_.recovery_attempts << "\n";
    report << "  Successful Recoveries: " << statistics_.successful_recoveries << "\n";
    report << "  Recovery Success Rate: " << (statistics_.recovery_attempts > 0 ?
        100.0 * statistics_.successful_recoveries / statistics_.recovery_attempts : 0.0) << "%\n\n";
    
    // Performance impact
    report << "PERFORMANCE IMPACT:\n";
    report << "  Total Safety Time: " << statistics_.total_safety_time_ns * 1e-6 << " ms\n";
    report << "  Average Safety Time: " << statistics_.avg_safety_time_ns << " ns/check\n";
    report << "  Maximum Safety Time: " << statistics_.max_safety_time_ns << " ns/check\n";
    report << "  Safety Overhead: " << (statistics_.total_safety_time_ns / 1e9) << " seconds\n\n";
    
    // Per-type statistics
    report << "PER-TYPE STATISTICS:\n";
    for (const auto& [type, count] : statistics_.violation_counts) {
        if (count > 0) {
            report << "  " << violationTypeToString(type) << ": " << count;
            auto time_it = statistics_.violation_times.find(type);
            if (time_it != statistics_.violation_times.end() && count > 0) {
                report << " (avg time: " << (time_it->second / count) << " ns)";
            }
            report << "\n";
        }
    }
    
    // Recent violation history (last 10)
    if (!violation_history_.empty()) {
        report << "\nRECENT VIOLATIONS (Last 10):\n";
        size_t start_idx = violation_history_.size() > 10 ? violation_history_.size() - 10 : 0;
        for (size_t i = start_idx; i < violation_history_.size(); ++i) {
            report << "  " << formatViolationEvent(violation_history_[i]) << "\n";
        }
    }
    
    report << "\n=== End Report ===\n";
    report.close();
    
    std::cout << "Numerical safety report generated: " << report_file << std::endl;
}

void NumericalSafetyLogger::flush() {
    std::lock_guard<std::mutex> lock(log_mutex_);
    writeBuffer();
}

SafetyStatistics NumericalSafetyLogger::getStatistics() const {
    std::lock_guard<std::mutex> lock(log_mutex_);
    return statistics_;
}

void NumericalSafetyLogger::reset() {
    std::lock_guard<std::mutex> lock(log_mutex_);
    statistics_ = SafetyStatistics{};
    violation_history_.clear();
    log_buffer_.clear();
    buffer_count_ = 0;
    
    if (log_file_.is_open()) {
        log_file_ << "[RESET] Safety statistics and history cleared\n";
        log_file_.flush();
    }
}

NumericalSafetyLogger::~NumericalSafetyLogger() {
    if (enabled_ && log_file_.is_open()) {
        // Final flush and statistics summary
        flush();
        
        log_file_ << "\n=== ICARION Numerical Safety Log Session Ended ===\n";
        log_file_ << "Total Checks: " << statistics_.total_checks << "\n";
        log_file_ << "Total Violations: " << statistics_.total_violations << "\n";
        log_file_ << "Total Safety Time: " << statistics_.total_safety_time_ns * 1e-6 << " ms\n";
        log_file_ << "===================================================\n\n";
        
        log_file_.close();
    }
}

void NumericalSafetyLogger::writeBuffer() {
    if (!enabled_ || !log_file_.is_open() || log_buffer_.empty()) return;
    
    for (const auto& entry : log_buffer_) {
        log_file_ << entry;
    }
    log_file_.flush();
    
    log_buffer_.clear();
    buffer_count_ = 0;
}

std::string NumericalSafetyLogger::formatViolationEvent(const ViolationEvent& event) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(9);
    
    // Timestamp
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        event.timestamp.time_since_epoch()).count();
    
    ss << "[" << time_ms << "ms] ";
    ss << violationTypeToString(event.type) << " - ";
    ss << "Ion " << event.ion_index << ", Step " << event.step_number << ", ";
    ss << "Time " << event.simulation_time << "s";
    
    if (verbose_mode_) {
        ss << "\n  Position: (" << event.position.x << ", " << event.position.y << ", " << event.position.z << ")";
        ss << "\n  Velocity: (" << event.velocity.x << ", " << event.velocity.y << ", " << event.velocity.z << ")";
        ss << "\n  Context: " << event.violation_context;
        if (event.violation_magnitude > 0.0) {
            ss << "\n  Magnitude: " << event.violation_magnitude;
        }
        if (event.recovery_attempted) {
            ss << "\n  Recovery: " << (event.recovery_successful ? "Success" : "Failed");
        }
        ss << "\n  Check Time: " << event.safety_check_time_ns << "ns";
    } else {
        ss << ", Context: " << event.violation_context;
    }
    
    return ss.str() + "\n";
}

std::string NumericalSafetyLogger::violationTypeToString(ViolationType type) {
    switch (type) {
        case ViolationType::NAN_POSITION: return "NaN_POSITION";
        case ViolationType::NAN_VELOCITY: return "NaN_VELOCITY";
        case ViolationType::INF_POSITION: return "INF_POSITION";
        case ViolationType::INF_VELOCITY: return "INF_VELOCITY";
        case ViolationType::BOUNDS_POSITION: return "BOUNDS_POSITION";
        case ViolationType::BOUNDS_VELOCITY: return "BOUNDS_VELOCITY";
        case ViolationType::BOUNDS_ACCELERATION: return "BOUNDS_ACCELERATION";
        case ViolationType::REJECTED_STEP: return "REJECTED_STEP";
        case ViolationType::TIMESTEP_TOO_SMALL: return "TIMESTEP_TOO_SMALL";
        default: return "UNKNOWN";
    }
}

void NumericalSafetyLogger::addToBuffer(const std::string& entry) {
    log_buffer_.push_back(entry);
    buffer_count_++;
    
    // Auto-flush when buffer is full
    if (buffer_count_ >= buffer_size_) {
        writeBuffer();
    }
}

} // namespace safety
} // namespace ICARION
