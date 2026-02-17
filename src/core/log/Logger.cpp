// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include "Logger.h"
#include <spdlog/pattern_formatter.h>
#include <algorithm>
#include <iostream>

namespace ICARION {
namespace log {

// === Static member initialization ===
std::unordered_map<std::string, std::shared_ptr<spdlog::logger>> Logger::loggers_;
spdlog::level::level_enum Logger::current_level_ = spdlog::level::info;
std::vector<spdlog::sink_ptr> Logger::sinks_;

// === Logger Implementation ===

void Logger::init(
    const std::string& level,
    const std::string& log_file,
    const std::string& format,
    size_t max_file_size_mb
) {
    // Parse log level
    current_level_ = parse_level(level);
    
    // Parse format (text or json)
    bool json_format = (format == "json");
    
    // Clear existing sinks and loggers
    sinks_.clear();
    loggers_.clear();
    
    // === Console Sink (colored output) ===
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    
    if (json_format) {
        // JSON format for machine processing
        // {"time":"2025-01-15T14:30:22.123","level":"info","cat":"config","msg":"Loading species"}
        console_sink->set_pattern(R"({"time":"%Y-%m-%dT%H:%M:%S.%e","level":"%l","cat":"%n","msg":"%v"})");
    } else {
        // Human-readable format
        // [2025-01-15 14:30:22.123] [config] [info] Loading species
        console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%n] [%^%l%$] %v");
    }
    
    console_sink->set_level(current_level_);
    sinks_.push_back(console_sink);
    
    // === File Sink (optional, rotating) ===
    if (!log_file.empty()) {
        size_t max_size = max_file_size_mb * 1024 * 1024;  // MB to bytes
        size_t max_files = 3;  // Keep 3 rotated files (log, log.1, log.2)
        
        try {
            auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
                log_file, max_size, max_files
            );
            
            if (json_format) {
                file_sink->set_pattern(R"({"time":"%Y-%m-%dT%H:%M:%S.%e","level":"%l","cat":"%n","msg":"%v"})");
            } else {
                // File format includes full timestamps (no color codes)
                file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%n] [%l] %v");
            }
            
            file_sink->set_level(current_level_);
            sinks_.push_back(file_sink);
            
            std::cout << "Logging to file: " << log_file 
                      << " (max size: " << max_file_size_mb << " MB, rotating)\n";
        } catch (const spdlog::spdlog_ex& ex) {
            std::cerr << "Failed to create log file sink: " << ex.what() << "\n";
        }
    }
    
    // Create default logger with shared sinks
    auto default_logger = std::make_shared<spdlog::logger>("default", sinks_.begin(), sinks_.end());
    default_logger->set_level(current_level_);
    default_logger->flush_on(spdlog::level::err);  // Auto-flush on errors
    spdlog::set_default_logger(default_logger);
    
    // Set global flush interval (every 3 seconds)
    spdlog::flush_every(std::chrono::seconds(3));
    
    // Log initialization message
    get("main")->info("ICARION Logging initialized (level: {}, format: {}{})", 
                      level,
                      json_format ? "json" : "text",
                      log_file.empty() ? "" : ", file: " + log_file);
}

std::shared_ptr<spdlog::logger> Logger::get(const std::string& category) {
    // Check if logger already exists
    auto it = loggers_.find(category);
    if (it != loggers_.end()) {
        return it->second;
    }
    
    // Create new logger with shared sinks
    auto logger = std::make_shared<spdlog::logger>(category, sinks_.begin(), sinks_.end());
    logger->set_level(current_level_);
    logger->flush_on(spdlog::level::err);  // Auto-flush on errors
    
    // Store for future access
    loggers_[category] = logger;
    
    return logger;
}

void Logger::shutdown() {
    // Flush all loggers
    spdlog::shutdown();
}

bool Logger::is_debug_enabled() {
    return current_level_ <= spdlog::level::debug;
}

spdlog::level::level_enum Logger::parse_level(const std::string& level) {
    // Convert to lowercase for comparison
    std::string lower = level;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    if (lower == "debug" || lower == "d") {
        return spdlog::level::debug;
    }
    if (lower == "info" || lower == "i") {
        return spdlog::level::info;
    }
    if (lower == "warn" || lower == "warning" || lower == "w") {
        return spdlog::level::warn;
    }
    if (lower == "error" || lower == "err" || lower == "e") {
        return spdlog::level::err;
    }
    if (lower == "critical" || lower == "crit" || lower == "c") {
        return spdlog::level::critical;
    }
    
    std::cerr << "Warning: Unknown log level '" << level << "', using INFO\n";
    return spdlog::level::info;
}

// === Timer Implementation ===

Timer::Timer(const std::string& name)
    : name_(name), start_(std::chrono::steady_clock::now()) {
}

Timer::~Timer() {
    auto end = std::chrono::steady_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_).count();
    
    // Log to "perf" category at DEBUG level
    Logger::perf()->debug("{} took {} ms", name_, duration_ms);
}

// Legacy debug_log function for backward compatibility
void debug_log(const std::string& msg) {
    auto logger = Logger::main();
    if (logger) {
        logger->debug("{}", msg);
    } else {
        // Fallback for early initialization (before Logger::init())
        std::cerr << "[DEBUG] " << msg << std::endl;
    }
}

}  // namespace log
}  // namespace ICARION
