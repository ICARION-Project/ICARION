// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file startupBanner.h
 * @brief Startup banner with basic system information (optional)
 * 
 * Prints an ASCII banner plus build/config/system info. Purely cosmetic; callers
 * should skip in batch/CI runs where stdout noise matters.
 */

#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <sys/utsname.h>
#include <unistd.h>
#include <fstream>

namespace ICARION {
namespace utils {

/**
 * @brief Print professional startup banner with system info
 * 
 * @param version ICARION version string
 * @param git_hash Git commit hash
 * @param config_file Configuration file path
 * @param log_level Log level string
 * @param log_file Optional log file path
 */
inline void print_startup_banner(
    const std::string& version,
    const std::string& git_hash,
    const std::string& config_file,
    const std::string& log_level,
    const std::string& log_file = ""
) {
    std::cout << "\n";
    // Lightweight banner; callers may disable in non-interactive contexts
    std::cout << "ICARION - Ion Collision And Reaction IntegratiON\n";
    std::cout << "High-Performance Trajectory Simulator\n";
    std::cout << "\n";
    
    // Version and build info
    std::cout << "   Version:      " << version << "\n";
    std::cout << "   Git Commit:   " << git_hash;
    
    // Get compile timestamp if available
    std::cout << " (" << __DATE__ << " " << __TIME__ << ")\n";
    
    // Build type and features
    std::cout << "   Build Type:   ";
    #ifdef NDEBUG
        std::cout << "Release";
    #else
        std::cout << "Debug";
    #endif
    
    #ifdef _OPENMP
        std::cout << " with OpenMP";
    #endif
    
    #ifdef ICARION_ENABLE_CUDA
        std::cout << " + CUDA";
    #endif
    
    std::cout << "\n";
    
    // Compiler info
    std::cout << "   Compiler:     ";
    #ifdef __GNUC__
        std::cout << "GCC " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
    #elif defined(__clang__)
        std::cout << "Clang " << __clang_major__ << "." << __clang_minor__;
    #elif defined(_MSC_VER)
        std::cout << "MSVC " << _MSC_VER;
    #else
        std::cout << "Unknown";
    #endif
    std::cout << "\n\n";
    
    std::cout << "   License:      GPL-3.0-only\n";
    std::cout << "   Support:      https://github.com/ICARION-Project/ICARION/issues\n";
    std::cout << "\n";
    
    std::cout << "────────────────────────────────────────────────────────────────────────────\n\n";
    
    // Configuration section
    std::cout << " Configuration\n\n";
    std::cout << "   File:         " << config_file << "\n";
    std::cout << "   Log Level:    " << log_level << "\n";
    if (!log_file.empty()) {
        std::cout << "   Log File:     " << log_file << "\n";
    }
    std::cout << "\n";
    
    std::cout << "────────────────────────────────────────────────────────────────────────────\n\n";
    
    // System information
    std::cout << " System Information\n\n";
    
    // Hostname
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    std::cout << "   Hostname:     " << hostname << "\n";
    
    // OS
    struct utsname sys_info;
    uname(&sys_info);
    std::cout << "   OS:           " << sys_info.sysname << " " << sys_info.release;
    
    // Try to get distro info
    std::ifstream os_release("/etc/os-release");
    if (os_release) {
        std::string line;
        while (std::getline(os_release, line)) {
            if (line.find("PRETTY_NAME=") == 0) {
                std::string distro = line.substr(13);
                distro.erase(std::remove(distro.begin(), distro.end(), '"'), distro.end());
                std::cout << " (" << distro << ")";
                break;
            }
        }
    }
    std::cout << "\n";
    
    // CPU info
    std::ifstream cpuinfo("/proc/cpuinfo");
    if (cpuinfo) {
        std::string line;
        std::string cpu_model;
        int cpu_cores = 0;
        while (std::getline(cpuinfo, line)) {
            if (line.find("model name") != std::string::npos && cpu_model.empty()) {
                cpu_model = line.substr(line.find(":") + 2);
            }
            if (line.find("processor") != std::string::npos) {
                cpu_cores++;
            }
        }
        if (!cpu_model.empty()) {
            std::cout << "   CPU:          " << cpu_model << " (" << cpu_cores << " cores)\n";
        }
    }
    
    // Memory
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo) {
        std::string line;
        long mem_total_kb = 0;
        while (std::getline(meminfo, line)) {
            if (line.find("MemTotal:") == 0) {
                std::istringstream(line.substr(10)) >> mem_total_kb;
                break;
            }
        }
        if (mem_total_kb > 0) {
            std::cout << "   Memory:       " << (mem_total_kb / 1024 / 1024) << " GB RAM\n";
        }
    }
    
    #ifdef ICARION_ENABLE_CUDA
        std::cout << "   GPU:          CUDA enabled\n";
    #endif
    
    std::cout << "\n";
    
    // Threading info
    #ifdef _OPENMP
        std::cout << "   Threads:      " << omp_get_max_threads() << " OpenMP threads\n";
    #else
        std::cout << "   Threads:      Single-threaded (OpenMP disabled)\n";
    #endif
    
    std::cout << "\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
    
    // Start timestamp
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm* tm_ptr = std::gmtime(&now_time);
    
    std::cout << " Starting simulation at " 
              << std::put_time(tm_ptr, "%Y-%m-%d %H:%M:%S UTC") << "\n";
    std::cout << "\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
}

/**
 * @brief Print completion summary with performance statistics
 * 
 * @param elapsed_s Wall time in seconds
 * @param num_ions Total number of ions
 * @param num_steps Total number of integration steps
 * @param output_file Output HDF5 file path
 * @param active_ions Number of ions still active
 * @param file_size_mb Approximate output file size in MB
 */
inline void print_completion_summary(
    double elapsed_s,
    size_t num_ions,
    size_t num_steps,
    const std::string& output_file,
    int active_ions = 0,
    double file_size_mb = 0.0
) {
    std::cout << "\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
    std::cout << " Simulation Completed Successfully\n";
    std::cout << "\n";
    std::cout << "────────────────────────────────────────────────────────────────────────────\n";
    std::cout << "\n";
    std::cout << " Performance\n\n";
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "   Wall Time:        " << elapsed_s << " seconds\n";
    
    if (num_steps > 0) {
        double steps_per_sec = num_steps / elapsed_s;
        double time_per_step_us = (elapsed_s / num_steps) * 1e6;
        
        std::cout << "   Steps/sec:        " << std::setprecision(0) << steps_per_sec << " steps/s\n";
        std::cout << "   Time per step:    " << std::setprecision(1) << time_per_step_us << " μs\n";
    }
    
    if (num_ions > 0 && num_steps > 0) {
        double ions_per_sec = (num_ions * num_steps) / elapsed_s;
        std::cout << "   Ion-steps/sec:    " << std::setprecision(0) << ions_per_sec << "\n";
    }
    
    std::cout << "\n";
    std::cout << "────────────────────────────────────────────────────────────────────────────\n";
    std::cout << "\n";
    std::cout << " Output\n\n";
    // Ensure .h5 extension is only shown once
    std::string output_display = output_file;
    if (output_display.size() < 3 || output_display.substr(output_display.size() - 3) != ".h5") {
        output_display += ".h5";
    }
    std::cout << "   Trajectory:       " << output_display;
    if (file_size_mb > 0) {
        std::cout << " (" << std::setprecision(1) << file_size_mb << " MB)";
    }
    std::cout << "\n";
    
    if (active_ions > 0) {
        std::cout << "\n";
        std::cout << "   Active ions:      " << active_ions << " ions still in simulation region\n";
    }
    
    std::cout << "\n";
    std::cout << "   HDF5 Structure:\n";
    std::cout << "     ├─ /metadata/      (build info, config, git hash)\n";
    std::cout << "     ├─ /trajectory/    (time, position, velocity)\n";
    std::cout << "     └─ /species/       (species properties)\n";
    
    std::cout << "\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
}

/**
 * @brief Print progress bar during simulation
 * 
 * @param current_step Current step number
 * @param total_steps Total number of steps
 * @param elapsed_s Elapsed time in seconds
 * @param active_ions Number of active ions
 * @param total_ions Total number of ions
 */
inline void print_progress_bar(
    size_t current_step,
    size_t total_steps,
    double elapsed_s,
    int active_ions = -1,
    int total_ions = -1
) {
    const int bar_width = 40;
    double progress = static_cast<double>(current_step) / total_steps;
    int filled = static_cast<int>(progress * bar_width);
    
    // Estimate time remaining
    double steps_per_sec = current_step / elapsed_s;
    double remaining_steps = total_steps - current_step;
    double eta_s = remaining_steps / steps_per_sec;
    
    // Print progress line (with carriage return to overwrite)
    std::cout << "\r";
    std::cout << " Step " << current_step << " / " << total_steps;
    std::cout << " (" << std::fixed << std::setprecision(1) << (progress * 100.0) << "%)  ";
    
    // Progress bar
    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < filled) {
            std::cout << "━";
        } else if (i == filled) {
            std::cout << "╸";
        } else {
            std::cout << "─";
        }
    }
    std::cout << "]  ";
    
    // ETA
    if (eta_s < 60) {
        std::cout << "ETA: " << std::setprecision(0) << eta_s << "s";
    } else if (eta_s < 3600) {
        std::cout << "ETA: " << std::setprecision(1) << (eta_s / 60) << "m";
    } else {
        std::cout << "ETA: " << std::setprecision(1) << (eta_s / 3600) << "h";
    }
    
    // Active ions (if provided)
    if (active_ions >= 0 && total_ions > 0) {
        std::cout << "  Active: " << active_ions << "/" << total_ions;
    }
    
    std::cout << std::flush;
}

} // namespace utils
} // namespace ICARION
