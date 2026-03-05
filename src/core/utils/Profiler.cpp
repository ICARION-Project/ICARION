// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "Profiler.h"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ICARION {
namespace profiling {

// === ProfileSection Implementation ===

ProfileSection::ProfileSection(const std::string& name)
    : name_(name), start_(std::chrono::steady_clock::now()) {
    Profiler::getInstance().startSection(name_);
}

ProfileSection::~ProfileSection() {
    if (Profiler::getInstance().isEnabled()) {
        Profiler::getInstance().endSection(name_);
    }
}

// === Profiler Implementation ===

Profiler& Profiler::getInstance() {
    static Profiler instance;
    return instance;
}

void Profiler::enable(bool enable) {
    std::lock_guard<std::mutex> lock(mutex_);
    enabled_ = enable;
    if (enable) {
        sections_.clear();
    }
}

void Profiler::startSection(const std::string& name) {
    if (!enabled_) return;
    
    std::lock_guard<std::mutex> lock(mutex_);
    auto& data = sections_[name];
    auto& starts = data.start_times[std::this_thread::get_id()];
    starts.push_back(std::chrono::steady_clock::now());
}

void Profiler::endSection(const std::string& name) {
    if (!enabled_) return;
    
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = sections_.find(name);
    if (it == sections_.end()) {
        return;  // Section not started
    }
    
    auto& data = it->second;
    auto thread_it = data.start_times.find(std::this_thread::get_id());
    if (thread_it == data.start_times.end() || thread_it->second.empty()) {
        return;  // Section not started on this thread
    }

    auto start = thread_it->second.back();
    thread_it->second.pop_back();
    if (thread_it->second.empty()) {
        data.start_times.erase(thread_it);
    }

    auto end = std::chrono::steady_clock::now();
    auto duration_ms = std::chrono::duration<double, std::milli>(
        end - start).count();
    
    data.total_time_ms += duration_ms;
    data.call_count++;
    data.min_ms = std::min(data.min_ms, duration_ms);
    data.max_ms = std::max(data.max_ms, duration_ms);
}

void Profiler::reset() {
    std::lock_guard<std::mutex> lock(mutex_);
    sections_.clear();
}

std::vector<TimingResult> Profiler::getResults() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::vector<TimingResult> results;
    results.reserve(sections_.size());
    
    for (const auto& [name, data] : sections_) {
        TimingResult result;
        result.name = name;
        result.duration_ms = data.total_time_ms;
        result.call_count = data.call_count;
        result.min_ms = data.min_ms;
        result.max_ms = data.max_ms;
        result.avg_ms = data.call_count > 0 ? 
            data.total_time_ms / data.call_count : 0.0;
        results.push_back(result);
    }
    
    // Sort by total time (descending)
    std::sort(results.begin(), results.end(),
        [](const TimingResult& a, const TimingResult& b) {
            return a.duration_ms > b.duration_ms;
        });
    
    return results;
}

void Profiler::printSummary(std::ostream& os) const {
    auto results = getResults();
    
    if (results.empty()) {
        os << "\nNo profiling data collected.\n";
        return;
    }
    
    // Calculate total time
    double total_ms = 0.0;
    for (const auto& r : results) {
        total_ms += r.duration_ms;
    }
    
    os << "\n";
    os << "ICARION Performance Profile (profiling must be enabled explicitly)\n\n";
    
    os << std::fixed << std::setprecision(3);
    os << std::left;
    
    // Header
    os << std::setw(40) << "Section"
       << std::right
       << std::setw(14) << "Total (ms)"
       << std::setw(12) << "Calls"
       << std::setw(14) << "Avg (ms)"
       << std::setw(14) << "Min (ms)"
       << std::setw(14) << "Max (ms)"
       << std::setw(10) << "%\n";
    os << std::string(82, '-') << "\n";
    
    os << std::left;
    
    // Results
    for (const auto& r : results) {
        double percent = (total_ms > 0.0) ? (r.duration_ms / total_ms) * 100.0 : 0.0;
        
        os << std::setw(40) << r.name
           << std::right
           << std::setw(14) << std::setprecision(3) << r.duration_ms
           << std::setw(12) << r.call_count
           << std::setw(14) << std::setprecision(3) << r.avg_ms
           << std::setw(14) << std::setprecision(3) << r.min_ms
           << std::setw(14) << std::setprecision(3) << r.max_ms
           << std::setw(9) << std::setprecision(1) << percent << "%\n";
        os << std::left;
    }
    
    os << std::string(82, '-') << "\n";
    os << std::setw(40) << "TOTAL"
       << std::right
       << std::setw(14) << std::setprecision(3) << total_ms << "\n\n";
}

void Profiler::exportJSON(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Failed to open profile output: " + filename);
    }
    
    auto results = getResults();
    
    file << "{\n";
    file << "  \"profiling_results\": [\n";
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        file << "    {\n";
        file << "      \"name\": \"" << r.name << "\",\n";
        file << "      \"total_ms\": " << r.duration_ms << ",\n";
        file << "      \"calls\": " << r.call_count << ",\n";
        file << "      \"avg_ms\": " << r.avg_ms << ",\n";
        file << "      \"min_ms\": " << r.min_ms << ",\n";
        file << "      \"max_ms\": " << r.max_ms << "\n";
        file << "    }" << (i + 1 < results.size() ? "," : "") << "\n";
    }
    
    file << "  ]\n";
    file << "}\n";
    
    file.close();
}

void Profiler::exportCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Failed to open profile output: " + filename);
    }
    
    auto results = getResults();
    
    file << "Section,Total_ms,Calls,Avg_ms,Min_ms,Max_ms\n";
    for (const auto& r : results) {
        file << r.name << ","
             << r.duration_ms << ","
             << r.call_count << ","
             << r.avg_ms << ","
             << r.min_ms << ","
             << r.max_ms << "\n";
    }
    
    file.close();
}

} // namespace profiling
} // namespace ICARION
