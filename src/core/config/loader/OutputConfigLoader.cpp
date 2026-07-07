// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "OutputConfigLoader.h"
#include "../conversion/EnumMapper.h"

namespace ICARION::config {

OutputConfig OutputConfigLoader::load(const Json::Value& json) {
    OutputConfig config;
    
    // Output folder
    if (json.isMember("folder") && json["folder"].isString()) {
        config.folder = json["folder"].asString();
    }
    
    // Trajectory file
    if (json.isMember("trajectory_file") && json["trajectory_file"].isString()) {
        config.trajectory_file = json["trajectory_file"].asString();
    }

    // Trajectory mode
    if (json.isMember("trajectory_mode") && json["trajectory_mode"].isString()) {
        config.trajectory_mode = json["trajectory_mode"].asString();
    }
    
    // Print progress
    if (json.isMember("print_progress") && json["print_progress"].isBool()) {
        config.print_progress = json["print_progress"].asBool();
    }

    // Optional RAM cap for trajectory buffer
    if (json.isMember("buffer_byte_cap") && json["buffer_byte_cap"].isNumeric()) {
        auto cap = json["buffer_byte_cap"].asInt64();
        if (cap < 0) {
            throw std::runtime_error("output.buffer_byte_cap must be >= 0");
        }
        config.buffer_byte_cap = static_cast<size_t>(cap);
    }

    if (json.isMember("deep_analysis") && json["deep_analysis"].isObject()) {
        const auto& da = json["deep_analysis"];
        if (da.isMember("mode") && da["mode"].isString()) {
            config.deep_analysis.mode_type = EnumMapper::parse_deep_analysis_mode(da["mode"].asString());
        }
        if (da.isMember("domain_filter_index") && da["domain_filter_index"].isInt()) {
            config.deep_analysis.domain_filter_index = da["domain_filter_index"].asInt();
        }
        if (da.isMember("sample_every_n") && da["sample_every_n"].isNumeric()) {
            auto n = da["sample_every_n"].asInt64();
            if (n <= 0) {
                throw std::runtime_error("output.deep_analysis.sample_every_n must be >= 1");
            }
            config.deep_analysis.sample_every_n = static_cast<size_t>(n);
        }
        if (da.isMember("max_events_per_ion") && da["max_events_per_ion"].isNumeric()) {
            auto n = da["max_events_per_ion"].asInt64();
            if (n < 0) {
                throw std::runtime_error("output.deep_analysis.max_events_per_ion must be >= 0");
            }
            config.deep_analysis.max_events_per_ion = static_cast<size_t>(n);
        }
    }

    if (json.isMember("deep_analysis_mode") && json["deep_analysis_mode"].isString()) {
        config.deep_analysis.mode_type = EnumMapper::parse_deep_analysis_mode(json["deep_analysis_mode"].asString());
    }
    if (json.isMember("deep_analysis_domain_filter_index") && json["deep_analysis_domain_filter_index"].isInt()) {
        config.deep_analysis.domain_filter_index = json["deep_analysis_domain_filter_index"].asInt();
    }
    if (json.isMember("deep_analysis_sample_every_n") && json["deep_analysis_sample_every_n"].isNumeric()) {
        auto n = json["deep_analysis_sample_every_n"].asInt64();
        if (n <= 0) {
            throw std::runtime_error("output.deep_analysis_sample_every_n must be >= 1");
        }
        config.deep_analysis.sample_every_n = static_cast<size_t>(n);
    }
    if (json.isMember("deep_analysis_max_events_per_ion") && json["deep_analysis_max_events_per_ion"].isNumeric()) {
        auto n = json["deep_analysis_max_events_per_ion"].asInt64();
        if (n < 0) {
            throw std::runtime_error("output.deep_analysis_max_events_per_ion must be >= 0");
        }
        config.deep_analysis.max_events_per_ion = static_cast<size_t>(n);
    }
    
    // Validate
    config.validate();

    config.deep_analysis.mode =
        EnumMapper::deep_analysis_mode_to_string(config.deep_analysis.mode_type);
    
    return config;
}

} // namespace ICARION::config
