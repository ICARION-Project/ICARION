// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "EHSSSamples.h"
#include <json/json.h>
#include <fstream>

namespace ICARION::physics {

namespace {

bool load_json_file(const std::filesystem::path& path, Json::Value& root, std::string& error_msg) {
    std::ifstream file(path);
    if (!file.is_open()) {
        error_msg = "Failed to open EHSS samples file: " + path.string();
        return false;
    }

    Json::CharReaderBuilder builder;
    std::string errs;
    if (!Json::parseFromStream(builder, file, &root, &errs)) {
        error_msg = "Failed to parse EHSS samples file: " + errs;
        return false;
    }

    return true;
}

}  // namespace

bool load_ehss_samples_file(
    const std::filesystem::path& path,
    EHSSOrientationSamples& out,
    std::string* error_msg
) {
    Json::Value root;
    std::string err;
    if (!load_json_file(path, root, err)) {
        if (error_msg) {
            *error_msg = err;
        }
        return false;
    }

    if (!root.isObject()) {
        if (error_msg) {
            *error_msg = "EHSS samples root must be a JSON object";
        }
        return false;
    }

    const Json::Value& orientations = root["orientations_quat"];
    if (!orientations.isArray() || orientations.empty()) {
        if (error_msg) {
            *error_msg = "EHSS samples must contain non-empty 'orientations_quat' array";
        }
        return false;
    }

    out.orientations_quat.clear();
    out.areas_by_gas_m2.clear();
    out.orientations_quat.reserve(orientations.size());

    for (const auto& q : orientations) {
        if (!q.isArray() || q.size() != 4) {
            if (error_msg) {
                *error_msg = "Each orientation must be an array of 4 numbers [w,x,y,z]";
            }
            return false;
        }
        std::array<double, 4> quat{};
        for (Json::ArrayIndex i = 0; i < 4; ++i) {
            if (!q[i].isNumeric()) {
                if (error_msg) {
                    *error_msg = "Orientation quaternion entries must be numeric";
                }
                return false;
            }
            quat[i] = q[i].asDouble();
        }
        out.orientations_quat.push_back(quat);
    }

    const Json::Value& areas = root["areas_by_gas_m2"];
    if (!areas.isObject()) {
        if (error_msg) {
            *error_msg = "EHSS samples must contain 'areas_by_gas_m2' object";
        }
        return false;
    }

    const auto gas_names = areas.getMemberNames();
    for (const auto& gas : gas_names) {
        const Json::Value& values = areas[gas];
        if (!values.isArray()) {
            if (error_msg) {
                *error_msg = "areas_by_gas_m2 entries must be arrays";
            }
            return false;
        }
        if (values.size() != orientations.size()) {
            if (error_msg) {
                *error_msg = "areas_by_gas_m2 arrays must match orientations_quat size";
            }
            return false;
        }
        std::vector<double> area_samples;
        area_samples.reserve(values.size());
        for (const auto& v : values) {
            if (!v.isNumeric()) {
                if (error_msg) {
                    *error_msg = "areas_by_gas_m2 entries must be numeric";
                }
                return false;
            }
            area_samples.push_back(v.asDouble());
        }
        out.areas_by_gas_m2.emplace(gas, std::move(area_samples));
    }

    if (root.isMember("n_orientations") && root["n_orientations"].isInt()) {
        const int n_orientations = root["n_orientations"].asInt();
        if (n_orientations != static_cast<int>(out.orientations_quat.size())) {
            if (error_msg) {
                *error_msg = "n_orientations does not match orientations_quat size";
            }
            return false;
        }
    }

    return true;
}

} // namespace ICARION::physics
