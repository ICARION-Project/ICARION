// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "EHSSOfflineSampleSet.h"
#include "core/io/hdf5Utils.h"
#include <json/json.h>
#include <fstream>
#include <cmath>
#include <H5Cpp.h>

namespace ICARION::physics {

namespace {

constexpr int EXPECTED_EHSS_OFFLINE_VERSION = 1;

bool validate_optional_string_field(
    const Json::Value& root,
    const char* key,
    const char* expected,
    std::string& error_msg
) {
    if (!root.isMember(key)) {
        return true;
    }
    if (!root[key].isString()) {
        error_msg = std::string("'") + key + "' must be a string when present";
        return false;
    }
    if (expected != nullptr && root[key].asString() != expected) {
        error_msg = std::string("'") + key + "' must be '" + expected + "'";
        return false;
    }
    return true;
}

bool validate_optional_int_field(
    const Json::Value& root,
    const char* key,
    int expected,
    std::string& error_msg
) {
    if (!root.isMember(key)) {
        return true;
    }
    if (!root[key].isInt()) {
        error_msg = std::string("'") + key + "' must be an integer when present";
        return false;
    }
    if (root[key].asInt() != expected) {
        error_msg = std::string("'") + key + "' does not match payload dimensions";
        return false;
    }
    return true;
}

bool validate_json_metadata(const Json::Value& root, const EHSSOfflineSampleSet& out, std::string& error_msg) {
    if (root.isMember("gas")) {
        if (!root["gas"].isString() || root["gas"].asString().empty()) {
            error_msg = "'gas' must be a non-empty string when present";
            return false;
        }
    }
    if (!validate_optional_string_field(root, "format", EHSS_OFFLINE_SAMPLE_SET_FORMAT, error_msg)) {
        return false;
    }
    if (!validate_optional_string_field(root, "units", EHSS_OFFLINE_SAMPLE_SET_UNITS, error_msg)) {
        return false;
    }
    if (!validate_optional_int_field(root, "version", EXPECTED_EHSS_OFFLINE_VERSION, error_msg)) {
        return false;
    }
    if (!validate_optional_int_field(root, "n_orientations", static_cast<int>(out.sigma_eff_m2.size()), error_msg)) {
        return false;
    }
    if (!validate_optional_int_field(root, "n_mu_samples", static_cast<int>(out.mu_stride), error_msg)) {
        return false;
    }
    return true;
}

bool validate_physical_values(const EHSSOfflineSampleSet& out, std::string& error_msg) {
    for (double sigma : out.sigma_eff_m2) {
        if (!std::isfinite(sigma)) {
            error_msg = "sigma_eff_m2 entries must be finite";
            return false;
        }
        if (sigma <= 0.0) {
            error_msg = "sigma_eff_m2 entries must be > 0";
            return false;
        }
    }

    for (double mu : out.mu_samples_flat) {
        if (!std::isfinite(mu)) {
            error_msg = "mu_samples entries must be finite";
            return false;
        }
        if (mu < -1.0 || mu > 1.0) {
            error_msg = "mu_samples entries must lie in [-1, 1]";
            return false;
        }
    }

    if (!out.valid()) {
        error_msg = "EHSS offline samples are internally inconsistent";
        return false;
    }

    return true;
}

bool load_json_file(const std::filesystem::path& path, Json::Value& root, std::string& error_msg) {
    std::ifstream file(path);
    if (!file.is_open()) {
        error_msg = "Failed to open EHSS offline samples file: " + path.string();
        return false;
    }

    Json::CharReaderBuilder builder;
    std::string errs;
    if (!Json::parseFromStream(builder, file, &root, &errs)) {
        error_msg = "Failed to parse EHSS offline samples file: " + errs;
        return false;
    }

    return true;
}

bool load_json_samples(const std::filesystem::path& path, EHSSOfflineSampleSet& out, std::string& error_msg) {
    Json::Value root;
    if (!load_json_file(path, root, error_msg)) {
        return false;
    }

    if (!root.isObject()) {
        error_msg = "EHSS offline samples root must be a JSON object";
        return false;
    }

    out.gas.clear();
    if (root.isMember("gas") && root["gas"].isString()) {
        out.gas = root["gas"].asString();
    }

    const Json::Value& sigma = root["sigma_eff_m2"];
    if (!sigma.isArray() || sigma.empty()) {
        error_msg = "EHSS offline samples must contain non-empty 'sigma_eff_m2' array";
        return false;
    }

    const Json::Value& mu = root["mu_samples"];
    if (!mu.isArray() || mu.empty()) {
        error_msg = "EHSS offline samples must contain non-empty 'mu_samples' array";
        return false;
    }

    if (mu.size() != sigma.size()) {
        error_msg = "mu_samples size must match sigma_eff_m2 size";
        return false;
    }

    out.sigma_eff_m2.clear();
    out.mu_samples_flat.clear();
    out.mu_stride = 0;

    out.sigma_eff_m2.reserve(sigma.size());
    for (const auto& v : sigma) {
        if (!v.isNumeric()) {
            error_msg = "sigma_eff_m2 entries must be numeric";
            return false;
        }
        out.sigma_eff_m2.push_back(v.asDouble());
    }

    const auto first_row = mu[0];
    if (!first_row.isArray() || first_row.empty()) {
        error_msg = "mu_samples rows must be non-empty arrays";
        return false;
    }
    out.mu_stride = first_row.size();
    out.mu_samples_flat.reserve(static_cast<size_t>(mu.size()) * out.mu_stride);

    for (const auto& row : mu) {
        if (!row.isArray() || row.size() != out.mu_stride) {
            error_msg = "mu_samples rows must all have the same length";
            return false;
        }
        for (const auto& entry : row) {
            if (!entry.isNumeric()) {
                error_msg = "mu_samples entries must be numeric";
                return false;
            }
            out.mu_samples_flat.push_back(entry.asDouble());
        }
    }

        return validate_json_metadata(root, out, error_msg)
                && validate_physical_values(out, error_msg);
}

bool load_hdf5_samples(const std::filesystem::path& path, EHSSOfflineSampleSet& out, std::string& error_msg) {
    try {
        H5::H5File file(path.string(), H5F_ACC_RDONLY);

        out.gas = ICARION::io::read_hdf5_attr_string(file, "gas");
        if (file.attrExists("gas") && out.gas.empty()) {
            error_msg = "'gas' attribute must be non-empty when present";
            return false;
        }
        if (file.attrExists("format") && ICARION::io::read_hdf5_attr_string(file, "format") != EHSS_OFFLINE_SAMPLE_SET_FORMAT) {
            error_msg = std::string("'format' attribute must be '") + EHSS_OFFLINE_SAMPLE_SET_FORMAT + "'";
            return false;
        }
        if (file.attrExists("units") && ICARION::io::read_hdf5_attr_string(file, "units") != EHSS_OFFLINE_SAMPLE_SET_UNITS) {
            error_msg = std::string("'units' attribute must be '") + EHSS_OFFLINE_SAMPLE_SET_UNITS + "'";
            return false;
        }

        long long attr_value = 0;
        if (ICARION::io::read_hdf5_attr_long_long(file, "version", attr_value) && attr_value != EXPECTED_EHSS_OFFLINE_VERSION) {
            error_msg = "'version' attribute is not supported";
            return false;
        }

        if (!file.nameExists("sigma_eff_m2") || !file.nameExists("mu_samples")) {
            error_msg = "HDF5 samples must contain datasets 'sigma_eff_m2' and 'mu_samples'";
            return false;
        }

        H5::DataSet sigma_dset = file.openDataSet("sigma_eff_m2");
        H5::DataSpace sigma_space = sigma_dset.getSpace();
        int sigma_rank = sigma_space.getSimpleExtentNdims();
        if (sigma_rank != 1) {
            error_msg = "sigma_eff_m2 must be 1D";
            return false;
        }
        hsize_t sigma_dims[1];
        sigma_space.getSimpleExtentDims(sigma_dims);
        const size_t n_orient = static_cast<size_t>(sigma_dims[0]);

        H5::DataSet mu_dset = file.openDataSet("mu_samples");
        H5::DataSpace mu_space = mu_dset.getSpace();
        int mu_rank = mu_space.getSimpleExtentNdims();
        if (mu_rank != 2) {
            error_msg = "mu_samples must be 2D (N x M)";
            return false;
        }
        hsize_t mu_dims[2];
        mu_space.getSimpleExtentDims(mu_dims);
        const size_t mu_orient = static_cast<size_t>(mu_dims[0]);
        const size_t mu_stride = static_cast<size_t>(mu_dims[1]);

        if (mu_orient != n_orient) {
            error_msg = "mu_samples first dimension must match sigma_eff_m2 size";
            return false;
        }

        out.sigma_eff_m2.assign(n_orient, 0.0);
        out.mu_samples_flat.assign(mu_orient * mu_stride, 0.0);
        out.mu_stride = mu_stride;

        if (ICARION::io::read_hdf5_attr_long_long(file, "n_orientations", attr_value)
            && attr_value != static_cast<long long>(n_orient)) {
            error_msg = "'n_orientations' attribute does not match payload dimensions";
            return false;
        }
        if (ICARION::io::read_hdf5_attr_long_long(file, "n_mu_samples", attr_value)
            && attr_value != static_cast<long long>(mu_stride)) {
            error_msg = "'n_mu_samples' attribute does not match payload dimensions";
            return false;
        }

        sigma_dset.read(out.sigma_eff_m2.data(), H5::PredType::NATIVE_DOUBLE);
        mu_dset.read(out.mu_samples_flat.data(), H5::PredType::NATIVE_DOUBLE);
        return validate_physical_values(out, error_msg);
    } catch (const H5::Exception& e) {
        error_msg = std::string("Failed to read HDF5 samples: ") + e.getCDetailMsg();
        return false;
    }
}

}  // namespace

bool load_ehss_offline_sample_set_file(
    const std::filesystem::path& path,
    EHSSOfflineSampleSet& out,
    std::string* error_msg
) {
    std::string err;
    const auto ext = path.extension().string();
    bool ok = false;
    if (ext == ".h5" || ext == ".hdf5") {
        ok = load_hdf5_samples(path, out, err);
    } else {
        ok = load_json_samples(path, out, err);
    }

    if (!ok && error_msg) {
        *error_msg = err;
    }
    return ok;
}

} // namespace ICARION::physics
