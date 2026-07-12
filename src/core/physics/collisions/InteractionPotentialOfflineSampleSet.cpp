// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "InteractionPotentialOfflineSampleSet.h"
#include "core/io/hdf5Utils.h"

#include <H5Cpp.h>
#include <cmath>

namespace ICARION::physics {

namespace {

bool validate_physical_values(const InteractionPotentialOfflineSampleSet& out, std::string& error_msg) {
    if (out.logv_bins.empty()) {
        error_msg = "logv_bins must be non-empty";
        return false;
    }
    for (size_t i = 0; i < out.logv_bins.size(); ++i) {
        const double v = out.logv_bins[i];
        if (!std::isfinite(v)) {
            error_msg = "logv_bins entries must be finite";
            return false;
        }
        if (i > 0 && !(out.logv_bins[i - 1] < v)) {
            error_msg = "logv_bins entries must be strictly increasing";
            return false;
        }
    }

    for (double sigma : out.sigma_mt_m2) {
        if (!std::isfinite(sigma)) {
            error_msg = "sigma_mt_m2 entries must be finite";
            return false;
        }
        if (sigma <= 0.0) {
            error_msg = "sigma_mt_m2 entries must be > 0";
            return false;
        }
    }

    for (double sigma : out.sigma_event_m2) {
        if (!std::isfinite(sigma)) {
            error_msg = "sigma_event_m2 entries must be finite";
            return false;
        }
        if (sigma <= 0.0) {
            error_msg = "sigma_event_m2 entries must be > 0";
            return false;
        }
    }

    for (double bmax : out.b_max_m) {
        if (!std::isfinite(bmax)) {
            error_msg = "b_max_m entries must be finite";
            return false;
        }
        if (bmax <= 0.0) {
            error_msg = "b_max_m entries must be > 0";
            return false;
        }
    }

    auto validate_finite = [&](const std::vector<double>& values, const char* name) -> bool {
        for (double x : values) {
            if (!std::isfinite(x)) {
                error_msg = std::string(name) + " entries must be finite";
                return false;
            }
        }
        return true;
    };

    if (!validate_finite(out.q1_m2, "q1_m2") ||
        !validate_finite(out.q2_m2, "q2_m2") ||
        !validate_finite(out.q3_m2, "q3_m2") ||
        !validate_finite(out.q12_m2, "q12_m2") ||
        !validate_finite(out.q13_m2, "q13_m2") ||
        !validate_finite(out.cdf_values, "cdf_values") ||
        !validate_finite(out.dp_samples, "dp_samples") ||
        !validate_finite(out.dp_stats, "dp_stats")) {
        return false;
    }

    for (size_t i = 2; i < out.dp_stats.size(); i += 4) {
        if (out.dp_stats[i] < 0.0 || out.dp_stats[i + 1] < 0.0) {
            error_msg = "dp_stats variances must be >= 0";
            return false;
        }
    }

    if (!out.valid()) {
        error_msg = "InteractionPotential sample content failed validation";
        return false;
    }
    return true;
}

bool read_double_cells(
    H5::H5File& file,
    const char* name,
    size_t n_orient,
    size_t n_bins,
    std::vector<double>& target,
    std::string& error_msg
) {
    H5::DataSet dset = file.openDataSet(name);
    H5::DataSpace space = dset.getSpace();
    if (space.getSimpleExtentNdims() != 2) {
        error_msg = std::string(name) + " must be 2D (N x K)";
        return false;
    }
    hsize_t dims[2];
    space.getSimpleExtentDims(dims);
    if (dims[0] != n_orient || dims[1] != n_bins) {
        error_msg = std::string(name) + " dimensions mismatch";
        return false;
    }
    target.assign(n_orient * n_bins, 0.0);
    dset.read(target.data(), H5::PredType::NATIVE_DOUBLE);
    return true;
}

bool read_long_long_cells(
    H5::H5File& file,
    const char* name,
    size_t n_orient,
    size_t n_bins,
    std::vector<long long>& target,
    std::string& error_msg
) {
    H5::DataSet dset = file.openDataSet(name);
    H5::DataSpace space = dset.getSpace();
    if (space.getSimpleExtentNdims() != 2) {
        error_msg = std::string(name) + " must be 2D (N x K)";
        return false;
    }
    hsize_t dims[2];
    space.getSimpleExtentDims(dims);
    if (dims[0] != n_orient || dims[1] != n_bins) {
        error_msg = std::string(name) + " dimensions mismatch";
        return false;
    }
    target.assign(n_orient * n_bins, 0);
    dset.read(target.data(), H5::PredType::NATIVE_LLONG);
    return true;
}

bool load_hdf5_samples(const std::filesystem::path& path, InteractionPotentialOfflineSampleSet& out, std::string& error_msg) {
    try {
        H5::H5File file(path.string(), H5F_ACC_RDONLY);

        if (!file.attrExists("format") ||
            ICARION::io::read_hdf5_attr_string(file, "format") != INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_FORMAT) {
            error_msg = std::string("'format' attribute is required and must be '")
                + INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_FORMAT + "'";
            return false;
        }
        if (!file.attrExists("units") ||
            ICARION::io::read_hdf5_attr_string(file, "units") != INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_UNITS) {
            error_msg = std::string("'units' attribute is required and must be '")
                + INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_UNITS + "'";
            return false;
        }
        out.gas = ICARION::io::read_hdf5_attr_string(file, "gas");
        if (!file.attrExists("gas") || out.gas.empty()) {
            error_msg = "'gas' attribute is required and must be non-empty";
            return false;
        }

        long long attr_value = 0;
        if (!ICARION::io::read_hdf5_attr_long_long(file, "version", attr_value)) {
            error_msg = "'version' attribute is required";
            return false;
        }
        if (attr_value != INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_VERSION) {
            error_msg = "'version' attribute is not supported";
            return false;
        }

        if (!file.nameExists("logv_bins") || !file.nameExists("orientations_quat") ||
            !file.nameExists("sigma_mt_m2") || !file.nameExists("sigma_event_m2") || !file.nameExists("b_max_m")) {
            error_msg = "HDF5 samples must contain logv_bins, orientations_quat, sigma_mt_m2, sigma_event_m2, b_max_m";
            return false;
        }

        H5::DataSet logv_dset = file.openDataSet("logv_bins");
        H5::DataSpace logv_space = logv_dset.getSpace();
        int logv_rank = logv_space.getSimpleExtentNdims();
        if (logv_rank != 1) {
            error_msg = "logv_bins must be 1D";
            return false;
        }
        hsize_t logv_dims[1];
        logv_space.getSimpleExtentDims(logv_dims);
        out.n_bins = static_cast<size_t>(logv_dims[0]);

        H5::DataSet sigma_dset = file.openDataSet("sigma_mt_m2");
        H5::DataSpace sigma_space = sigma_dset.getSpace();
        int sigma_rank = sigma_space.getSimpleExtentNdims();
        if (sigma_rank != 2) {
            error_msg = "sigma_mt_m2 must be 2D (N x K)";
            return false;
        }
        hsize_t sigma_dims[2];
        sigma_space.getSimpleExtentDims(sigma_dims);
        out.n_orient = static_cast<size_t>(sigma_dims[0]);
        if (out.n_orient == 0) {
            error_msg = "sigma_mt_m2 must have at least one orientation row";
            return false;
        }
        if (static_cast<size_t>(sigma_dims[1]) != out.n_bins) {
            error_msg = "sigma_mt_m2 second dimension must match logv_bins";
            return false;
        }

        if (ICARION::io::read_hdf5_attr_long_long(file, "n_orientations", attr_value)
            && attr_value != static_cast<long long>(out.n_orient)) {
            error_msg = "'n_orientations' attribute does not match payload dimensions";
            return false;
        }
        if (ICARION::io::read_hdf5_attr_long_long(file, "v_bins", attr_value)
            && attr_value != static_cast<long long>(out.n_bins)) {
            error_msg = "'v_bins' attribute does not match payload dimensions";
            return false;
        }

        out.logv_bins.assign(out.n_bins, 0.0);
        logv_dset.read(out.logv_bins.data(), H5::PredType::NATIVE_DOUBLE);

        out.sigma_mt_m2.assign(out.n_orient * out.n_bins, 0.0);
        sigma_dset.read(out.sigma_mt_m2.data(), H5::PredType::NATIVE_DOUBLE);

        {
            H5::DataSet orient_dset = file.openDataSet("orientations_quat");
            H5::DataSpace orient_space = orient_dset.getSpace();
            if (orient_space.getSimpleExtentNdims() != 2) {
                error_msg = "orientations_quat must be 2D (N x 4)";
                return false;
            }
            hsize_t orient_dims[2];
            orient_space.getSimpleExtentDims(orient_dims);
            if (orient_dims[0] != out.n_orient || orient_dims[1] != 4) {
                error_msg = "orientations_quat dimensions mismatch";
                return false;
            }
        }

        if (!read_double_cells(file, "sigma_event_m2", out.n_orient, out.n_bins, out.sigma_event_m2, error_msg)) {
            return false;
        }

        if (file.nameExists("q1_m2")) {
            if (!read_double_cells(file, "q1_m2", out.n_orient, out.n_bins, out.q1_m2, error_msg)) return false;
        }
        if (file.nameExists("q2_m2")) {
            if (!read_double_cells(file, "q2_m2", out.n_orient, out.n_bins, out.q2_m2, error_msg)) return false;
        }
        if (file.nameExists("q3_m2")) {
            if (!read_double_cells(file, "q3_m2", out.n_orient, out.n_bins, out.q3_m2, error_msg)) return false;
        }
        if (file.nameExists("q12_m2")) {
            if (!read_double_cells(file, "q12_m2", out.n_orient, out.n_bins, out.q12_m2, error_msg)) return false;
        }
        if (file.nameExists("q13_m2")) {
            if (!read_double_cells(file, "q13_m2", out.n_orient, out.n_bins, out.q13_m2, error_msg)) return false;
        }

        if (!read_double_cells(file, "b_max_m", out.n_orient, out.n_bins, out.b_max_m, error_msg)) {
            return false;
        }

        if (file.nameExists("cdf_offsets")) {
            if (!read_long_long_cells(file, "cdf_offsets", out.n_orient, out.n_bins, out.cdf_offsets, error_msg)) {
                return false;
            }
        }
        if (file.nameExists("cdf_counts")) {
            if (!read_long_long_cells(file, "cdf_counts", out.n_orient, out.n_bins, out.cdf_counts, error_msg)) {
                return false;
            }
        }
        if (file.nameExists("cdf_values")) {
            H5::DataSet vals = file.openDataSet("cdf_values");
            H5::DataSpace vals_space = vals.getSpace();
            if (vals_space.getSimpleExtentNdims() != 1) {
                error_msg = "cdf_values must be 1D";
                return false;
            }
            hsize_t dims[1];
            vals_space.getSimpleExtentDims(dims);
            out.cdf_values.assign(static_cast<size_t>(dims[0]), 0.0);
            vals.read(out.cdf_values.data(), H5::PredType::NATIVE_DOUBLE);
        }
        if (file.nameExists("dp_samples")) {
            H5::DataSet dps = file.openDataSet("dp_samples");
            H5::DataSpace dps_space = dps.getSpace();
            int rank = dps_space.getSimpleExtentNdims();
            if (rank != 2) {
                error_msg = "dp_samples must be 2D (N x 3)";
                return false;
            }
            hsize_t dims[2];
            dps_space.getSimpleExtentDims(dims);
            if (dims[1] != 3) {
                error_msg = "dp_samples second dimension must be 3";
                return false;
            }
            out.dp_samples.assign(static_cast<size_t>(dims[0] * 3), 0.0);
            dps.read(out.dp_samples.data(), H5::PredType::NATIVE_DOUBLE);
        }
        if (file.nameExists("dp_stats")) {
            H5::DataSet stats = file.openDataSet("dp_stats");
            H5::DataSpace stats_space = stats.getSpace();
            int rank = stats_space.getSimpleExtentNdims();
            if (rank != 3) {
                error_msg = "dp_stats must be 3D (N x K x 4)";
                return false;
            }
            hsize_t dims[3];
            stats_space.getSimpleExtentDims(dims);
            if (dims[0] != out.n_orient || dims[1] != out.n_bins || dims[2] != 4) {
                error_msg = "dp_stats dimensions mismatch";
                return false;
            }
            out.dp_stats.assign(out.n_orient * out.n_bins * 4, 0.0);
            stats.read(out.dp_stats.data(), H5::PredType::NATIVE_DOUBLE);
        }

        return validate_physical_values(out, error_msg);
    } catch (const H5::Exception& e) {
        error_msg = std::string("Failed to read InteractionPotential HDF5 samples: ") + e.getCDetailMsg();
        return false;
    }
}

}  // namespace

bool load_interaction_potential_offline_sample_set_file(
    const std::filesystem::path& path,
    InteractionPotentialOfflineSampleSet& out,
    std::string* error_msg
) {
    std::string err;
    const auto ext = path.extension().string();
    bool ok = false;
    if (ext == ".h5" || ext == ".hdf5") {
        ok = load_hdf5_samples(path, out, err);
    } else {
        err = "Only HDF5 InteractionPotential samples supported";
        ok = false;
    }
    if (!ok && error_msg) {
        *error_msg = err;
    }
    return ok;
}

} // namespace ICARION::physics
