// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "core/config/types/EnvironmentConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/physics/collisions/InteractionPotentialCollisionHandler.h"
#include "core/physics/collisions/InteractionPotentialOfflineSampleSet.h"
#include "core/types/IonEnsemble.h"
#include "utils/constants.h"

#include <H5Cpp.h>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <vector>

using namespace ICARION;

namespace {

std::filesystem::path make_temp_h5_path(const std::string& tag) {
    const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("icarion_ipm_" + tag + "_" + std::to_string(stamp) + ".h5");
    std::error_code ec;
    std::filesystem::remove(path, ec);
    return path;
}

void write_common_ipm_metadata(H5::H5File& file, const std::string& gas_tag, size_t n_orient = 1) {
    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);

    const long long version = physics::INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_VERSION;
    H5::Attribute version_attr = file.createAttribute("version", H5::PredType::NATIVE_LLONG, scalar_space);
    version_attr.write(H5::PredType::NATIVE_LLONG, &version);

    const char* format = physics::INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_FORMAT;
    H5::Attribute format_attr = file.createAttribute("format", str_type, scalar_space);
    format_attr.write(str_type, &format);

    const char* units = physics::INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_UNITS;
    H5::Attribute units_attr = file.createAttribute("units", str_type, scalar_space);
    units_attr.write(str_type, &units);

    const char* gas_c = gas_tag.c_str();
    H5::Attribute gas_attr = file.createAttribute("gas", str_type, scalar_space);
    gas_attr.write(str_type, &gas_c);

    const hsize_t dims[2] = {static_cast<hsize_t>(n_orient), 4};
    H5::DataSpace space(2, dims);
    H5::DataSet orientations = file.createDataSet("orientations_quat", H5::PredType::NATIVE_DOUBLE, space);
    std::vector<double> values(n_orient * 4, 0.0);
    for (size_t i = 0; i < n_orient; ++i) {
        values[i * 4] = 1.0;
    }
    orientations.write(values.data(), H5::PredType::NATIVE_DOUBLE);
}

void write_samples_file(const std::filesystem::path& path, const std::string& gas_tag) {
    H5::H5File file(path.string(), H5F_ACC_TRUNC);

    write_common_ipm_metadata(file, gas_tag);

    {
        const hsize_t dims[1] = {1};
        H5::DataSpace space(1, dims);
        H5::DataSet ds = file.createDataSet("logv_bins", H5::PredType::NATIVE_DOUBLE, space);
        const double values[1] = {std::log(200.0)};
        ds.write(values, H5::PredType::NATIVE_DOUBLE);
    }
    {
        const hsize_t dims[2] = {1, 1};
        H5::DataSpace space(2, dims);
        H5::DataSet sigma = file.createDataSet("sigma_mt_m2", H5::PredType::NATIVE_DOUBLE, space);
        const double sigma_values[1] = {1.0e-18};
        sigma.write(sigma_values, H5::PredType::NATIVE_DOUBLE);

        H5::DataSet sigma_event = file.createDataSet("sigma_event_m2", H5::PredType::NATIVE_DOUBLE, space);
        sigma_event.write(sigma_values, H5::PredType::NATIVE_DOUBLE);

        H5::DataSet bmax = file.createDataSet("b_max_m", H5::PredType::NATIVE_DOUBLE, space);
        const double bmax_values[1] = {1.0e-9};
        bmax.write(bmax_values, H5::PredType::NATIVE_DOUBLE);
    }
    {
        const hsize_t cell_dims[2] = {1, 1};
        H5::DataSpace cell_space(2, cell_dims);
        H5::DataSet offsets = file.createDataSet("cdf_offsets", H5::PredType::NATIVE_LLONG, cell_space);
        const long long offsets_values[1] = {0};
        offsets.write(offsets_values, H5::PredType::NATIVE_LLONG);

        H5::DataSet counts = file.createDataSet("cdf_counts", H5::PredType::NATIVE_LLONG, cell_space);
        const long long counts_values[1] = {1};
        counts.write(counts_values, H5::PredType::NATIVE_LLONG);

        const hsize_t sample_dims[1] = {1};
        H5::DataSpace space(1, sample_dims);
        H5::DataSet cdf = file.createDataSet("cdf_values", H5::PredType::NATIVE_DOUBLE, space);
        const double cdf_values[1] = {1.0};
        cdf.write(cdf_values, H5::PredType::NATIVE_DOUBLE);
    }
    {
        const hsize_t dims[2] = {1, 3};
        H5::DataSpace space(2, dims);
        H5::DataSet dp = file.createDataSet("dp_samples", H5::PredType::NATIVE_DOUBLE, space);
        const double dp_values[3] = {0.0, 0.0, 1.0e-22};
        dp.write(dp_values, H5::PredType::NATIVE_DOUBLE);
    }
}

void write_dp_stats_only_samples_file(
    const std::filesystem::path& path,
    const std::string& gas_tag,
    const std::vector<double>& logv_bins,
    const std::vector<double>& sigma_mt_m2,
    const std::vector<double>& dp_stats,
    size_t n_orient = 1
) {
    H5::H5File file(path.string(), H5F_ACC_TRUNC);

    write_common_ipm_metadata(file, gas_tag, n_orient);

    {
        const hsize_t dims[1] = {static_cast<hsize_t>(logv_bins.size())};
        H5::DataSpace space(1, dims);
        H5::DataSet ds = file.createDataSet("logv_bins", H5::PredType::NATIVE_DOUBLE, space);
        ds.write(logv_bins.data(), H5::PredType::NATIVE_DOUBLE);
    }
    {
        const hsize_t dims[2] = {
            static_cast<hsize_t>(n_orient),
            static_cast<hsize_t>(logv_bins.size())
        };
        H5::DataSpace space(2, dims);

        H5::DataSet sigma = file.createDataSet("sigma_mt_m2", H5::PredType::NATIVE_DOUBLE, space);
        sigma.write(sigma_mt_m2.data(), H5::PredType::NATIVE_DOUBLE);

        H5::DataSet sigma_event = file.createDataSet("sigma_event_m2", H5::PredType::NATIVE_DOUBLE, space);
        sigma_event.write(sigma_mt_m2.data(), H5::PredType::NATIVE_DOUBLE);

        H5::DataSet bmax = file.createDataSet("b_max_m", H5::PredType::NATIVE_DOUBLE, space);
        std::vector<double> bmax_values(n_orient * logv_bins.size(), 1.0e-9);
        bmax.write(bmax_values.data(), H5::PredType::NATIVE_DOUBLE);
    }
    {
        const hsize_t dims[3] = {
            static_cast<hsize_t>(n_orient),
            static_cast<hsize_t>(logv_bins.size()),
            4
        };
        H5::DataSpace space(3, dims);
        H5::DataSet stats = file.createDataSet("dp_stats", H5::PredType::NATIVE_DOUBLE, space);
        stats.write(dp_stats.data(), H5::PredType::NATIVE_DOUBLE);
    }
}

config::SpeciesDatabase make_species_db(const std::filesystem::path& sample_path) {
    config::SpeciesDatabase species_db;
    config::SpeciesProperties sp;
    sp.id = "X+";
    sp.mass_amu = 28.0;
    sp.charge = 1;
    sp.ipm_samples_file = sample_path.string();
    species_db.species.emplace(sp.id, sp);
    return species_db;
}

IonState make_test_ion(double vx_m_s) {
    IonState ion;
    ion.species_id = "X+";
    ion.mass_kg = 28.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1.0e-18;
    ion.vel = Vec3{vx_m_s, 0.0, 0.0};
    return ion;
}

bool near_rate(double actual, double expected) {
    return std::abs(actual - expected) <= std::max(1.0e3, 1.0e-12 * std::abs(expected));
}

}  // namespace

TEST_CASE("InteractionPotentialOfflineSampleSet validates CDF offset/count bounds", "[collision][ipm][validation]") {
    physics::InteractionPotentialOfflineSampleSet s;
    s.n_orient = 1;
    s.n_bins = 1;
    s.logv_bins = {0.0};
    s.sigma_mt_m2 = {1.0e-18};
    s.sigma_event_m2 = {1.0e-18};
    s.b_max_m = {1.0e-9};
    s.cdf_offsets = {3};
    s.cdf_counts = {1};
    s.cdf_values = {0.1, 0.5};
    s.dp_samples = {0.0, 0.0, 1.0e-22};

    REQUIRE_FALSE(s.valid());
}

TEST_CASE("InteractionPotentialOfflineSampleSet validates cdf_values to dp_samples consistency", "[collision][ipm][validation]") {
    physics::InteractionPotentialOfflineSampleSet s;
    s.n_orient = 1;
    s.n_bins = 1;
    s.logv_bins = {0.0};
    s.sigma_mt_m2 = {1.0e-18};
    s.sigma_event_m2 = {1.0e-18};
    s.b_max_m = {1.0e-9};
    s.cdf_offsets = {0};
    s.cdf_counts = {2};
    s.cdf_values = {0.2, 1.0};
    s.dp_samples = {0.0, 0.0, 1.0e-22};  // only one sample

    REQUIRE_FALSE(s.valid());
}

TEST_CASE("InteractionPotentialOfflineSampleSet validates CDF monotonicity", "[collision][ipm][validation]") {
    physics::InteractionPotentialOfflineSampleSet s;
    s.n_orient = 1;
    s.n_bins = 1;
    s.logv_bins = {0.0};
    s.sigma_mt_m2 = {1.0e-18};
    s.sigma_event_m2 = {1.0e-18};
    s.b_max_m = {1.0e-9};
    s.cdf_offsets = {0};
    s.cdf_counts = {2};
    s.cdf_values = {0.9, 0.3};  // non-monotonic
    s.dp_samples = {0.0, 0.0, 1.0e-22, 0.0, 0.0, 2.0e-22};

    REQUIRE_FALSE(s.valid());
}

TEST_CASE("InteractionPotentialOfflineSampleSet loader rejects inconsistent HDF5 metadata", "[collision][ipm][loader]") {
    const auto sample_path = make_temp_h5_path("meta_mismatch");
    write_samples_file(sample_path, "He");

    {
        H5::H5File file(sample_path.string(), H5F_ACC_RDWR);
        H5::DataSpace scalar_space(H5S_SCALAR);
        const long long bad_orient = 2;
        H5::Attribute attr = file.createAttribute("n_orientations", H5::PredType::NATIVE_LLONG, scalar_space);
        attr.write(H5::PredType::NATIVE_LLONG, &bad_orient);
    }

    physics::InteractionPotentialOfflineSampleSet loaded;
    std::string error;
    REQUIRE_FALSE(physics::load_interaction_potential_offline_sample_set_file(sample_path, loaded, &error));
    REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("n_orientations"));

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialOfflineSampleSet loader rejects non-physical sigma", "[collision][ipm][loader]") {
    const auto sample_path = make_temp_h5_path("sigma_nonphysical");

    {
        H5::H5File file(sample_path.string(), H5F_ACC_TRUNC);
        write_common_ipm_metadata(file, "He");

        {
            const hsize_t dims[1] = {1};
            H5::DataSpace space(1, dims);
            H5::DataSet ds = file.createDataSet("logv_bins", H5::PredType::NATIVE_DOUBLE, space);
            const double values[1] = {std::log(200.0)};
            ds.write(values, H5::PredType::NATIVE_DOUBLE);
        }
        {
            const hsize_t dims[2] = {1, 1};
            H5::DataSpace space(2, dims);
            H5::DataSet sigma = file.createDataSet("sigma_mt_m2", H5::PredType::NATIVE_DOUBLE, space);
            const double sigma_values[1] = {0.0};
            sigma.write(sigma_values, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet sigma_event = file.createDataSet("sigma_event_m2", H5::PredType::NATIVE_DOUBLE, space);
            const double sigma_event_values[1] = {1.0e-18};
            sigma_event.write(sigma_event_values, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet bmax = file.createDataSet("b_max_m", H5::PredType::NATIVE_DOUBLE, space);
            const double bmax_values[1] = {1.0e-9};
            bmax.write(bmax_values, H5::PredType::NATIVE_DOUBLE);
        }
    }

    physics::InteractionPotentialOfflineSampleSet loaded;
    std::string error;
    REQUIRE_FALSE(physics::load_interaction_potential_offline_sample_set_file(sample_path, loaded, &error));
    REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("sigma_mt_m2"));

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialOfflineSampleSet loader rejects unsupported HDF5 versions", "[collision][ipm][loader]") {
    const auto sample_path = make_temp_h5_path("old_version");
    write_samples_file(sample_path, "He");

    {
        H5::H5File file(sample_path.string(), H5F_ACC_RDWR);
        H5::Attribute attr = file.openAttribute("version");
        const long long unsupported_version = 2;
        attr.write(H5::PredType::NATIVE_LLONG, &unsupported_version);
    }

    physics::InteractionPotentialOfflineSampleSet loaded;
    std::string error;
    REQUIRE_FALSE(physics::load_interaction_potential_offline_sample_set_file(sample_path, loaded, &error));
    REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("version"));

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialOfflineSampleSet loader requires core HDF5 attributes", "[collision][ipm][loader]") {
    const std::vector<std::string> required_attrs = {"format", "units", "gas"};

    for (const auto& attr_name : required_attrs) {
        const auto sample_path = make_temp_h5_path("missing_" + attr_name);
        write_samples_file(sample_path, "He");

        {
            H5::H5File file(sample_path.string(), H5F_ACC_RDWR);
            H5Adelete(file.getId(), attr_name.c_str());
        }

        physics::InteractionPotentialOfflineSampleSet loaded;
        std::string error;
        REQUIRE_FALSE(physics::load_interaction_potential_offline_sample_set_file(sample_path, loaded, &error));
        REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring(attr_name));

        std::error_code ec;
        std::filesystem::remove(sample_path, ec);
    }
}

TEST_CASE("InteractionPotentialOfflineSampleSet loader validates HDF5 cell-table shapes", "[collision][ipm][loader]") {
    SECTION("b_max_m must be a 2D cell table") {
        const auto sample_path = make_temp_h5_path("bmax_rank");
        {
            H5::H5File file(sample_path.string(), H5F_ACC_TRUNC);
            write_common_ipm_metadata(file, "He");

            const hsize_t one_dim[1] = {1};
            H5::DataSpace vec_space(1, one_dim);
            const double logv[1] = {std::log(200.0)};
            file.createDataSet("logv_bins", H5::PredType::NATIVE_DOUBLE, vec_space).write(logv, H5::PredType::NATIVE_DOUBLE);

            const hsize_t cell_dims[2] = {1, 1};
            H5::DataSpace cell_space(2, cell_dims);
            const double sigma[1] = {1.0e-18};
            file.createDataSet("sigma_mt_m2", H5::PredType::NATIVE_DOUBLE, cell_space).write(sigma, H5::PredType::NATIVE_DOUBLE);
            file.createDataSet("sigma_event_m2", H5::PredType::NATIVE_DOUBLE, cell_space).write(sigma, H5::PredType::NATIVE_DOUBLE);
            file.createDataSet("b_max_m", H5::PredType::NATIVE_DOUBLE, vec_space).write(sigma, H5::PredType::NATIVE_DOUBLE);
        }

        physics::InteractionPotentialOfflineSampleSet loaded;
        std::string error;
        REQUIRE_FALSE(physics::load_interaction_potential_offline_sample_set_file(sample_path, loaded, &error));
        REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("b_max_m"));

        std::error_code ec;
        std::filesystem::remove(sample_path, ec);
    }

    SECTION("CDF lookup tables must match cell-table shape") {
        const auto sample_path = make_temp_h5_path("cdf_shape");
        write_samples_file(sample_path, "He");
        {
            H5::H5File file(sample_path.string(), H5F_ACC_RDWR);
            H5Ldelete(file.getId(), "cdf_offsets", H5P_DEFAULT);
            const hsize_t bad_dims[1] = {1};
            H5::DataSpace bad_space(1, bad_dims);
            const long long values[1] = {0};
            file.createDataSet("cdf_offsets", H5::PredType::NATIVE_LLONG, bad_space).write(values, H5::PredType::NATIVE_LLONG);
        }

        physics::InteractionPotentialOfflineSampleSet loaded;
        std::string error;
        REQUIRE_FALSE(physics::load_interaction_potential_offline_sample_set_file(sample_path, loaded, &error));
        REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("cdf_offsets"));

        std::error_code ec;
        std::filesystem::remove(sample_path, ec);
    }

    SECTION("cdf_values must be one-dimensional") {
        const auto sample_path = make_temp_h5_path("cdf_values_rank");
        write_samples_file(sample_path, "He");
        {
            H5::H5File file(sample_path.string(), H5F_ACC_RDWR);
            H5Ldelete(file.getId(), "cdf_values", H5P_DEFAULT);
            const hsize_t bad_dims[2] = {1, 1};
            H5::DataSpace bad_space(2, bad_dims);
            const double values[1] = {1.0};
            file.createDataSet("cdf_values", H5::PredType::NATIVE_DOUBLE, bad_space).write(values, H5::PredType::NATIVE_DOUBLE);
        }

        physics::InteractionPotentialOfflineSampleSet loaded;
        std::string error;
        REQUIRE_FALSE(physics::load_interaction_potential_offline_sample_set_file(sample_path, loaded, &error));
        REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("cdf_values"));

        std::error_code ec;
        std::filesystem::remove(sample_path, ec);
    }
}

TEST_CASE("InteractionPotentialCollisionHandler throws on gas mismatch between samples and environment", "[collision][ipm][safety]") {
    const auto sample_path = make_temp_h5_path("gas_mismatch");
    write_samples_file(sample_path, "He");

    config::SpeciesDatabase species_db;
    config::SpeciesProperties sp;
    sp.id = "X+";
    sp.mass_amu = 28.0;
    sp.charge = 1;
    sp.ipm_samples_file = sample_path.string();
    species_db.species.emplace(sp.id, sp);

    physics::InteractionPotentialCollisionHandler handler(false, &species_db);

    IonState ion;
    ion.species_id = "X+";
    ion.mass_kg = 28.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1.0e-18;
    ion.vel = Vec3{400.0, 0.0, 0.0};

    auto ensemble = core::IonEnsemble::from_legacy({ion});
    auto view = ensemble.collision_data(0);

    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.gas_species = "N2";
    env.gas_mass_kg = 28.0 * AMU_TO_KG;
    env.particle_density_m_3 = 1.0e25;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};

    physics::PhysicsRng rng(1234);

    REQUIRE_THROWS_AS(handler.handle_collision(view, 1.0, rng, env), std::runtime_error);

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler destructor flushes logs without collected samples", "[collision][ipm][logging]") {
    const auto sample_path = make_temp_h5_path("flush_logs");
    write_samples_file(sample_path, "He");

    const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const auto log_prefix = (std::filesystem::temp_directory_path() /
        ("icarion_ipm_flush_" + std::to_string(stamp) + ".csv")).string();

    config::SpeciesDatabase species_db;
    config::SpeciesProperties sp;
    sp.id = "X+";
    sp.mass_amu = 28.0;
    sp.charge = 1;
    sp.ipm_samples_file = sample_path.string();
    species_db.species.emplace(sp.id, sp);

    REQUIRE_NOTHROW([&]() {
        physics::InteractionPotentialCollisionHandler handler(
            false,
            &species_db,
            "random",
            0,
            log_prefix,
            log_prefix
        );
    }());

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
    std::filesystem::remove(log_prefix, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler ignores inactive mixture components", "[collision][ipm][multigas]") {
    const auto sample_path = make_temp_h5_path("inactive_component");
    write_samples_file(sample_path, "He");

    config::SpeciesDatabase species_db;
    config::SpeciesProperties sp;
    sp.id = "X+";
    sp.mass_amu = 28.0;
    sp.charge = 1;
    sp.ipm_samples_file = sample_path.string();
    species_db.species.emplace(sp.id, sp);

    physics::InteractionPotentialCollisionHandler handler(false, &species_db);

    IonState ion;
    ion.species_id = "X+";
    ion.mass_kg = 28.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1.0e-18;
    ion.vel = Vec3{400.0, 0.0, 0.0};

    auto ensemble = core::IonEnsemble::from_legacy({ion});
    auto view = ensemble.collision_data(0);

    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 1000.0;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.gas_mixture = {
        {"N2", 0.5, -1.0, -1.0},
        {"He", 0.5, -1.0, -1.0}
    };
    env.gas_mixture[0].participates_in_collisions = false;
    env.gas_mixture[1].participates_in_collisions = true;
    env.compute_derived_properties();

    physics::PhysicsRng rng(42);
    REQUIRE_NOTHROW(handler.handle_collision(view, 1e-6, rng, env));

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler rejects multiple active mixture components", "[collision][ipm][multigas]") {
    const auto sample_path = make_temp_h5_path("multiple_active_components");
    write_samples_file(sample_path, "He");

    auto species_db = make_species_db(sample_path);
    physics::InteractionPotentialCollisionHandler handler(false, &species_db);

    auto ensemble = core::IonEnsemble::from_legacy({make_test_ion(400.0)});
    auto view = ensemble.collision_data(0);

    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 1000.0;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.gas_mixture = {
        {"He", 0.5, -1.0, -1.0},
        {"Ar", 0.5, -1.0, -1.0}
    };
    env.gas_mixture[0].participates_in_collisions = true;
    env.gas_mixture[1].participates_in_collisions = true;
    env.compute_derived_properties();

    physics::PhysicsRng rng(11);
    REQUIRE_THROWS_WITH(
        handler.handle_collision(view, 1.0e-6, rng, env),
        Catch::Matchers::ContainsSubstring("multiple active components")
    );

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler returns false when mixture has no active components", "[collision][ipm][multigas]") {
    const auto sample_path = make_temp_h5_path("no_active_components");
    write_samples_file(sample_path, "He");

    auto species_db = make_species_db(sample_path);
    physics::InteractionPotentialCollisionHandler handler(false, &species_db);

    auto ensemble = core::IonEnsemble::from_legacy({make_test_ion(400.0)});
    auto view = ensemble.collision_data(0);

    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 1000.0;
    env.gas_species = "He";
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.gas_mixture = {
        {"He", 0.5, -1.0, -1.0},
        {"Ar", 0.5, -1.0, -1.0}
    };
    env.gas_mixture[0].participates_in_collisions = false;
    env.gas_mixture[1].participates_in_collisions = false;
    env.compute_derived_properties();

    physics::PhysicsRng rng(2032);
    REQUIRE_FALSE(handler.handle_collision(view, 1.0e-6, rng, env));

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler clamps out-of-range relative velocities to edge bins", "[collision][ipm][coverage]") {
    const auto sample_path = make_temp_h5_path("coverage_clamp");
    write_dp_stats_only_samples_file(
        sample_path,
        "He",
        {std::log(200.0)},
        {1.0e-18},
        {1.0e-22, 0.0, 0.0, 0.0}
    );

    auto species_db = make_species_db(sample_path);
    physics::InteractionPotentialCollisionHandler handler(false, &species_db);

    config::EnvironmentConfig env;
    env.temperature_K = 0.0;
    env.gas_species = "He";
    env.gas_mass_kg = 4.0 * AMU_TO_KG;
    env.particle_density_m_3 = 1.0e25;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};

    physics::PhysicsRng low_rng(123);
    auto low_ensemble = core::IonEnsemble::from_legacy({make_test_ion(100.0)});
    auto low_view = low_ensemble.collision_data(0);
    REQUIRE(handler.handle_collision(low_view, 1.0e-6, low_rng, env));
    CHECK(near_rate(handler.get_stats().average_collision_rate, 1.0e9));
    CHECK(low_view.kin.vel().x > 100.0);

    handler.reset_stats();

    physics::PhysicsRng high_rng(456);
    auto high_ensemble = core::IonEnsemble::from_legacy({make_test_ion(400.0)});
    auto high_view = high_ensemble.collision_data(0);
    REQUIRE(handler.handle_collision(high_view, 1.0e-6, high_rng, env));
    CHECK(near_rate(handler.get_stats().average_collision_rate, 4.0e9));
    CHECK(high_view.kin.vel().x > 400.0);

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler falls back to dp_stats when CDF samples are absent", "[collision][ipm][sampling]") {
    const auto sample_path = make_temp_h5_path("dp_stats_fallback");
    write_dp_stats_only_samples_file(
        sample_path,
        "He",
        {std::log(250.0)},
        {2.0e-18},
        {3.0e-22, 0.0, 0.0, 0.0}
    );

    auto species_db = make_species_db(sample_path);
    physics::InteractionPotentialCollisionHandler handler(false, &species_db, "fixed", 0);

    auto ensemble = core::IonEnsemble::from_legacy({make_test_ion(250.0)});
    auto view = ensemble.collision_data(0);

    config::EnvironmentConfig env;
    env.temperature_K = 0.0;
    env.gas_species = "He";
    env.gas_mass_kg = 4.0 * AMU_TO_KG;
    env.particle_density_m_3 = 1.0e25;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};

    physics::PhysicsRng rng(999);
    REQUIRE(handler.handle_collision(view, 1.0e-6, rng, env));
    CHECK(near_rate(handler.get_stats().average_collision_rate, 5.0e9));
    CHECK(view.kin.vel().x > 250.0);

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler fixed orientation selects the requested offline sample", "[collision][ipm][orientation]") {
    const auto sample_path = make_temp_h5_path("fixed_orientation");
    write_dp_stats_only_samples_file(
        sample_path,
        "He",
        {std::log(300.0)},
        {1.0e-18, 5.0e-18},
        {
            1.0e-22, 0.0, 0.0, 0.0,
            5.0e-22, 0.0, 0.0, 0.0
        },
        2
    );

    auto species_db = make_species_db(sample_path);
    physics::InteractionPotentialCollisionHandler handler(false, &species_db, "fixed", 3);

    auto ensemble = core::IonEnsemble::from_legacy({make_test_ion(300.0)});
    auto view = ensemble.collision_data(0);

    config::EnvironmentConfig env;
    env.temperature_K = 0.0;
    env.gas_species = "He";
    env.gas_mass_kg = 4.0 * AMU_TO_KG;
    env.particle_density_m_3 = 1.0e25;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};

    physics::PhysicsRng rng(2026);
    REQUIRE(handler.handle_collision(view, 1.0e-6, rng, env));
    CHECK(near_rate(handler.get_stats().average_collision_rate, 1.5e10));
    CHECK(view.kin.vel().x > 300.0);

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler random orientation samples multiple offline orientations", "[collision][ipm][orientation]") {
    const auto sample_path = make_temp_h5_path("random_orientation");
    write_dp_stats_only_samples_file(
        sample_path,
        "He",
        {std::log(300.0)},
        {1.0e-18, 5.0e-18},
        {
            1.0e-22, 0.0, 0.0, 0.0,
            5.0e-22, 0.0, 0.0, 0.0
        },
        2
    );

    auto species_db = make_species_db(sample_path);
    physics::InteractionPotentialCollisionHandler handler(false, &species_db, "random", 0);

    config::EnvironmentConfig env;
    env.temperature_K = 0.0;
    env.gas_species = "He";
    env.gas_mass_kg = 4.0 * AMU_TO_KG;
    env.particle_density_m_3 = 1.0e25;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};

    bool saw_first_orientation = false;
    bool saw_second_orientation = false;
    physics::PhysicsRng rng(7);

    for (int attempt = 0; attempt < 32; ++attempt) {
        auto ensemble = core::IonEnsemble::from_legacy({make_test_ion(300.0)});
        auto view = ensemble.collision_data(0);

        handler.reset_stats();
        REQUIRE(handler.handle_collision(view, 1.0e-6, rng, env));
        const double rate = handler.get_stats().average_collision_rate;
        if (near_rate(rate, 3.0e9)) {
            saw_first_orientation = true;
        } else if (near_rate(rate, 1.5e10)) {
            saw_second_orientation = true;
        } else {
            FAIL("Unexpected sigma_mt selected by random orientation path");
        }

        if (saw_first_orientation && saw_second_orientation) {
            break;
        }
    }

    CHECK(saw_first_orientation);
    CHECK(saw_second_orientation);

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}

TEST_CASE("InteractionPotentialCollisionHandler parses orientation_mode case-insensitively", "[collision][ipm][orientation]") {
    const auto sample_path = make_temp_h5_path("orientation_case_insensitive");
    write_dp_stats_only_samples_file(
        sample_path,
        "He",
        {std::log(300.0)},
        {1.0e-18, 5.0e-18},
        {
            1.0e-22, 0.0, 0.0, 0.0,
            5.0e-22, 0.0, 0.0, 0.0
        },
        2
    );

    auto species_db = make_species_db(sample_path);

    config::EnvironmentConfig env;
    env.temperature_K = 0.0;
    env.gas_species = "He";
    env.gas_mass_kg = 4.0 * AMU_TO_KG;
    env.particle_density_m_3 = 1.0e25;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};

    SECTION("Fixed is treated like fixed") {
        physics::InteractionPotentialCollisionHandler handler(false, &species_db, "Fixed", 3);
        auto ensemble = core::IonEnsemble::from_legacy({make_test_ion(300.0)});
        auto view = ensemble.collision_data(0);
        physics::PhysicsRng rng(2040);

        REQUIRE(handler.handle_collision(view, 1.0e-6, rng, env));
        CHECK(near_rate(handler.get_stats().average_collision_rate, 1.5e10));
    }

    SECTION("RANDOM is treated like random") {
        physics::InteractionPotentialCollisionHandler handler(false, &species_db, "RANDOM", 0);
        bool saw_first_orientation = false;
        bool saw_second_orientation = false;
        physics::PhysicsRng rng(2041);

        for (int attempt = 0; attempt < 32; ++attempt) {
            auto ensemble = core::IonEnsemble::from_legacy({make_test_ion(300.0)});
            auto view = ensemble.collision_data(0);

            handler.reset_stats();
            REQUIRE(handler.handle_collision(view, 1.0e-6, rng, env));
            const double rate = handler.get_stats().average_collision_rate;
            if (near_rate(rate, 3.0e9)) {
                saw_first_orientation = true;
            } else if (near_rate(rate, 1.5e10)) {
                saw_second_orientation = true;
            } else {
                FAIL("Unexpected sigma_mt selected by uppercase random orientation path");
            }

            if (saw_first_orientation && saw_second_orientation) {
                break;
            }
        }

        CHECK(saw_first_orientation);
        CHECK(saw_second_orientation);
    }

    std::error_code ec;
    std::filesystem::remove(sample_path, ec);
}
