// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/catch_test_macros.hpp>
#include <H5Cpp.h>
#include <chrono>
#include <filesystem>
#include <fstream>

#include "core/physics/collisions/EHSSOfflineSampleSet.h"
#include "core/physics/collisions/EHSSSamples.h"

using ICARION::physics::EHSSOfflineSampleSet;
using ICARION::physics::EHSSOrientationSamples;
using ICARION::physics::load_ehss_offline_sample_set_file;
using ICARION::physics::load_ehss_samples_file;

namespace {

std::filesystem::path make_temp_h5_path(const std::string& tag) {
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path()
        / ("ehss_offline_samples_" + tag + "_" + std::to_string(stamp) + ".h5");
}

void write_ehss_offline_h5_fixture(
    const std::filesystem::path& path,
    long long n_orientations_attr,
    long long n_mu_samples_attr,
    const std::vector<double>& sigma,
    const std::vector<double>& mu_flat,
    size_t mu_stride
) {
    H5::H5File file(path.string(), H5F_ACC_TRUNC);

    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);

    const char* format = "ehss_offline_samples";
    H5::Attribute format_attr = file.createAttribute("format", str_type, scalar_space);
    format_attr.write(str_type, &format);

    const char* units = "sigma_eff_m2,mu";
    H5::Attribute units_attr = file.createAttribute("units", str_type, scalar_space);
    units_attr.write(str_type, &units);

    const char* gas = "He";
    H5::Attribute gas_attr = file.createAttribute("gas", str_type, scalar_space);
    gas_attr.write(str_type, &gas);

    const long long version = 1;
    H5::Attribute version_attr = file.createAttribute("version", H5::PredType::NATIVE_LLONG, scalar_space);
    version_attr.write(H5::PredType::NATIVE_LLONG, &version);

    H5::Attribute orient_attr =
        file.createAttribute("n_orientations", H5::PredType::NATIVE_LLONG, scalar_space);
    orient_attr.write(H5::PredType::NATIVE_LLONG, &n_orientations_attr);

    H5::Attribute mu_count_attr =
        file.createAttribute("n_mu_samples", H5::PredType::NATIVE_LLONG, scalar_space);
    mu_count_attr.write(H5::PredType::NATIVE_LLONG, &n_mu_samples_attr);

    {
        const hsize_t dims[1] = {static_cast<hsize_t>(sigma.size())};
        H5::DataSpace space(1, dims);
        H5::DataSet dset = file.createDataSet("sigma_eff_m2", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(sigma.data(), H5::PredType::NATIVE_DOUBLE);
    }

    {
        const hsize_t dims[2] = {
            static_cast<hsize_t>(sigma.size()),
            static_cast<hsize_t>(mu_stride)
        };
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("mu_samples", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(mu_flat.data(), H5::PredType::NATIVE_DOUBLE);
    }
}

} // namespace

TEST_CASE("EHSS samples file loads", "[collision][ehss][samples]") {
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("ehss_samples_" + std::to_string(stamp) + ".json");

    std::ofstream ofs(path);
    REQUIRE(ofs.is_open());
    ofs << "{\n"
        << "  \"version\": 1,\n"
        << "  \"species_id\": \"X+\",\n"
        << "  \"n_orientations\": 2,\n"
        << "  \"orientations_quat\": [[1,0,0,0],[0,1,0,0]],\n"
        << "  \"areas_by_gas_m2\": {\n"
        << "    \"He\": [1.0e-18, 2.0e-18]\n"
        << "  }\n"
        << "}\n";
    ofs.close();

    EHSSOrientationSamples samples;
    std::string error;
    REQUIRE(load_ehss_samples_file(path, samples, &error));
    REQUIRE(samples.orientations_quat.size() == 2);
    REQUIRE(samples.areas_by_gas_m2.count("He") == 1);
    REQUIRE(samples.areas_by_gas_m2.at("He").size() == 2);

    std::filesystem::remove(path);
}

TEST_CASE("EHSS offline samples file loads", "[collision][ehss][samples][offline]") {
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("ehss_offline_samples_" + std::to_string(stamp) + ".json");

    std::ofstream ofs(path);
    REQUIRE(ofs.is_open());
    ofs << "{\n"
        << "  \"gas\": \"He\",\n"
        << "  \"sigma_eff_m2\": [1.0e-18, 2.0e-18],\n"
        << "  \"mu_samples\": [\n"
        << "    [-0.5, 0.0, 0.5],\n"
        << "    [-0.2, 0.2, 0.9]\n"
        << "  ]\n"
        << "}\n";
    ofs.close();

    EHSSOfflineSampleSet samples;
    std::string error;
    REQUIRE(load_ehss_offline_sample_set_file(path, samples, &error));
    REQUIRE(samples.gas == "He");
    REQUIRE(samples.sigma_eff_m2.size() == 2);
    REQUIRE(samples.mu_stride == 3);
    REQUIRE(samples.mu_at(1, 2) == Catch::Approx(0.9));

    std::filesystem::remove(path);
}

TEST_CASE("EHSS offline samples reject non-physical values", "[collision][ehss][samples][offline]") {
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("ehss_offline_samples_bad_" + std::to_string(stamp) + ".json");

    std::ofstream ofs(path);
    REQUIRE(ofs.is_open());
    ofs << "{\n"
        << "  \"gas\": \"He\",\n"
        << "  \"sigma_eff_m2\": [0.0, 2.0e-18],\n"
        << "  \"mu_samples\": [\n"
        << "    [-0.5, 0.0, 1.5],\n"
        << "    [-0.2, 0.2, 0.9]\n"
        << "  ]\n"
        << "}\n";
    ofs.close();

    EHSSOfflineSampleSet samples;
    std::string error;
    REQUIRE_FALSE(load_ehss_offline_sample_set_file(path, samples, &error));
    REQUIRE_FALSE(error.empty());

    std::filesystem::remove(path);
}

TEST_CASE("EHSS offline samples reject inconsistent metadata", "[collision][ehss][samples][offline]") {
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("ehss_offline_samples_meta_bad_" + std::to_string(stamp) + ".json");

    std::ofstream ofs(path);
    REQUIRE(ofs.is_open());
    ofs << "{\n"
        << "  \"version\": 2,\n"
        << "  \"format\": \"wrong_format\",\n"
        << "  \"gas\": \"He\",\n"
        << "  \"n_orientations\": 3,\n"
        << "  \"n_mu_samples\": 4,\n"
        << "  \"sigma_eff_m2\": [1.0e-18, 2.0e-18],\n"
        << "  \"mu_samples\": [\n"
        << "    [-0.5, 0.0, 0.5],\n"
        << "    [-0.2, 0.2, 0.9]\n"
        << "  ]\n"
        << "}\n";
    ofs.close();

    EHSSOfflineSampleSet samples;
    std::string error;
    REQUIRE_FALSE(load_ehss_offline_sample_set_file(path, samples, &error));
    REQUIRE_FALSE(error.empty());

    std::filesystem::remove(path);
}

TEST_CASE("EHSS offline HDF5 samples reject inconsistent metadata", "[collision][ehss][samples][offline]") {
    std::filesystem::path path = make_temp_h5_path("meta_bad");

    const std::vector<double> sigma = {1.0e-18, 2.0e-18};
    const std::vector<double> mu = {
        -0.5, 0.0, 0.5,
        -0.2, 0.2, 0.9
    };
    write_ehss_offline_h5_fixture(path, 3, 3, sigma, mu, 3);

    EHSSOfflineSampleSet samples;
    std::string error;
    REQUIRE_FALSE(load_ehss_offline_sample_set_file(path, samples, &error));
    REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("n_orientations"));

    std::filesystem::remove(path);
}

TEST_CASE("EHSS offline HDF5 samples reject non-physical values", "[collision][ehss][samples][offline]") {
    std::filesystem::path path = make_temp_h5_path("nonphysical");

    const std::vector<double> sigma = {1.0e-18, 0.0};
    const std::vector<double> mu = {
        -0.5, 0.0, 0.5,
        -0.2, 0.2, 0.9
    };
    write_ehss_offline_h5_fixture(path, 2, 3, sigma, mu, 3);

    EHSSOfflineSampleSet samples;
    std::string error;
    REQUIRE_FALSE(load_ehss_offline_sample_set_file(path, samples, &error));
    REQUIRE_THAT(error, Catch::Matchers::ContainsSubstring("sigma_eff_m2"));

    std::filesystem::remove(path);
}
