// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
#include "ReproducibilityMetadata.h"

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <filesystem>
#include <iomanip>
#include <openssl/sha.h>

#ifndef ICARION_BUILD_CONFIGURATION
#error "ICARION_BUILD_CONFIGURATION must be supplied by CMake"
#endif
#ifndef ICARION_GIT_STATE
#error "ICARION_GIT_STATE must be supplied by CMake"
#endif
#ifndef ICARION_GIT_STATE_CAPTURE
#error "ICARION_GIT_STATE_CAPTURE must be supplied by CMake"
#endif
#ifndef ICARION_CMAKE_GPU_ENABLED
#error "ICARION_CMAKE_GPU_ENABLED must be supplied by CMake"
#endif
#ifndef ICARION_CMAKE_CORE_ONLY_ENABLED
#error "ICARION_CMAKE_CORE_ONLY_ENABLED must be supplied by CMake"
#endif
#ifdef ICARION_USE_GPU
static_assert(ICARION_CMAKE_GPU_ENABLED == 1, "ICARION_USE_GPU disagrees with CMake GPU configuration");
#else
static_assert(ICARION_CMAKE_GPU_ENABLED == 0, "ICARION_USE_GPU missing for a GPU build");
#endif
#ifdef ICARION_BUILD_CORE_ONLY
static_assert(ICARION_CMAKE_CORE_ONLY_ENABLED == 1, "ICARION_BUILD_CORE_ONLY disagrees with CMake");
#else
static_assert(ICARION_CMAKE_CORE_ONLY_ENABLED == 0, "ICARION_BUILD_CORE_ONLY missing for core-only build");
#endif

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/utsname.h>
#include <unistd.h>
#endif

namespace ICARION::io::reproducibility {
namespace {
void string(H5::Group& group, const std::string& name, const std::string& value) {
    H5::StrType type(H5::PredType::C_S1, value.size() + 1);
    H5::DataSpace space(H5S_SCALAR);
    auto dataset = group.createDataSet(name, type, space);
    dataset.write(value.c_str(), type);
}
template<class T> void scalar(H5::Group& group, const std::string& name, const T& value, const H5::PredType& type) {
    H5::DataSpace space(H5S_SCALAR);
    auto dataset = group.createDataSet(name, type, space);
    dataset.write(&value, type);
}
void boolean(H5::Group& group, const std::string& name, bool value) {
    hbool_t stored = value ? 1 : 0;
    scalar(group, name, stored, H5::PredType::NATIVE_HBOOL);
}
void integer(H5::Group& g, const std::string& n, int v) { scalar(g, n, v, H5::PredType::NATIVE_INT); }
void uinteger(H5::Group& g, const std::string& n, std::uint64_t v) { scalar(g, n, v, H5::PredType::NATIVE_UINT64); }
void real(H5::Group& g, const std::string& n, double v) { scalar(g, n, v, H5::PredType::NATIVE_DOUBLE); }

void os_kernel(std::string& os, std::string& kernel) {
#ifdef _WIN32
    os = "Windows"; kernel = "Windows";
#else
    struct utsname info{};
    if (uname(&info) == 0) { os = info.sysname; kernel = std::string(info.release) + " " + info.version; }
    else { os = kernel = "unknown"; }
#endif
}
std::string compiler_name() {
#if defined(__clang__)
    return "Clang";
#elif defined(__GNUC__)
    return "GNU";
#elif defined(_MSC_VER)
    return "MSVC";
#else
    return "unknown";
#endif
}
std::string compiler_version() {
#if defined(__clang__) || defined(__GNUC__)
    return __VERSION__;
#elif defined(_MSC_VER)
    return std::to_string(_MSC_VER);
#else
    return "unknown";
#endif
}
}

std::string utc_timestamp_now() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t value = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &value);
#else
    gmtime_r(&value, &utc);
#endif
    char buffer[32]{};
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &utc);
    return buffer;
}

InputFile capture_input_file(const std::string& key, const std::string& path,
                             const std::string& role, bool used, bool required) {
    InputFile result;
    result.key = key;
    result.filename = path.empty() ? "" : std::filesystem::path(path).filename().string();
    result.role = role;
    result.used = used;
    result.required = required;
    if (!used || path.empty()) return result;
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        if (required) throw std::runtime_error("Cannot snapshot required input: " + path);
        return result;
    }
    std::ostringstream bytes;
    bytes << input.rdbuf();
    if (input.bad()) throw std::runtime_error("Cannot read input snapshot: " + path);
    result.content = bytes.str();
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(result.content.data()), result.content.size(), digest);
    std::ostringstream hex;
    for (unsigned char value : digest) hex << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(value);
    result.sha256 = hex.str();
    return result;
}

static void write_ipm_system_metadata(H5::Group& system) {
    std::string os, kernel;
    os_kernel(os, kernel);
    string(system, "operating_system", os);
    string(system, "kernel", kernel);
    string(system, "timestamp_utc", utc_timestamp_now());
    const int cpu_count = static_cast<int>(std::thread::hardware_concurrency());
    integer(system, "logical_cpu_count", cpu_count);
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line, cpu_model;
    while (std::getline(cpuinfo, line)) {
        if (cpu_model.empty() && line.rfind("model name", 0) == 0) cpu_model = line.substr(line.find(':') + 2);
    }
    string(system, "cpu_model", cpu_model.empty() ? "unknown" : cpu_model);
    std::ifstream meminfo("/proc/meminfo");
    std::uint64_t memory_bytes = 0;
    while (std::getline(meminfo, line)) {
        if (line.rfind("MemTotal:", 0) == 0) {
            std::istringstream(line.substr(9)) >> memory_bytes;
            memory_bytes *= 1024;
            break;
        }
    }
    uinteger(system, "total_memory_bytes", memory_bytes);
}

void write_ipm_metadata(H5::H5File& file, const IpmMetadata& m, int format_version) {
    auto metadata = file.createGroup("/metadata");
    auto schema = metadata.createGroup("schema");
    string(schema, "metadata_schema_version", "1.0.0");
    string(schema, "data_format", "ipm_offline_samples");
    integer(schema, "data_format_version", format_version);
    string(schema, "description", "ICARION interaction-potential offline collision samples with reproducibility metadata");

    auto software = metadata.createGroup("software");
#ifdef ICARION_VERSION
    string(software, "icarion_version", ICARION_VERSION);
#else
    string(software, "icarion_version", "unknown");
#endif
#ifdef GIT_HASH
    string(software, "git_hash", GIT_HASH);
#else
    string(software, "git_hash", "unknown");
#endif
    string(software, "git_state", ICARION_GIT_STATE);
    string(software, "git_state_capture", ICARION_GIT_STATE_CAPTURE);
    string(software, "build_type", ICARION_BUILD_CONFIGURATION);
    string(software, "compiler", compiler_name());
    string(software, "compiler_version", compiler_version());
    integer(software, "cxx_standard", static_cast<int>(__cplusplus));
#ifdef _OPENMP
    boolean(software, "openmp_enabled", true);
    integer(software, "openmp_thread_count", m.openmp_threads);
#else
    boolean(software, "openmp_enabled", false);
    integer(software, "openmp_thread_count", 1);
#endif
#ifdef ICARION_USE_GPU
    boolean(software, "cuda_enabled", true);
#else
    boolean(software, "cuda_enabled", false);
#endif
#ifdef ICARION_BUILD_CORE_ONLY
    string(software, "build_mode", "core-only");
#else
    string(software, "build_mode", "full");
#endif

    auto system = metadata.createGroup("system");
    write_ipm_system_metadata(system);

    auto rng = metadata.createGroup("rng");
    uinteger(rng, "base_seed", m.seed); string(rng, "algorithm", "std::mt19937_64");
    string(rng, "orientation_sampling_mode", m.orientation_sampling_mode);
    string(rng, "cell_seed_scheme", "splitmix64(base_seed XOR (0x9E3779B97F4A7C15 * (orientation_index + 1)) XOR (0xBF58476D1CE4E5B9 * (velocity_bin_index + 1)))");
    boolean(rng, "seed_supplied_explicitly", m.seed_explicit);
    string(rng, "seed_origin", m.seed_explicit ? "command_line" : "std::hash(species_id + '|' + gas_id + '|ipm_v1'); implementation-dependent std::hash");

    auto species = metadata.createGroup("species");
    string(species, "species_id", m.species_id); real(species, "charge_e", m.charge_e);
    real(species, "ion_mass_u", m.ion_mass_u); real(species, "ion_mass_kg", m.ion_mass_kg);
    real(species, "reference_temperature_K", m.reference_temperature_K);
    string(species, "database_entry_json", m.species_entry_json);
    string(species, "molecular_geometry_filename", m.geometry_filename);
    integer(species, "number_of_atoms", static_cast<int>(m.atom_count));
    real(species, "geometry_total_mass_u", m.geometry_total_mass_u);

    auto neutral = metadata.createGroup("neutral");
    string(neutral, "gas_id", m.gas_id); string(neutral, "gas_model", m.gas_model);
    real(neutral, "gas_mass_u", m.gas_mass_u); real(neutral, "gas_mass_kg", m.gas_mass_kg);
    real(neutral, "sigma_m", m.sigma_m); real(neutral, "epsilon_J", m.epsilon_J);
    real(neutral, "polarizability_m3", m.polarizability_m3);
    real(neutral, "polarizability_parallel_m3", m.polarizability_parallel_m3);
    real(neutral, "polarizability_perpendicular_m3", m.polarizability_perpendicular_m3);
    real(neutral, "n2_bond_m", m.n2_bond_m); boolean(neutral, "n2_quadrupole_enabled", m.n2_quadrupole);
    real(neutral, "n2_q_site_e", m.n2_q_site_e); real(neutral, "n2_q_center_e", m.n2_q_center_e);
    boolean(neutral, "n2_average_lj", m.n2_average_lj);
    boolean(neutral, "n2_anisotropic_polarizability", m.n2_anisotropic_polarizability);
    string(neutral, "sigma_source", m.sigma_source);
    string(neutral, "epsilon_source", m.epsilon_source);
    string(neutral, "polarizability_source", m.polarizability_source);
    string(neutral, "polarizability_parallel_source", m.polarizability_parallel_source);
    string(neutral, "polarizability_perpendicular_source", m.polarizability_perpendicular_source);

    auto pre = metadata.createGroup("precompute");
    string(pre, "potential_model", m.potential_model); string(pre, "parameter_model", m.parameter_model);
    string(pre, "mixing_rule", m.mixing_rule); string(pre, "polarization_model", m.polarization_model);
    string(pre, "orientation_grid_type", m.orientation_sampling_mode);
    integer(pre, "number_of_orientations", m.orientations); integer(pre, "monte_carlo_trials_per_cell", m.trials_per_cell);
    integer(pre, "number_of_velocity_bins", m.velocity_bins); real(pre, "temperature_K", m.temperature_K);
    real(pre, "minimum_relative_velocity_m_s", m.velocity_min_m_s); real(pre, "maximum_relative_velocity_m_s", m.velocity_max_m_s);
    real(pre, "deflection_epsilon", m.deflection_epsilon); real(pre, "impact_parameter_initial_guess_m", m.impact_parameter_guess_m);
    real(pre, "impact_parameter_growth_factor", m.impact_parameter_growth_factor);
    real(pre, "impact_parameter_relative_tolerance", m.impact_parameter_relative_tolerance);
    real(pre, "maximum_non_asymptotic_fraction", m.maximum_non_asymptotic_fraction);
    real(pre, "integration_eta_dt", m.integration_eta_dt); integer(pre, "maximum_trajectory_steps", m.maximum_trajectory_steps);
    real(pre, "sigma_scale_factor", m.sigma_scale_factor); real(pre, "epsilon_scale_factor", m.epsilon_scale_factor);
    real(pre, "polarization_damping_radius_m", m.polarization_damping_radius_m);
    real(pre, "mmff_energy_scale_factor", m.mmff_energy_scale_factor); real(pre, "mmff_distance_scale_factor", m.mmff_distance_scale_factor);
    boolean(pre, "full_cdf_mode", m.full_cdf); integer(pre, "checkpoint_interval_cells", m.checkpoint_interval_cells);
    boolean(pre, "resume_mode_requested", m.resume_requested); integer(pre, "effective_openmp_thread_count", m.openmp_threads);
    string(pre, "command_line", m.command_line); string(pre, "resolved_options_json", m.resolved_options_json);

    auto inputs = metadata.createGroup("inputs"); auto hashes = inputs.createGroup("hashes"); auto blobs = inputs.createGroup("blobs");
    std::unordered_map<std::string, std::string> embedded_by_hash;
    for (const auto& input : m.inputs) {
        string(hashes, input.key + "_filename", input.filename);
        string(hashes, input.key + "_role", input.role);
        boolean(hashes, input.key + "_used", input.used);
        if (!input.used) { string(hashes, input.key + "_sha256", "not_used"); continue; }
        if (input.sha256.empty()) {
            if (input.required) throw std::runtime_error("Required input snapshot has no hash: " + input.key);
            string(hashes, input.key + "_sha256", "unavailable");
        } else {
            string(hashes, input.key + "_sha256", input.sha256);
            const auto found = embedded_by_hash.find(input.sha256);
            if (found == embedded_by_hash.end()) {
                string(blobs, input.key, input.content); embedded_by_hash.emplace(input.sha256, input.key);
            } else {
                string(blobs, input.key + "_duplicate_of", found->second);
            }
        }
    }
    // Convenient species aliases make the two primary provenance values discoverable.
    for (const auto& input : m.inputs) if (input.key == "molecular_geometry" && input.used) {
        string(species, "molecular_geometry_sha256", input.sha256);
        string(species, "molecular_geometry_content", input.content);
    }

    auto completion = metadata.createGroup("completion");
    boolean(completion, "success", m.success); uinteger(completion, "completed_cells", m.completed_cells);
    uinteger(completion, "total_cells", m.total_cells); string(completion, "start_timestamp_utc", m.start_timestamp_utc);
    string(completion, "completion_timestamp_utc", m.completion_timestamp_utc);
    real(completion, "wall_clock_runtime_s", m.wall_clock_runtime_s); boolean(completion, "is_checkpoint", m.checkpoint);
    real(completion, "accumulated_wall_clock_runtime_s", m.accumulated_wall_clock_runtime_s);
    boolean(completion, "resume_mode_used", m.resume_used);
}
} // namespace ICARION::io::reproducibility
