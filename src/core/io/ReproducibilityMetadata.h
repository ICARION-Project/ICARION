// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <H5Cpp.h>
#include <cstdint>
#include <string>
#include <vector>
#include "UserAnnotation.h"

namespace ICARION::io::reproducibility {

struct InputFile {
    std::string key;
    std::string filename;
    std::string role;
    std::string sha256;
    std::string content;
    bool used = false;
    bool required = false;
};

struct IpmMetadata {
    std::string species_id;
    double charge_e = 0.0;
    double ion_mass_u = 0.0;
    double ion_mass_kg = 0.0;
    double reference_temperature_K = 0.0;
    std::string species_entry_json;
    std::string geometry_filename;
    std::size_t atom_count = 0;
    double geometry_total_mass_u = 0.0;

    std::string gas_id;
    std::string gas_model;
    double gas_mass_u = 0.0;
    double gas_mass_kg = 0.0;
    double sigma_m = 0.0;
    double epsilon_J = 0.0;
    double polarizability_m3 = 0.0;
    double polarizability_parallel_m3 = 0.0;
    double polarizability_perpendicular_m3 = 0.0;
    double n2_bond_m = 0.0;
    bool n2_quadrupole = false;
    double n2_q_site_e = 0.0;
    double n2_q_center_e = 0.0;
    bool n2_average_lj = true;
    bool n2_anisotropic_polarizability = false;
    std::string sigma_source;
    std::string epsilon_source;
    std::string polarizability_source;
    std::string polarizability_parallel_source;
    std::string polarizability_perpendicular_source;

    std::uint64_t seed = 0;
    bool seed_explicit = false;
    std::string orientation_sampling_mode;

    std::string potential_model;
    std::string parameter_model;
    std::string mixing_rule;
    std::string polarization_model;
    int orientations = 0;
    int trials_per_cell = 0;
    int velocity_bins = 0;
    double temperature_K = 0.0;
    double velocity_min_m_s = 0.0;
    double velocity_max_m_s = 0.0;
    double deflection_epsilon = 0.0;
    double impact_parameter_guess_m = 0.0;
    double impact_parameter_growth_factor = 0.0;
    double impact_parameter_relative_tolerance = 0.0;
    double maximum_non_asymptotic_fraction = 0.0;
    double integration_eta_dt = 0.0;
    int maximum_trajectory_steps = 0;
    double sigma_scale_factor = 0.0;
    double epsilon_scale_factor = 0.0;
    double polarization_damping_radius_m = 0.0;
    double mmff_energy_scale_factor = 0.0;
    double mmff_distance_scale_factor = 0.0;
    bool full_cdf = false;
    int checkpoint_interval_cells = 0;
    bool resume_requested = false;
    bool resume_used = false;
    int openmp_threads = 1;
    std::string command_line;
    std::string resolved_options_json;
    std::vector<InputFile> inputs;

    std::string start_timestamp_utc;
    std::string completion_timestamp_utc;
    double wall_clock_runtime_s = 0.0;
    double accumulated_wall_clock_runtime_s = 0.0;
    std::uint64_t completed_cells = 0;
    std::uint64_t total_cells = 0;
    bool success = false;
    bool checkpoint = true;
    UserAnnotation user_annotation;
};

std::string utc_timestamp_now();
InputFile capture_input_file(const std::string& key, const std::string& path,
                             const std::string& role, bool used, bool required);
void write_ipm_metadata(H5::H5File& file, const IpmMetadata& metadata, int data_format_version);

} // namespace ICARION::io::reproducibility
