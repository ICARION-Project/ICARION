// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "LegacyAdapter.h"
#include <stdexcept>

namespace ICARION::config {

// ====================================================================
// Enum conversion helpers
// ====================================================================

// NOTE: convert_collision_model() REMOVED in Phase 2E
// SSOT: core::CollisionModel is now an alias to config::CollisionModel
// No conversion needed - both types are identical!

static core::SolverType convert_solver_type(config::SolverType solver) {
    switch (solver) {
        case config::SolverType::RK4:
            return core::SolverType::RK4;
        case config::SolverType::RK45:
            return core::SolverType::RK45;
        case config::SolverType::Boris:
            return core::SolverType::Boris;
        default:
            return core::SolverType::RK45;
    }
}

core::GlobalParams LegacyAdapter::to_global_params(const FullConfig& config) {
    core::GlobalParams g{};
    
    // === Simulation parameters ===
    g.dt_s = config.simulation.dt_s;
    g.sim_time_steps = config.simulation.total_steps;
    g.write_interval = config.simulation.write_interval;
    g.t_eval = config.simulation.t_eval;
    
    // === Execution flags ===
    g.enable_gpu = config.simulation.enable_gpu;
    g.parallelization = config.simulation.enable_openmp;
    g.rng_seed = config.simulation.rng_seed;
    
    // === Physics ===
    g.collisionModel = config.physics.collision_model;  // SSOT: Direct assignment (Phase 2E)
    g.enable_reactions = config.physics.enable_reactions;
    g.enable_space_charge = config.physics.enable_space_charge;
    g.enable_ou_thermalization = config.physics.enable_ou_thermalization;
    g.force_ou_for_stochastic_models = config.physics.force_ou_for_stochastic;
    
    // === Output ===
    g.output_file = config.output.folder + "/" + config.output.trajectory_file;
    g.print_results = config.output.print_progress;
    
    // === Database/file paths ===
    g.species_database_file = config.species_database_path;
    g.reaction_file = config.reaction_database_path;
    g.ion_cloud_file = config.ion_cloud_path;
    
    // === Continue mode ===
    g.continue_from = config.simulation.continue_from;
    g.continue_time_s = config.simulation.continue_time_s;
    g.auto_continue_if_active = config.simulation.auto_continue_if_active;
    
    // === Derived: num_ions will be set by init_ions() ===
    g.num_ions = 0;  // Filled later by ion loader
    
    return g;
}

std::vector<core::InstrumentDomain> LegacyAdapter::to_instrument_domains(const FullConfig& config) {
    std::vector<core::InstrumentDomain> domains;
    domains.reserve(config.domains.size());
    
    for (const auto& domain_cfg : config.domains) {
        domains.push_back(convert_domain(domain_cfg));
    }
    
    return domains;
}

core::InstrumentDomain LegacyAdapter::convert_domain(const DomainConfig& domain) {
    core::InstrumentDomain dom{};
    
    // === Identification ===
    dom.index = domain.domain_index;
    dom.instrument = domain.instrument;
    dom.solver_type = convert_solver_type(domain.solver);
    
    // === Geometry ===
    dom.geom = convert_geometry(domain.geometry);
    
    // === Environment ===
    dom.env = convert_environment(domain.environment);
    
    // === Fields ===
    convert_fields(domain.fields, dom);
    
    // === Field arrays (if any) ===
    // Note: Field array loading happens separately via load_field_array()
    // We just need to pass the file paths for now
    dom.FA_file = "";  // Legacy single file (unused in new system)
    dom.fieldArrayLoaded = false;
    
    // Convert field array terms
    for (const auto& term : domain.fields.field_array_terms) {
        core::InstrumentDomain::FieldArrayTerm legacy_term;
        legacy_term.file = term.file;
        legacy_term.constant = term.constant;
        legacy_term.phase_rad = term.phase_rad;
        legacy_term.frequency_Hz = term.frequency_Hz;
        legacy_term.loaded = false;
        
        // Convert scale kind enum
        switch (term.kind) {
            case FieldsConfig::FieldArrayTerm::ScaleKind::Constant:
                legacy_term.kind = core::InstrumentDomain::FAScaleKind::Constant;
                break;
            case FieldsConfig::FieldArrayTerm::ScaleKind::DC_Axial:
                legacy_term.kind = core::InstrumentDomain::FAScaleKind::DC_Axial;
                break;
            case FieldsConfig::FieldArrayTerm::ScaleKind::DC_Quad:
                legacy_term.kind = core::InstrumentDomain::FAScaleKind::DC_Quad;
                break;
            case FieldsConfig::FieldArrayTerm::ScaleKind::DC_Radial:
                legacy_term.kind = core::InstrumentDomain::FAScaleKind::DC_Radial;
                break;
            case FieldsConfig::FieldArrayTerm::ScaleKind::RF:
                legacy_term.kind = core::InstrumentDomain::FAScaleKind::RF;
                break;
        }
        
        dom.FA_terms.push_back(legacy_term);
    }
    
    return dom;
}

core::Geometry LegacyAdapter::convert_geometry(const GeometryConfig& geom) {
    core::Geometry g{};
    
    g.length_m = geom.length_m;
    g.radius_m = geom.radius_m;
    g.radius_in_m = geom.radius_in_m;
    g.radius_out_m = geom.radius_out_m;
    g.radius_char_m = geom.radius_char_m;
    g.acc_length_m = geom.acc_length_m;
    g.end_aperture_m = geom.end_aperture_m;
    g.origin_m = geom.origin_m;
    
    // Bounds (GeometryConfig stores as Vec3)
    // Note: GeometryConfig::compute_bounds() must have been called
    g.min_bound = geom.min_bound;
    g.max_bound = geom.max_bound;
    
    return g;
}

core::Environment LegacyAdapter::convert_environment(const EnvironmentConfig& env) {
    core::Environment e{};
    
    // Input parameters
    e.pressure_Pa = env.pressure_Pa;
    e.temperature_K = env.temperature_K;
    e.neutral_species_id = env.gas_species;
    e.gas_velocity_m_s = env.gas_velocity_m_s;
    
    // Derived quantities (computed by EnvironmentConfig::compute_derived_properties())
    e.neutral_mass_kg = env.gas_mass_kg;
    e.particle_density_m_3 = env.particle_density_m_3;
    e.neutral_polarizability_m3 = env.gas_polarizability_m3;
    e.neutral_radius_m = env.gas_radius_m;
    e.mean_thermal_velocity_m_s = env.mean_thermal_velocity_m_s;
    
    return e;
}

void LegacyAdapter::convert_fields(const FieldsConfig& fields, core::InstrumentDomain& dom) {
    // === DC fields ===
    dom.DC.axial_V = fields.dc.axial_V;
    dom.DC.EN_Td = fields.dc.EN_Td;
    dom.DC.EN_Vm2 = fields.dc.EN_Vm2;
    dom.DC.quad_V = fields.dc.quad_V;
    dom.DC.radial_V = fields.dc.radial_V;
    // Note: DC voltage sweeps not yet in new config system
    dom.DC.enable_radial_voltage_sweep = false;
    
    // === RF fields ===
    dom.RF.voltage_V = fields.rf.voltage_V;
    dom.RF.frequency_Hz = fields.rf.frequency_Hz;
    dom.RF.angular_frequency_rad_s = fields.rf.angular_frequency_rad_s;
    dom.RF.phase_rad = fields.rf.phase_rad;
    
    // === AC fields ===
    dom.AC.voltage_V = fields.ac.voltage_V;
    dom.AC.frequency_Hz = fields.ac.frequency_Hz;
    dom.AC.angular_frequency_rad_s = fields.ac.angular_frequency_rad_s;
    
    // Voltage sweep
    dom.AC.enable_voltage_sweep = fields.ac.enable_voltage_sweep;
    dom.AC.amplitude_slope_V_s = fields.ac.amplitude_slope_V_s;
    dom.AC.start_time_s = fields.ac.start_time_s;
    dom.AC.rise_time_s = fields.ac.rise_time_s;
    
    // Frequency sweep
    dom.AC.enable_frequency_sweep = fields.ac.enable_frequency_sweep;
    dom.AC.ac_start_freq_Hz = fields.ac.frequency_start_Hz;
    dom.AC.ac_sweep_slope_Hz_per_s = fields.ac.frequency_sweep_slope_Hz_s;
    
    // LQIT lock
    dom.AC.lqit_lock_enable = fields.ac.lqit_lock_enable;
    dom.AC.lqit_lock_phase_rad = fields.ac.lqit_lock_phase_rad;
    dom.AC.lqit_lock_bandwidth_Hz = fields.ac.lqit_lock_bandwidth_Hz;
    
    // Voltage time table
    dom.AC.voltage_time_table = fields.ac.voltage_time_table;
    
    // === Magnetic fields ===
    dom.B.enabled = fields.magnetic.enabled;
    dom.B.field_strength_T = fields.magnetic.field_strength_T;
    dom.B.field_gradient_T_m = fields.magnetic.field_gradient_T_m;
}

} // namespace ICARION::config
