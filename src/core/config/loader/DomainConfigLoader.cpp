// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "DomainConfigLoader.h"
#include "../conversion/UnitConverter.h"
#include <stdexcept>

namespace ICARION::config {

DomainConfig DomainConfigLoader::load(const Json::Value& json, const std::string& default_integrator) {
    DomainConfig config;
    
    // === Identification ===
    if (json.isMember("name") && json["name"].isString()) {
        config.name = json["name"].asString();
    } else {
        throw std::runtime_error("Domain missing required 'name' field");
    }
    
    if (json.isMember("instrument") && json["instrument"].isString()) {
        std::string instr_str = json["instrument"].asString();
        config.instrument = EnumMapper::parse_instrument(instr_str);
    } else {
        throw std::runtime_error("Domain '" + config.name + "' missing required 'instrument' field");
    }
    
    // === Geometry ===
    if (json.isMember("geometry")) {
        config.geometry = load_geometry(json["geometry"]);
    } else {
        throw std::runtime_error("Domain '" + config.name + "' missing required 'geometry' section");
    }
    
    // === Environment ===
    if (json.isMember("env") || json.isMember("environment")) {
        // Support both "env" (new) and "environment" (legacy)
        const Json::Value& env_json = json.isMember("env") ? json["env"] : json["environment"];
        config.environment = load_environment(env_json);
    } else {
        throw std::runtime_error("Domain '" + config.name + "' missing required 'env' section");
    }
    
    // === Fields ===
    if (json.isMember("fields")) {
        config.fields = load_fields(json["fields"]);
    }
    // Fields are optional (e.g., pure drift with only gas flow)
    
    // === Solver ===
    // Use domain-specific integrator if provided, otherwise fall back to global simulation.integrator
    std::string solver_str = default_integrator;  // Start with global default
    if (json.isMember("integrator") && json["integrator"].isString()) {
        solver_str = json["integrator"].asString();  // Override with domain-specific
    }
    config.solver = EnumMapper::parse_solver(solver_str);
    
    // === Coordinate transforms (future/multi-domain) ===
    // For now, keep identity transforms
    
    return config;
}

GeometryConfig DomainConfigLoader::load_geometry(const Json::Value& json) {
    GeometryConfig geom;
    
    // Basic dimensions
    if (json.isMember("length_m") && json["length_m"].isNumeric()) {
        geom.length_m = json["length_m"].asDouble();
    }
    
    if (json.isMember("radius_m") && json["radius_m"].isNumeric()) {
        geom.radius_m = json["radius_m"].asDouble();
    }
    
    // Orbitrap-specific
    if (json.isMember("radius_in_m") && json["radius_in_m"].isNumeric()) {
        geom.radius_in_m = json["radius_in_m"].asDouble();
    }
    
    if (json.isMember("radius_out_m") && json["radius_out_m"].isNumeric()) {
        geom.radius_out_m = json["radius_out_m"].asDouble();
    }
    
    if (json.isMember("radius_char_m") && json["radius_char_m"].isNumeric()) {
        geom.radius_char_m = json["radius_char_m"].asDouble();
    }
    
    // Compute hyperbolic boundary constants for Orbitrap
    // Hyperboloid equation: z² - r²/2 = C, where C = -0.5 * R²
    geom.orbitrap_C_in  = -0.5 * geom.radius_in_m  * geom.radius_in_m;
    geom.orbitrap_C_out = -0.5 * geom.radius_out_m * geom.radius_out_m;
    
    // Multi-domain specific
    if (json.isMember("acc_length_m") && json["acc_length_m"].isNumeric()) {
        geom.acc_length_m = json["acc_length_m"].asDouble();
    }
    
    if (json.isMember("end_aperture_m") && json["end_aperture_m"].isNumeric()) {
        geom.end_aperture_m = json["end_aperture_m"].asDouble();
    }
    
    // Origin
    if (json.isMember("origin_m")) {
        geom.origin_m = load_vec3(json["origin_m"]);
    }
    
    return geom;
}

EnvironmentConfig DomainConfigLoader::load_environment(const Json::Value& json) {
    EnvironmentConfig env;
    
    // Pressure
    if (json.isMember("pressure_Pa") && json["pressure_Pa"].isNumeric()) {
        env.pressure_Pa = json["pressure_Pa"].asDouble();
    }
    
    // Temperature
    if (json.isMember("temperature_K") && json["temperature_K"].isNumeric()) {
        env.temperature_K = json["temperature_K"].asDouble();
    }
    
    // Gas species
    if (json.isMember("gas_species") && json["gas_species"].isString()) {
        env.gas_species = json["gas_species"].asString();
    }
    
    // Gas velocity
    if (json.isMember("gas_velocity_m_s")) {
        env.gas_velocity_m_s = load_vec3(json["gas_velocity_m_s"]);
    }
    
    return env;
}

FieldsConfig DomainConfigLoader::load_fields(const Json::Value& json) {
    FieldsConfig fields;
    
    // DC fields
    if (json.isMember("DC")) {
        fields.dc = load_dc_fields(json["DC"]);
    }
    
    // RF fields
    if (json.isMember("RF")) {
        fields.rf = load_rf_fields(json["RF"]);
    }
    
    // AC fields
    if (json.isMember("AC")) {
        fields.ac = load_ac_fields(json["AC"]);
    }
    
    // Magnetic fields
    if (json.isMember("B") || json.isMember("magnetic")) {
        const Json::Value& mag_json = json.isMember("B") ? json["B"] : json["magnetic"];
        fields.magnetic = load_magnetic_fields(mag_json);
    }
    
    // Field arrays
    if (json.isMember("field_array_terms") && json["field_array_terms"].isArray()) {
        for (const auto& term_json : json["field_array_terms"]) {
            FieldsConfig::FieldArrayTerm term;
            
            // Required: file path
            if (term_json.isMember("file") && term_json["file"].isString()) {
                term.file = term_json["file"].asString();
            } else {
                throw std::runtime_error("field_array_term missing required 'file' field");
            }
            
            // Scale type (default: Constant)
            if (term_json.isMember("scale_type") && term_json["scale_type"].isString()) {
                std::string kind_str = term_json["scale_type"].asString();
                if (kind_str == "Constant") {
                    term.kind = FieldsConfig::FieldArrayTerm::ScaleKind::Constant;
                } else if (kind_str == "DC_Axial") {
                    term.kind = FieldsConfig::FieldArrayTerm::ScaleKind::DC_Axial;
                } else if (kind_str == "DC_Quad") {
                    term.kind = FieldsConfig::FieldArrayTerm::ScaleKind::DC_Quad;
                } else if (kind_str == "DC_Radial") {
                    term.kind = FieldsConfig::FieldArrayTerm::ScaleKind::DC_Radial;
                } else if (kind_str == "RF") {
                    term.kind = FieldsConfig::FieldArrayTerm::ScaleKind::RF;
                } else {
                    throw std::runtime_error("Unknown scale_type: " + kind_str);
                }
            }
            
            // Optional: scaling parameters
            if (term_json.isMember("constant_V") && term_json["constant_V"].isNumeric()) {
                term.constant = term_json["constant_V"].asDouble();
            }
            if (term_json.isMember("phase_rad") && term_json["phase_rad"].isNumeric()) {
                term.phase_rad = term_json["phase_rad"].asDouble();
            }
            if (term_json.isMember("frequency_Hz") && term_json["frequency_Hz"].isNumeric()) {
                term.frequency_Hz = term_json["frequency_Hz"].asDouble();
            }
            
            fields.field_array_terms.push_back(term);
        }
    }
    
    return fields;
}

DCFieldConfig DomainConfigLoader::load_dc_fields(const Json::Value& json) {
    DCFieldConfig dc;
    
    // Voltage specification
    if (json.isMember("axial_V") && json["axial_V"].isNumeric()) {
        dc.axial_V = json["axial_V"].asDouble();
    }
    
    if (json.isMember("quad_V") && json["quad_V"].isNumeric()) {
        dc.quad_V = json["quad_V"].asDouble();
    }
    
    if (json.isMember("radial_V") && json["radial_V"].isNumeric()) {
        dc.radial_V = json["radial_V"].asDouble();
    }
    
    // Field strength specification (alternative)
    if (json.isMember("EN_Td") && json["EN_Td"].isNumeric()) {
        dc.EN_Td = json["EN_Td"].asDouble();
        dc.EN_Vm2 = UnitConverter::townsend_to_Vm2(dc.EN_Td);
    }
    
    return dc;
}

RFFieldConfig DomainConfigLoader::load_rf_fields(const Json::Value& json) {
    RFFieldConfig rf;
    
    if (json.isMember("voltage_V") && json["voltage_V"].isNumeric()) {
        rf.voltage_V = json["voltage_V"].asDouble();
    }
    
    if (json.isMember("frequency_Hz") && json["frequency_Hz"].isNumeric()) {
        rf.frequency_Hz = json["frequency_Hz"].asDouble();
    }
    
    if (json.isMember("phase_rad") && json["phase_rad"].isNumeric()) {
        rf.phase_rad = json["phase_rad"].asDouble();
    }
    
    return rf;
}

ACFieldConfig DomainConfigLoader::load_ac_fields(const Json::Value& json) {
    ACFieldConfig ac;
    
    if (json.isMember("voltage_V") && json["voltage_V"].isNumeric()) {
        ac.voltage_V = json["voltage_V"].asDouble();
    }
    
    if (json.isMember("frequency_Hz") && json["frequency_Hz"].isNumeric()) {
        ac.frequency_Hz = json["frequency_Hz"].asDouble();
    }
    
    // Voltage sweep
    if (json.isMember("enable_voltage_sweep") && json["enable_voltage_sweep"].isBool()) {
        ac.enable_voltage_sweep = json["enable_voltage_sweep"].asBool();
        
        if (ac.enable_voltage_sweep) {
            if (json.isMember("amplitude_slope_V_s") && json["amplitude_slope_V_s"].isNumeric()) {
                ac.amplitude_slope_V_s = json["amplitude_slope_V_s"].asDouble();
            }
            if (json.isMember("start_time_s") && json["start_time_s"].isNumeric()) {
                ac.start_time_s = json["start_time_s"].asDouble();
            }
            if (json.isMember("rise_time_s") && json["rise_time_s"].isNumeric()) {
                ac.rise_time_s = json["rise_time_s"].asDouble();
            }
        }
    }
    
    // Frequency sweep
    if (json.isMember("enable_frequency_sweep") && json["enable_frequency_sweep"].isBool()) {
        ac.enable_frequency_sweep = json["enable_frequency_sweep"].asBool();
        
        if (ac.enable_frequency_sweep) {
            if (json.isMember("frequency_start_Hz") && json["frequency_start_Hz"].isNumeric()) {
                ac.frequency_start_Hz = json["frequency_start_Hz"].asDouble();
            }
            if (json.isMember("frequency_sweep_slope_Hz_s") && json["frequency_sweep_slope_Hz_s"].isNumeric()) {
                ac.frequency_sweep_slope_Hz_s = json["frequency_sweep_slope_Hz_s"].asDouble();
            }
        }
    }
    
    // LQIT phase lock
    if (json.isMember("lqit_lock_enable") && json["lqit_lock_enable"].isBool()) {
        ac.lqit_lock_enable = json["lqit_lock_enable"].asBool();
        
        if (ac.lqit_lock_enable) {
            if (json.isMember("lqit_lock_phase_rad") && json["lqit_lock_phase_rad"].isNumeric()) {
                ac.lqit_lock_phase_rad = json["lqit_lock_phase_rad"].asDouble();
            }
            if (json.isMember("lqit_lock_bandwidth_Hz") && json["lqit_lock_bandwidth_Hz"].isNumeric()) {
                ac.lqit_lock_bandwidth_Hz = json["lqit_lock_bandwidth_Hz"].asDouble();
            }
        }
    }
    
    // Voltage time table (future waveform support)
    if (json.isMember("voltage_time_table") && json["voltage_time_table"].isArray()) {
        for (const auto& entry : json["voltage_time_table"]) {
            if (entry.isArray() && entry.size() >= 2) {
                double t = entry[0].asDouble();
                double v = entry[1].asDouble();
                ac.voltage_time_table.emplace_back(t, v);
            }
        }
    }
    
    return ac;
}

MagneticFieldConfig DomainConfigLoader::load_magnetic_fields(const Json::Value& json) {
    MagneticFieldConfig mag;
    
    if (json.isMember("enabled") && json["enabled"].isBool()) {
        mag.enabled = json["enabled"].asBool();
    }
    
    if (json.isMember("field_strength_T")) {
        mag.field_strength_T = load_vec3(json["field_strength_T"]);
        if (!mag.enabled) {
            mag.enabled = true;  // Auto-enable if field specified
        }
    }
    
    if (json.isMember("field_gradient_T_m")) {
        mag.field_gradient_T_m = load_vec3(json["field_gradient_T_m"]);
    }
    
    return mag;
}

Vec3 DomainConfigLoader::load_vec3(const Json::Value& json, const Vec3& default_val) {
    if (!json.isArray() || json.size() < 3) {
        return default_val;
    }
    
    return Vec3{
        json[0].asDouble(),
        json[1].asDouble(),
        json[2].asDouble()
    };
}

} // namespace ICARION::config
