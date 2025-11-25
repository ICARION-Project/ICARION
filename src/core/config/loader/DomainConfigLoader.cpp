// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "DomainConfigLoader.h"
#include "WaveformLoader.h"
#include "../conversion/UnitConverter.h"
#include <stdexcept>
#include <iostream>

namespace ICARION::config {

DomainConfig DomainConfigLoader::load(
    const Json::Value& json, 
    const std::string& default_integrator,
    const std::map<std::string, Waveform>& global_waveforms
) {
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
        config.fields = load_fields(json["fields"], global_waveforms);
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

    // Gas mixture
    if (json.isMember("gas_mixture") && json["gas_mixture"].isArray()) {
        for (const auto& comp_json : json["gas_mixture"]) {
            GasMixtureComponent comp;
            if (comp_json.isMember("species") && comp_json["species"].isString()) {
                comp.species = comp_json["species"].asString();
            }
            if (comp_json.isMember("mole_fraction") && comp_json["mole_fraction"].isNumeric()) {
                comp.mole_fraction = comp_json["mole_fraction"].asDouble();
            }
            if (comp_json.isMember("cross_section_m2") && comp_json["cross_section_m2"].isNumeric()) {
                comp.cross_section_m2 = comp_json["cross_section_m2"].asDouble();
            }
            if (comp_json.isMember("polarizability_m3") && comp_json["polarizability_m3"].isNumeric()) {
                comp.polarizability_m3 = comp_json["polarizability_m3"].asDouble();
            }
            if (comp_json.isMember("participates_in_collisions") && comp_json["participates_in_collisions"].isBool()) {
                comp.participates_in_collisions = comp_json["participates_in_collisions"].asBool();
            }
            if (comp_json.isMember("participates_in_reactions") && comp_json["participates_in_reactions"].isBool()) {
                comp.participates_in_reactions = comp_json["participates_in_reactions"].asBool();
            }
            env.gas_mixture.push_back(comp);
        }
    }
    
    // Gas velocity
    if (json.isMember("gas_velocity_m_s")) {
        env.gas_velocity_m_s = load_vec3(json["gas_velocity_m_s"]);
    }
    
    return env;
}

FieldsConfig DomainConfigLoader::load_fields(const Json::Value& json, const std::map<std::string, Waveform>& global_waveforms) {
    FieldsConfig fields;
    
    // v1.1: Load domain-local waveform library first (for @reference resolution)
    if (json.isMember("waveforms") && json["waveforms"].isObject()) {
        try {
            fields.waveform_library = WaveformLoader::load_library(json["waveforms"]);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load waveform library: ") + e.what());
        }
    }
    
    // DC fields (pass both local and global libraries)
    if (json.isMember("DC")) {
        fields.dc = load_dc_fields(json["DC"], fields.waveform_library, global_waveforms);
    }
    
    // RF fields (pass both local and global libraries)
    if (json.isMember("RF")) {
        fields.rf = load_rf_fields(json["RF"], fields.waveform_library, global_waveforms);
    }
    
    // AC fields (pass both local and global libraries)
    if (json.isMember("AC")) {
        fields.ac = load_ac_fields(json["AC"], fields.waveform_library, global_waveforms);
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

DCFieldConfig DomainConfigLoader::load_dc_fields(const Json::Value& json, const std::map<std::string, Waveform>& local_library, const std::map<std::string, Waveform>& global_library) {
    DCFieldConfig dc;
    
    // v1.1: Voltage specification (static or waveform)
    if (json.isMember("axial_V")) {
        try {
            dc.axial_V = WaveformLoader::load_value_or_waveform(json["axial_V"], local_library, global_library);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load DC axial_V: ") + e.what());
        }
    }
    
    if (json.isMember("quad_V")) {
        try {
            dc.quad_V = WaveformLoader::load_value_or_waveform(json["quad_V"], local_library, global_library);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load DC quad_V: ") + e.what());
        }
    }
    
    if (json.isMember("radial_V")) {
        try {
            dc.radial_V = WaveformLoader::load_value_or_waveform(json["radial_V"], local_library, global_library);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load DC radial_V: ") + e.what());
        }
    }
    
    // Field strength specification (alternative) - v1.1: can also be waveform
    if (json.isMember("EN_Td")) {
        try {
            dc.EN_Td = WaveformLoader::load_value_or_waveform(json["EN_Td"], local_library, global_library);
            // If static, compute EN_Vm2
            if (dc.EN_Td.constant_value.has_value()) {
                dc.EN_Vm2 = UnitConverter::townsend_to_Vm2(dc.EN_Td.constant_value.value());
            }
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load DC EN_Td: ") + e.what());
        }
    }
    
    return dc;
}

RFFieldConfig DomainConfigLoader::load_rf_fields(const Json::Value& json, const std::map<std::string, Waveform>& local_library, const std::map<std::string, Waveform>& global_library) {
    RFFieldConfig rf;
    
    // v1.1: voltage_V (static or waveform)
    if (json.isMember("voltage_V")) {
        try {
            rf.voltage_V = WaveformLoader::load_value_or_waveform(json["voltage_V"], local_library, global_library);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load RF voltage_V: ") + e.what());
        }
    }
    
    // v1.1: frequency_Hz (static or waveform for chirps)
    if (json.isMember("frequency_Hz")) {
        try {
            rf.frequency_Hz = WaveformLoader::load_value_or_waveform(json["frequency_Hz"], local_library, global_library);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load RF frequency_Hz: ") + e.what());
        }
    }
    
    if (json.isMember("phase_rad") && json["phase_rad"].isNumeric()) {
        rf.phase_rad = json["phase_rad"].asDouble();
    }
    
    return rf;
}

ACFieldConfig DomainConfigLoader::load_ac_fields(const Json::Value& json, const std::map<std::string, Waveform>& local_library, const std::map<std::string, Waveform>& global_library) {
    ACFieldConfig ac;
    
    // v1.1: voltage_V (static or waveform)
    if (json.isMember("voltage_V")) {
        try {
            ac.voltage_V = WaveformLoader::load_value_or_waveform(json["voltage_V"], local_library, global_library);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load AC voltage_V: ") + e.what());
        }
    }
    
    // v1.1: frequency_Hz (static or waveform)
    if (json.isMember("frequency_Hz")) {
        try {
            ac.frequency_Hz = WaveformLoader::load_value_or_waveform(json["frequency_Hz"], local_library, global_library);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Failed to load AC frequency_Hz: ") + e.what());
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
