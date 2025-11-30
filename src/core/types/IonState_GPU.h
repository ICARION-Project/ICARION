// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once
#include "core/types/Vec3.h"

/**
 * @struct IonStateGPU
 * @brief Represents the physical and chemical state of a single ion for GPU computations.
 * @details
 * Stores position, velocity, species identity (as integer), mass, charge, mobility,
 * collision cross-section, and optional metadata for reactions.
 * Designed for efficient transfer and computation on GPU devices.
 */
struct IonStateGPU {
    Vec3 pos;
    Vec3 vel;

    double mass_kg;
    double reduced_mobility_cm2_Vs;
    double ion_charge_C;
    double CCS_m2;
    double birth_time_s;

    int history_index;
    int active;
    int born;
    int current_domain_index;
    int species_id_int;   // numeric ID for GPU lookup

    double domain_neutral_mass_kg;
    double domain_temperature_K;
    double domain_particle_density_m3;
    double domain_neutral_polarizability_m3;
    Vec3 domain_gas_velocity_m_s;

    double t;
    double dt;

    __host__ __device__
    IonStateGPU() :
        mass_kg(0.0), reduced_mobility_cm2_Vs(0.0), ion_charge_C(0.0), CCS_m2(0.0),
        birth_time_s(0.0), history_index(-1), active(1), born(1),
        current_domain_index(0), species_id_int(0),
        domain_neutral_mass_kg(0.0), domain_temperature_K(0.0),
        domain_particle_density_m3(0.0), domain_neutral_polarizability_m3(0.0),
        domain_gas_velocity_m_s(0.0, 0.0, 0.0), t(0.0), dt(0.0) {}

    __host__ __device__
    IonStateGPU(const Vec3& pos_in, const Vec3& vel_in)
        : pos(pos_in), vel(vel_in),
        mass_kg(0.0), reduced_mobility_cm2_Vs(0.0), ion_charge_C(0.0), CCS_m2(0.0),
        birth_time_s(0.0), history_index(-1), active(1), born(1),
        current_domain_index(0), species_id_int(0),
        domain_neutral_mass_kg(0.0), domain_temperature_K(0.0),
        domain_particle_density_m3(0.0), domain_neutral_polarizability_m3(0.0),
        domain_gas_velocity_m_s(0.0, 0.0, 0.0), t(0.0), dt(0.0) {}

};
