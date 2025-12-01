// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <array>
#include <stdexcept>
#include <string>

#include "core/utils/mathUtils.h"

namespace ICARION {
namespace core {

/**
 * @struct IonState
 * @brief Represents the physical and chemical state of a single ion
 *
 * @details
 * Complete ion state for trajectory integration and chemical tracking:
 * - Kinematic: position, velocity (6D phase space)
 * - Physical: mass, charge, mobility, collision cross-section
 * - Chemical: species identity, reaction history
 * - Simulation: activity status, domain location, birth time
 * 
 * Supports vector arithmetic for Runge-Kutta integrators (RK4, RK45).
 * CUDA-compatible via __host__ __device__ annotations.
 * Units are SI unless otherwise noted (exceptions: mobility in cm²/(V·s), CCS in m²).
 * 
 * @note Cache line aligned (64 bytes) to prevent false sharing in OpenMP parallel loops.
 *       Each IonState occupies its own cache line(s), eliminating cache coherency traffic
 *       between threads processing adjacent ions.
 */
struct alignas(64) IonState {
    // Kinematic parameters
    Vec3 pos;  ///< Position vector [m]
    Vec3 vel;  ///< Velocity vector [m/s]

    // Chemical / physical identity
    std::string species_id;               ///< Species name (e.g., "H3O+", "NO2+")
    double      mass_kg;                  ///< Ion mass [kg]
    double      reduced_mobility_cm2_Vs;  ///< Reduced mobility [cm²/(V·s)] at STP
    double      ion_charge_C;             ///< Ion charge [C] (typically ±1.602e-19 C)
    double      CCS_m2;                   ///< Collision cross-section [m²]
    double      birth_time_s;             ///< Time when ion was created [s] (0 for initial ions)
    double      death_time_s = -1.0;      ///< Time when ion was deactivated [s] (-1 = still alive)
    int         history_index = -1;       ///< Index in trajectory history buffer
    bool active = true;                   ///< false if ion was lost (boundary/detector)
    bool born = false;                    ///< true if created by reaction during simulation
    int current_domain_index = 0;         ///< Index of current instrument domain
    
    double t = 0.0;   ///< Current simulation time [s] (for adaptive integrators)
    double dt = 0.0;  ///< Individual timestep [s] per ion for adaptive RK45 solver

    /** @brief Default constructor - initializes ion to inactive state */
    __host__ __device__
    IonState() : mass_kg(0), ion_charge_C(0), reduced_mobility_cm2_Vs(0),
                 CCS_m2(0), active(true), born(true), t(0), dt(0) {}

    __host__ __device__
    IonState(const Vec3& pos_in, const Vec3& vel_in)
        : pos(pos_in), vel(vel_in),
        mass_kg(0), ion_charge_C(0), reduced_mobility_cm2_Vs(0),
        CCS_m2(0), active(true), born(true), t(0), dt(0) {}
    /**
     * @brief Add two ion states (component-wise for pos/vel).
     *
     * Preserves species identity and physical properties from the left-hand side.
     */

    inline IonState operator+(const IonState& other) const noexcept {
        IonState result(*this);
        result.pos                     = pos + other.pos;
        result.vel                     = vel + other.vel;
        return result;
    }

    /**
     * @brief Scale an ion state by a scalar factor.
     *
     * Multiplies position and velocity by @p s, preserves other properties.
     */

    inline IonState operator*(double s) const noexcept {
        IonState result(*this);
        result.pos                     = pos * s;
        result.vel                     = vel * s;
        return result;
    }

    /**
     * @brief In-place addition of another ion state (component-wise).
     *
     * Adds position and velocity of @p other to this state.
     */

    IonState& operator+=(const IonState& other) {
        this->pos += other.pos;
        this->vel += other.vel;

        return *this;
    }

    /**
     * @brief Checks for unphysical values.
     *
     * Checks ion parameters for unphysical values (e.g. negative masses and CCS values)
     * or undefined starting conditions (position and velocity). Avoids division by zero.
     *
     */
    void validate() const {
        if (!(mass_kg > 0.0) || std::isnan(mass_kg)) {
            throw std::runtime_error("IonState: mass must be positive and finite.");
        }
        if (reduced_mobility_cm2_Vs == 0.0 || std::isnan(reduced_mobility_cm2_Vs)) {
            throw std::runtime_error("IonState: ion mobility must be non-zero and finite.");
        }
        if (ion_charge_C == 0.0 || std::isnan(ion_charge_C)) {
            throw std::runtime_error("IonState: charge must be non-zero and finite.");
        }
        if (CCS_m2 <= 0.0 || std::isnan(CCS_m2)) {
            throw std::runtime_error("IonState: CCS must be > 0 and finite.");
        }
        if (species_id.empty()) {
            throw std::runtime_error("IonState: species_id is empty.");
        }
        if (std::isnan(pos.x) || std::isnan(pos.y) || std::isnan(pos.z)) {
            throw std::runtime_error("IonState: position contains NaN.");
        }
        if (std::isnan(vel.x) || std::isnan(vel.y) || std::isnan(vel.z)) {
            throw std::runtime_error("IonState: velocity contains NaN.");
        }
    }
};

}  // namespace core
}  // namespace ICARION

// Bring IonState into global namespace for backward compatibility
using ICARION::core::IonState;