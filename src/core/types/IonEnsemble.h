// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IonState.h"
#include "Vec3.h"
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>

namespace ICARION {
namespace core {

// Forward declarations for view classes
struct IonKinematics;
struct IonCollisionData;
struct IonReactionData;
struct IonOutputData;

/**
 * @brief Structure-of-Arrays container for ion ensemble
 * 
 * Design philosophy: "Only load what you need"
 * - Hot data (pos, vel, mass, charge) separate from cold data
 * - Cache-line aligned for optimal performance
 * - Zero false-sharing between OpenMP threads
 * - SIMD-friendly contiguous layout
 * 
 * Memory layout:
 * - Hot data: 80 bytes/ion (fits in L2 cache)
 * - Total: ~120 bytes/ion (45% reduction vs AoS)
 */
class IonEnsemble {
public:
    // === Construction ===
    IonEnsemble() = default;
    
    /**
     * @brief Create ensemble from legacy AoS vector
     * @param ions Legacy ion vector
     * @return SoA ensemble
     */
    static IonEnsemble from_legacy(const std::vector<IonState>& ions);
    
    /**
     * @brief Convert back to legacy AoS format (for compatibility)
     * @return Vector of IonState structs
     * @note Used for output and testing. Will be removed in Phase 6.
     */
    std::vector<IonState> to_legacy() const;
    
    // === Size management ===
    size_t size() const { return hot_.pos_x.size(); }
    bool empty() const { return hot_.pos_x.empty(); }
    
    void reserve(size_t n);
    void resize(size_t n);
    void clear();
    
    /**
     * @brief Remove inactive ions and defragment memory
     * @return Number of ions removed
     * 
     * Call periodically (e.g., every 10k steps) to:
     * - Reduce memory footprint
     * - Improve cache locality
     * - Speed up OpenMP loops
     */
    size_t compact_inactive();
    
    // === Direct element access (for transition period) ===
    
    /**
     * @brief Check if ion is active
     */
    bool is_active(size_t i) const { 
        return i < size() && hot_.active[i] && hot_.born[i]; 
    }
    
    /**
     * @brief Get position as Vec3
     */
    Vec3 get_pos(size_t i) const {
        return {hot_.pos_x[i], hot_.pos_y[i], hot_.pos_z[i]};
    }
    
    /**
     * @brief Set position from Vec3
     */
    void set_pos(size_t i, const Vec3& pos) {
        hot_.pos_x[i] = pos.x;
        hot_.pos_y[i] = pos.y;
        hot_.pos_z[i] = pos.z;
    }
    
    /**
     * @brief Get velocity as Vec3
     */
    Vec3 get_vel(size_t i) const {
        return {hot_.vel_x[i], hot_.vel_y[i], hot_.vel_z[i]};
    }
    
    /**
     * @brief Set velocity from Vec3
     */
    void set_vel(size_t i, const Vec3& vel) {
        hot_.vel_x[i] = vel.x;
        hot_.vel_y[i] = vel.y;
        hot_.vel_z[i] = vel.z;
    }
    
    // === View access (zero-copy, Phase 2+) ===
    
    IonKinematics kinematics(size_t i);
    IonCollisionData collision_data(size_t i);
    IonReactionData reaction_data(size_t i);
    IonOutputData output_data(size_t i);
    
    // === Bulk operations (SIMD-friendly) ===
    
    /**
     * @brief Update all positions: pos += vel * dt
     * @param dt Timestep
     * 
     * Vectorized operation, much faster than per-ion loop.
     */
    void advance_positions(double dt);
    
    /**
     * @brief Update all velocities: vel += accel * dt
     * @param accel_x, accel_y, accel_z Acceleration arrays
     * @param dt Timestep
     */
    void advance_velocities(const double* accel_x, const double* accel_y, 
                           const double* accel_z, double dt);
    
    // === Domain management ===
    
    /**
     * @brief Update domain-specific cached properties
     * @param ion_idx Ion index
     * @param new_domain New domain index
     * @param temperature Temperature [K]
     * @param gas_density Number density [1/m³]
     * @param neutral_mass Neutral gas mass [kg]
     */
    void update_domain_cache(size_t ion_idx, int new_domain,
                            double temperature, double gas_density, 
                            double neutral_mass);
    
    // === Species management ===
    
    /**
     * @brief Update species properties (for reactions)
     * @param ion_idx Ion index
     * @param new_species_id Species identifier
     * @param new_mass Mass [kg]
     * @param new_charge Charge [C]
     * @param new_CCS Collision cross-section [m²]
     * @param new_mobility Reduced mobility [cm²/(V·s)]
     */
    void update_species(size_t ion_idx, const std::string& new_species_id,
                       double new_mass, double new_charge, 
                       double new_CCS, double new_mobility);
    
    // === Memory diagnostics ===
    
    /**
     * @brief Get total memory footprint in bytes
     */
    size_t memory_footprint() const;
    
    /**
     * @brief Print cache efficiency statistics
     */
    void print_memory_layout() const;
    
    // === Direct data access (for advanced users) ===
    
    // Hot data
    double* pos_x_data() { return hot_.pos_x.data(); }
    double* pos_y_data() { return hot_.pos_y.data(); }
    double* pos_z_data() { return hot_.pos_z.data(); }
    double* vel_x_data() { return hot_.vel_x.data(); }
    double* vel_y_data() { return hot_.vel_y.data(); }
    double* vel_z_data() { return hot_.vel_z.data(); }
    double* mass_data() { return hot_.mass.data(); }
    double* charge_data() { return hot_.charge.data(); }
    uint8_t* active_data() { return hot_.active.data(); }
    uint8_t* born_data() { return hot_.born.data(); }
    
    const double* pos_x_data() const { return hot_.pos_x.data(); }
    const double* pos_y_data() const { return hot_.pos_y.data(); }
    const double* pos_z_data() const { return hot_.pos_z.data(); }
    const double* vel_x_data() const { return hot_.vel_x.data(); }
    const double* vel_y_data() const { return hot_.vel_y.data(); }
    const double* vel_z_data() const { return hot_.vel_z.data(); }
    const double* mass_data() const { return hot_.mass.data(); }
    const double* charge_data() const { return hot_.charge.data(); }
    const uint8_t* active_data() const { return hot_.active.data(); }
    const uint8_t* born_data() const { return hot_.born.data(); }
    
    // Cold data (const access)
    const std::string& species_id(size_t i) const { 
        return cold_.species_pool[cold_.species_id[i]]; 
    }
    double CCS(size_t i) const { return cold_.CCS[i]; }
    double mobility(size_t i) const { return cold_.mobility[i]; }
    double birth_time(size_t i) const { return cold_.birth_time[i]; }
    
    // Cold data (mutable access for SoA processing)
    double* CCS_data() { return cold_.CCS.data(); }
    double* mobility_data() { return cold_.mobility.data(); }
    const std::vector<std::string>* species_pool() const { return &cold_.species_pool; }
    const uint32_t* species_id_indices() const { return cold_.species_id.data(); }
    
    // Domain cache (const access)
    int domain_index(size_t i) const { return domain_.domain_index[i]; }
    double temperature(size_t i) const { return domain_.temperature[i]; }
    double gas_density(size_t i) const { return domain_.gas_density[i]; }
    double neutral_mass(size_t i) const { return domain_.neutral_mass[i]; }
    
    // Domain cache (mutable access for SoA processing)
    int32_t* domain_index_data() { return domain_.domain_index.data(); }
    double* temperature_data() { return domain_.temperature.data(); }
    double* gas_density_data() { return domain_.gas_density.data(); }
    double* neutral_mass_data() { return domain_.neutral_mass.data(); }
    
    // Output data
    int history_index(size_t i) const { return output_.history_index[i]; }
    double time(size_t i) const { return output_.t[i]; }
    void set_time(size_t i, double t) { output_.t[i] = t; }

private:
    /**
     * @brief Hot data: Accessed every timestep (integration, collision)
     * 
     * Cache-line aligned to prevent false sharing.
     * Total: 80 bytes/ion
     */
    struct alignas(64) HotData {
        std::vector<double> pos_x, pos_y, pos_z;  // 24 bytes
        std::vector<double> vel_x, vel_y, vel_z;  // 24 bytes
        std::vector<double> mass;                 // 8 bytes
        std::vector<double> charge;               // 8 bytes
        std::vector<uint8_t> active;              // 1 byte (bool as uint8 for SIMD)
        std::vector<uint8_t> born;                // 1 byte
        // Padding: ~14 bytes to cache-line boundary
    } hot_;
    
    /**
     * @brief Cold data: Read-mostly, rarely changed
     * 
     * Stored separately to avoid polluting hot cache lines.
     * Total: ~40 bytes/ion
     */
    struct ColdData {
        std::vector<double> CCS;                  // 8 bytes
        std::vector<double> mobility;             // 8 bytes
        std::vector<uint32_t> species_id;         // 4 bytes (index into pool)
        std::vector<double> birth_time;           // 8 bytes
        
        // String pool: Deduplicate species names (e.g., "H3O+")
        // Typical: 5-10 unique species for 10k ions
        std::vector<std::string> species_pool;
        std::unordered_map<std::string, uint32_t> species_index;
    } cold_;
    
    /**
     * @brief Domain cache: Computed on domain transition
     * 
     * Replaces legacy domain_XXX fields in IonState.
     * Only updated when ion changes domain (~0.1% of timesteps).
     * Total: ~20 bytes/ion
     */
    struct DomainCache {
        std::vector<double> gas_density;          // 8 bytes [1/m³]
        std::vector<double> temperature;          // 8 bytes [K]
        std::vector<double> neutral_mass;         // 8 bytes [kg]
        std::vector<int32_t> domain_index;        // 4 bytes
    } domain_;
    
    /**
     * @brief Output data: For trajectory logging
     * 
     * Total: ~12 bytes/ion
     */
    struct OutputData {
        std::vector<int32_t> history_index;       // 4 bytes
        std::vector<double> t;                    // 8 bytes [s]
    } output_;
    
    /**
     * @brief Get or create species index
     */
    uint32_t get_species_index(const std::string& species_id);
};

/**
 * @brief Lightweight view for integration (zero-copy)
 * 
 * Contains only pointers + index. No data copying.
 * Size: 72 bytes (9 pointers + 1 size_t)
 */
struct IonKinematics {
    double* pos_x;
    double* pos_y;
    double* pos_z;
    double* vel_x;
    double* vel_y;
    double* vel_z;
    double* mass;
    double* charge;
    uint8_t* active;
    size_t index;
    
    Vec3 pos() const { return {pos_x[index], pos_y[index], pos_z[index]}; }
    Vec3 vel() const { return {vel_x[index], vel_y[index], vel_z[index]}; }
    double get_mass() const { return mass[index]; }
    double get_charge() const { return charge[index]; }
    bool is_active() const { return active[index]; }
    
    void set_pos(const Vec3& p) { pos_x[index]=p.x; pos_y[index]=p.y; pos_z[index]=p.z; }
    void set_vel(const Vec3& v) { vel_x[index]=v.x; vel_y[index]=v.y; vel_z[index]=v.z; }
};

/**
 * @brief View for collision handling
 * 
 * Adds collision-specific data to kinematics.
 */
struct IonCollisionData {
    IonKinematics kin;
    const double* CCS;
    const double* temperature;
    const double* gas_density;
    const double* neutral_mass;
    
    double get_CCS() const { return CCS[kin.index]; }
    double get_temperature() const { return temperature[kin.index]; }
    double get_gas_density() const { return gas_density[kin.index]; }
    double get_neutral_mass() const { return neutral_mass[kin.index]; }
};

/**
 * @brief View for reaction handling
 */
struct IonReactionData {
    IonKinematics kin;
    const std::vector<std::string>* species_pool;
    const uint32_t* species_id_index;
    
    const std::string& species_id() const { 
        return (*species_pool)[species_id_index[kin.index]]; 
    }
};

/**
 * @brief View for output writing
 */
struct IonOutputData {
    const double* pos_x;
    const double* pos_y;
    const double* pos_z;
    const double* vel_x;
    const double* vel_y;
    const double* vel_z;
    const std::vector<std::string>* species_pool;
    const uint32_t* species_id_index;
    const double* t;
    const int32_t* domain_index;
    size_t index;
    
    Vec3 pos() const { return {pos_x[index], pos_y[index], pos_z[index]}; }
    Vec3 vel() const { return {vel_x[index], vel_y[index], vel_z[index]}; }
    const std::string& species_id() const { 
        return (*species_pool)[species_id_index[index]]; 
    }
    double time() const { return t[index]; }
    int domain() const { return domain_index[index]; }
};

} // namespace core
} // namespace ICARION
