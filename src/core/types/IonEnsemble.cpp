// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "IonEnsemble.h"
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace ICARION {
namespace core {

// === Construction ===

IonEnsemble IonEnsemble::from_legacy(const std::vector<IonState>& ions) {
    IonEnsemble ensemble;
    const size_t n = ions.size();
    
    if (n == 0) return ensemble;
    
    ensemble.reserve(n);
    ensemble.resize(n);
    
    for (size_t i = 0; i < n; ++i) {
        const auto& ion = ions[i];
        
        // Hot data
        ensemble.hot_.pos_x[i] = ion.pos.x;
        ensemble.hot_.pos_y[i] = ion.pos.y;
        ensemble.hot_.pos_z[i] = ion.pos.z;
        ensemble.hot_.vel_x[i] = ion.vel.x;
        ensemble.hot_.vel_y[i] = ion.vel.y;
        ensemble.hot_.vel_z[i] = ion.vel.z;
        ensemble.hot_.mass[i] = ion.mass_kg;
        ensemble.hot_.charge[i] = ion.ion_charge_C;
        ensemble.hot_.active[i] = ion.active ? 1 : 0;
        ensemble.hot_.born[i] = ion.born ? 1 : 0;
        
        // Cold data
        ensemble.cold_.CCS[i] = ion.CCS_m2;
        ensemble.cold_.mobility[i] = ion.reduced_mobility_cm2_Vs;
        ensemble.cold_.birth_time[i] = ion.birth_time_s;
        ensemble.cold_.species_id[i] = ensemble.get_species_index(ion.species_id);
        
        // Domain cache - initialize to defaults (updated by DomainManager)
        ensemble.domain_.temperature[i] = 0.0;
        ensemble.domain_.gas_density[i] = 0.0;
        ensemble.domain_.neutral_mass[i] = 0.0;
        ensemble.domain_.domain_index[i] = ion.current_domain_index;
        
        // Output data
        ensemble.output_.history_index[i] = ion.history_index;
        ensemble.output_.t[i] = ion.t;
    }
    
    return ensemble;
}

std::vector<IonState> IonEnsemble::to_legacy() const {
    std::vector<IonState> ions;
    const size_t n = size();
    ions.reserve(n);
    
    for (size_t i = 0; i < n; ++i) {
        IonState ion;
        
        // Hot data
        ion.pos = {hot_.pos_x[i], hot_.pos_y[i], hot_.pos_z[i]};
        ion.vel = {hot_.vel_x[i], hot_.vel_y[i], hot_.vel_z[i]};
        ion.mass_kg = hot_.mass[i];
        ion.ion_charge_C = hot_.charge[i];
        ion.active = (hot_.active[i] != 0);
        ion.born = (hot_.born[i] != 0);
        
        // Cold data
        ion.CCS_m2 = cold_.CCS[i];
        ion.reduced_mobility_cm2_Vs = cold_.mobility[i];
        ion.birth_time_s = cold_.birth_time[i];
        ion.species_id = cold_.species_pool[cold_.species_id[i]];
        
        // Domain cache
        ion.current_domain_index = domain_.domain_index[i];
        
        // Output data
        ion.history_index = output_.history_index[i];
        ion.t = output_.t[i];
        ion.dt = 0.0;  // Not stored in SoA
        
        ions.push_back(ion);
    }
    
    return ions;
}

// === Size management ===

void IonEnsemble::reserve(size_t n) {
    // Hot data
    hot_.pos_x.reserve(n);
    hot_.pos_y.reserve(n);
    hot_.pos_z.reserve(n);
    hot_.vel_x.reserve(n);
    hot_.vel_y.reserve(n);
    hot_.vel_z.reserve(n);
    hot_.mass.reserve(n);
    hot_.charge.reserve(n);
    hot_.active.reserve(n);
    hot_.born.reserve(n);
    
    // Cold data
    cold_.CCS.reserve(n);
    cold_.mobility.reserve(n);
    cold_.species_id.reserve(n);
    cold_.birth_time.reserve(n);
    
    // Domain cache
    domain_.gas_density.reserve(n);
    domain_.temperature.reserve(n);
    domain_.neutral_mass.reserve(n);
    domain_.domain_index.reserve(n);
    
    // Output data
    output_.history_index.reserve(n);
    output_.t.reserve(n);
}

void IonEnsemble::resize(size_t n) {
    // Hot data
    hot_.pos_x.resize(n);
    hot_.pos_y.resize(n);
    hot_.pos_z.resize(n);
    hot_.vel_x.resize(n);
    hot_.vel_y.resize(n);
    hot_.vel_z.resize(n);
    hot_.mass.resize(n);
    hot_.charge.resize(n);
    hot_.active.resize(n);
    hot_.born.resize(n);
    
    // Cold data
    cold_.CCS.resize(n);
    cold_.mobility.resize(n);
    cold_.species_id.resize(n);
    cold_.birth_time.resize(n);
    
    // Domain cache
    domain_.gas_density.resize(n);
    domain_.temperature.resize(n);
    domain_.neutral_mass.resize(n);
    domain_.domain_index.resize(n);
    
    // Output data
    output_.history_index.resize(n);
    output_.t.resize(n);
}

void IonEnsemble::clear() {
    // Hot data
    hot_.pos_x.clear();
    hot_.pos_y.clear();
    hot_.pos_z.clear();
    hot_.vel_x.clear();
    hot_.vel_y.clear();
    hot_.vel_z.clear();
    hot_.mass.clear();
    hot_.charge.clear();
    hot_.active.clear();
    hot_.born.clear();
    
    // Cold data
    cold_.CCS.clear();
    cold_.mobility.clear();
    cold_.species_id.clear();
    cold_.birth_time.clear();
    cold_.species_pool.clear();
    cold_.species_index.clear();
    
    // Domain cache
    domain_.gas_density.clear();
    domain_.temperature.clear();
    domain_.neutral_mass.clear();
    domain_.domain_index.clear();
    
    // Output data
    output_.history_index.clear();
    output_.t.clear();
}

size_t IonEnsemble::compact_inactive() {
    const size_t n = size();
    size_t write_idx = 0;
    
    for (size_t read_idx = 0; read_idx < n; ++read_idx) {
        if (is_active(read_idx)) {
            if (write_idx != read_idx) {
                // Move active ion to write position
                hot_.pos_x[write_idx] = hot_.pos_x[read_idx];
                hot_.pos_y[write_idx] = hot_.pos_y[read_idx];
                hot_.pos_z[write_idx] = hot_.pos_z[read_idx];
                hot_.vel_x[write_idx] = hot_.vel_x[read_idx];
                hot_.vel_y[write_idx] = hot_.vel_y[read_idx];
                hot_.vel_z[write_idx] = hot_.vel_z[read_idx];
                hot_.mass[write_idx] = hot_.mass[read_idx];
                hot_.charge[write_idx] = hot_.charge[read_idx];
                hot_.active[write_idx] = hot_.active[read_idx];
                hot_.born[write_idx] = hot_.born[read_idx];
                
                cold_.CCS[write_idx] = cold_.CCS[read_idx];
                cold_.mobility[write_idx] = cold_.mobility[read_idx];
                cold_.species_id[write_idx] = cold_.species_id[read_idx];
                cold_.birth_time[write_idx] = cold_.birth_time[read_idx];
                
                domain_.gas_density[write_idx] = domain_.gas_density[read_idx];
                domain_.temperature[write_idx] = domain_.temperature[read_idx];
                domain_.neutral_mass[write_idx] = domain_.neutral_mass[read_idx];
                domain_.domain_index[write_idx] = domain_.domain_index[read_idx];
                
                output_.history_index[write_idx] = output_.history_index[read_idx];
                output_.t[write_idx] = output_.t[read_idx];
            }
            ++write_idx;
        }
    }
    
    size_t removed = n - write_idx;
    resize(write_idx);
    return removed;
}

// === View access ===

IonKinematics IonEnsemble::kinematics(size_t i) {
    return {
        hot_.pos_x.data(),
        hot_.pos_y.data(),
        hot_.pos_z.data(),
        hot_.vel_x.data(),
        hot_.vel_y.data(),
        hot_.vel_z.data(),
        hot_.mass.data(),
        hot_.charge.data(),
        hot_.active.data(),
        i
    };
}

IonCollisionData IonEnsemble::collision_data(size_t i) {
    return {
        kinematics(i),
        cold_.CCS.data(),
        domain_.temperature.data(),
        domain_.gas_density.data(),
        domain_.neutral_mass.data()
    };
}

IonReactionData IonEnsemble::reaction_data(size_t i) {
    return {
        kinematics(i),
        &cold_.species_pool,
        cold_.species_id.data()
    };
}

IonOutputData IonEnsemble::output_data(size_t i) {
    return {
        hot_.pos_x.data(),
        hot_.pos_y.data(),
        hot_.pos_z.data(),
        hot_.vel_x.data(),
        hot_.vel_y.data(),
        hot_.vel_z.data(),
        &cold_.species_pool,
        cold_.species_id.data(),
        output_.t.data(),
        domain_.domain_index.data(),
        i
    };
}

// === Bulk operations ===

void IonEnsemble::advance_positions(double dt) {
    const size_t n = size();
    
    // SIMD-friendly loop (compiler can auto-vectorize)
    for (size_t i = 0; i < n; ++i) {
        if (is_active(i)) {
            hot_.pos_x[i] += hot_.vel_x[i] * dt;
            hot_.pos_y[i] += hot_.vel_y[i] * dt;
            hot_.pos_z[i] += hot_.vel_z[i] * dt;
        }
    }
}

void IonEnsemble::advance_velocities(const double* accel_x, const double* accel_y, 
                                     const double* accel_z, double dt) {
    const size_t n = size();
    
    for (size_t i = 0; i < n; ++i) {
        if (is_active(i)) {
            hot_.vel_x[i] += accel_x[i] * dt;
            hot_.vel_y[i] += accel_y[i] * dt;
            hot_.vel_z[i] += accel_z[i] * dt;
        }
    }
}

// === Domain management ===

void IonEnsemble::update_domain_cache(size_t ion_idx, int new_domain,
                                      double temperature, double gas_density, 
                                      double neutral_mass) {
    domain_.domain_index[ion_idx] = new_domain;
    domain_.temperature[ion_idx] = temperature;
    domain_.gas_density[ion_idx] = gas_density;
    domain_.neutral_mass[ion_idx] = neutral_mass;
}

// === Species management ===

void IonEnsemble::update_species(size_t ion_idx, const std::string& new_species_id,
                                 double new_mass, double new_charge, 
                                 double new_CCS, double new_mobility) {
    hot_.mass[ion_idx] = new_mass;
    hot_.charge[ion_idx] = new_charge;
    cold_.CCS[ion_idx] = new_CCS;
    cold_.mobility[ion_idx] = new_mobility;
    cold_.species_id[ion_idx] = get_species_index(new_species_id);
}

// === Memory diagnostics ===

size_t IonEnsemble::memory_footprint() const {
    size_t hot = size() * (
        sizeof(double) * 8 +  // pos, vel, mass, charge
        sizeof(uint8_t) * 2   // active, born
    );
    
    size_t cold = size() * (
        sizeof(double) * 2 +  // CCS, mobility
        sizeof(uint32_t) +    // species_id index
        sizeof(double)        // birth_time
    ) + cold_.species_pool.size() * 32;  // Approximate string overhead
    
    size_t domain = size() * (
        sizeof(double) * 3 +  // gas_density, temperature, neutral_mass
        sizeof(int32_t)       // domain_index
    );
    
    size_t output = size() * (
        sizeof(int32_t) +     // history_index
        sizeof(double)        // t
    );
    
    return hot + cold + domain + output;
}

void IonEnsemble::print_memory_layout() const {
    std::cout << "\n=== IonEnsemble Memory Layout ===\n";
    std::cout << "Total ions: " << size() << "\n";
    std::cout << "Active ions: " << std::count(hot_.active.begin(), hot_.active.end(), 1) << "\n";
    std::cout << "\nMemory usage:\n";
    std::cout << "  Hot data:    " << std::setw(8) << (size() * 80) / 1024 << " KB  (80 bytes/ion)\n";
    std::cout << "  Cold data:   " << std::setw(8) << (size() * 24 + cold_.species_pool.size() * 32) / 1024 << " KB  (~24 bytes/ion + strings)\n";
    std::cout << "  Domain cache:" << std::setw(8) << (size() * 28) / 1024 << " KB  (28 bytes/ion)\n";
    std::cout << "  Output data: " << std::setw(8) << (size() * 12) / 1024 << " KB  (12 bytes/ion)\n";
    std::cout << "  TOTAL:       " << std::setw(8) << memory_footprint() / 1024 << " KB  (~" 
              << memory_footprint() / size() << " bytes/ion)\n";
    std::cout << "\nUnique species: " << cold_.species_pool.size() << "\n";
    for (size_t i = 0; i < cold_.species_pool.size(); ++i) {
        size_t count = std::count(cold_.species_id.begin(), cold_.species_id.end(), i);
        std::cout << "  " << cold_.species_pool[i] << ": " << count << " ions\n";
    }
    std::cout << "===============================\n\n";
}

// === Private helpers ===

uint32_t IonEnsemble::get_species_index(const std::string& species_id) {
    auto it = cold_.species_index.find(species_id);
    if (it != cold_.species_index.end()) {
        return it->second;
    }
    
    // Add new species to pool
    uint32_t idx = static_cast<uint32_t>(cold_.species_pool.size());
    cold_.species_pool.push_back(species_id);
    cold_.species_index[species_id] = idx;
    return idx;
}

} // namespace core
} // namespace ICARION
