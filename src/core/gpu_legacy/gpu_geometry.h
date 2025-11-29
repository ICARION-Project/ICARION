// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef GPU_GEOMETRY_H
#define GPU_GEOMETRY_H

#include "core/types/Vec3.h"
#include <cuda_runtime.h>
#include <iostream>

/**
 * @brief Structure to hold molecular geometry data on GPU.
 *
 * For each species, stores atom positions and radii in device memory.
 * Used for EHSS collision calculations on GPU.
 */
struct GeometryDataGPU {
    Vec3* centers;      ///< Atom positions in molecule frame [m] (device pointer)
    double* radii;      ///< Atom radii [m] (device pointer)
    int num_atoms;      ///< Number of atoms in molecule
    double Rn;          ///< Neutral radius [m] (for impact parameter calculation)
};

/**
 * @brief Structure to hold multiple species geometries on GPU.
 *
 * Supports multi-species simulations by storing all geometries in flat arrays
 * with offset information for each species.
 */
struct MultiSpeciesGeometryGPU {
    Vec3* all_centers;          ///< All atom centers concatenated (device pointer)
    double* all_radii;          ///< All atom radii concatenated (device pointer)
    int* atom_offsets;          ///< Start index for each species in all_centers/all_radii (device pointer)
    int* num_atoms_per_species; ///< Number of atoms for each species (device pointer)
    double* Rn_per_species;     ///< Neutral radius for each species (device pointer)
    int* species_indices;       ///< Species index for each ion (device pointer)
    double* precomputed_CCS;    ///< Precomputed CCS for each species (device pointer) - PERFORMANCE OPTIMIZATION
    int num_species;            ///< Total number of species
    int total_atoms;            ///< Total atoms across all species
};

/**
 * @brief Allocate and copy geometry data to GPU.
 *
 * @param[in] h_centers Host array of atom centers
 * @param[in] h_radii Host array of atom radii
 * @param[in] num_atoms Number of atoms
 * @param[in] Rn Neutral radius
 * @param[out] d_geom GPU geometry structure (device pointers will be allocated)
 */
inline void allocate_geometry_gpu(const std::vector<Vec3>& h_centers,
                                  const std::vector<double>& h_radii,
                                  int num_atoms,
                                  double Rn,
                                  GeometryDataGPU& d_geom) {
    d_geom.num_atoms = num_atoms;
    d_geom.Rn = Rn;
    
    // Allocate device memory
    cudaMalloc(&d_geom.centers, num_atoms * sizeof(Vec3));
    cudaMalloc(&d_geom.radii, num_atoms * sizeof(double));
    
    // Copy data to device
    cudaMemcpy(d_geom.centers, h_centers.data(), num_atoms * sizeof(Vec3), 
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_geom.radii, h_radii.data(), num_atoms * sizeof(double),
               cudaMemcpyHostToDevice);
}

/**
 * @brief Free geometry data on GPU.
 *
 * @param[in,out] d_geom GPU geometry structure (device pointers will be freed)
 */
inline void free_geometry_gpu(GeometryDataGPU& d_geom) {
    if (d_geom.centers != nullptr) {
        cudaFree(d_geom.centers);
        d_geom.centers = nullptr;
    }
    if (d_geom.radii != nullptr) {
        cudaFree(d_geom.radii);
        d_geom.radii = nullptr;
    }
    d_geom.num_atoms = 0;
}

/**
 * @brief Allocate and copy multi-species geometry data to GPU.
 *
 * @param[in] geometry_map Map from species name to (centers, radii) vectors
 * @param[in] ions Vector of ions (to determine species order and indices)
 * @param[in] neutral_radius_m Radius of neutral gas molecule (e.g., He or N2)
 * @param[out] d_multi_geom GPU multi-species geometry structure
 */
inline void allocate_multi_species_geometry_gpu(
    const std::unordered_map<std::string, std::pair<std::vector<Vec3>, std::vector<double>>>& geometry_map,
    const std::vector<IonState>& ions,
    double neutral_radius_m,
    MultiSpeciesGeometryGPU& d_multi_geom,
    const std::unordered_map<std::string, double>& mobcal_ccs_map = {}) {
    
    // Build species list and ion->species mapping
    std::vector<std::string> species_list;
    std::unordered_map<std::string, int> species_to_index;
    
    for (const auto& kv : geometry_map) {
        species_to_index[kv.first] = species_list.size();
        species_list.push_back(kv.first);
    }
    
    d_multi_geom.num_species = species_list.size();
    
    // Prepare host data
    std::vector<Vec3> all_centers_host;
    std::vector<double> all_radii_host;
    std::vector<int> atom_offsets_host;
    std::vector<int> num_atoms_host;
    std::vector<double> Rn_host;
    std::vector<double> precomputed_CCS_host;  // PERFORMANCE: Precompute CCS on CPU
    std::vector<int> species_indices_host(ions.size());
    
    int current_offset = 0;
    for (const auto& sp_name : species_list) {
        const auto& geometry = geometry_map.at(sp_name);
        const auto& centers = geometry.first;
        const auto& radii = geometry.second;
        
        atom_offsets_host.push_back(current_offset);
        num_atoms_host.push_back(centers.size());
        
        // Use neutral gas radius from domain parameters
        Rn_host.push_back(neutral_radius_m);
        
        // PERFORMANCE OPTIMIZATION: Precompute CCS from geometry on CPU (once per species)
        // If a MobCal CCS value is provided for this species, prefer it and use that
        // as the precomputed CCS (explicit MobCal override).
        double ccs = 0.0;
        auto it_mob = mobcal_ccs_map.find(sp_name);
        if (it_mob != mobcal_ccs_map.end() && it_mob->second > 0.0) {
            ccs = it_mob->second; // Use MobCal override
        } else if (!centers.empty()) {
            // Compute center of mass
            Vec3 com = {0.0, 0.0, 0.0};
            for (const auto& c : centers) {
                com.x += c.x;
                com.y += c.y;
                com.z += c.z;
            }
            com.x /= centers.size();
            com.y /= centers.size();
            com.z /= centers.size();
            
            // Find maximum radial extent
            double R_max = 0.0;
            for (size_t a = 0; a < centers.size(); ++a) {
                double dx = centers[a].x - com.x;
                double dy = centers[a].y - com.y;
                double dz = centers[a].z - com.z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz) + radii[a];
                if (dist > R_max) R_max = dist;
            }
            
            // Projection approximation
            double R_eff = R_max + neutral_radius_m;
            ccs = 3.14159265358979323846 * R_eff * R_eff;
        }
    precomputed_CCS_host.push_back(ccs);
        
        // Concatenate geometry data
        all_centers_host.insert(all_centers_host.end(), centers.begin(), centers.end());
        all_radii_host.insert(all_radii_host.end(), radii.begin(), radii.end());
        
        current_offset += centers.size();
    }
    
    d_multi_geom.total_atoms = all_centers_host.size();

    // Debug: print precomputed per-species CCS values (host-side)
    if (!precomputed_CCS_host.empty()) {
        std::cout << "Precomputed CCS per species (host-side, Å^2):" << std::endl;
        for (size_t si = 0; si < species_list.size(); ++si) {
            std::cout << "  " << species_list[si] << " -> " << precomputed_CCS_host[si] * 1e20 << " Ų" << std::endl;
        }
    }
    
    // Map each ion to its species index
    for (size_t i = 0; i < ions.size(); ++i) {
        auto it = species_to_index.find(ions[i].species_id);
        species_indices_host[i] = (it != species_to_index.end()) ? it->second : 0;
    }
    
    // Allocate device memory
    cudaMalloc(&d_multi_geom.all_centers, d_multi_geom.total_atoms * sizeof(Vec3));
    cudaMalloc(&d_multi_geom.all_radii, d_multi_geom.total_atoms * sizeof(double));
    cudaMalloc(&d_multi_geom.atom_offsets, d_multi_geom.num_species * sizeof(int));
    cudaMalloc(&d_multi_geom.num_atoms_per_species, d_multi_geom.num_species * sizeof(int));
    cudaMalloc(&d_multi_geom.Rn_per_species, d_multi_geom.num_species * sizeof(double));
    cudaMalloc(&d_multi_geom.precomputed_CCS, d_multi_geom.num_species * sizeof(double));  // NEW
    cudaMalloc(&d_multi_geom.species_indices, ions.size() * sizeof(int));
    
    // Copy data to device
    cudaMemcpy(d_multi_geom.all_centers, all_centers_host.data(), 
               d_multi_geom.total_atoms * sizeof(Vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_multi_geom.all_radii, all_radii_host.data(),
               d_multi_geom.total_atoms * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_multi_geom.atom_offsets, atom_offsets_host.data(),
               d_multi_geom.num_species * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_multi_geom.num_atoms_per_species, num_atoms_host.data(),
               d_multi_geom.num_species * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_multi_geom.Rn_per_species, Rn_host.data(),
               d_multi_geom.num_species * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_multi_geom.precomputed_CCS, precomputed_CCS_host.data(),  // NEW
               d_multi_geom.num_species * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_multi_geom.species_indices, species_indices_host.data(),
               ions.size() * sizeof(int), cudaMemcpyHostToDevice);
}

/**
 * @brief Free multi-species geometry data on GPU.
 *
 * @param[in,out] d_multi_geom GPU multi-species geometry structure
 */
inline void free_multi_species_geometry_gpu(MultiSpeciesGeometryGPU& d_multi_geom) {
    if (d_multi_geom.all_centers != nullptr) cudaFree(d_multi_geom.all_centers);
    if (d_multi_geom.all_radii != nullptr) cudaFree(d_multi_geom.all_radii);
    if (d_multi_geom.atom_offsets != nullptr) cudaFree(d_multi_geom.atom_offsets);
    if (d_multi_geom.num_atoms_per_species != nullptr) cudaFree(d_multi_geom.num_atoms_per_species);
    if (d_multi_geom.Rn_per_species != nullptr) cudaFree(d_multi_geom.Rn_per_species);
    if (d_multi_geom.precomputed_CCS != nullptr) cudaFree(d_multi_geom.precomputed_CCS);  // NEW
    if (d_multi_geom.species_indices != nullptr) cudaFree(d_multi_geom.species_indices);
    
    d_multi_geom.all_centers = nullptr;
    d_multi_geom.all_radii = nullptr;
    d_multi_geom.atom_offsets = nullptr;
    d_multi_geom.num_atoms_per_species = nullptr;
    d_multi_geom.Rn_per_species = nullptr;
    d_multi_geom.precomputed_CCS = nullptr;  // NEW
    d_multi_geom.species_indices = nullptr;
    d_multi_geom.num_species = 0;
    d_multi_geom.total_atoms = 0;
}

#endif // GPU_GEOMETRY_H
