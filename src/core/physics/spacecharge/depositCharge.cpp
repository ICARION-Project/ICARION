/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        depositCharge.h
 *   @brief       Utility functions for charge deposition onto a 3D grid in ICARION.
 *
 * @details
 * Provides:
 * - Function to deposit ion charges onto a 3D grid.
 * - Support for different deposition methods (NGP, CIC, TSC).
 * - Charge density computation from ion positions and charges.
 *
 *   @date        2025-10-18
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#include "core/physics/spacecharge/depositCharge.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Deposit ion charges onto a 3D grid to compute charge density ρ [C/m³].
 *
 * Supports simple nearest-grid-point (NGP) deposition.
 * Optionally extendable to cloud-in-cell (CIC) or TSC weighting.
 *
 * @param ions  Vector of all active ions.
 * @param grid  Simulation grid (geometry & spacing).
 * @return std::vector<double> Charge density array (size Nx*Ny*Nz).
 */
std::vector<double> deposit_charge(const std::vector<IonState>& ions,
                                   const Grid3D& grid)
{
    const int Nx = grid.Nx, Ny = grid.Ny, Nz = grid.Nz;
    std::vector<double> rho(Nx * Ny * Nz, 0.0);

    const double inv_dx = 1.0 / grid.dx;
    const double inv_dy = 1.0 / grid.dy;
    const double inv_dz = 1.0 / grid.dz;
    
    // Performance optimization for many ions
    const size_t num_ions = ions.size();
    if (num_ions > 100000) {
        std::cout << "[deposit_charge] High-performance mode for " << num_ions << " ions..." << std::endl;
    }

    // Pre-calculate cell volume (constant for uniform grid)
    const double cell_volume = grid.dx * grid.dy * grid.dz;
    const double inv_cell_volume = 1.0 / cell_volume;
    
    if (num_ions > 1000000) {
        // Ultra-high performance mode for millions of ions
        // Use thread-local storage to reduce atomic contention
#ifdef _OPENMP
        const int num_threads = omp_get_max_threads();
#else
        const int num_threads = 1;
#endif
        std::vector<std::vector<double>> thread_rho(num_threads, std::vector<double>(rho.size(), 0.0));
        
        #pragma omp parallel
        {
#ifdef _OPENMP
            const int thread_id = omp_get_thread_num();
#else
            const int thread_id = 0;
#endif
            auto& local_rho = thread_rho[thread_id];
            
            #pragma omp for nowait
            for (size_t ion_idx = 0; ion_idx < ions.size(); ++ion_idx) {
                const auto& ion = ions[ion_idx];
                if (!ion.active || !ion.born) continue;

                // Fast grid mapping
                int i = static_cast<int>((ion.pos.x - grid.origin_m.x) * inv_dx);
                int j = static_cast<int>((ion.pos.y - grid.origin_m.y) * inv_dy);
                int k = static_cast<int>((ion.pos.z - grid.origin_m.z) * inv_dz);

                if (i >= 0 && j >= 0 && k >= 0 && i < Nx && j < Ny && k < Nz) {
                    int idx = grid.index(i,j,k);
                    local_rho[idx] += ion.ion_charge_C * inv_cell_volume;
                }
            }
        }
        
        // Combine thread-local results
        #pragma omp parallel for
        for (size_t idx = 0; idx < rho.size(); ++idx) {
            for (int t = 0; t < num_threads; ++t) {
                rho[idx] += thread_rho[t][idx];
            }
        }
        
    } else {
        // Standard parallel deposition
        #pragma omp parallel for
        for (size_t ion_idx = 0; ion_idx < ions.size(); ++ion_idx) {
            const auto& ion = ions[ion_idx];
            if (!ion.active || !ion.born)
                continue;

            // --- lokale Koordinaten relativ zum Ursprung des Grids ---
            double x_rel = (ion.pos.x - grid.origin_m.x);
            double y_rel = (ion.pos.y - grid.origin_m.y);
            double z_rel = (ion.pos.z - grid.origin_m.z);

            // --- Index im Gitter ---
            int i = static_cast<int>(std::floor(x_rel * inv_dx));
            int j = static_cast<int>(std::floor(y_rel * inv_dy));
            int k = static_cast<int>(std::floor(z_rel * inv_dz));

            if (i < 0 || j < 0 || k < 0 || i >= Nx || j >= Ny || k >= Nz)
                continue; // Ion außerhalb des Grids

            int idx = grid.index(i,j,k);

            // --- Atomic add for thread safety ---
            #pragma omp atomic
            rho[idx] += ion.ion_charge_C * inv_cell_volume;
        }
    }
    double total_charge = 0.0;
    for (double r : rho) total_charge += r * (grid.dx * grid.dy * grid.dz);
    std::cout << "[deposit_charge] Total deposited charge = "
            << total_charge << " C" << std::endl;
    return rho;
}
