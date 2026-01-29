// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "core/physics/spacecharge/depositCharge.h"
#include "core/log/Logger.h"
#include "core/utils/safety/numericalSafetyGuards.h"
#include <cmath>
#include <algorithm>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

// Performance thresholds
namespace {
    constexpr size_t HIGH_PERFORMANCE_THRESHOLD = 100000;    ///< Enable high-perf mode for >100k ions
    constexpr size_t ULTRA_HIGH_PERFORMANCE_THRESHOLD = 1000000; ///< Use thread-local storage for >1M ions
}

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
                                   const Grid3D& grid,
                                   const ICARION::config::IDomainGeometry* geometry)
{
    const int Nx = grid.Nx, Ny = grid.Ny, Nz = grid.Nz;
    std::vector<double> rho(Nx * Ny * Nz, 0.0);

    const double inv_dx = 1.0 / grid.dx;
    const double inv_dy = 1.0 / grid.dy;
    const double inv_dz = 1.0 / grid.dz;
    
    // Performance optimization for many ions
    const size_t num_ions = ions.size();
    const bool debug_enabled = ICARION::log::Logger::is_debug_enabled();
    if (debug_enabled && num_ions > HIGH_PERFORMANCE_THRESHOLD) {
        ICARION::log::Logger::physics()->debug(
            "[deposit_charge] High-performance mode for {} ions...", num_ions);
    }

    // Pre-calculate cell volume (constant for uniform grid)
    const double cell_volume = grid.dx * grid.dy * grid.dz;
    
    // Safety check: Prevent division by zero
    if (!ICARION::safety::is_finite_value(cell_volume) || cell_volume <= 0.0) {
        throw std::runtime_error("[deposit_charge] Invalid cell volume: " + std::to_string(cell_volume) + " (grid spacing may be zero or negative)");
    }
    
    const double inv_cell_volume = 1.0 / cell_volume;
    
    // Validate grid resolution vs. ion distribution
    if (debug_enabled && num_ions > 0) {
        // Estimate ion cloud size from first active ion's typical scale
        double grid_domain_size = std::max({grid.dx * Nx, grid.dy * Ny, grid.dz * Nz});
        double max_cell_size = std::max({grid.dx, grid.dy, grid.dz});
        
        // Warning: Grid too coarse (cell size > 10% of domain)
        if (max_cell_size > 0.1 * grid_domain_size && Nx < 32 && Ny < 32 && Nz < 32) {
            ICARION::log::Logger::physics()->debug(
                "[deposit_charge] WARNING: Coarse grid detected ({}x{}x{}), cell size = {} μm. "
                "Consider finer resolution for better accuracy.",
                Nx, Ny, Nz, max_cell_size * 1e6);
        }
        
        // Critical: Very few grid points with many ions
        size_t total_cells = static_cast<size_t>(Nx) * Ny * Nz;
        if (num_ions > total_cells * 10) {
            ICARION::log::Logger::physics()->debug(
                "[deposit_charge] WARNING: High ion density! {} ions on {} grid cells "
                "({} ions/cell avg). Grid may be under-resolved. Consider increasing resolution.",
                num_ions, total_cells, (num_ions / static_cast<double>(total_cells)));
        }
    }
    
    if (num_ions > ULTRA_HIGH_PERFORMANCE_THRESHOLD) {
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
                if (geometry && !geometry->contains(ion.pos)) continue;

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
        // Standard CIC (Cloud-In-Cell) deposition with trilinear weighting
        // Distributes charge over 8 surrounding grid nodes for smooth fields
        #pragma omp parallel for
        for (size_t ion_idx = 0; ion_idx < ions.size(); ++ion_idx) {
            const auto& ion = ions[ion_idx];
            if (!ion.active || !ion.born)
                continue;
            if (geometry && !geometry->contains(ion.pos))
                continue;

            // --- Position relative to grid origin ---
            double x_rel = (ion.pos.x - grid.origin_m.x) * inv_dx;
            double y_rel = (ion.pos.y - grid.origin_m.y) * inv_dy;
            double z_rel = (ion.pos.z - grid.origin_m.z) * inv_dz;

            // --- Lower-left grid node (base indices) ---
            int i0 = static_cast<int>(std::floor(x_rel));
            int j0 = static_cast<int>(std::floor(y_rel));
            int k0 = static_cast<int>(std::floor(z_rel));

            // --- Bounds check: Need i0, i0+1, j0, j0+1, k0, k0+1 all in bounds ---
            if (i0 < 0 || j0 < 0 || k0 < 0 || i0 >= Nx-1 || j0 >= Ny-1 || k0 >= Nz-1)
                continue; // Ion too close to boundary for CIC (would need 8 nodes)

            // --- Fractional position within cell [0,1] ---
            double fx = x_rel - i0;
            double fy = y_rel - j0;
            double fz = z_rel - k0;

            // --- Trilinear weights (sum = 1.0) ---
            // w[i][j][k] = weight for node (i0+i, j0+j, k0+k)
            double w[2][2][2];
            w[0][0][0] = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
            w[1][0][0] = fx * (1.0 - fy) * (1.0 - fz);
            w[0][1][0] = (1.0 - fx) * fy * (1.0 - fz);
            w[1][1][0] = fx * fy * (1.0 - fz);
            w[0][0][1] = (1.0 - fx) * (1.0 - fy) * fz;
            w[1][0][1] = fx * (1.0 - fy) * fz;
            w[0][1][1] = (1.0 - fx) * fy * fz;
            w[1][1][1] = fx * fy * fz;

            // --- Distribute charge to 8 surrounding nodes ---
            const double charge_density = ion.ion_charge_C * inv_cell_volume;
            
            for (int di = 0; di <= 1; ++di) {
                for (int dj = 0; dj <= 1; ++dj) {
                    for (int dk = 0; dk <= 1; ++dk) {
                        int idx = grid.index(i0 + di, j0 + dj, k0 + dk);
                        double weighted_charge = charge_density * w[di][dj][dk];
                        
                        // Atomic add for thread safety
                        #pragma omp atomic
                        rho[idx] += weighted_charge;
                    }
                }
            }
        }
    }
    if (debug_enabled) {
        double total_charge = 0.0;
        for (double r : rho) total_charge += r * (grid.dx * grid.dy * grid.dz);
        ICARION::log::Logger::physics()->debug(
            "[deposit_charge] Total deposited charge = {} C", total_charge);
    }
    return rho;
}

std::vector<double> deposit_charge(const ICARION::core::IonEnsemble& ions,
                                   const Grid3D& grid,
                                   const ICARION::config::IDomainGeometry* geometry)
{
    const int Nx = grid.Nx, Ny = grid.Ny, Nz = grid.Nz;
    std::vector<double> rho(Nx * Ny * Nz, 0.0);

    const double inv_dx = 1.0 / grid.dx;
    const double inv_dy = 1.0 / grid.dy;
    const double inv_dz = 1.0 / grid.dz;
    const double cell_volume = grid.dx * grid.dy * grid.dz;

    if (!ICARION::safety::is_finite_value(cell_volume) || cell_volume <= 0.0) {
        throw std::runtime_error("[deposit_charge] Invalid cell volume");
    }
    const double inv_cell_volume = 1.0 / cell_volume;

    const size_t num_ions = ions.size();
    const bool debug_enabled = ICARION::log::Logger::is_debug_enabled();

    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();
    const auto* charge = ions.charge_data();
    const auto* active = ions.active_data();
    const auto* born = ions.born_data();

    if (num_ions > ULTRA_HIGH_PERFORMANCE_THRESHOLD) {
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
            for (size_t ion_idx = 0; ion_idx < num_ions; ++ion_idx) {
                if (active[ion_idx] == 0 || born[ion_idx] == 0) continue;
                if (geometry) {
                    Vec3 global_pos{pos_x[ion_idx], pos_y[ion_idx], pos_z[ion_idx]};
                    if (!geometry->contains(global_pos)) continue;
                }

                int i = static_cast<int>((pos_x[ion_idx] - grid.origin_m.x) * inv_dx);
                int j = static_cast<int>((pos_y[ion_idx] - grid.origin_m.y) * inv_dy);
                int k = static_cast<int>((pos_z[ion_idx] - grid.origin_m.z) * inv_dz);

                if (i >= 0 && j >= 0 && k >= 0 && i < Nx && j < Ny && k < Nz) {
                    int idx = grid.index(i,j,k);
                    local_rho[idx] += charge[ion_idx] * inv_cell_volume;
                }
            }
        }

        #pragma omp parallel for
        for (size_t idx = 0; idx < rho.size(); ++idx) {
            for (int t = 0; t < num_threads; ++t) {
                rho[idx] += thread_rho[t][idx];
            }
        }
    } else {
        #pragma omp parallel for
        for (size_t ion_idx = 0; ion_idx < num_ions; ++ion_idx) {
            if (active[ion_idx] == 0 || born[ion_idx] == 0) continue;
            if (geometry) {
                Vec3 global_pos{pos_x[ion_idx], pos_y[ion_idx], pos_z[ion_idx]};
                if (!geometry->contains(global_pos)) continue;
            }

            double x_rel = (pos_x[ion_idx] - grid.origin_m.x) * inv_dx;
            double y_rel = (pos_y[ion_idx] - grid.origin_m.y) * inv_dy;
            double z_rel = (pos_z[ion_idx] - grid.origin_m.z) * inv_dz;

            int i0 = static_cast<int>(std::floor(x_rel));
            int j0 = static_cast<int>(std::floor(y_rel));
            int k0 = static_cast<int>(std::floor(z_rel));

            if (i0 < 0 || j0 < 0 || k0 < 0 || i0 >= Nx-1 || j0 >= Ny-1 || k0 >= Nz-1)
                continue;

            double fx = x_rel - i0;
            double fy = y_rel - j0;
            double fz = z_rel - k0;

            double w[2][2][2];
            w[0][0][0] = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
            w[1][0][0] = fx * (1.0 - fy) * (1.0 - fz);
            w[0][1][0] = (1.0 - fx) * fy * (1.0 - fz);
            w[1][1][0] = fx * fy * (1.0 - fz);
            w[0][0][1] = (1.0 - fx) * (1.0 - fy) * fz;
            w[1][0][1] = fx * (1.0 - fy) * fz;
            w[0][1][1] = (1.0 - fx) * fy * fz;
            w[1][1][1] = fx * fy * fz;

            const double charge_density = charge[ion_idx] * inv_cell_volume;

            for (int di = 0; di <= 1; ++di) {
                for (int dj = 0; dj <= 1; ++dj) {
                    for (int dk = 0; dk <= 1; ++dk) {
                        int idx = grid.index(i0 + di, j0 + dj, k0 + dk);
                        double weighted_charge = charge_density * w[di][dj][dk];
                        #pragma omp atomic
                        rho[idx] += weighted_charge;
                    }
                }
            }
        }
    }

    if (debug_enabled) {
        double total_charge = 0.0;
        for (double r : rho) total_charge += r * (grid.dx * grid.dy * grid.dz);
        ICARION::log::Logger::physics()->debug(
            "[deposit_charge] Total deposited charge = {} C", total_charge);
    }
    return rho;
}
