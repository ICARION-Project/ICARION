/**
 * =====================================================================
 *  _______ _____            _____ ______
 * |__   __|  __ \     /\   / ____|  ____|
 *    | |  | |__) |   /  \ | |    | |__
 *    | |  |  _  /   / /\ \| |    |  __|
 *    | |  | | \ \  / ____ \ |____| |____
 *    |_|  |_|  \_\/_/    \_\_____|______|
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A C++ framework for simulating ion trajectories in electric fields
 *   under the influence of gas.
 *
 *   @file        integrate_rk4_kernel.cu
 *   @brief       GPU Kernel: RK4 Integration Step
 *
 *   @details
 *   Integrates ion trajectories using only the RK4 method on the GPU.
 *   The RK4 method uses a fixed time step `dt_s` for each integration step.
 *   This kernel computes one RK4 step for each ion in parallel on the GPU.
 *
 *
 *   @date        <2025-10-17>
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License

 *
 * =====================================================================
 */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "core/types/Vec3.h"
#include "core/types/IonState_GPU.h"
#include "core/param/GPUParams.h"
#include "core/gpu/gpu_geometry.h"
#include "core/gpu/computeAccelerationDevice.cuh"
#include "core/gpu/collision_hs_gpu.cuh"
#include "core/gpu/collision_ehss_gpu.cuh"
#include <curand_kernel.h>

// Simple derivative container: same shape wie dein CPU-Return {vel_total, acc_total}
struct Deriv { Vec3 posdot; Vec3 veldot; };

/**
 * @brief Check if position is inside the domain (local coordinates assumed)
 */
// Check whether a global position is inside the given domain (no rotation assumed)
__device__ inline bool is_inside_domain_gpu_globalpos(const Vec3& global_pos, const DomainGPU& dom) {
    // Transform global position into local domain coordinates using precomputed
    // rotation (GLOBAL -> LOCAL) stored in the domain geometry. local = R * (global - origin)
    Vec3 diff;
    diff.x = global_pos.x - dom.geom.origin_m.x;
    diff.y = global_pos.y - dom.geom.origin_m.y;
    diff.z = global_pos.z - dom.geom.origin_m.z;

    // Apply GLOBAL->LOCAL rotation rows
    Vec3 local_pos;
    local_pos.x = dom.geom.rot_row0.x * diff.x + dom.geom.rot_row0.y * diff.y + dom.geom.rot_row0.z * diff.z;
    local_pos.y = dom.geom.rot_row1.x * diff.x + dom.geom.rot_row1.y * diff.y + dom.geom.rot_row1.z * diff.z;
    local_pos.z = dom.geom.rot_row2.x * diff.x + dom.geom.rot_row2.y * diff.y + dom.geom.rot_row2.z * diff.z;

    double r = sqrt(local_pos.x * local_pos.x + local_pos.y * local_pos.y);
    const double tol = 0.0;

    if (dom.instrument != 3) { // non-Orbitrap
        return (local_pos.z >= tol && local_pos.z < dom.geom.length_m) && (r < dom.geom.radius_m);
    }

    // Simplified Orbitrap check (rotation doesn't affect radial bounds check)
    const double Rin = dom.geom.radius_in_m;
    const double Rout = dom.geom.radius_out_m;
    return (r >= Rin) && (r <= Rout);
}

// f(t, y): berechne posdot (= vel_total) und veldot (= acc_total)
__device__ inline Deriv f_eval(double t,
                               const IonStateGPU& y,
                               const GlobalParamsGPU& g,
                               const DomainGPU& dom)
{
    Deriv d;
    compute_accelerations_device(t, y, g, dom, d.posdot, d.veldot);
    return d;
}

__global__ void integrate_rk4_step_kernel(IonStateGPU* ions, int n,
                                          GlobalParamsGPU g, const DomainGPU* doms, int n_domains,
                                          double t_start, double dt_local, double t_target,
                                          curandState* rng_states,
                                          const MultiSpeciesGeometryGPU* multi_geom)
{
    // PERFORMANCE: Use grid-stride loop to process multiple ions per thread
    // This ensures GPU stays busy even with small n_ions
    int stride = blockDim.x * gridDim.x;
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += stride) {
    IonStateGPU y = ions[i];
    // determine current domain index for this ion
    int initial_dom = y.current_domain_index;
    if (initial_dom < 0) initial_dom = 0;
    if (initial_dom >= n_domains) initial_dom = 0;
    DomainGPU curDom = doms[initial_dom];
        if (!y.active) continue;
        
        // PERFORMANCE: Multi-step batching - integrate from t_start to t_target
        // This dramatically reduces kernel launch overhead
        double t_current = t_start;
        while (t_current < t_target && y.active) {
            double dt = (t_current + dt_local <= t_target) ? dt_local : (t_target - t_current);
            if (dt <= 0.0) break;
            
            // k1
            Deriv k1 = f_eval(t_current, y, g, curDom);

            // k2
            IonStateGPU y2 = y;
            y2.pos += k1.posdot * (0.5 * dt);
            y2.vel += k1.veldot * (0.5 * dt);
            Deriv k2 = f_eval(t_current + 0.5*dt, y2, g, curDom);

            // k3
            IonStateGPU y3 = y;
            y3.pos += k2.posdot * (0.5 * dt);
            y3.vel += k2.veldot * (0.5 * dt);
            Deriv k3 = f_eval(t_current + 0.5*dt, y3, g, curDom);

            // k4
            IonStateGPU y4 = y;
            y4.pos += k3.posdot * dt;
            y4.vel += k3.veldot * dt;
            Deriv k4 = f_eval(t_current + dt, y4, g, curDom);

            // RK4-Kombination
            const double w = dt / 6.0;
            const Vec3 dpos = (k1.posdot + (k2.posdot + k3.posdot)*2.0 + k4.posdot) * w;
            const Vec3 dvel = (k1.veldot + (k2.veldot + k3.veldot)*2.0 + k4.veldot) * w;

            y.pos += dpos;
            y.vel += dvel;
            t_current += dt;
            y.t = t_current;

            // Handle stochastic collisions if enabled and RNG states provided
            if (rng_states != nullptr) {
                // Build collision parameters from current domain environment
                EHSSParamsGPU collision_params;
                collision_params.n = curDom.env.particle_density_m_3;
                collision_params.dt = dt;
                collision_params.mi = y.mass_kg;
                collision_params.mn = curDom.env.neutral_mass_kg;
                collision_params.kB = 1.380649e-23;  // Boltzmann constant
                collision_params.Tn = curDom.env.temperature_K;
                collision_params.sigma_eff = y.CCS_m2;
                collision_params.ubx = curDom.env.flow_velocity_m_s.x;
                collision_params.uby = curDom.env.flow_velocity_m_s.y;
                collision_params.ubz = curDom.env.flow_velocity_m_s.z;

                if (g.collisionModel == 4) {  // CollisionModel::HSMC = 4
                    // Hard-sphere Monte Carlo (isotropic scattering)
                    handle_collision_hs_gpu(y.vel, collision_params, &rng_states[i]);
                }
                else if (g.collisionModel == 3 && multi_geom != nullptr && multi_geom->all_centers != nullptr) {  // CollisionModel::EHSS = 3
                    // EHSS with multi-species molecular geometry
                    
                    // Get species index for this ion
                    int species_idx = multi_geom->species_indices[i];
                    
                    // Bounds check on species_idx
                    if (species_idx < 0 || species_idx >= multi_geom->num_species) {
                        species_idx = 0;  // Fallback to first species if invalid
                    }
                    
                    // Lookup geometry parameters (still needed for EHSS collision handler)
                    int atom_offset = multi_geom->atom_offsets[species_idx];
                    int num_atoms = multi_geom->num_atoms_per_species[species_idx];
                    double Rn = multi_geom->Rn_per_species[species_idx];
                    
                    // Pointers to this species' geometry data
                    Vec3* species_centers = multi_geom->all_centers + atom_offset;
                    double* species_radii = multi_geom->all_radii + atom_offset;
                    
                    // PERFORMANCE OPTIMIZATION: Use precomputed CCS (computed once on CPU at initialization)
                    // This eliminates the expensive per-ion center-of-mass and R_max loops
                    double ccs_from_geometry = multi_geom->precomputed_CCS[species_idx];
                    
                    // Sample neutral velocity
                    Vec3 v_neutral = sample_neutral_velocity_gpu(collision_params, &rng_states[i]);
                    
                    // Compute collision probability using precomputed CCS
                    Vec3 vrel;
                    vrel.x = y.vel.x - v_neutral.x;
                    vrel.y = y.vel.y - v_neutral.y;
                    vrel.z = y.vel.z - v_neutral.z;
                    double v_rel = sqrt(vrel.x * vrel.x + vrel.y * vrel.y + vrel.z * vrel.z);
                    
                    // Use geometry-based CCS if available, otherwise fallback
                    double ccs_effective = (ccs_from_geometry > 0.0) ? ccs_from_geometry : collision_params.sigma_eff;
                    double P = 1.0 - exp(-collision_params.n * ccs_effective * v_rel * collision_params.dt);
                    
                    // Check if collision occurs
                    if (curand_uniform_double(&rng_states[i]) < P) {
                        // Build geometry parameters
                        EHSSGeometryParamsGPU geom_params;
                        geom_params.n = collision_params.n;
                        geom_params.dt = collision_params.dt;
                        geom_params.mi = collision_params.mi;
                        geom_params.mn = collision_params.mn;
                        geom_params.kB = collision_params.kB;
                        geom_params.Tn = collision_params.Tn;
                        geom_params.sigma_eff = ccs_effective;  // Use geometry-based CCS
                        geom_params.ubx = collision_params.ubx;
                        geom_params.uby = collision_params.uby;
                        geom_params.ubz = collision_params.ubz;
                        geom_params.Rn = Rn;
                        geom_params.num_atoms = num_atoms;
                        geom_params.max_attempts = 256;
                        
                        // Perform EHSS collision with species-specific geometry
                        y.vel = collide_ehss_gpu_geometry(y.vel, v_neutral, geom_params,
                                                          species_centers, species_radii, &rng_states[i]);
                    }
                }
            }

            // Check if ion has left the current domain; if so, try to find another domain
            if (!is_inside_domain_gpu_globalpos(y.pos, curDom)) {
                int new_dom = -1;
                for (int di = 0; di < n_domains; ++di) {
                    if (is_inside_domain_gpu_globalpos(y.pos, doms[di])) { new_dom = di; break; }
                }
                if (new_dom >= 0) {
                    y.current_domain_index = new_dom;
                    // update curDom for subsequent steps
                    curDom = doms[new_dom];
                    y.domain_neutral_mass_kg = curDom.env.neutral_mass_kg;
                    y.domain_temperature_K = curDom.env.temperature_K;
                    y.domain_particle_density_m3 = curDom.env.particle_density_m_3;
                    y.domain_gas_velocity_m_s = curDom.env.flow_velocity_m_s;
                } else {
                    y.active = 0;  // Deactivate ion that left all domains
                }
            }
        }  // End of while (t_current < t_target) loop

        ions[i] = y;
    }  // End of grid-stride loop
}

// ------------------------
// Host-Wrapper
// ------------------------
static inline const char* _cudaErrStr(cudaError_t e) { return cudaGetErrorString(e); }

extern "C" void integrate_rk4_step_gpu(IonStateGPU* d_ions, int n,
                                       const GlobalParamsGPU& g, const DomainGPU* d_doms, int n_domains,
                                       double t_local, double dt_local, double t_target,
                                       curandState* d_rng_states,
                                       const MultiSpeciesGeometryGPU* d_multi_geom)
{
    if (n == 0) return;
    
    // PERFORMANCE: Launch many blocks to saturate GPU, regardless of n_ions
    // Grid-stride loop in kernel ensures correctness even with n < total_threads
    const int threads = 256;
    const int min_blocks = 256;  // Launch 256 blocks = 65,536 threads to saturate RTX 5070 Ti
    const int blocks = std::max((n + threads - 1) / threads, min_blocks);

    integrate_rk4_step_kernel<<<blocks, threads>>>(d_ions, n, g, d_doms, n_domains, t_local, dt_local, t_target, 
                                                   d_rng_states, d_multi_geom);
    // NOTE: No sync here - caller will sync when needed (e.g., before memcpy)
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "[CUDA] integrate_rk4_step_gpu launch failed: %s\n", _cudaErrStr(err));
    }
}

// Kernel to count active ions using reduction
__global__ void count_active_kernel(IonStateGPU* ions, int n, int* d_count) {
    __shared__ int local_count[128];  // Shared memory for thread block
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
    
    // Each thread checks its ion
    local_count[tid] = (i < n && ions[i].active) ? 1 : 0;
    __syncthreads();
    
    // Reduction within block
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            local_count[tid] += local_count[tid + stride];
        }
        __syncthreads();
    }
    
    // First thread in block writes result
    if (tid == 0) {
        atomicAdd(d_count, local_count[0]);
    }
}

extern "C" int count_active_ions_gpu(IonStateGPU* d_ions, int n) {
    int* d_count;
    cudaMalloc(&d_count, sizeof(int));
    cudaMemset(d_count, 0, sizeof(int));
    
    const int threads = 256;
    const int blocks = (n + threads - 1) / threads;
    
    count_active_kernel<<<blocks, threads>>>(d_ions, n, d_count);
    cudaDeviceSynchronize();
    
    int h_count = 0;
    cudaMemcpy(&h_count, d_count, sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_count);
    
    return h_count;
}
