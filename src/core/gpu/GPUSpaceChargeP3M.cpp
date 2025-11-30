// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file GPUSpaceChargeP3M.cpp
 * @brief Host-side C++ wrapper for P³M GPU kernels
 */

#ifdef ICARION_USE_GPU

#include "GPUSpaceChargeP3M.h"
#include <cstdio>
#include <vector>

namespace ICARION {
namespace gpu {

// Forward declare CUDA kernels (implemented in .cu file)
extern void launch_p2g_cic_kernel(
    const Vec3* d_positions,
    const double* d_charges,
    float* d_charge_density,
    int n_ions,
    Vec3 domain_min,
    Vec3 inv_spacing,
    int nx, int ny, int nz,
    cudaStream_t stream
);

extern void launch_poisson_solve_fourier_kernel(
    cufftComplex* d_rho_fft,
    cufftComplex* d_phi_fft,
    int nx, int ny, int nz_half,
    float dx, float dy, float dz,
    float epsilon_0,
    cudaStream_t stream
);

extern void launch_compute_E_field_kernel(
    const float* d_potential,
    float* d_Ex, float* d_Ey, float* d_Ez,
    int nx, int ny, int nz,
    float dx, float dy, float dz,
    cudaStream_t stream
);

extern void launch_g2p_cic_kernel(
    const Vec3* d_positions,
    Vec3* d_E_fields,
    const float* d_Ex_grid,
    const float* d_Ey_grid,
    const float* d_Ez_grid,
    int n_ions,
    Vec3 domain_min,
    Vec3 inv_spacing,
    int nx, int ny, int nz,
    cudaStream_t stream
);

// Implementation of main compute function (continued from .cu file)
bool GPUSpaceChargeP3M::compute_space_charge_field(
    const std::vector<IonState>& ions,
    std::vector<Vec3>& E_field_out
) {
    if (!initialized_) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Solver not initialized\n");
        return false;
    }
    
    size_t n_ions = ions.size();
    if (n_ions == 0) {
        E_field_out.clear();
        return true;
    }
    
    E_field_out.resize(n_ions);
    
    // Grid parameters
    int nx = config_.grid_nx;
    int ny = config_.grid_ny;
    int nz = config_.grid_nz;
    int grid_size = nx * ny * nz;
    
    Vec3 domain_size = config_.domain_max - config_.domain_min;
    Vec3 spacing{domain_size.x / nx, domain_size.y / ny, domain_size.z / nz};
    Vec3 inv_spacing{1.0 / spacing.x, 1.0 / spacing.y, 1.0 / spacing.z};
    
    // Timing
    cudaEvent_t start, stop, p2g_end, fft_end, poisson_end, g2p_end;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventCreate(&p2g_end);
    cudaEventCreate(&fft_end);
    cudaEventCreate(&poisson_end);
    cudaEventCreate(&g2p_end);
    
    cudaEventRecord(start);
    
    // =========================================================================
    // Step 1: Upload ion data to GPU
    // =========================================================================
    
    // Allocate ion buffers (if not already sized correctly)
    if (d_ion_positions_ == nullptr) {
        cudaMalloc(&d_ion_positions_, n_ions * sizeof(Vec3));
        cudaMalloc(&d_ion_charges_, n_ions * sizeof(double));
        cudaMalloc(&d_ion_E_fields_, n_ions * sizeof(Vec3));
    }
    
    // Extract positions and charges from ions
    std::vector<Vec3> positions(n_ions);
    std::vector<double> charges(n_ions);
    for (size_t i = 0; i < n_ions; ++i) {
        positions[i] = ions[i].pos;
        charges[i] = ions[i].ion_charge_C;
    }
    
    cudaMemcpy(d_ion_positions_, positions.data(), n_ions * sizeof(Vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ion_charges_, charges.data(), n_ions * sizeof(double), cudaMemcpyHostToDevice);
    
    // =========================================================================
    // Step 2: Clear charge density grid
    // =========================================================================
    cudaMemset(d_charge_density_, 0, grid_size * sizeof(float));
    
    // =========================================================================
    // Step 3: Particle-to-Grid (P²G) - scatter charges
    // =========================================================================
    launch_p2g_cic_kernel(
        d_ion_positions_, d_ion_charges_, d_charge_density_,
        n_ions, config_.domain_min, inv_spacing,
        nx, ny, nz, nullptr
    );
    cudaEventRecord(p2g_end);
    
    // =========================================================================
    // Step 4: FFT forward (charge density → Fourier space)
    // =========================================================================
    cufftResult fft_result = cufftExecR2C(fft_plan_forward_, d_charge_density_, d_potential_fft_);
    if (fft_result != CUFFT_SUCCESS) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Forward FFT failed: %d\n", fft_result);
        cudaEventDestroy(start); cudaEventDestroy(stop);
        cudaEventDestroy(p2g_end); cudaEventDestroy(fft_end); 
        cudaEventDestroy(poisson_end); cudaEventDestroy(g2p_end);
        return false;
    }
    cudaEventRecord(fft_end);
    
    // =========================================================================
    // Step 5: Solve Poisson equation in Fourier space
    // =========================================================================
    int nz_half = nz/2 + 1;
    launch_poisson_solve_fourier_kernel(
        d_potential_fft_, d_potential_fft_,  // In-place solve
        nx, ny, nz_half,
        spacing.x, spacing.y, spacing.z,
        config_.epsilon_0, nullptr
    );
    cudaEventRecord(poisson_end);
    
    // =========================================================================
    // Step 6: FFT inverse (potential in Fourier space → real space)
    // =========================================================================
    fft_result = cufftExecC2R(fft_plan_inverse_, d_potential_fft_, d_potential_);
    if (fft_result != CUFFT_SUCCESS) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Inverse FFT failed: %d\n", fft_result);
        cudaEventDestroy(start); cudaEventDestroy(stop);
        cudaEventDestroy(p2g_end); cudaEventDestroy(fft_end); 
        cudaEventDestroy(poisson_end); cudaEventDestroy(g2p_end);
        return false;
    }
    
    // Normalize FFT output (cuFFT doesn't normalize)
    float scale = 1.0f / grid_size;
    // TODO: Add kernel to scale potential or incorporate into E-field computation
    
    // =========================================================================
    // Step 7: Compute E-field from potential: E = -∇φ
    // =========================================================================
    launch_compute_E_field_kernel(
        d_potential_, d_Ex_, d_Ey_, d_Ez_,
        nx, ny, nz,
        spacing.x, spacing.y, spacing.z, nullptr
    );
    
    // =========================================================================
    // Step 8: Grid-to-Particle (G²P) - interpolate E-field to ions
    // =========================================================================
    launch_g2p_cic_kernel(
        d_ion_positions_, d_ion_E_fields_,
        d_Ex_, d_Ey_, d_Ez_,
        n_ions, config_.domain_min, inv_spacing,
        nx, ny, nz, nullptr
    );
    cudaEventRecord(g2p_end);
    
    // =========================================================================
    // Step 9: Download E-fields to CPU
    // =========================================================================
    cudaMemcpy(E_field_out.data(), d_ion_E_fields_, n_ions * sizeof(Vec3), cudaMemcpyDeviceToHost);
    
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    
    // Extract timing
    float total_ms, p2g_ms, fft_ms, poisson_ms, g2p_ms;
    cudaEventElapsedTime(&total_ms, start, stop);
    cudaEventElapsedTime(&p2g_ms, start, p2g_end);
    cudaEventElapsedTime(&fft_ms, p2g_end, fft_end);
    cudaEventElapsedTime(&poisson_ms, fft_end, poisson_end);
    cudaEventElapsedTime(&g2p_ms, poisson_end, g2p_end);
    
    // Update statistics
    stats_.total_calls++;
    stats_.total_time_ms += total_ms;
    stats_.avg_time_ms = stats_.total_time_ms / stats_.total_calls;
    stats_.last_p2g_time_ms = p2g_ms;
    stats_.last_fft_time_ms = fft_ms;
    stats_.last_poisson_time_ms = poisson_ms;
    stats_.last_g2p_time_ms = g2p_ms;
    stats_.last_n_ions = n_ions;
    
    // Cleanup events
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaEventDestroy(p2g_end);
    cudaEventDestroy(fft_end);
    cudaEventDestroy(poisson_end);
    cudaEventDestroy(g2p_end);
    
    return true;
}

} // namespace gpu
} // namespace ICARION

#endif // ICARION_USE_GPU
