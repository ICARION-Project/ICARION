// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file GPUSpaceChargeP3M.cpp
 * @brief Host-side C++ wrapper for P³M GPU kernels
 */

#ifdef ICARION_USE_GPU

#include "GPUSpaceChargeP3M.h"
#include <cstdio>
#include <vector>

namespace icarion {
namespace gpu {

// Kernel launch wrappers are implemented as member functions in .cu file
// Implementation of main compute function
bool GPUSpaceChargeP3M::compute_space_charge_field_raw(
    const std::vector<Vec3>& positions,
    const std::vector<double>& charges,
    std::vector<Vec3>& E_field_out
) {
    if (!initialized_) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Solver not initialized\n");
        return false;
    }
    
    size_t n_ions = positions.size();
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
    
    cudaMemcpy(d_ion_positions_, positions.data(), n_ions * sizeof(Vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ion_charges_, charges.data(), n_ions * sizeof(double), cudaMemcpyHostToDevice);
    
    // =========================================================================
    // Step 2: Clear charge density grid
    // =========================================================================
    cudaMemset(d_charge_density_, 0, grid_size * sizeof(double));
    
    // =========================================================================
    // Step 3: Particle-to-Grid (P²G) - scatter charges
    // =========================================================================
    launch_p2g_cic_kernel(
        d_ion_positions_, d_ion_charges_,
        n_ions, d_charge_density_
    );
    cudaDeviceSynchronize();  // Wait for P2G to complete
    
    cudaEventRecord(p2g_end);
    
    // =========================================================================
    // Step 4: FFT forward (charge density → Fourier space)
    // =========================================================================
    cufftResult fft_result = cufftExecD2Z(fft_plan_forward_, d_charge_density_, d_potential_fft_);
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
    launch_poisson_solve_fourier_kernel(
        d_potential_fft_, d_potential_fft_,  // In-place solve
        config_.epsilon_0
    );
    cudaDeviceSynchronize();  // Wait for Poisson solve
    cudaEventRecord(poisson_end);
    
    // =========================================================================
    // Step 6: FFT inverse (potential in Fourier space → real space)
    // =========================================================================
    fft_result = cufftExecZ2D(fft_plan_inverse_, d_potential_fft_, d_potential_);
    if (fft_result != CUFFT_SUCCESS) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Inverse FFT failed: %d\n", fft_result);
        cudaEventDestroy(start); cudaEventDestroy(stop);
        cudaEventDestroy(p2g_end); cudaEventDestroy(fft_end); 
        cudaEventDestroy(poisson_end); cudaEventDestroy(g2p_end);
        return false;
    }
    
    // Normalize FFT output (cuFFT doesn't normalize)
    double scale = 1.0 / grid_size;
    launch_scale_potential_kernel(d_potential_, scale);
    cudaDeviceSynchronize();  // Wait for scaling
    
    // =========================================================================
    // Step 7: Compute E-field from potential: E = -∇φ
    // =========================================================================
    launch_compute_E_field_kernel(
        d_potential_, d_Ex_, d_Ey_, d_Ez_
    );
    cudaDeviceSynchronize();  // Wait for E-field computation
    
    // =========================================================================
    // Step 8: Grid-to-Particle (G²P) - interpolate E-field to ions
    // =========================================================================
    launch_g2p_cic_kernel(
        d_ion_positions_, n_ions,
        d_Ex_, d_Ey_, d_Ez_,
        d_ion_E_fields_
    );
    cudaDeviceSynchronize();  // Wait for G2P
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

bool GPUSpaceChargeP3M::compute_space_charge_field(
    const std::vector<IonState>& ions,
    std::vector<Vec3>& E_field_out
) {
    std::vector<Vec3> positions;
    std::vector<double> charges;
    positions.reserve(ions.size());
    charges.reserve(ions.size());
    for (const auto& ion : ions) {
        positions.push_back(ion.pos);
        charges.push_back(ion.ion_charge_C);
    }
    return compute_space_charge_field_raw(positions, charges, E_field_out);
}

bool GPUSpaceChargeP3M::compute_space_charge_field(
    const core::IonEnsemble& ions,
    std::vector<Vec3>& E_field_out
) {
    const size_t n = ions.size();
    std::vector<Vec3> positions;
    std::vector<double> charges;
    positions.reserve(n);
    charges.reserve(n);

    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();
    const auto* charge = ions.charge_data();
    for (size_t i = 0; i < n; ++i) {
        positions.push_back(Vec3{pos_x[i], pos_y[i], pos_z[i]});
        charges.push_back(charge[i]);
    }

    return compute_space_charge_field_raw(positions, charges, E_field_out);
}

} // namespace gpu
} // namespace icarion

#endif // ICARION_USE_GPU
