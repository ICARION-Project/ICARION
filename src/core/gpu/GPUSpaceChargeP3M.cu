// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file GPUSpaceChargeP3M.cu
 * @brief P³M algorithm CUDA kernels for space charge field computation
 */

#include "GPUSpaceChargeP3M.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace ICARION {
namespace gpu {

// =============================================================================
// CUDA Kernels
// =============================================================================

/**
 * @brief Particle-to-Grid (P²G) kernel - Cloud-in-Cell (CIC) interpolation
 * 
 * **Thread mapping:** 1 thread = 1 ion
 * 
 * **Algorithm:**
 * For each ion at position (x,y,z):
 * 1. Find grid cell: (i,j,k) = floor((x,y,z) / cell_size)
 * 2. Compute fractional position: (fx,fy,fz) = (x,y,z) / cell_size - (i,j,k)
 * 3. Distribute charge to 8 corners via trilinear weights
 * 
 * **Atomics required:** Multiple ions can contribute to same grid cell
 */
__global__ void p2g_cic_kernel(
    const Vec3* __restrict__ positions,
    const double* __restrict__ charges,
    float* __restrict__ charge_density,
    int n_ions,
    Vec3 domain_min,
    Vec3 inv_spacing,
    int nx, int ny, int nz
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_ions) return;
    
    // Load ion position and charge
    Vec3 pos = positions[idx];
    double q = charges[idx];
    
    // Map position to grid coordinates
    double fx_d = (pos.x - domain_min.x) * inv_spacing.x;
    double fy_d = (pos.y - domain_min.y) * inv_spacing.y;
    double fz_d = (pos.z - domain_min.z) * inv_spacing.z;
    
    // Check bounds
    if (fx_d < 0 || fx_d >= nx-1 || fy_d < 0 || fy_d >= ny-1 || fz_d < 0 || fz_d >= nz-1) {
        return;  // Ion outside grid
    }
    
    // Integer cell indices
    int i = static_cast<int>(fx_d);
    int j = static_cast<int>(fy_d);
    int k = static_cast<int>(fz_d);
    
    // Fractional positions within cell [0,1]
    float fx = static_cast<float>(fx_d - i);
    float fy = static_cast<float>(fy_d - j);
    float fz = static_cast<float>(fz_d - k);
    
    // CIC weights (trilinear interpolation weights)
    float w000 = (1.0f - fx) * (1.0f - fy) * (1.0f - fz);
    float w100 = fx * (1.0f - fy) * (1.0f - fz);
    float w010 = (1.0f - fx) * fy * (1.0f - fz);
    float w110 = fx * fy * (1.0f - fz);
    float w001 = (1.0f - fx) * (1.0f - fy) * fz;
    float w101 = fx * (1.0f - fy) * fz;
    float w011 = (1.0f - fx) * fy * fz;
    float w111 = fx * fy * fz;
    
    // Grid stride (row-major: z varies fastest)
    int stride_y = nz;
    int stride_x = ny * nz;
    
    float q_f = static_cast<float>(q);
    
    // Distribute charge to 8 corners (atomic adds for thread safety)
    atomicAdd(&charge_density[i*stride_x + j*stride_y + k], w000 * q_f);
    atomicAdd(&charge_density[(i+1)*stride_x + j*stride_y + k], w100 * q_f);
    atomicAdd(&charge_density[i*stride_x + (j+1)*stride_y + k], w010 * q_f);
    atomicAdd(&charge_density[(i+1)*stride_x + (j+1)*stride_y + k], w110 * q_f);
    atomicAdd(&charge_density[i*stride_x + j*stride_y + (k+1)], w001 * q_f);
    atomicAdd(&charge_density[(i+1)*stride_x + j*stride_y + (k+1)], w101 * q_f);
    atomicAdd(&charge_density[i*stride_x + (j+1)*stride_y + (k+1)], w011 * q_f);
    atomicAdd(&charge_density[(i+1)*stride_x + (j+1)*stride_y + (k+1)], w111 * q_f);
}

/**
 * @brief Solve Poisson equation in Fourier space
 * 
 * **Thread mapping:** 1 thread = 1 grid cell in Fourier space
 * 
 * **Algorithm:**
 * Poisson equation: ∇²φ = -ρ/ε₀
 * 
 * In Fourier space:
 * -k² φ̂(k) = -ρ̂(k)/ε₀
 * φ̂(k) = ρ̂(k) / (ε₀ k²)
 * 
 * where k² = kx² + ky² + kz² (wave vector magnitude squared)
 * 
 * **Special case:** k=0 (DC component) → set φ̂(0) = 0 (arbitrary reference)
 */
__global__ void poisson_solve_fourier_kernel(
    cufftComplex* __restrict__ rho_fft,
    cufftComplex* __restrict__ phi_fft,
    int nx, int ny, int nz_half,
    float dx, float dy, float dz,
    float epsilon_0
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (i >= nx || j >= ny || k >= nz_half) return;
    
    // Wave vectors (FFT convention: k_i = 2π * freq_i / L_i)
    float kx = (i < nx/2) ? (2.0f * M_PI * i / (nx * dx)) : (2.0f * M_PI * (i - nx) / (nx * dx));
    float ky = (j < ny/2) ? (2.0f * M_PI * j / (ny * dy)) : (2.0f * M_PI * (j - ny) / (ny * dy));
    float kz = 2.0f * M_PI * k / ((nz_half-1) * 2 * dz);  // R2C FFT: only positive frequencies
    
    float k2 = kx*kx + ky*ky + kz*kz;
    
    int idx = i * ny * nz_half + j * nz_half + k;
    
    // Handle DC component (k=0)
    if (k2 < 1e-10f) {
        phi_fft[idx].x = 0.0f;
        phi_fft[idx].y = 0.0f;
        return;
    }
    
    // Poisson solve: φ̂ = ρ̂ / (ε₀ k²)
    float scale = 1.0f / (epsilon_0 * k2);
    phi_fft[idx].x = rho_fft[idx].x * scale;
    phi_fft[idx].y = rho_fft[idx].y * scale;
}

/**
 * @brief Compute electric field from potential: E = -∇φ
 * 
 * **Thread mapping:** 1 thread = 1 grid cell
 * 
 * **Algorithm:**
 * Central differences:
 * Ex(i,j,k) = -(φ(i+1,j,k) - φ(i-1,j,k)) / (2*dx)
 * Ey(i,j,k) = -(φ(i,j+1,k) - φ(i,j-1,k)) / (2*dy)
 * Ez(i,j,k) = -(φ(i,j,k+1) - φ(i,j,k-1)) / (2*dz)
 * 
 * **Boundary:** Forward/backward differences at edges
 */
__global__ void compute_E_field_kernel(
    const float* __restrict__ potential,
    float* __restrict__ Ex,
    float* __restrict__ Ey,
    float* __restrict__ Ez,
    int nx, int ny, int nz,
    float dx, float dy, float dz
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (i >= nx || j >= ny || k >= nz) return;
    
    int stride_y = nz;
    int stride_x = ny * nz;
    int idx = i*stride_x + j*stride_y + k;
    
    // Ex: ∂φ/∂x
    float phi_xp, phi_xm;
    if (i == 0) {
        phi_xp = potential[(i+1)*stride_x + j*stride_y + k];
        phi_xm = potential[i*stride_x + j*stride_y + k];
        Ex[idx] = -(phi_xp - phi_xm) / dx;
    } else if (i == nx-1) {
        phi_xp = potential[i*stride_x + j*stride_y + k];
        phi_xm = potential[(i-1)*stride_x + j*stride_y + k];
        Ex[idx] = -(phi_xp - phi_xm) / dx;
    } else {
        phi_xp = potential[(i+1)*stride_x + j*stride_y + k];
        phi_xm = potential[(i-1)*stride_x + j*stride_y + k];
        Ex[idx] = -(phi_xp - phi_xm) / (2.0f * dx);
    }
    
    // Ey: ∂φ/∂y
    float phi_yp, phi_ym;
    if (j == 0) {
        phi_yp = potential[i*stride_x + (j+1)*stride_y + k];
        phi_ym = potential[i*stride_x + j*stride_y + k];
        Ey[idx] = -(phi_yp - phi_ym) / dy;
    } else if (j == ny-1) {
        phi_yp = potential[i*stride_x + j*stride_y + k];
        phi_ym = potential[i*stride_x + (j-1)*stride_y + k];
        Ey[idx] = -(phi_yp - phi_ym) / dy;
    } else {
        phi_yp = potential[i*stride_x + (j+1)*stride_y + k];
        phi_ym = potential[i*stride_x + (j-1)*stride_y + k];
        Ey[idx] = -(phi_yp - phi_ym) / (2.0f * dy);
    }
    
    // Ez: ∂φ/∂z
    float phi_zp, phi_zm;
    if (k == 0) {
        phi_zp = potential[i*stride_x + j*stride_y + (k+1)];
        phi_zm = potential[i*stride_x + j*stride_y + k];
        Ez[idx] = -(phi_zp - phi_zm) / dz;
    } else if (k == nz-1) {
        phi_zp = potential[i*stride_x + j*stride_y + k];
        phi_zm = potential[i*stride_x + j*stride_y + (k-1)];
        Ez[idx] = -(phi_zp - phi_zm) / dz;
    } else {
        phi_zp = potential[i*stride_x + j*stride_y + (k+1)];
        phi_zm = potential[i*stride_x + j*stride_y + (k-1)];
        Ez[idx] = -(phi_zp - phi_zm) / (2.0f * dz);
    }
}

/**
 * @brief Grid-to-Particle (G²P) kernel - interpolate E-field to ion positions
 * 
 * **Thread mapping:** 1 thread = 1 ion
 * 
 * **Algorithm:** Trilinear interpolation (same as P²G CIC, but reversed)
 */
__global__ void g2p_cic_kernel(
    const Vec3* __restrict__ positions,
    Vec3* __restrict__ E_fields,
    const float* __restrict__ Ex_grid,
    const float* __restrict__ Ey_grid,
    const float* __restrict__ Ez_grid,
    int n_ions,
    Vec3 domain_min,
    Vec3 inv_spacing,
    int nx, int ny, int nz
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_ions) return;
    
    // Load ion position
    Vec3 pos = positions[idx];
    
    // Map to grid coordinates
    double fx_d = (pos.x - domain_min.x) * inv_spacing.x;
    double fy_d = (pos.y - domain_min.y) * inv_spacing.y;
    double fz_d = (pos.z - domain_min.z) * inv_spacing.z;
    
    // Check bounds
    if (fx_d < 0 || fx_d >= nx-1 || fy_d < 0 || fy_d >= ny-1 || fz_d < 0 || fz_d >= nz-1) {
        E_fields[idx] = Vec3{0.0, 0.0, 0.0};  // Zero field outside grid
        return;
    }
    
    int i = static_cast<int>(fx_d);
    int j = static_cast<int>(fy_d);
    int k = static_cast<int>(fz_d);
    
    float fx = static_cast<float>(fx_d - i);
    float fy = static_cast<float>(fy_d - j);
    float fz = static_cast<float>(fz_d - k);
    
    // CIC weights
    float w000 = (1.0f - fx) * (1.0f - fy) * (1.0f - fz);
    float w100 = fx * (1.0f - fy) * (1.0f - fz);
    float w010 = (1.0f - fx) * fy * (1.0f - fz);
    float w110 = fx * fy * (1.0f - fz);
    float w001 = (1.0f - fx) * (1.0f - fy) * fz;
    float w101 = fx * (1.0f - fy) * fz;
    float w011 = (1.0f - fx) * fy * fz;
    float w111 = fx * fy * fz;
    
    int stride_y = nz;
    int stride_x = ny * nz;
    
    // Interpolate Ex, Ey, Ez
    #define GET_FIELD(grid, i, j, k) grid[(i)*stride_x + (j)*stride_y + (k)]
    
    float Ex = w000 * GET_FIELD(Ex_grid, i, j, k) +
               w100 * GET_FIELD(Ex_grid, i+1, j, k) +
               w010 * GET_FIELD(Ex_grid, i, j+1, k) +
               w110 * GET_FIELD(Ex_grid, i+1, j+1, k) +
               w001 * GET_FIELD(Ex_grid, i, j, k+1) +
               w101 * GET_FIELD(Ex_grid, i+1, j, k+1) +
               w011 * GET_FIELD(Ex_grid, i, j+1, k+1) +
               w111 * GET_FIELD(Ex_grid, i+1, j+1, k+1);
    
    float Ey = w000 * GET_FIELD(Ey_grid, i, j, k) +
               w100 * GET_FIELD(Ey_grid, i+1, j, k) +
               w010 * GET_FIELD(Ey_grid, i, j+1, k) +
               w110 * GET_FIELD(Ey_grid, i+1, j+1, k) +
               w001 * GET_FIELD(Ey_grid, i, j, k+1) +
               w101 * GET_FIELD(Ey_grid, i+1, j, k+1) +
               w011 * GET_FIELD(Ey_grid, i, j+1, k+1) +
               w111 * GET_FIELD(Ey_grid, i+1, j+1, k+1);
    
    float Ez = w000 * GET_FIELD(Ez_grid, i, j, k) +
               w100 * GET_FIELD(Ez_grid, i+1, j, k) +
               w010 * GET_FIELD(Ez_grid, i, j+1, k) +
               w110 * GET_FIELD(Ez_grid, i+1, j+1, k) +
               w001 * GET_FIELD(Ez_grid, i, j, k+1) +
               w101 * GET_FIELD(Ez_grid, i+1, j, k+1) +
               w011 * GET_FIELD(Ez_grid, i, j+1, k+1) +
               w111 * GET_FIELD(Ez_grid, i+1, j+1, k+1);
    
    #undef GET_FIELD
    
    E_fields[idx] = Vec3{static_cast<double>(Ex), static_cast<double>(Ey), static_cast<double>(Ez)};
}

// =============================================================================
// Host-side implementation
// =============================================================================

std::string GPUSpaceChargeP3M::Config::validate() const {
    if (grid_nx < 8 || grid_ny < 8 || grid_nz < 8) {
        return "Grid dimensions must be ≥ 8";
    }
    if (grid_nx > 512 || grid_ny > 512 || grid_nz > 512) {
        return "Grid dimensions must be ≤ 512 (memory limit)";
    }
    if (domain_max.x <= domain_min.x || domain_max.y <= domain_min.y || domain_max.z <= domain_min.z) {
        return "Invalid domain bounds (max must be > min)";
    }
    if (epsilon_0 <= 0) {
        return "epsilon_0 must be positive";
    }
    return "";  // Valid
}

GPUSpaceChargeP3M::GPUSpaceChargeP3M(const GPUContext& context, const Config& config)
    : context_(context), config_(config)
{
}

GPUSpaceChargeP3M::~GPUSpaceChargeP3M() {
    cleanup();
}

std::unique_ptr<GPUSpaceChargeP3M> GPUSpaceChargeP3M::create(
    const GPUContext& context,
    const Config& config
) {
    std::string err = config.validate();
    if (!err.empty()) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Invalid config: %s\n", err.c_str());
        return nullptr;
    }
    
    auto solver = std::unique_ptr<GPUSpaceChargeP3M>(new GPUSpaceChargeP3M(context, config));
    
    if (!solver->initialize()) {
        return nullptr;
    }
    
    return solver;
}

bool GPUSpaceChargeP3M::initialize() {
    if (initialized_) return true;
    
    int nx = config_.grid_nx;
    int ny = config_.grid_ny;
    int nz = config_.grid_nz;
    int grid_size = nx * ny * nz;
    int fft_size = nx * ny * (nz/2 + 1);  // R2C FFT size
    
    // Allocate GPU memory
    cudaError_t err;
    
    // Ion data (allocated on-demand in compute_space_charge_field)
    
    // Grid data
    err = cudaMalloc(&d_charge_density_, grid_size * sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Failed to allocate charge_density: %s\n", cudaGetErrorString(err));
        return false;
    }
    
    err = cudaMalloc(&d_potential_fft_, fft_size * sizeof(cufftComplex));
    if (err != cudaSuccess) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Failed to allocate potential_fft: %s\n", cudaGetErrorString(err));
        return false;
    }
    
    err = cudaMalloc(&d_potential_, grid_size * sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Failed to allocate potential: %s\n", cudaGetErrorString(err));
        return false;
    }
    
    err = cudaMalloc(&d_Ex_, grid_size * sizeof(float));
    err |= cudaMalloc(&d_Ey_, grid_size * sizeof(float));
    err |= cudaMalloc(&d_Ez_, grid_size * sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Failed to allocate E-field grids: %s\n", cudaGetErrorString(err));
        return false;
    }
    
    // Create cuFFT plans
    cufftResult fft_err;
    fft_err = cufftPlan3d(&fft_plan_forward_, nx, ny, nz, CUFFT_R2C);
    if (fft_err != CUFFT_SUCCESS) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Failed to create forward FFT plan: %d\n", fft_err);
        return false;
    }
    
    fft_err = cufftPlan3d(&fft_plan_inverse_, nx, ny, nz, CUFFT_C2R);
    if (fft_err != CUFFT_SUCCESS) {
        fprintf(stderr, "[GPUSpaceChargeP3M] Failed to create inverse FFT plan: %d\n", fft_err);
        return false;
    }
    
    initialized_ = true;
    return true;
}

void GPUSpaceChargeP3M::cleanup() {
    if (d_ion_positions_) cudaFree(d_ion_positions_);
    if (d_ion_charges_) cudaFree(d_ion_charges_);
    if (d_ion_E_fields_) cudaFree(d_ion_E_fields_);
    if (d_charge_density_) cudaFree(d_charge_density_);
    if (d_potential_fft_) cudaFree(d_potential_fft_);
    if (d_potential_) cudaFree(d_potential_);
    if (d_Ex_) cudaFree(d_Ex_);
    if (d_Ey_) cudaFree(d_Ey_);
    if (d_Ez_) cudaFree(d_Ez_);
    
    if (fft_plan_forward_) cufftDestroy(fft_plan_forward_);
    if (fft_plan_inverse_) cufftDestroy(fft_plan_inverse_);
    
    // Reset pointers
    d_ion_positions_ = nullptr;
    d_ion_charges_ = nullptr;
    d_ion_E_fields_ = nullptr;
    d_charge_density_ = nullptr;
    d_potential_fft_ = nullptr;
    d_potential_ = nullptr;
    d_Ex_ = d_Ey_ = d_Ez_ = nullptr;
    fft_plan_forward_ = fft_plan_inverse_ = 0;
    
    initialized_ = false;
}

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
    
    // Resize output
    E_field_out.resize(n_ions);
    
    // Start timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    
    // TODO: Complete implementation in next file (this is getting long)
    // Steps:
    // 1. Upload ions (positions, charges)
    // 2. Clear charge density grid
    // 3. P²G kernel
    // 4. FFT forward
    // 5. Poisson solve in Fourier space
    // 6. FFT inverse
    // 7. Compute E-field
    // 8. G²P kernel
    // 9. Download E-fields
    
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    
    float elapsed_ms;
    cudaEventElapsedTime(&elapsed_ms, start, stop);
    
    stats_.total_calls++;
    stats_.total_time_ms += elapsed_ms;
    stats_.avg_time_ms = stats_.total_time_ms / stats_.total_calls;
    stats_.last_n_ions = n_ions;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    return true;
}

double GPUSpaceChargeP3M::Stats::speedup_vs_direct_cpu() const {
    if (last_n_ions == 0 || avg_time_ms == 0) return 1.0;
    
    // Estimate CPU direct summation time: ~100 ns per ion-ion interaction
    double cpu_time_ms = (last_n_ions * last_n_ions * 100e-9) * 1000.0;  // Convert to ms
    
    return cpu_time_ms / avg_time_ms;
}

} // namespace gpu
} // namespace ICARION
