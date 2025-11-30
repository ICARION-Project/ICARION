// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "IonState_GPU.h"
#include "core/gpu/core/GPUContext.h"
#include <cuda_runtime.h>
#include <vector>

namespace icarion {
namespace gpu {

void IonStateGPU::allocate(size_t N) {
    if (is_allocated()) {
        free();
    }
    
    count = N;
    
    if (N == 0) {
        return;
    }
    
    // Allocate device memory for all arrays
    CUDA_CHECK(cudaMalloc(&x, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&y, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&z, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&vx, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&vy, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&vz, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&mass, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&charge, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&species_id, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&domain_id, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&active, N * sizeof(bool)));
}

void IonStateGPU::free() {
    if (!is_allocated()) {
        return;
    }
    
    cudaFree(x);
    cudaFree(y);
    cudaFree(z);
    cudaFree(vx);
    cudaFree(vy);
    cudaFree(vz);
    cudaFree(mass);
    cudaFree(charge);
    cudaFree(species_id);
    cudaFree(domain_id);
    cudaFree(active);
    
    x = y = z = nullptr;
    vx = vy = vz = nullptr;
    mass = charge = nullptr;
    species_id = domain_id = nullptr;
    active = nullptr;
    count = 0;
}

namespace ion_state_conversion {

void upload_ions(
    const std::vector<IonState>& ions_cpu,
    IonStateGPU& ions_gpu,
    cudaStream_t stream
) {
    size_t N = ions_cpu.size();
    if (N == 0) {
        return;
    }
    
    if (!ions_gpu.is_allocated() || ions_gpu.count != N) {
        throw std::runtime_error("IonStateGPU not properly allocated");
    }
    
    // Create temporary host arrays (SoA format)
    std::vector<double> x_h(N), y_h(N), z_h(N);
    std::vector<double> vx_h(N), vy_h(N), vz_h(N);
    std::vector<double> mass_h(N), charge_h(N);
    std::vector<int> species_id_h(N), domain_id_h(N);
    std::vector<uint8_t> active_h(N);  // Use uint8_t instead of bool for .data()
    
    // Convert AoS -> SoA
    for (size_t i = 0; i < N; i++) {
        const auto& ion = ions_cpu[i];
        x_h[i] = ion.pos.x;
        y_h[i] = ion.pos.y;
        z_h[i] = ion.pos.z;
        vx_h[i] = ion.vel.x;
        vy_h[i] = ion.vel.y;
        vz_h[i] = ion.vel.z;
        mass_h[i] = ion.mass_kg;
        charge_h[i] = ion.ion_charge_C;
        // species_id is string in IonState - we'll use history_index for now
        species_id_h[i] = ion.history_index;
        domain_id_h[i] = ion.current_domain_index;
        active_h[i] = ion.active;
    }
    
    // Upload to device
    if (stream == 0) {
        // Synchronous
        CUDA_CHECK(cudaMemcpy(ions_gpu.x, x_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.y, y_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.z, z_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.vx, vx_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.vy, vy_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.vz, vz_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.mass, mass_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.charge, charge_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.species_id, species_id_h.data(), N * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.domain_id, domain_id_h.data(), N * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.active, active_h.data(), N * sizeof(bool), cudaMemcpyHostToDevice));
    } else {
        // Asynchronous
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.x, x_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.y, y_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.z, z_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.vx, vx_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.vy, vy_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.vz, vz_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.mass, mass_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.charge, charge_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.species_id, species_id_h.data(), N * sizeof(int), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.domain_id, domain_id_h.data(), N * sizeof(int), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.active, active_h.data(), N * sizeof(bool), cudaMemcpyHostToDevice, stream));
    }
}

void download_ions(
    const IonStateGPU& ions_gpu,
    std::vector<IonState>& ions_cpu,
    cudaStream_t stream
) {
    size_t N = ions_gpu.count;
    if (N == 0) {
        ions_cpu.clear();
        return;
    }
    
    if (!ions_gpu.is_allocated()) {
        throw std::runtime_error("IonStateGPU not allocated");
    }
    
    // Create temporary host arrays (SoA format)
    std::vector<double> x_h(N), y_h(N), z_h(N);
    std::vector<double> vx_h(N), vy_h(N), vz_h(N);
    std::vector<double> mass_h(N), charge_h(N);
    std::vector<int> species_id_h(N), domain_id_h(N);
    std::vector<uint8_t> active_h(N);  // Use uint8_t instead of bool for .data()
    
    // Download from device
    if (stream == 0) {
        // Synchronous
        CUDA_CHECK(cudaMemcpy(x_h.data(), ions_gpu.x, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(y_h.data(), ions_gpu.y, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(z_h.data(), ions_gpu.z, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(vx_h.data(), ions_gpu.vx, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(vy_h.data(), ions_gpu.vy, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(vz_h.data(), ions_gpu.vz, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(mass_h.data(), ions_gpu.mass, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(charge_h.data(), ions_gpu.charge, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(species_id_h.data(), ions_gpu.species_id, N * sizeof(int), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(domain_id_h.data(), ions_gpu.domain_id, N * sizeof(int), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(active_h.data(), ions_gpu.active, N * sizeof(bool), cudaMemcpyDeviceToHost));
    } else {
        // Asynchronous
        CUDA_CHECK(cudaMemcpyAsync(x_h.data(), ions_gpu.x, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(y_h.data(), ions_gpu.y, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(z_h.data(), ions_gpu.z, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(vx_h.data(), ions_gpu.vx, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(vy_h.data(), ions_gpu.vy, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(vz_h.data(), ions_gpu.vz, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(mass_h.data(), ions_gpu.mass, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(charge_h.data(), ions_gpu.charge, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(species_id_h.data(), ions_gpu.species_id, N * sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(domain_id_h.data(), ions_gpu.domain_id, N * sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(active_h.data(), ions_gpu.active, N * sizeof(bool), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
    }
    
    // Convert SoA -> AoS
    ions_cpu.resize(N);
    for (size_t i = 0; i < N; i++) {
        auto& ion = ions_cpu[i];
        ion.pos = {x_h[i], y_h[i], z_h[i]};
        ion.vel = {vx_h[i], vy_h[i], vz_h[i]};
        ion.mass_kg = mass_h[i];
        ion.ion_charge_C = charge_h[i];
        ion.history_index = species_id_h[i];
        ion.current_domain_index = domain_id_h[i];
        ion.active = active_h[i];
    }
}

void upload_positions_velocities(
    const std::vector<IonState>& ions_cpu,
    IonStateGPU& ions_gpu,
    cudaStream_t stream
) {
    size_t N = ions_cpu.size();
    if (N == 0) {
        return;
    }
    
    if (!ions_gpu.is_allocated() || ions_gpu.count != N) {
        throw std::runtime_error("IonStateGPU not properly allocated");
    }
    
    // Create temporary host arrays
    std::vector<double> x_h(N), y_h(N), z_h(N);
    std::vector<double> vx_h(N), vy_h(N), vz_h(N);
    
    // Convert AoS -> SoA (positions and velocities only)
    for (size_t i = 0; i < N; i++) {
        const auto& ion = ions_cpu[i];
        x_h[i] = ion.pos.x;
        y_h[i] = ion.pos.y;
        z_h[i] = ion.pos.z;
        vx_h[i] = ion.vel.x;
        vy_h[i] = ion.vel.y;
        vz_h[i] = ion.vel.z;
    }
    
    // Upload to device
    if (stream == 0) {
        CUDA_CHECK(cudaMemcpy(ions_gpu.x, x_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.y, y_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.z, z_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.vx, vx_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.vy, vy_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ions_gpu.vz, vz_h.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    } else {
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.x, x_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.y, y_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.z, z_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.vx, vx_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.vy, vy_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu.vz, vz_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, stream));
    }
}

void download_positions_velocities(
    const IonStateGPU& ions_gpu,
    std::vector<IonState>& ions_cpu,
    cudaStream_t stream
) {
    size_t N = ions_gpu.count;
    if (N == 0 || ions_cpu.size() != N) {
        throw std::runtime_error("Size mismatch in download_positions_velocities");
    }
    
    if (!ions_gpu.is_allocated()) {
        throw std::runtime_error("IonStateGPU not allocated");
    }
    
    // Create temporary host arrays
    std::vector<double> x_h(N), y_h(N), z_h(N);
    std::vector<double> vx_h(N), vy_h(N), vz_h(N);
    
    // Download from device
    if (stream == 0) {
        CUDA_CHECK(cudaMemcpy(x_h.data(), ions_gpu.x, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(y_h.data(), ions_gpu.y, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(z_h.data(), ions_gpu.z, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(vx_h.data(), ions_gpu.vx, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(vy_h.data(), ions_gpu.vy, N * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(vz_h.data(), ions_gpu.vz, N * sizeof(double), cudaMemcpyDeviceToHost));
    } else {
        CUDA_CHECK(cudaMemcpyAsync(x_h.data(), ions_gpu.x, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(y_h.data(), ions_gpu.y, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(z_h.data(), ions_gpu.z, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(vx_h.data(), ions_gpu.vx, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(vy_h.data(), ions_gpu.vy, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(vz_h.data(), ions_gpu.vz, N * sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
    }
    
    // Convert SoA -> AoS (positions and velocities only)
    for (size_t i = 0; i < N; i++) {
        auto& ion = ions_cpu[i];
        ion.pos = {x_h[i], y_h[i], z_h[i]};
        ion.vel = {vx_h[i], vy_h[i], vz_h[i]};
    }
}

} // namespace ion_state_conversion

} // namespace gpu
} // namespace icarion
