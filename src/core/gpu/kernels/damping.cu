// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "damping.cuh"
#include <cmath>
#include <stdexcept>

namespace icarion {
namespace gpu {

__global__ void damping_kernel(
    double* __restrict__ vx,
    double* __restrict__ vy,
    double* __restrict__ vz,
    const bool* __restrict__ active,
    int N,
    DeviceDamping damping,
    double dt
) {
    if (!damping.enabled) {
        return;
    }
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += gridDim.x * blockDim.x) {
        if (!active[i]) {
            continue;
        }
        float nu = damping.nu_per_ion ? damping.nu_per_ion[i] : damping.nu_const;
        if (nu <= 0.0f) {
            continue;
        }
        double factor = exp(-static_cast<double>(nu) * dt);
        vx[i] *= factor;
        vy[i] *= factor;
        vz[i] *= factor;
    }
}

void launch_damping_kernel(
    const IonStateGPU& ions,
    const DeviceDamping& damping,
    double dt,
    cudaStream_t stream
) {
    if (!damping.enabled || ions.count == 0) {
        return;
    }

    int N = static_cast<int>(ions.count);
    int threads = 256;
    int blocks = (N + threads - 1) / threads;
    blocks = min(blocks, 2048);

    damping_kernel<<<blocks, threads, 0, stream>>>(
        ions.vx,
        ions.vy,
        ions.vz,
        ions.active,
        N,
        damping,
        dt
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("Damping kernel launch failed: ") + cudaGetErrorString(err));
    }
}

}  // namespace gpu
}  // namespace icarion
