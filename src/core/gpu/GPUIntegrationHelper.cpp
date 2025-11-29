#include "GPUIntegrationHelper.h"
#include "integrate_rk4_batch.cuh"
#include <chrono>
#include <stdexcept>

namespace icarion {
namespace gpu {

std::unique_ptr<GPUIntegrationHelper> GPUIntegrationHelper::create(
    const GPUContext& context,
    size_t threshold
) {
    if (!context.is_valid()) {
        return nullptr;
    }
    
    auto memory_pool = std::make_unique<GPUMemoryPool>(context);
    
    return std::unique_ptr<GPUIntegrationHelper>(
        new GPUIntegrationHelper(context, std::move(memory_pool), threshold)
    );
}

GPUIntegrationHelper::GPUIntegrationHelper(
    const GPUContext& context,
    std::unique_ptr<GPUMemoryPool> memory_pool,
    size_t threshold
)
    : context_(context)
    , memory_pool_(std::move(memory_pool))
    , allocated_capacity_(0)
    , threshold_(threshold)
    , enabled_(true)
    , stats_{}
{
}

GPUIntegrationHelper::~GPUIntegrationHelper() {
    if (ions_gpu_in_.is_allocated()) {
        ions_gpu_in_.free();
    }
    if (ions_gpu_out_.is_allocated()) {
        ions_gpu_out_.free();
    }
}

bool GPUIntegrationHelper::integrate_batch_rk4(
    std::vector<ICARION::core::IonState>& ions,
    double dt,
    double t
) {
    if (!enabled_ || ions.size() < threshold_) {
        return false;  // Caller should use CPU path
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        size_t N = ions.size();
        
        // Allocate/reallocate GPU buffers if needed
        if (N > allocated_capacity_) {
            if (ions_gpu_in_.is_allocated()) {
                ions_gpu_in_.free();
            }
            if (ions_gpu_out_.is_allocated()) {
                ions_gpu_out_.free();
            }
            
            // Add 20% headroom
            size_t new_capacity = N + N / 5;
            ions_gpu_in_.allocate(new_capacity);
            ions_gpu_out_.allocate(new_capacity);
            allocated_capacity_ = new_capacity;
        }
        
        ions_gpu_in_.count = N;
        ions_gpu_out_.count = N;
        
        // TODO Phase 3: Extract fields from ForceRegistry
        // For now, use zero fields (testing phase)
        Vec3 E_field = {0.0, 0.0, 0.0};
        Vec3 B_field = {0.0, 0.0, 0.0};
        
        // Upload (async)
        ion_state_conversion::upload_ions(ions, ions_gpu_in_, context_.get_stream());
        
        // Initialize output buffer with input state (preserves mass, charge, species_id, domain_id)
        // Kernel will only update pos/vel/active
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.mass, ions_gpu_in_.mass, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.charge, ions_gpu_in_.charge, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.species_id, ions_gpu_in_.species_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.domain_id, ions_gpu_in_.domain_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        
        // Integrate (async)
        integrate_rk4_batch(
            ions_gpu_in_,
            ions_gpu_out_,
            E_field,
            B_field,
            dt,
            context_.get_stream()
        );
        
        // Download (async + sync)
        ion_state_conversion::download_ions(ions_gpu_out_, ions, context_.get_stream());
        
        context_.synchronize();
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        
        stats_.gpu_integrations++;
        stats_.total_ions_gpu += N;
        stats_.total_time_ms += elapsed_ms;
        
        return true;
    }
    catch (const std::exception& e) {
        fprintf(stderr, "GPU integration failed: %s\n", e.what());
        enabled_ = false;
        return false;
    }
}

} // namespace gpu
} // namespace icarion
