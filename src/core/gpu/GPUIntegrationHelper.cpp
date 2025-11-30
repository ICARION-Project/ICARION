#include "GPUIntegrationHelper.h"
#include "integrate_rk4_batch.cuh"
#include "integrate_rk45_batch.cuh"
#include "integrate_boris_batch.cuh"
#include "FieldArrayGPU_conversion.h"
#include "fieldsolver/utils/GridFieldProvider.h"
#include "core/io/fieldArrayLoader.h"
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
    , field_array_gpu_{}
    , has_gpu_fields_(false)
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
    if (has_gpu_fields_) {
        free_field_array_gpu(field_array_gpu_);
    }
}

const FieldArray* GPUIntegrationHelper::try_extract_field_array(
    const IFieldProvider* provider
) const {
    if (!provider) {
        return nullptr;
    }
    
    // Try to downcast to GridFieldProvider
    const GridFieldProvider* grid_provider = dynamic_cast<const GridFieldProvider*>(provider);
    if (!grid_provider) {
        return nullptr;  // Not a grid-based provider (e.g., BEM, analytical)
    }
    
    // Extract FieldArray pointer (returns nullptr if using FieldSnapshot)
    return grid_provider->get_field_array();
}

bool GPUIntegrationHelper::integrate_batch_rk4(
    std::vector<ICARION::core::IonState>& ions,
    double dt,
    double t,
    const IFieldProvider* field_provider
) {
    (void)t;  // Time not currently used (fields are time-independent for now)
    
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
        
        // Try to extract field array from provider
        const FieldArray* field_array = try_extract_field_array(field_provider);
        
        // Upload fields to GPU if available
        bool use_field_interpolation = false;
        if (field_array && is_field_valid_for_gpu(*field_array)) {
            // Free previous GPU fields if any
            if (has_gpu_fields_) {
                free_field_array_gpu(field_array_gpu_);
                has_gpu_fields_ = false;
            }
            
            // Upload new fields
            try {
                upload_field_array_to_gpu(*field_array, field_array_gpu_);
                has_gpu_fields_ = true;
                use_field_interpolation = true;
            }
            catch (const std::exception& e) {
                fprintf(stderr, "Warning: Failed to upload fields to GPU: %s\n", e.what());
                fprintf(stderr, "         Falling back to zero fields.\n");
                use_field_interpolation = false;
            }
        }
        
        // Upload ions (async)
        ion_state_conversion::upload_ions(ions, ions_gpu_in_, context_.get_stream());
        
        // Initialize output buffer with input state (preserves mass, charge, species_id, domain_id)
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.mass, ions_gpu_in_.mass, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.charge, ions_gpu_in_.charge, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.species_id, ions_gpu_in_.species_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.domain_id, ions_gpu_in_.domain_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        
        // Integrate with fields or zero fields
        if (use_field_interpolation) {
            integrate_rk4_batch_with_fields(
                ions_gpu_in_,
                ions_gpu_out_,
                field_array_gpu_,
                dt,
                context_.get_stream()
            );
        } else {
            // Zero fields (free-particle motion)
            Vec3 E_zero = {0.0, 0.0, 0.0};
            Vec3 B_zero = {0.0, 0.0, 0.0};
            integrate_rk4_batch(
                ions_gpu_in_,
                ions_gpu_out_,
                E_zero,
                B_zero,
                dt,
                context_.get_stream()
            );
        }
        
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
        fprintf(stderr, "GPU RK4 integration failed: %s\n", e.what());
        enabled_ = false;
        return false;
    }
}

bool GPUIntegrationHelper::integrate_batch_rk45(
    std::vector<ICARION::core::IonState>& ions,
    double dt,
    double t,
    const IFieldProvider* field_provider,
    double atol,
    double rtol
) {
    (void)t;  // Time not currently used
    
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
            
            size_t new_capacity = N + N / 5;
            ions_gpu_in_.allocate(new_capacity);
            ions_gpu_out_.allocate(new_capacity);
            allocated_capacity_ = new_capacity;
        }
        
        ions_gpu_in_.count = N;
        ions_gpu_out_.count = N;
        
        // Try to extract field array from provider
        const FieldArray* field_array = try_extract_field_array(field_provider);
        
        // Upload fields to GPU if available
        bool use_field_interpolation = false;
        if (field_array && is_field_valid_for_gpu(*field_array)) {
            if (has_gpu_fields_) {
                free_field_array_gpu(field_array_gpu_);
                has_gpu_fields_ = false;
            }
            
            try {
                upload_field_array_to_gpu(*field_array, field_array_gpu_);
                has_gpu_fields_ = true;
                use_field_interpolation = true;
            }
            catch (const std::exception& e) {
                fprintf(stderr, "Warning: Failed to upload fields to GPU: %s\n", e.what());
                use_field_interpolation = false;
            }
        }
        
        // Upload ions (async)
        ion_state_conversion::upload_ions(ions, ions_gpu_in_, context_.get_stream());
        
        // Initialize output buffer
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.mass, ions_gpu_in_.mass, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.charge, ions_gpu_in_.charge, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.species_id, ions_gpu_in_.species_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.domain_id, ions_gpu_in_.domain_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        
        // RK45 parameters
        RK45Params params;
        params.atol = atol;
        params.rtol = rtol;
        params.safety_factor = 0.9;
        params.min_step_factor = 0.2;
        params.max_step_factor = 5.0;
        params.max_substeps = 1000;
        
        // Integrate with fields or zero fields
        if (use_field_interpolation) {
            integrate_rk45_batch_with_fields(
                ions_gpu_in_,
                ions_gpu_out_,
                field_array_gpu_,
                dt,
                params,
                context_.get_stream()
            );
        } else {
            Vec3 E_zero = {0.0, 0.0, 0.0};
            Vec3 B_zero = {0.0, 0.0, 0.0};
            integrate_rk45_batch(
                ions_gpu_in_,
                ions_gpu_out_,
                E_zero,
                B_zero,
                dt,
                params,
                context_.get_stream()
            );
        }
        
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
        fprintf(stderr, "GPU RK45 integration failed: %s\n", e.what());
        enabled_ = false;
        return false;
    }
}

bool GPUIntegrationHelper::integrate_batch_boris(
    std::vector<ICARION::core::IonState>& ions,
    double dt,
    double t,
    const IFieldProvider* field_provider
) {
    (void)t;  // Time not currently used
    
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
            
            size_t new_capacity = N + N / 5;
            ions_gpu_in_.allocate(new_capacity);
            ions_gpu_out_.allocate(new_capacity);
            allocated_capacity_ = new_capacity;
        }
        
        ions_gpu_in_.count = N;
        ions_gpu_out_.count = N;
        
        // Try to extract field array from provider
        const FieldArray* field_array = try_extract_field_array(field_provider);
        
        // Upload fields to GPU if available
        bool use_field_interpolation = false;
        if (field_array && is_field_valid_for_gpu(*field_array)) {
            if (has_gpu_fields_) {
                free_field_array_gpu(field_array_gpu_);
                has_gpu_fields_ = false;
            }
            
            try {
                upload_field_array_to_gpu(*field_array, field_array_gpu_);
                has_gpu_fields_ = true;
                use_field_interpolation = true;
            }
            catch (const std::exception& e) {
                fprintf(stderr, "Warning: Failed to upload fields to GPU: %s\n", e.what());
                use_field_interpolation = false;
            }
        }
        
        // Upload ions (async)
        ion_state_conversion::upload_ions(ions, ions_gpu_in_, context_.get_stream());
        
        // Initialize output buffer
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.mass, ions_gpu_in_.mass, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.charge, ions_gpu_in_.charge, N * sizeof(double), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.species_id, ions_gpu_in_.species_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        CUDA_CHECK(cudaMemcpyAsync(ions_gpu_out_.domain_id, ions_gpu_in_.domain_id, N * sizeof(int), 
                                  cudaMemcpyDeviceToDevice, context_.get_stream()));
        
        // Integrate with fields or zero fields
        if (use_field_interpolation) {
            integrate_boris_batch_with_fields(
                ions_gpu_in_,
                ions_gpu_out_,
                field_array_gpu_,
                dt,
                context_.get_stream()
            );
        } else {
            Vec3 E_zero = {0.0, 0.0, 0.0};
            Vec3 B_zero = {0.0, 0.0, 0.0};
            integrate_boris_batch(
                ions_gpu_in_,
                ions_gpu_out_,
                E_zero,
                B_zero,
                dt,
                context_.get_stream()
            );
        }
        
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
        fprintf(stderr, "GPU Boris integration failed: %s\n", e.what());
        enabled_ = false;
        return false;
    }
}

} // namespace gpu
} // namespace icarion
