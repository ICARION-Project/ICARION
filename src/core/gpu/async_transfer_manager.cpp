/**
 * GPU Async Transfer Manager Implementation
 * 
 * Date: 2025-01-23
 */

#ifdef USE_CUDA

#include "async_transfer_manager.h"
#include <iostream>
#include <chrono>
#include <algorithm>

namespace ICARION {
namespace gpu {

void AsyncTransferManager::initialize(size_t max_ions) {
    if (initialized_) return;
    
    std::cout << "[AsyncTransfer] Initializing for " << max_ions << " ions..." << std::endl;
    
    // Create CUDA streams
    for (int i = 0; i < NUM_STREAMS; ++i) {
        cudaStreamCreate(&streams_[i]);
    }
    
    // Allocate pinned memory for faster transfers
    max_buffer_size_ = max_ions;
    
    // Allocate pinned host memory
    cudaHostAlloc(&pinned_cpu_buffer_, max_ions * sizeof(IonState), cudaHostAllocDefault);
    cudaHostAlloc(&pinned_gpu_buffer_, max_ions * sizeof(IonStateGPU), cudaHostAllocDefault);
    
    // Create event pool for synchronization
    event_pool_.reserve(10);
    for (int i = 0; i < 10; ++i) {
        cudaEvent_t event;
        cudaEventCreate(&event);
        event_pool_.push_back(event);
    }
    
    initialized_ = true;
    
    std::cout << "[AsyncTransfer] Initialized with " << NUM_STREAMS 
              << " streams and pinned memory buffers" << std::endl;
}

cudaEvent_t AsyncTransferManager::startHostToDeviceAsync(const std::vector<IonState>& cpu_ions,
                                                         IonStateGPU* gpu_ions, 
                                                         size_t count,
                                                         int stream_id) {
    if (!initialized_) {
        initialize(count);
    }
    
    total_transfers_++;
    async_transfers_++;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Convert CPU ions to GPU format in pinned memory
    convertToGPU(cpu_ions, pinned_gpu_buffer_, count);
    
    // Get event for this transfer
    cudaEvent_t transfer_event = getEvent();
    
    // Async copy to device
    stream_id = std::min(stream_id, NUM_STREAMS - 1);
    cudaMemcpyAsync(gpu_ions, pinned_gpu_buffer_, 
                    count * sizeof(IonStateGPU), 
                    cudaMemcpyHostToDevice, 
                    streams_[stream_id]);
    
    // Record event
    cudaEventRecord(transfer_event, streams_[stream_id]);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    total_transfer_time_ms_ += duration.count() / 1000.0;
    
    std::cout << "[AsyncTransfer] H2D started: " << count << " ions, stream " 
              << stream_id << ", " << duration.count() << " μs" << std::endl;
    
    return transfer_event;
}

cudaEvent_t AsyncTransferManager::startDeviceToHostAsync(const IonStateGPU* gpu_ions,
                                                         std::vector<IonState>& cpu_ions,
                                                         size_t count,
                                                         int stream_id) {
    if (!initialized_) {
        initialize(count);
    }
    
    total_transfers_++;
    async_transfers_++;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Get event for this transfer
    cudaEvent_t transfer_event = getEvent();
    
    // Async copy from device
    stream_id = std::min(stream_id, NUM_STREAMS - 1);
    cudaMemcpyAsync(pinned_gpu_buffer_, gpu_ions, 
                    count * sizeof(IonStateGPU), 
                    cudaMemcpyDeviceToHost, 
                    streams_[stream_id]);
    
    // Record event
    cudaEventRecord(transfer_event, streams_[stream_id]);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    total_transfer_time_ms_ += duration.count() / 1000.0;
    
    std::cout << "[AsyncTransfer] D2H started: " << count << " ions, stream " 
              << stream_id << ", " << duration.count() << " μs" << std::endl;
    
    return transfer_event;
}

void AsyncTransferManager::waitForTransfer(cudaEvent_t event) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Wait for transfer completion
    cudaEventSynchronize(event);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    
    std::cout << "[AsyncTransfer] Transfer completed, wait time: " 
              << duration.count() << " μs" << std::endl;
    
    // Return event to pool
    returnEvent(event);
}

cudaStream_t AsyncTransferManager::getComputeStream(int stream_id) {
    if (!initialized_) {
        initialize(1000); // Default size
    }
    
    stream_id = std::min(stream_id, NUM_STREAMS - 1);
    return streams_[stream_id];
}

void AsyncTransferManager::synchronizeAll() {
    if (!initialized_) return;
    
    for (int i = 0; i < NUM_STREAMS; ++i) {
        cudaStreamSynchronize(streams_[i]);
    }
}

void AsyncTransferManager::convertToGPU(const std::vector<IonState>& cpu_ions, 
                                        IonStateGPU* gpu_buffer, 
                                        size_t count) {
    for (size_t i = 0; i < count && i < cpu_ions.size(); ++i) {
        const auto& cpu_ion = cpu_ions[i];
        auto& gpu_ion = gpu_buffer[i];
        
        gpu_ion.pos = cpu_ion.pos;
        gpu_ion.vel = cpu_ion.vel;
        gpu_ion.mass_kg = cpu_ion.mass_kg;
        gpu_ion.reduced_mobility_cm2_Vs = cpu_ion.reduced_mobility_cm2_Vs;
        gpu_ion.ion_charge_C = cpu_ion.ion_charge_C;
        gpu_ion.CCS_m2 = cpu_ion.CCS_m2;
        gpu_ion.birth_time_s = cpu_ion.birth_time_s;
        gpu_ion.history_index = cpu_ion.history_index;
        gpu_ion.active = cpu_ion.active ? 1 : 0;
        gpu_ion.born = cpu_ion.born ? 1 : 0;
        gpu_ion.current_domain_index = cpu_ion.current_domain_index;
        // Note: species_id_int would need mapping from string to int
        gpu_ion.species_id_int = 0; // Simplified for now
    }
}

void AsyncTransferManager::convertToCPU(const IonStateGPU* gpu_buffer, 
                                        std::vector<IonState>& cpu_ions, 
                                        size_t count) {
    if (cpu_ions.size() < count) {
        cpu_ions.resize(count);
    }
    
    for (size_t i = 0; i < count; ++i) {
        const auto& gpu_ion = gpu_buffer[i];
        auto& cpu_ion = cpu_ions[i];
        
        cpu_ion.pos = gpu_ion.pos;
        cpu_ion.vel = gpu_ion.vel;
        cpu_ion.mass_kg = gpu_ion.mass_kg;
        cpu_ion.reduced_mobility_cm2_Vs = gpu_ion.reduced_mobility_cm2_Vs;
        cpu_ion.ion_charge_C = gpu_ion.ion_charge_C;
        cpu_ion.CCS_m2 = gpu_ion.CCS_m2;
        cpu_ion.birth_time_s = gpu_ion.birth_time_s;
        cpu_ion.history_index = gpu_ion.history_index;
        cpu_ion.active = (gpu_ion.active != 0);
        cpu_ion.born = (gpu_ion.born != 0);
        cpu_ion.current_domain_index = gpu_ion.current_domain_index;
        // Note: would need reverse mapping from int to string for species_id
        cpu_ion.species_id = "unknown"; // Simplified for now
    }
}

cudaEvent_t AsyncTransferManager::getEvent() {
    if (event_pool_.empty()) {
        cudaEvent_t event;
        cudaEventCreate(&event);
        return event;
    }
    
    cudaEvent_t event = event_pool_.back();
    event_pool_.pop_back();
    return event;
}

void AsyncTransferManager::returnEvent(cudaEvent_t event) {
    event_pool_.push_back(event);
}

void AsyncTransferManager::printStats() const {
    std::cout << "\n=== Async Transfer Statistics ===" << std::endl;
    std::cout << "Total transfers: " << total_transfers_ << std::endl;
    std::cout << "Async transfers: " << async_transfers_ << std::endl;
    std::cout << "Total transfer time: " << total_transfer_time_ms_ << " ms" << std::endl;
    
    if (total_transfers_ > 0) {
        std::cout << "Average transfer time: " 
                  << (total_transfer_time_ms_ / total_transfers_) << " ms" << std::endl;
    }
    
    std::cout << "Active streams: " << NUM_STREAMS << std::endl;
    std::cout << "Event pool size: " << event_pool_.size() << std::endl;
    std::cout << "===================================" << std::endl;
}

AsyncTransferManager::~AsyncTransferManager() {
    if (!initialized_) return;
    
    std::cout << "[AsyncTransfer] Cleaning up..." << std::endl;
    
    // Synchronize all streams before cleanup
    synchronizeAll();
    
    // Destroy streams
    for (int i = 0; i < NUM_STREAMS; ++i) {
        cudaStreamDestroy(streams_[i]);
    }
    
    // Free pinned memory
    if (pinned_cpu_buffer_) {
        cudaFreeHost(pinned_cpu_buffer_);
    }
    if (pinned_gpu_buffer_) {
        cudaFreeHost(pinned_gpu_buffer_);
    }
    
    // Destroy events
    for (auto event : event_pool_) {
        cudaEventDestroy(event);
    }
    
    std::cout << "[AsyncTransfer] Cleanup complete" << std::endl;
}

} // namespace gpu
} // namespace ICARION

#endif // USE_CUDA