/**
 * GPU Async Transfer Manager
 * 
 * Implements overlapped computation and data movement using CUDA streams.
 * Expected 1.1-1.3x speedup from Phase 2 Section 9.2 optimization.
 * 
 * Key benefits:
 * - Overlap CPU-GPU data transfers with GPU computation
 * - Use pinned memory for faster transfers
 * - Multiple streams for pipeline parallelism
 * 
 * Date: 2025-01-23
 */

#pragma once

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <vector>
#include <memory>
#include "utils/IonState_GPU.h"
#include "core/types/IonState.h"

namespace ICARION {
namespace gpu {

/**
 * @brief GPU Async Transfer Manager
 * 
 * Manages asynchronous transfers between CPU and GPU using CUDA streams
 * to overlap computation with data movement.
 */
class AsyncTransferManager {
public:
    static AsyncTransferManager& getInstance() {
        static AsyncTransferManager instance;
        return instance;
    }

    /**
     * @brief Initialize async transfer system
     * @param max_ions Maximum number of ions to support
     */
    void initialize(size_t max_ions);

    /**
     * @brief Start async transfer of ions from CPU to GPU
     * @param cpu_ions Source CPU ion data
     * @param gpu_ions Destination GPU ion data
     * @param count Number of ions to transfer
     * @param stream_id Stream to use for transfer (0-2)
     * @return Transfer handle for synchronization
     */
    cudaEvent_t startHostToDeviceAsync(const std::vector<IonState>& cpu_ions,
                                       IonStateGPU* gpu_ions, 
                                       size_t count,
                                       int stream_id = 0);

    /**
     * @brief Start async transfer of ions from GPU to CPU
     * @param gpu_ions Source GPU ion data
     * @param cpu_ions Destination CPU ion data
     * @param count Number of ions to transfer
     * @param stream_id Stream to use for transfer (0-2)
     * @return Transfer handle for synchronization
     */
    cudaEvent_t startDeviceToHostAsync(const IonStateGPU* gpu_ions,
                                       std::vector<IonState>& cpu_ions,
                                       size_t count,
                                       int stream_id = 0);

    /**
     * @brief Wait for a specific transfer to complete
     * @param event Event handle from start*Async call
     */
    void waitForTransfer(cudaEvent_t event);

    /**
     * @brief Get CUDA stream for computation overlap
     * @param stream_id Stream ID (0-2)
     * @return CUDA stream handle
     */
    cudaStream_t getComputeStream(int stream_id = 1);

    /**
     * @brief Synchronize all streams
     */
    void synchronizeAll();

    /**
     * @brief Get performance statistics
     */
    void printStats() const;

    /**
     * @brief Cleanup streams and events
     */
    ~AsyncTransferManager();

private:
    AsyncTransferManager() = default;
    AsyncTransferManager(const AsyncTransferManager&) = delete;
    AsyncTransferManager& operator=(const AsyncTransferManager&) = delete;

    // CUDA streams for async operations
    static const int NUM_STREAMS = 3;
    cudaStream_t streams_[NUM_STREAMS];
    
    // Events for synchronization
    std::vector<cudaEvent_t> event_pool_;
    
    // Pinned memory buffers for faster transfers
    IonState* pinned_cpu_buffer_ = nullptr;
    IonStateGPU* pinned_gpu_buffer_ = nullptr;
    size_t max_buffer_size_ = 0;
    
    // Performance tracking
    mutable size_t total_transfers_ = 0;
    mutable size_t async_transfers_ = 0;
    mutable double total_transfer_time_ms_ = 0;
    
    bool initialized_ = false;
    
    /**
     * @brief Convert IonState to IonStateGPU
     */
    void convertToGPU(const std::vector<IonState>& cpu_ions, 
                      IonStateGPU* gpu_buffer, 
                      size_t count);
    
    /**
     * @brief Convert IonStateGPU to IonState
     */
    void convertToCPU(const IonStateGPU* gpu_buffer, 
                      std::vector<IonState>& cpu_ions, 
                      size_t count);
    
    /**
     * @brief Get next available event from pool
     */
    cudaEvent_t getEvent();
    
    /**
     * @brief Return event to pool
     */
    void returnEvent(cudaEvent_t event);
};

#endif // USE_CUDA

} // namespace gpu
} // namespace ICARION