// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifndef ICARION_GPU_CONTEXT_H
#define ICARION_GPU_CONTEXT_H

#include <cuda_runtime.h>
#include <memory>
#include <string>
#include <vector>

namespace icarion {
namespace gpu {

/**
 * @brief GPU device properties and capabilities
 */
struct DeviceProperties {
    int device_id;
    std::string name;
    size_t total_memory;
    size_t free_memory;
    int compute_capability_major;
    int compute_capability_minor;
    int multiprocessor_count;
    int max_threads_per_block;
    int max_threads_per_multiprocessor;
    int warp_size;
    bool concurrent_kernels;
    bool unified_addressing;
};

/**
 * @brief RAII wrapper for CUDA device context and stream management
 * 
 * Manages CUDA device initialization, stream creation, and error handling.
 * Provides automatic cleanup on destruction. Used by GPU helpers; not all
 * helpers are fully integrated or validated.
 * 
 * Usage:
 *   auto context = GPUContext::create(0);  // Device 0
 *   cudaStream_t stream = context->get_stream();
 *   // ... use stream ...
 *   context->synchronize();
 */
class GPUContext {
public:
    /**
     * @brief Create GPU context for specified device
     * @param device_id CUDA device ID (default: 0)
     * @return Unique pointer to GPUContext, or nullptr if GPU unavailable
     */
    static std::unique_ptr<GPUContext> create(int device_id = 0);

    /**
     * @brief Check if CUDA is available on this system
     * @return true if at least one CUDA device is available
     */
    static bool is_cuda_available();

    /**
     * @brief Get number of available CUDA devices
     * @return Number of CUDA devices, or 0 if none
     */
    static int get_device_count();

    /**
     * @brief Destructor - cleans up CUDA resources
     */
    ~GPUContext();

    // Disable copy
    GPUContext(const GPUContext&) = delete;
    GPUContext& operator=(const GPUContext&) = delete;

    // Enable move
    GPUContext(GPUContext&& other) noexcept;
    GPUContext& operator=(GPUContext&& other) noexcept;

    /**
     * @brief Get the CUDA stream for async operations
     * @return cudaStream_t handle
     */
    cudaStream_t get_stream() const { return stream_; }

    /**
     * @brief Get device properties
     * @return DeviceProperties struct with device info
     */
    const DeviceProperties& get_properties() const { return properties_; }

    /**
     * @brief Get device ID
     * @return CUDA device ID
     */
    int get_device_id() const { return properties_.device_id; }

    /**
     * @brief Synchronize the stream (wait for all operations to complete)
     */
    void synchronize() const;

    /**
     * @brief Check if context is valid and operational
     * @return true if context can be used
     */
    bool is_valid() const { return initialized_; }

    /**
     * @brief Get last CUDA error as string
     * @return Error message, or empty string if no error
     */
    std::string get_last_error() const;

    /**
     * @brief Reset device and clear all allocations (for testing)
     * WARNING: Invalidates all CUDA resources!
     */
    void reset_device();

private:
    /**
     * @brief Private constructor - use create() factory method
     */
    GPUContext(int device_id);

    /**
     * @brief Initialize CUDA device and streams
     * @return true if initialization successful
     */
    bool initialize();

    /**
     * @brief Query and store device properties
     */
    void query_device_properties();

    /**
     * @brief Cleanup CUDA resources
     */
    void cleanup();

    int device_id_;                 ///< CUDA device ID
    cudaStream_t stream_;           ///< CUDA stream for async ops
    DeviceProperties properties_;   ///< Device capabilities
    bool initialized_;              ///< Initialization status
};

/**
 * @brief RAII wrapper for CUDA error checking
 * 
 * Usage:
 *   CUDA_CHECK(cudaMalloc(&ptr, size));
 */
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            throw std::runtime_error( \
                std::string("CUDA error at ") + __FILE__ + ":" + \
                std::to_string(__LINE__) + " - " + \
                cudaGetErrorString(err)); \
        } \
    } while (0)

/**
 * @brief Check CUDA error without throwing
 * @return true if no error
 */
inline bool cuda_check_error(cudaError_t err, const char* file, int line) {
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA error at %s:%d - %s\n", 
                file, line, cudaGetErrorString(err));
        return false;
    }
    return true;
}

#define CUDA_CHECK_NOTHROW(call) \
    cuda_check_error(call, __FILE__, __LINE__)

} // namespace gpu
} // namespace icarion

#endif // ICARION_GPU_CONTEXT_H
