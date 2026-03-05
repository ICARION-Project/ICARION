// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_GPU_MEMORY_POOL_H
#define ICARION_GPU_MEMORY_POOL_H

#include "core/gpu/core/GPUContext.h"
#include <cuda_runtime.h>
#include <memory>
#include <vector>
#include <unordered_map>
#include <cstddef>

namespace icarion {
namespace gpu {

/**
 * @brief GPU device buffer with automatic cleanup
 */
template<typename T>
class DeviceBuffer {
public:
    DeviceBuffer() : ptr_(nullptr), size_(0) {}
    
    explicit DeviceBuffer(size_t count)
        : ptr_(nullptr), size_(count)
    {
        if (count > 0) {
            CUDA_CHECK(cudaMalloc(&ptr_, count * sizeof(T)));
        }
    }
    
    ~DeviceBuffer() {
        if (ptr_ != nullptr) {
            cudaFree(ptr_);
        }
    }
    
    // Disable copy
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    
    // Enable move
    DeviceBuffer(DeviceBuffer&& other) noexcept
        : ptr_(other.ptr_), size_(other.size_)
    {
        other.ptr_ = nullptr;
        other.size_ = 0;
    }
    
    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            if (ptr_ != nullptr) {
                cudaFree(ptr_);
            }
            ptr_ = other.ptr_;
            size_ = other.size_;
            other.ptr_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }
    
    T* get() { return ptr_; }
    const T* get() const { return ptr_; }
    size_t size() const { return size_; }
    
    /**
     * @brief Upload data from host to device (sync)
     */
    void upload(const T* host_data, size_t count) {
        if (count > size_) {
            throw std::runtime_error("DeviceBuffer::upload: count exceeds buffer size");
        }
        CUDA_CHECK(cudaMemcpy(ptr_, host_data, count * sizeof(T), cudaMemcpyHostToDevice));
    }
    
    /**
     * @brief Download data from device to host (sync)
     */
    void download(T* host_data, size_t count) const {
        if (count > size_) {
            throw std::runtime_error("DeviceBuffer::download: count exceeds buffer size");
        }
        CUDA_CHECK(cudaMemcpy(host_data, ptr_, count * sizeof(T), cudaMemcpyDeviceToHost));
    }
    
    /**
     * @brief Upload data asynchronously
     */
    void upload_async(const T* host_data, size_t count, cudaStream_t stream) {
        if (count > size_) {
            throw std::runtime_error("DeviceBuffer::upload_async: count exceeds buffer size");
        }
        CUDA_CHECK(cudaMemcpyAsync(ptr_, host_data, count * sizeof(T), 
                                   cudaMemcpyHostToDevice, stream));
    }
    
    /**
     * @brief Download data asynchronously
     */
    void download_async(T* host_data, size_t count, cudaStream_t stream) const {
        if (count > size_) {
            throw std::runtime_error("DeviceBuffer::download_async: count exceeds buffer size");
        }
        CUDA_CHECK(cudaMemcpyAsync(host_data, ptr_, count * sizeof(T), 
                                   cudaMemcpyDeviceToHost, stream));
    }

private:
    T* ptr_;
    size_t size_;
};

/**
 * @brief Pinned host buffer for fast GPU transfers
 */
template<typename T>
class PinnedBuffer {
public:
    PinnedBuffer() : ptr_(nullptr), size_(0) {}
    
    explicit PinnedBuffer(size_t count)
        : ptr_(nullptr), size_(count)
    {
        if (count > 0) {
            CUDA_CHECK(cudaMallocHost(&ptr_, count * sizeof(T)));
        }
    }
    
    ~PinnedBuffer() {
        if (ptr_ != nullptr) {
            cudaFreeHost(ptr_);
        }
    }
    
    // Disable copy
    PinnedBuffer(const PinnedBuffer&) = delete;
    PinnedBuffer& operator=(const PinnedBuffer&) = delete;
    
    // Enable move
    PinnedBuffer(PinnedBuffer&& other) noexcept
        : ptr_(other.ptr_), size_(other.size_)
    {
        other.ptr_ = nullptr;
        other.size_ = 0;
    }
    
    PinnedBuffer& operator=(PinnedBuffer&& other) noexcept {
        if (this != &other) {
            if (ptr_ != nullptr) {
                cudaFreeHost(ptr_);
            }
            ptr_ = other.ptr_;
            size_ = other.size_;
            other.ptr_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }
    
    T* get() { return ptr_; }
    const T* get() const { return ptr_; }
    size_t size() const { return size_; }
    
    T& operator[](size_t i) { return ptr_[i]; }
    const T& operator[](size_t i) const { return ptr_[i]; }

private:
    T* ptr_;
    size_t size_;
};

/**
 * @brief Memory pool with buffer reuse to reduce allocation overhead
 * 
 * Manages device and pinned host memory with automatic reuse.
 * Reduces cudaMalloc/cudaFree overhead by keeping buffers alive.
 * 
 * Usage:
 *   GPUMemoryPool pool(context);
 *   auto dev_buf = pool.get_device_buffer<double>(1000);
 *   auto host_buf = pool.get_pinned_buffer<double>(1000);
 *   // ... use buffers ...
 *   pool.release_device_buffer(std::move(dev_buf));
 *   pool.release_pinned_buffer(std::move(host_buf));
 */
class GPUMemoryPool {
public:
    /**
     * @brief Constructor
     * @param context GPU context (must outlive the pool)
     */
    explicit GPUMemoryPool(const GPUContext& context);
    
    /**
     * @brief Destructor - frees all cached buffers
     */
    ~GPUMemoryPool();
    
    // Disable copy
    GPUMemoryPool(const GPUMemoryPool&) = delete;
    GPUMemoryPool& operator=(const GPUMemoryPool&) = delete;
    
    /**
     * @brief Get device buffer (allocates or reuses from pool)
     * @param count Number of elements
     * @return Device buffer
     */
    template<typename T>
    DeviceBuffer<T> get_device_buffer(size_t count);
    
    /**
     * @brief Get pinned host buffer (allocates or reuses from pool)
     * @param count Number of elements
     * @return Pinned host buffer
     */
    template<typename T>
    PinnedBuffer<T> get_pinned_buffer(size_t count);
    
    /**
     * @brief Return device buffer to pool for reuse
     * @param buffer Buffer to return
     */
    template<typename T>
    void release_device_buffer(DeviceBuffer<T>&& buffer);
    
    /**
     * @brief Return pinned buffer to pool for reuse
     * @param buffer Buffer to return
     */
    template<typename T>
    void release_pinned_buffer(PinnedBuffer<T>&& buffer);
    
    /**
     * @brief Clear all cached buffers (frees memory)
     */
    void clear();
    
    /**
     * @brief Get memory pool statistics
     */
    struct Stats {
        size_t device_buffers_allocated;
        size_t device_buffers_cached;
        size_t device_bytes_allocated;
        size_t device_bytes_cached;
        size_t pinned_buffers_allocated;
        size_t pinned_buffers_cached;
        size_t pinned_bytes_allocated;
        size_t pinned_bytes_cached;
    };
    
    Stats get_stats() const;

private:
    struct BufferKey {
        size_t size;        // Size in bytes
        size_t type_hash;   // Hash of type
        
        bool operator==(const BufferKey& other) const {
            return size == other.size && type_hash == other.type_hash;
        }
    };
    
    struct BufferKeyHash {
        size_t operator()(const BufferKey& key) const {
            return key.size ^ (key.type_hash << 1);
        }
    };
    
    const GPUContext& context_;
    
    // Free lists for buffer reuse
    std::unordered_map<BufferKey, std::vector<void*>, BufferKeyHash> device_free_list_;
    std::unordered_map<BufferKey, std::vector<void*>, BufferKeyHash> pinned_free_list_;
    
    // Statistics
    mutable Stats stats_;
};

// Template implementations

template<typename T>
DeviceBuffer<T> GPUMemoryPool::get_device_buffer(size_t count) {
    BufferKey key{count * sizeof(T), typeid(T).hash_code()};
    
    // Check free list
    auto it = device_free_list_.find(key);
    if (it != device_free_list_.end() && !it->second.empty()) {
        void* ptr = it->second.back();
        it->second.pop_back();
        stats_.device_bytes_cached -= key.size;
        stats_.device_buffers_cached--;
        
        // Create buffer from cached pointer
        DeviceBuffer<T> buffer;
        // Hackish but safe: we know the buffer is valid
        buffer = DeviceBuffer<T>(0);
        *reinterpret_cast<T**>(&buffer) = static_cast<T*>(ptr);
        *reinterpret_cast<size_t*>(reinterpret_cast<char*>(&buffer) + sizeof(T*)) = count;
        return buffer;
    }
    
    // Allocate new buffer
    stats_.device_buffers_allocated++;
    stats_.device_bytes_allocated += key.size;
    return DeviceBuffer<T>(count);
}

template<typename T>
PinnedBuffer<T> GPUMemoryPool::get_pinned_buffer(size_t count) {
    BufferKey key{count * sizeof(T), typeid(T).hash_code()};
    
    // Check free list
    auto it = pinned_free_list_.find(key);
    if (it != pinned_free_list_.end() && !it->second.empty()) {
        void* ptr = it->second.back();
        it->second.pop_back();
        stats_.pinned_bytes_cached -= key.size;
        stats_.pinned_buffers_cached--;
        
        // Create buffer from cached pointer
        PinnedBuffer<T> buffer;
        *reinterpret_cast<T**>(&buffer) = static_cast<T*>(ptr);
        *reinterpret_cast<size_t*>(reinterpret_cast<char*>(&buffer) + sizeof(T*)) = count;
        return buffer;
    }
    
    // Allocate new buffer
    stats_.pinned_buffers_allocated++;
    stats_.pinned_bytes_allocated += key.size;
    return PinnedBuffer<T>(count);
}

template<typename T>
void GPUMemoryPool::release_device_buffer(DeviceBuffer<T>&& buffer) {
    if (buffer.get() == nullptr || buffer.size() == 0) {
        return;
    }
    
    BufferKey key{buffer.size() * sizeof(T), typeid(T).hash_code()};
    
    // Add to free list
    device_free_list_[key].push_back(buffer.get());
    stats_.device_bytes_cached += key.size;
    stats_.device_buffers_cached++;
    
    // Prevent buffer from freeing memory in destructor
    *reinterpret_cast<T**>(&buffer) = nullptr;
    *reinterpret_cast<size_t*>(reinterpret_cast<char*>(&buffer) + sizeof(T*)) = 0;
}

template<typename T>
void GPUMemoryPool::release_pinned_buffer(PinnedBuffer<T>&& buffer) {
    if (buffer.get() == nullptr || buffer.size() == 0) {
        return;
    }
    
    BufferKey key{buffer.size() * sizeof(T), typeid(T).hash_code()};
    
    // Add to free list
    pinned_free_list_[key].push_back(buffer.get());
    stats_.pinned_bytes_cached += key.size;
    stats_.pinned_buffers_cached++;
    
    // Prevent buffer from freeing memory in destructor
    *reinterpret_cast<T**>(&buffer) = nullptr;
    *reinterpret_cast<size_t*>(reinterpret_cast<char*>(&buffer) + sizeof(T*)) = 0;
}

} // namespace gpu
} // namespace icarion

#endif // ICARION_GPU_MEMORY_POOL_H
