// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "GPUContext.h"
#include <stdexcept>
#include <sstream>

namespace icarion {
namespace gpu {

// Static factory method
std::unique_ptr<GPUContext> GPUContext::create(int device_id) {
    if (!is_cuda_available()) {
        return nullptr;
    }

    int device_count = get_device_count();
    if (device_id < 0 || device_id >= device_count) {
        return nullptr;
    }

    auto context = std::unique_ptr<GPUContext>(new GPUContext(device_id));
    if (!context->initialize()) {
        return nullptr;
    }

    return context;
}

bool GPUContext::is_cuda_available() {
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    
    // cudaErrorNoDevice or cudaErrorInsufficientDriver are expected when no GPU
    if (err == cudaErrorNoDevice || err == cudaErrorInsufficientDriver) {
        return false;
    }
    
    if (err != cudaSuccess) {
        return false;
    }
    
    return device_count > 0;
}

int GPUContext::get_device_count() {
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    
    if (err != cudaSuccess) {
        return 0;
    }
    
    return device_count;
}

// Private constructor
GPUContext::GPUContext(int device_id)
    : device_id_(device_id)
    , stream_(nullptr)
    , properties_{}
    , initialized_(false)
{
}

// Destructor
GPUContext::~GPUContext() {
    cleanup();
}

// Move constructor
GPUContext::GPUContext(GPUContext&& other) noexcept
    : device_id_(other.device_id_)
    , stream_(other.stream_)
    , properties_(other.properties_)
    , initialized_(other.initialized_)
{
    other.stream_ = nullptr;
    other.initialized_ = false;
}

// Move assignment
GPUContext& GPUContext::operator=(GPUContext&& other) noexcept {
    if (this != &other) {
        cleanup();
        
        device_id_ = other.device_id_;
        stream_ = other.stream_;
        properties_ = other.properties_;
        initialized_ = other.initialized_;
        
        other.stream_ = nullptr;
        other.initialized_ = false;
    }
    return *this;
}

bool GPUContext::initialize() {
    try {
        // Set device
        CUDA_CHECK(cudaSetDevice(device_id_));
        
        // Query device properties
        query_device_properties();
        
        // Create stream
        CUDA_CHECK(cudaStreamCreate(&stream_));
        
        initialized_ = true;
        return true;
    }
    catch (const std::exception& e) {
        cleanup();
        return false;
    }
}

void GPUContext::query_device_properties() {
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device_id_));
    
    properties_.device_id = device_id_;
    properties_.name = prop.name;
    properties_.compute_capability_major = prop.major;
    properties_.compute_capability_minor = prop.minor;
    properties_.multiprocessor_count = prop.multiProcessorCount;
    properties_.max_threads_per_block = prop.maxThreadsPerBlock;
    properties_.max_threads_per_multiprocessor = prop.maxThreadsPerMultiProcessor;
    properties_.warp_size = prop.warpSize;
    properties_.concurrent_kernels = prop.concurrentKernels;
    properties_.unified_addressing = prop.unifiedAddressing;
    
    // Query memory info
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    properties_.free_memory = free_mem;
    properties_.total_memory = total_mem;
}

void GPUContext::cleanup() {
    if (stream_ != nullptr) {
        cudaStreamDestroy(stream_);
        stream_ = nullptr;
    }
    initialized_ = false;
}

void GPUContext::synchronize() const {
    if (!initialized_) {
        throw std::runtime_error("GPUContext not initialized");
    }
    
    CUDA_CHECK(cudaStreamSynchronize(stream_));
}

std::string GPUContext::get_last_error() const {
    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess) {
        return "";
    }
    
    std::stringstream ss;
    ss << "CUDA error: " << cudaGetErrorString(err) 
       << " (code " << err << ")";
    return ss.str();
}

void GPUContext::reset_device() {
    cleanup();
    CUDA_CHECK(cudaDeviceReset());
}

} // namespace gpu
} // namespace icarion
