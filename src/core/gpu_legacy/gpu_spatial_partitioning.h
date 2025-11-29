/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   GPU-accelerated spatial partitioning system for O(N) collision detection
 *
 *   @file        gpu_spatial_partitioning.h  
 *   @brief       GPU spatial hashing for efficient neighbor searching
 *
 *   @details
 *   Transforms collision detection and spatial queries from O(N²) to O(N) using:
 *   - 3D uniform spatial grid with GPU hash tables
 *   - Spatial binning on GPU with CUDA kernels
 *   - Efficient neighbor cell traversal with memory coalescing
 *   - Memory pool integration for allocation efficiency
 *
 *   Target Performance: 5-20× speedup for collision detection in dense systems
 *   Compatibility: Integrates with existing memory pool and async transfer systems
 *
 *   @date        2025-01-27
 *   @version     1.0
 *   @author      AI Assistant
 *   @license     MIT License
 *
 * =====================================================================
 */

#ifndef GPU_SPATIAL_PARTITIONING_H
#define GPU_SPATIAL_PARTITIONING_H

#include "core/types/IonState.h"
#include "core/utils/mathUtils.h"
#include "gpu_memory_pool.h"
#include <cuda_runtime.h>
#include <vector>
#include <memory>

using namespace ICARION::core;

/**
 * @class GPUSpatialPartitioningSystem
 * @brief GPU-accelerated spatial partitioning for O(N) collision detection
 *
 * Implements spatial hashing on GPU to transform collision detection
 * from O(N²) pairwise checks to O(N) spatial grid lookups.
 *
 * Algorithm Overview:
 * 1. Create 3D uniform grid based on interaction radius
 * 2. Hash ion positions to grid cells on GPU
 * 3. Use cell-based neighbor lookup for interactions
 * 4. Support dynamic grid resizing based on ion distribution
 *
 * Performance Benefits:
 * - O(N) collision detection vs O(N²) brute force
 * - GPU memory coalescing for spatial queries  
 * - Efficient memory pool integration
 * - Adaptive grid sizing for optimal performance
 */
class GPUSpatialPartitioningSystem {
public:
    /**
     * @brief Initialize spatial partitioning system
     * @param initial_capacity Maximum number of ions to support
     * @param cell_size Spatial resolution of grid cells [m]
     * @param bounds_min Minimum coordinates of simulation domain [m]  
     * @param bounds_max Maximum coordinates of simulation domain [m]
     */
    explicit GPUSpatialPartitioningSystem(int initial_capacity = 10000,
                                          double cell_size = 1e-3,
                                          Vec3 bounds_min = {-0.1, -0.1, -0.1},
                                          Vec3 bounds_max = {0.1, 0.1, 0.1});
                                          
    /**
     * @brief Destructor - cleanup GPU memory
     */
    ~GPUSpatialPartitioningSystem();

    /**
     * @brief Update spatial grid with current ion positions
     * @param ions Ion states to partition spatially
     * @param stream CUDA stream for async execution
     * @return True if partitioning was successful
     */
    bool updateGrid(const std::vector<IonState>& ions, cudaStream_t stream = nullptr);
    
    /**
     * @brief Find all ions within interaction radius of target position
     * @param target_pos Position to search around [m]
     * @param radius Search radius [m]
     * @param neighbor_indices Output buffer for neighbor ion indices
     * @param max_neighbors Maximum neighbors to return
     * @param stream CUDA stream for async execution
     * @return Number of neighbors found (capped at max_neighbors)
     */
    int findNeighbors(const Vec3& target_pos, double radius, 
                      int* neighbor_indices, int max_neighbors,
                      cudaStream_t stream = nullptr);
    
    /**
     * @brief Get collision pairs using spatial partitioning (O(N) complexity)
     * @param collision_radius Interaction radius for collision detection [m]
     * @param collision_pairs Output buffer for ion index pairs (i,j)
     * @param max_pairs Maximum number of pairs to return
     * @param stream CUDA stream for async execution
     * @return Number of collision pairs found
     */
    int getCollisionPairs(double collision_radius, 
                          int* collision_pairs, int max_pairs,
                          cudaStream_t stream = nullptr);
    
    /**
     * @brief Compute space charge field using spatial partitioning
     * @param target_ion Ion for which to compute field
     * @param field_radius Maximum distance for field contributions [m]
     * @param stream CUDA stream for async execution
     * @return Electric field [V/m] from nearby ions
     */
    Vec3 computeSpaceChargeField(const IonState& target_ion, double field_radius,
                                 cudaStream_t stream = nullptr);
    
    /**
     * @brief Get partitioning performance statistics
     * @return Struct containing timing and efficiency metrics
     */
    struct PartitioningStats {
        double grid_build_time_ms;      ///< Time to build spatial grid
        double neighbor_search_time_ms; ///< Average neighbor search time  
        int total_cells;                ///< Number of active grid cells
        int max_ions_per_cell;          ///< Load balancing metric
        double memory_usage_mb;         ///< GPU memory consumption
        double effective_speedup;       ///< Measured speedup vs O(N²)
    };
    PartitioningStats getStats() const { return stats_; }
    
    /**
     * @brief Resize grid to optimize for current ion distribution
     * @param ions Current ion distribution for analysis
     * @return True if resize was beneficial
     */
    bool optimizeGrid(const std::vector<IonState>& ions);

private:
    // GPU memory management
    ICARION::gpu::GPUMemoryPool* memory_pool_;      ///< Shared memory pool for efficiency
    
    // Spatial grid parameters
    Vec3 bounds_min_, bounds_max_;    ///< Simulation domain bounds [m]
    double cell_size_;                ///< Grid cell size [m]
    int grid_nx_, grid_ny_, grid_nz_; ///< Grid dimensions
    int total_cells_;                 ///< Total number of cells
    
    // GPU data structures
    int* d_cell_ion_count_;           ///< Number of ions per cell
    int* d_cell_ion_indices_;         ///< Ion indices sorted by cell
    int* d_ion_cell_ids_;             ///< Cell ID for each ion
    int* d_cell_offsets_;             ///< Starting index for each cell in sorted array
    
    // Device memory buffers
    Vec3* d_ion_positions_;           ///< Ion positions on GPU
    int* d_temp_storage_;             ///< Temporary storage for sorting
    int capacity_;                    ///< Maximum ions supported
    int current_ion_count_;           ///< Current number of ions
    
    // Performance tracking
    mutable PartitioningStats stats_;
    cudaEvent_t start_event_, stop_event_;
    
    // Internal helper methods
    
    /**
     * @brief Convert 3D position to grid cell coordinates
     */
    __device__ __host__ int3 positionToCell(const Vec3& pos) const;
    
    /**
     * @brief Convert cell coordinates to linear cell ID
     */
    __device__ __host__ int cellToIndex(int cx, int cy, int cz) const;
    
    /**
     * @brief Hash ion positions to grid cells on GPU
     */
    void hashIonsToGrid(const std::vector<IonState>& ions, cudaStream_t stream);
    
    /**
     * @brief Sort ions by cell ID for memory coalescing
     */
    void sortIonsByCell(cudaStream_t stream);
    
    /**
     * @brief Build cell offset array for efficient neighbor access
     */  
    void buildCellOffsets(cudaStream_t stream);
    
    /**
     * @brief Check if current grid size is optimal for ion distribution
     */
    bool shouldResizeGrid(const std::vector<IonState>& ions) const;
    
    /**
     * @brief Allocate GPU memory for spatial grid structures
     */
    bool allocateDeviceMemory();
    
    /**
     * @brief Free all GPU memory
     */
    void deallocateDeviceMemory();
    
    /**
     * @brief Update performance statistics
     */
    void updateStats(double operation_time_ms);
};

// CUDA kernel function declarations
extern "C" {
    /**
     * @brief GPU kernel to hash ion positions to grid cells
     */
    __global__ void hashIonsKernel(const Vec3* ion_positions, int* ion_cell_ids, 
                                   int num_ions, Vec3 bounds_min, double inv_cell_size,
                                   int grid_nx, int grid_ny, int grid_nz);
    
    /**
     * @brief GPU kernel to find neighbors in spatial grid
     */
    __global__ void findNeighborsKernel(const Vec3& target_pos, double radius,
                                        const Vec3* ion_positions, const int* cell_offsets,
                                        const int* cell_ion_indices, int* neighbor_indices,
                                        int* neighbor_count, int max_neighbors,
                                        Vec3 bounds_min, double inv_cell_size,
                                        int grid_nx, int grid_ny, int grid_nz, int total_ions);
    
    /**
     * @brief GPU kernel to detect collision pairs using spatial partitioning
     */
    __global__ void findCollisionPairsKernel(const Vec3* ion_positions, 
                                             const int* cell_offsets, const int* cell_ion_indices,
                                             double collision_radius, int* collision_pairs,
                                             int* pair_count, int max_pairs,
                                             Vec3 bounds_min, double inv_cell_size,
                                             int grid_nx, int grid_ny, int grid_nz, int total_ions);
}

#endif // GPU_SPATIAL_PARTITIONING_H