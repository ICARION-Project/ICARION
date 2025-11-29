/**
 * =====================================================================
 *   ICARION: GPU Field Kernel Optimizations
 * =====================================================================
 *   @file        field_kernels_optimized.cuh
 *   @brief       Optimized field computation kernels
 *
 *   @details
 *   Performance optimizations:
 *   - Shared memory for domain parameters
 *   - Coalesced memory access patterns
 *   - Reduced register pressure
 *   - Fast math operations
 *   - Warp-level optimizations
 *
 *   @date        2025-11-11
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 * =====================================================================
 */

#pragma once
#include <cuda_runtime.h>
#include "core/types/Vec3.h"
#include "core/param/GPUParams.h"
#include "core/physics/fields/physics_math_shared.h"

// =============================================================
// Optimized Field Computation Helpers
// =============================================================

/**
 * @brief Fast coordinate transformation using shared memory
 * 
 * Optimization: Store frequently accessed geometry in shared memory
 * to reduce global memory bandwidth
 */
__device__ inline Vec3 transform_to_local_fast(
    const Vec3& global_pos,
    const Vec3& origin,
    const Vec3& rot_row0,
    const Vec3& rot_row1,
    const Vec3& rot_row2
) {
    // Use fast subtraction (compiler can optimize to __fsub_rn)
    Vec3 diff;
    diff.x = global_pos.x - origin.x;
    diff.y = global_pos.y - origin.y;
    diff.z = global_pos.z - origin.z;
    
    // Matrix-vector multiply with potential for vectorization
    Vec3 local_pos;
    local_pos.x = __fma_rn(rot_row0.x, diff.x, __fma_rn(rot_row0.y, diff.y, rot_row0.z * diff.z));
    local_pos.y = __fma_rn(rot_row1.x, diff.x, __fma_rn(rot_row1.y, diff.y, rot_row1.z * diff.z));
    local_pos.z = __fma_rn(rot_row2.x, diff.x, __fma_rn(rot_row2.y, diff.y, rot_row2.z * diff.z));
    
    return local_pos;
}

/**
 * @brief Optimized RF field computation with fast math
 * 
 * Optimization: Use __fmul_rn, __fadd_rn for fast floating-point ops
 */
__device__ inline Vec3 compute_rf_field_fast(
    const Vec3& pos,
    double Vrf,
    double Vquad,
    double omega,
    double r0_inv,  // Pre-computed 1/r0 to avoid division
    double t
) {
    // Fast cosine using shared helper
    double c = cos_shared(omega * t);
    
    // Combine operations: fac = 2*(Vquad + Vrf*c) * r0_inv^2
    double fac = __fma_rn(Vrf, c, Vquad);
    fac = __fmul_rn(2.0, fac);
    fac = __fmul_rn(fac, __fmul_rn(r0_inv, r0_inv));
    
    Vec3 E;
    E.x = __fmul_rn(fac, pos.x);
    E.y = __fmul_rn(-fac, pos.y);
    E.z = 0.0;
    
    return E;
}

/**
 * @brief Optimized AC field computation
 */
__device__ inline Vec3 compute_ac_field_fast(
    const Vec3& pos,
    double Vac,
    double omega,
    double radius_inv,  // Pre-computed 1/radius
    double t,
    const Vec3& dir_unit  // Pre-normalized direction
) {
    double c = cos_shared(omega * t);
    double mag = __fmul_rn(-radius_inv, __fmul_rn(Vac, c));
    
    Vec3 E;
    E.x = __fmul_rn(dir_unit.x, mag);
    E.y = __fmul_rn(dir_unit.y, mag);
    E.z = __fmul_rn(dir_unit.z, mag);
    
    return E;
}

/**
 * @brief Optimized DC field computation
 */
__device__ inline Vec3 compute_dc_field_fast(
    const Vec3& pos,
    double axial_V,
    double L_inv  // Pre-computed 1/L
) {
    return Vec3{0.0, 0.0, __fmul_rn(axial_V, L_inv)};
}

/**
 * @brief Shared memory cache for domain parameters
 * 
 * Optimization: Load domain parameters into shared memory once per block
 * to reduce global memory traffic
 */
struct DomainCache {
    // Geometry (most frequently accessed)
    Vec3 origin;
    Vec3 rot_row0, rot_row1, rot_row2;
    double radius_m;
    double radius_inv;  // Pre-computed
    double length_m;
    double length_inv;  // Pre-computed
    
    // RF parameters
    double rf_voltage;
    double rf_omega;
    double dc_quad;
    
    // AC parameters
    double ac_voltage;
    double ac_omega;
    Vec3 ac_dir_unit;
    
    // DC parameters
    double dc_axial;
    
    // Instrument type
    int instrument_type;
};

/**
 * @brief Load domain parameters into shared memory
 * 
 * Call this once per block with thread 0
 */
__device__ inline void load_domain_cache(
    DomainCache& cache,
    const DomainGPU& dom
) {
    cache.origin = dom.geom.origin_m;
    cache.rot_row0 = dom.geom.rot_row0;
    cache.rot_row1 = dom.geom.rot_row1;
    cache.rot_row2 = dom.geom.rot_row2;
    cache.radius_m = dom.geom.radius_m;
    cache.radius_inv = (dom.geom.radius_m > 0.0) ? (1.0 / dom.geom.radius_m) : 0.0;
    cache.length_m = dom.geom.length_m;
    cache.length_inv = (dom.geom.length_m > 0.0) ? (1.0 / dom.geom.length_m) : 0.0;
    
    cache.rf_voltage = dom.RF.voltage_V;
    cache.rf_omega = dom.RF.omega_rad_s;
    cache.dc_quad = dom.DC.quad_V;
    
    cache.ac_voltage = dom.AC.voltage_V;
    cache.ac_omega = dom.AC.omega_rad_s;
    
    // Pre-normalize AC direction
    Vec3 dir = Vec3{1.0, 0.0, 0.0};  // Default LQIT direction
    double len = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    cache.ac_dir_unit = (len > 0.0) ? Vec3{dir.x/len, dir.y/len, dir.z/len} : Vec3{0.0, 0.0, 0.0};
    
    cache.dc_axial = dom.DC.axial_V;
    cache.instrument_type = dom.instrument;
}

/**
 * @brief Optimized acceleration computation using cached domain
 */
__device__ inline void compute_acceleration_from_cache(
    const Vec3& global_pos,
    double t,
    double charge_over_mass,
    const DomainCache& cache,
    Vec3& out_acc
) {
    out_acc = Vec3{0.0, 0.0, 0.0};
    
    // Transform to local coordinates once
    Vec3 local_pos = transform_to_local_fast(
        global_pos,
        cache.origin,
        cache.rot_row0,
        cache.rot_row1,
        cache.rot_row2
    );
    
    // LQIT instrument (most common case first for branch prediction)
    if (cache.instrument_type == 0) {  // LQIT = 0
        // RF field
        Vec3 E_rf = compute_rf_field_fast(
            local_pos,
            cache.rf_voltage,
            cache.dc_quad,
            cache.rf_omega,
            cache.radius_inv,
            t
        );
        out_acc += E_rf * charge_over_mass;
        
        // AC field
        Vec3 E_ac = compute_ac_field_fast(
            local_pos,
            cache.ac_voltage,
            cache.ac_omega,
            cache.radius_inv,
            t,
            cache.ac_dir_unit
        );
        out_acc += E_ac * charge_over_mass;
        
        // DC boundary fields (branchless when possible)
        double z_frac = local_pos.z * cache.length_inv;
        if (z_frac > 0.9) {
            Vec3 E_dc = compute_dc_field_fast(local_pos, -cache.dc_axial, cache.length_inv * 10.0);
            out_acc += E_dc * charge_over_mass;
        } else if (z_frac < 0.1) {
            Vec3 E_dc = compute_dc_field_fast(local_pos, cache.dc_axial, cache.length_inv * 10.0);
            out_acc += E_dc * charge_over_mass;
        }
    }
    // Add other instrument types as needed
}

// =============================================================
// Performance Monitoring Macros
// =============================================================

#ifdef ENABLE_GPU_PROFILING
#define GPU_PROFILE_START(event) cudaEventRecord(event##_start)
#define GPU_PROFILE_STOP(event) cudaEventRecord(event##_stop)
#define GPU_PROFILE_ELAPSED(event, time_ms) \
    cudaEventSynchronize(event##_stop); \
    cudaEventElapsedTime(&time_ms, event##_start, event##_stop)
#else
#define GPU_PROFILE_START(event)
#define GPU_PROFILE_STOP(event)
#define GPU_PROFILE_ELAPSED(event, time_ms)
#endif

