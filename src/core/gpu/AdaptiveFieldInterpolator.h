// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#pragma once

#include "core/types/Vec3.h"
#include <vector>
#include <memory>

namespace icarion {
namespace gpu {

/**
 * @brief Adaptive field interpolation with gradient-based LOD selection
 * 
 * **Purpose:** Improve field interpolation accuracy by adaptively selecting
 * interpolation order based on local field gradient magnitude.
 * 
 * **Algorithm:**
 * - **High gradient** (|∇E| > threshold): Use cubic/quintic interpolation
 * - **Low gradient** (|∇E| ≤ threshold): Use linear/trilinear interpolation
 * 
 * **Benefits:**
 * - Better accuracy near electrodes (steep gradients)
 * - Lower cost in field-free regions (smooth gradients)
 * - Reduces particle loss due to poor interpolation
 * 
 * **Typical Use Case:**
 * Ion traps with sharp potential wells where linear interpolation causes
 * artificial heating and particle loss.
 * 
 * **Example:**
 * @code
 * AdaptiveFieldInterpolator interp(field_grid, gradient_threshold);
 * 
 * for (auto& ion : ions) {
 *     Vec3 E = interp.interpolate(ion.pos);  // Automatic LOD selection
 *     // interp uses cubic near electrodes, linear far away
 * }
 * @endcode
 * 
 * **Performance:**
 * - Gradient computation: O(1) per query (cached gradients)
 * - Interpolation: 2-10× slower than linear, but only where needed
 * - Overall: ~20% overhead for typical ion trap simulations
 * 
 * **Limitations:**
 * - Requires precomputed gradient field (extra memory: 3× field size)
 * - Not beneficial for smooth fields (e.g., uniform drift tubes)
 * - CPU-only for now (GPU version in Phase 13)
 */
class AdaptiveFieldInterpolator {
public:
    /**
     * @brief Interpolation order (Level of Detail)
     */
    enum class InterpolationOrder {
        Linear,      ///< 1st-order: trilinear (8 samples, C⁰ continuous)
        Cubic,       ///< 3rd-order: tricubic (64 samples, C¹ continuous)
        Quintic      ///< 5th-order: triquintic (216 samples, C² continuous)
    };
    
    /**
     * @brief Configuration for adaptive interpolation
     */
    struct Config {
        double gradient_threshold_low = 1e3;   ///< Below this: linear [V/m²]
        double gradient_threshold_high = 1e5;  ///< Above this: quintic [V/m²]
        InterpolationOrder max_order = InterpolationOrder::Cubic;  ///< Maximum order to use
        bool cache_gradients = true;           ///< Precompute and cache gradients
        bool enable_gpu = false;               ///< Use GPU for interpolation (Phase 13)
    };
    
    /**
     * @brief Construct adaptive interpolator
     * 
     * @param field_grid Electric field grid [V/m] (nx×ny×nz×3)
     * @param nx, ny, nz Grid dimensions
     * @param origin Grid origin [m]
     * @param spacing Cell spacing [m]
     * @param config Adaptive interpolation configuration
     * 
     * **Memory Usage:**
     * - Field grid: nx×ny×nz×3 × 8 bytes (Ex, Ey, Ez as double)
     * - Gradient grid: nx×ny×nz×9 × 8 bytes (∂Ex/∂x, ∂Ex/∂y, ...) if cached
     * - Total: ~4× field size for cached gradients
     * 
     * **Initialization Time:**
     * - Gradient computation: O(nx×ny×nz) via central differences
     * - Typical: ~10 ms for 128³ grid
     */
    // Constructor with custom config
    AdaptiveFieldInterpolator(
        const std::vector<double>& Ex_grid,
        const std::vector<double>& Ey_grid,
        const std::vector<double>& Ez_grid,
        int nx, int ny, int nz,
        Vec3 origin, Vec3 spacing,
        const Config& config
    );
    
    // Constructor with default config
    AdaptiveFieldInterpolator(
        const std::vector<double>& Ex_grid,
        const std::vector<double>& Ey_grid,
        const std::vector<double>& Ez_grid,
        int nx, int ny, int nz,
        Vec3 origin, Vec3 spacing
    );
    
    /**
     * @brief Interpolate electric field at position with adaptive LOD
     * 
     * @param pos Position in global coordinates [m]
     * @return Electric field [V/m]
     * 
     * **Algorithm:**
     * 1. Compute fractional grid position: (x-x0)/dx
     * 2. Estimate local gradient: |∇E| from cached gradient grid
     * 3. Select interpolation order based on gradient threshold
     * 4. Perform interpolation (linear/cubic/quintic)
     * 
     * **Complexity:**
     * - Linear: 8 samples, ~10 flops
     * - Cubic: 64 samples, ~200 flops
     * - Quintic: 216 samples, ~600 flops
     * 
     * **Boundary Handling:**
     * - Inside grid: Interpolate normally
     * - Outside grid: Return {0, 0, 0} (or extrapolate if enabled)
     */
    Vec3 interpolate(const Vec3& pos) const;
    
    /**
     * @brief Get interpolation order at position (for debugging)
     * 
     * @param pos Position in global coordinates [m]
     * @return Selected interpolation order
     */
    InterpolationOrder get_order_at(const Vec3& pos) const;
    
    /**
     * @brief Get statistics on interpolation order usage
     * 
     * @return {n_linear, n_cubic, n_quintic} calls since last reset
     */
    struct Stats {
        size_t n_linear = 0;
        size_t n_cubic = 0;
        size_t n_quintic = 0;
        
        double fraction_linear() const {
            size_t total = n_linear + n_cubic + n_quintic;
            return total > 0 ? static_cast<double>(n_linear) / total : 0.0;
        }
        
        double fraction_cubic() const {
            size_t total = n_linear + n_cubic + n_quintic;
            return total > 0 ? static_cast<double>(n_cubic) / total : 0.0;
        }
        
        double fraction_quintic() const {
            size_t total = n_linear + n_cubic + n_quintic;
            return total > 0 ? static_cast<double>(n_quintic) / total : 0.0;
        }
    };
    
    const Stats& get_stats() const { return stats_; }
    void reset_stats() { stats_ = Stats{}; }
    
private:
    // Field grid data
    std::vector<double> Ex_, Ey_, Ez_;  ///< Electric field components [V/m]
    int nx_, ny_, nz_;                  ///< Grid dimensions
    Vec3 origin_;                       ///< Grid origin [m]
    Vec3 spacing_;                      ///< Cell spacing [m]
    Vec3 inv_spacing_;                  ///< 1/spacing for fast division
    
    // Gradient cache (optional)
    std::vector<double> grad_Ex_x_, grad_Ex_y_, grad_Ex_z_;  ///< ∂Ex/∂x, ∂Ex/∂y, ∂Ex/∂z
    std::vector<double> grad_Ey_x_, grad_Ey_y_, grad_Ey_z_;  ///< ∂Ey/∂x, ∂Ey/∂y, ∂Ey/∂z
    std::vector<double> grad_Ez_x_, grad_Ez_y_, grad_Ez_z_;  ///< ∂Ez/∂x, ∂Ez/∂y, ∂Ez/∂z
    bool has_gradients_ = false;
    
    // Configuration
    Config config_;
    
    // Statistics (mutable for const interpolate())
    mutable Stats stats_;
    
    /**
     * @brief Compute field gradients via central differences
     * 
     * For each field component (Ex, Ey, Ez), compute:
     * - ∂E/∂x = (E(i+1,j,k) - E(i-1,j,k)) / (2*dx)
     * - ∂E/∂y = (E(i,j+1,k) - E(i,j-1,k)) / (2*dy)
     * - ∂E/∂z = (E(i,j,k+1) - E(i,j,k-1)) / (2*dz)
     * 
     * Forward/backward differences at boundaries.
     */
    void compute_gradients();
    
    /**
     * @brief Estimate gradient magnitude at position
     * 
     * @param pos Position in global coordinates [m]
     * @return |∇E| ≈ √(Σ (∂Ei/∂xj)²) [V/m²]
     * 
     * Uses trilinear interpolation of cached gradient field.
     */
    double estimate_gradient_magnitude(const Vec3& pos) const;
    
    /**
     * @brief Select interpolation order based on gradient
     * 
     * @param grad_magnitude Local gradient magnitude [V/m²]
     * @return Interpolation order to use
     * 
     * Decision logic:
     * - grad < threshold_low  → Linear
     * - grad > threshold_high → max_order (Cubic or Quintic)
     * - In between            → Cubic (smooth transition)
     */
    InterpolationOrder select_order(double grad_magnitude) const;
    
    /**
     * @brief Trilinear interpolation (8 samples)
     * 
     * Standard GPU-style interpolation:
     * E(x,y,z) = Σ w_ijk * E[i,j,k]
     * 
     * where w_ijk are trilinear weights based on fractional position.
     */
    Vec3 interpolate_linear(const Vec3& pos) const;
    
    /**
     * @brief Tricubic interpolation (64 samples)
     * 
     * Uses Catmull-Rom splines for C¹ continuity:
     * - 4×4×4 stencil around point
     * - Cubic polynomial interpolation in each dimension
     * - Guarantees smooth first derivatives
     */
    Vec3 interpolate_cubic(const Vec3& pos) const;
    
    /**
     * @brief Triquintic interpolation (216 samples)
     * 
     * Uses 5th-order polynomials for C² continuity:
     * - 6×6×6 stencil around point
     * - Guarantees smooth second derivatives
     * - Best accuracy, highest cost
     */
    Vec3 interpolate_quintic(const Vec3& pos) const;
    
    /**
     * @brief Convert position to fractional grid coordinates
     * 
     * @param pos Position in global coordinates [m]
     * @param[out] i, j, k Integer cell indices
     * @param[out] fx, fy, fz Fractional positions within cell [0, 1]
     * @return true if inside grid, false if outside
     */
    bool pos_to_grid(const Vec3& pos, int& i, int& j, int& k,
                     double& fx, double& fy, double& fz) const;
    
    /**
     * @brief Get field component at grid point (bounds-checked)
     * 
     * @param field Field component grid (Ex, Ey, or Ez)
     * @param i, j, k Grid indices
     * @return Field value, or 0.0 if out of bounds
     */
    double get_field_at(const std::vector<double>& field, int i, int j, int k) const;
};

}  // namespace gpu
}  // namespace icarion
