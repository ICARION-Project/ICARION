// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
//
// Experimental: Only linear interpolation is implemented/validated in v1.0.0;
// higher-order modes fall back to linear.

#include "AdaptiveFieldInterpolator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace icarion {
namespace gpu {

// Constructor with custom config
AdaptiveFieldInterpolator::AdaptiveFieldInterpolator(
    const std::vector<double>& Ex_grid,
    const std::vector<double>& Ey_grid,
    const std::vector<double>& Ez_grid,
    int nx, int ny, int nz,
    Vec3 origin, Vec3 spacing,
    const Config& config
) : Ex_(Ex_grid), Ey_(Ey_grid), Ez_(Ez_grid),
    nx_(nx), ny_(ny), nz_(nz),
    origin_(origin), spacing_(spacing),
    config_(config)
{
    // Validate inputs
    size_t expected_size = static_cast<size_t>(nx) * ny * nz;
    if (Ex_.size() != expected_size || Ey_.size() != expected_size || Ez_.size() != expected_size) {
        throw std::invalid_argument(
            "AdaptiveFieldInterpolator: Grid size mismatch. "
            "Expected " + std::to_string(expected_size) + " elements, got " +
            std::to_string(Ex_.size()) + " (Ex), " + std::to_string(Ey_.size()) + " (Ey), " +
            std::to_string(Ez_.size()) + " (Ez)"
        );
    }
    
    if (nx < 4 || ny < 4 || nz < 4) {
        throw std::invalid_argument(
            "AdaptiveFieldInterpolator: Grid too small for adaptive interpolation. "
            "Minimum: 4×4×4, got " + std::to_string(nx) + "×" + std::to_string(ny) + "×" + std::to_string(nz)
        );
    }
    
    // Precompute inverse spacing for fast division
    inv_spacing_ = Vec3{1.0 / spacing.x, 1.0 / spacing.y, 1.0 / spacing.z};
    
    // Compute gradients if enabled
    if (config_.cache_gradients) {
        compute_gradients();
    }
}

// Constructor with default config
AdaptiveFieldInterpolator::AdaptiveFieldInterpolator(
    const std::vector<double>& Ex_grid,
    const std::vector<double>& Ey_grid,
    const std::vector<double>& Ez_grid,
    int nx, int ny, int nz,
    Vec3 origin, Vec3 spacing
) : AdaptiveFieldInterpolator(Ex_grid, Ey_grid, Ez_grid, nx, ny, nz, origin, spacing,
                               Config{1e3, 1e5, InterpolationOrder::Cubic, true, false})
{}

void AdaptiveFieldInterpolator::compute_gradients() {
    // Allocate gradient storage
    size_t grid_size = static_cast<size_t>(nx_) * ny_ * nz_;
    grad_Ex_x_.resize(grid_size);
    grad_Ex_y_.resize(grid_size);
    grad_Ex_z_.resize(grid_size);
    grad_Ey_x_.resize(grid_size);
    grad_Ey_y_.resize(grid_size);
    grad_Ey_z_.resize(grid_size);
    grad_Ez_x_.resize(grid_size);
    grad_Ez_y_.resize(grid_size);
    grad_Ez_z_.resize(grid_size);
    
    // Compute gradients via central differences
    for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
            for (int k = 0; k < nz_; ++k) {
                int idx = i * (ny_ * nz_) + j * nz_ + k;
                
                // ∂Ex/∂x
                if (i == 0) {
                    // Forward difference
                    double Ex_p1 = get_field_at(Ex_, i+1, j, k);
                    double Ex_0 = get_field_at(Ex_, i, j, k);
                    grad_Ex_x_[idx] = (Ex_p1 - Ex_0) / spacing_.x;
                } else if (i == nx_-1) {
                    // Backward difference
                    double Ex_0 = get_field_at(Ex_, i, j, k);
                    double Ex_m1 = get_field_at(Ex_, i-1, j, k);
                    grad_Ex_x_[idx] = (Ex_0 - Ex_m1) / spacing_.x;
                } else {
                    // Central difference
                    double Ex_p1 = get_field_at(Ex_, i+1, j, k);
                    double Ex_m1 = get_field_at(Ex_, i-1, j, k);
                    grad_Ex_x_[idx] = (Ex_p1 - Ex_m1) / (2.0 * spacing_.x);
                }
                
                // ∂Ex/∂y (similar logic)
                if (j == 0) {
                    double Ex_p1 = get_field_at(Ex_, i, j+1, k);
                    double Ex_0 = get_field_at(Ex_, i, j, k);
                    grad_Ex_y_[idx] = (Ex_p1 - Ex_0) / spacing_.y;
                } else if (j == ny_-1) {
                    double Ex_0 = get_field_at(Ex_, i, j, k);
                    double Ex_m1 = get_field_at(Ex_, i, j-1, k);
                    grad_Ex_y_[idx] = (Ex_0 - Ex_m1) / spacing_.y;
                } else {
                    double Ex_p1 = get_field_at(Ex_, i, j+1, k);
                    double Ex_m1 = get_field_at(Ex_, i, j-1, k);
                    grad_Ex_y_[idx] = (Ex_p1 - Ex_m1) / (2.0 * spacing_.y);
                }
                
                // ∂Ex/∂z
                if (k == 0) {
                    double Ex_p1 = get_field_at(Ex_, i, j, k+1);
                    double Ex_0 = get_field_at(Ex_, i, j, k);
                    grad_Ex_z_[idx] = (Ex_p1 - Ex_0) / spacing_.z;
                } else if (k == nz_-1) {
                    double Ex_0 = get_field_at(Ex_, i, j, k);
                    double Ex_m1 = get_field_at(Ex_, i, j, k-1);
                    grad_Ex_z_[idx] = (Ex_0 - Ex_m1) / spacing_.z;
                } else {
                    double Ex_p1 = get_field_at(Ex_, i, j, k+1);
                    double Ex_m1 = get_field_at(Ex_, i, j, k-1);
                    grad_Ex_z_[idx] = (Ex_p1 - Ex_m1) / (2.0 * spacing_.z);
                }
                
                // Repeat for Ey gradients
                // ∂Ey/∂x
                if (i == 0) {
                    grad_Ey_x_[idx] = (get_field_at(Ey_, i+1, j, k) - get_field_at(Ey_, i, j, k)) / spacing_.x;
                } else if (i == nx_-1) {
                    grad_Ey_x_[idx] = (get_field_at(Ey_, i, j, k) - get_field_at(Ey_, i-1, j, k)) / spacing_.x;
                } else {
                    grad_Ey_x_[idx] = (get_field_at(Ey_, i+1, j, k) - get_field_at(Ey_, i-1, j, k)) / (2.0 * spacing_.x);
                }
                
                // ∂Ey/∂y
                if (j == 0) {
                    grad_Ey_y_[idx] = (get_field_at(Ey_, i, j+1, k) - get_field_at(Ey_, i, j, k)) / spacing_.y;
                } else if (j == ny_-1) {
                    grad_Ey_y_[idx] = (get_field_at(Ey_, i, j, k) - get_field_at(Ey_, i, j-1, k)) / spacing_.y;
                } else {
                    grad_Ey_y_[idx] = (get_field_at(Ey_, i, j+1, k) - get_field_at(Ey_, i, j-1, k)) / (2.0 * spacing_.y);
                }
                
                // ∂Ey/∂z
                if (k == 0) {
                    grad_Ey_z_[idx] = (get_field_at(Ey_, i, j, k+1) - get_field_at(Ey_, i, j, k)) / spacing_.z;
                } else if (k == nz_-1) {
                    grad_Ey_z_[idx] = (get_field_at(Ey_, i, j, k) - get_field_at(Ey_, i, j, k-1)) / spacing_.z;
                } else {
                    grad_Ey_z_[idx] = (get_field_at(Ey_, i, j, k+1) - get_field_at(Ey_, i, j, k-1)) / (2.0 * spacing_.z);
                }
                
                // Repeat for Ez gradients
                // ∂Ez/∂x
                if (i == 0) {
                    grad_Ez_x_[idx] = (get_field_at(Ez_, i+1, j, k) - get_field_at(Ez_, i, j, k)) / spacing_.x;
                } else if (i == nx_-1) {
                    grad_Ez_x_[idx] = (get_field_at(Ez_, i, j, k) - get_field_at(Ez_, i-1, j, k)) / spacing_.x;
                } else {
                    grad_Ez_x_[idx] = (get_field_at(Ez_, i+1, j, k) - get_field_at(Ez_, i-1, j, k)) / (2.0 * spacing_.x);
                }
                
                // ∂Ez/∂y
                if (j == 0) {
                    grad_Ez_y_[idx] = (get_field_at(Ez_, i, j+1, k) - get_field_at(Ez_, i, j, k)) / spacing_.y;
                } else if (j == ny_-1) {
                    grad_Ez_y_[idx] = (get_field_at(Ez_, i, j, k) - get_field_at(Ez_, i, j-1, k)) / spacing_.y;
                } else {
                    grad_Ez_y_[idx] = (get_field_at(Ez_, i, j+1, k) - get_field_at(Ez_, i, j-1, k)) / (2.0 * spacing_.y);
                }
                
                // ∂Ez/∂z
                if (k == 0) {
                    grad_Ez_z_[idx] = (get_field_at(Ez_, i, j, k+1) - get_field_at(Ez_, i, j, k)) / spacing_.z;
                } else if (k == nz_-1) {
                    grad_Ez_z_[idx] = (get_field_at(Ez_, i, j, k) - get_field_at(Ez_, i, j, k-1)) / spacing_.z;
                } else {
                    grad_Ez_z_[idx] = (get_field_at(Ez_, i, j, k+1) - get_field_at(Ez_, i, j, k-1)) / (2.0 * spacing_.z);
                }
            }
        }
    }
    
    has_gradients_ = true;
}

Vec3 AdaptiveFieldInterpolator::interpolate(const Vec3& pos) const {
    // Select interpolation order based on gradient
    InterpolationOrder order = InterpolationOrder::Linear;
    
    if (has_gradients_) {
        double grad_mag = estimate_gradient_magnitude(pos);
        order = select_order(grad_mag);
    }
    
    // Dispatch to appropriate interpolation method
    Vec3 result;
    switch (order) {
        case InterpolationOrder::Linear:
            result = interpolate_linear(pos);
            stats_.n_linear++;
            break;
        case InterpolationOrder::Cubic:
            result = interpolate_cubic(pos);
            stats_.n_cubic++;
            break;
        case InterpolationOrder::Quintic:
            result = interpolate_quintic(pos);
            stats_.n_quintic++;
            break;
    }
    
    return result;
}

AdaptiveFieldInterpolator::InterpolationOrder 
AdaptiveFieldInterpolator::get_order_at(const Vec3& pos) const {
    if (!has_gradients_) {
        return InterpolationOrder::Linear;
    }
    
    double grad_mag = estimate_gradient_magnitude(pos);
    return select_order(grad_mag);
}

double AdaptiveFieldInterpolator::estimate_gradient_magnitude(const Vec3& pos) const {
    // Use trilinear interpolation of cached gradient field
    // |∇E| ≈ √(Σ (∂Ei/∂xj)²)
    
    int i, j, k;
    double fx, fy, fz;
    if (!pos_to_grid(pos, i, j, k, fx, fy, fz)) {
        return 0.0;  // Outside grid → assume zero gradient
    }
    
    // Trilinear weights
    double w000 = (1 - fx) * (1 - fy) * (1 - fz);
    double w100 = fx * (1 - fy) * (1 - fz);
    double w010 = (1 - fx) * fy * (1 - fz);
    double w110 = fx * fy * (1 - fz);
    double w001 = (1 - fx) * (1 - fy) * fz;
    double w101 = fx * (1 - fy) * fz;
    double w011 = (1 - fx) * fy * fz;
    double w111 = fx * fy * fz;
    
    // Interpolate all 9 gradient components
    auto interp_component = [&](const std::vector<double>& grad) -> double {
        int idx000 = i * (ny_ * nz_) + j * nz_ + k;
        int idx100 = (i+1) * (ny_ * nz_) + j * nz_ + k;
        int idx010 = i * (ny_ * nz_) + (j+1) * nz_ + k;
        int idx110 = (i+1) * (ny_ * nz_) + (j+1) * nz_ + k;
        int idx001 = i * (ny_ * nz_) + j * nz_ + (k+1);
        int idx101 = (i+1) * (ny_ * nz_) + j * nz_ + (k+1);
        int idx011 = i * (ny_ * nz_) + (j+1) * nz_ + (k+1);
        int idx111 = (i+1) * (ny_ * nz_) + (j+1) * nz_ + (k+1);
        
        return w000 * grad[idx000] + w100 * grad[idx100] +
               w010 * grad[idx010] + w110 * grad[idx110] +
               w001 * grad[idx001] + w101 * grad[idx101] +
               w011 * grad[idx011] + w111 * grad[idx111];
    };
    
    double dEx_dx = interp_component(grad_Ex_x_);
    double dEx_dy = interp_component(grad_Ex_y_);
    double dEx_dz = interp_component(grad_Ex_z_);
    double dEy_dx = interp_component(grad_Ey_x_);
    double dEy_dy = interp_component(grad_Ey_y_);
    double dEy_dz = interp_component(grad_Ey_z_);
    double dEz_dx = interp_component(grad_Ez_x_);
    double dEz_dy = interp_component(grad_Ez_y_);
    double dEz_dz = interp_component(grad_Ez_z_);
    
    // Frobenius norm of gradient tensor
    double grad_sq = dEx_dx*dEx_dx + dEx_dy*dEx_dy + dEx_dz*dEx_dz +
                     dEy_dx*dEy_dx + dEy_dy*dEy_dy + dEy_dz*dEy_dz +
                     dEz_dx*dEz_dx + dEz_dy*dEz_dy + dEz_dz*dEz_dz;
    
    return std::sqrt(grad_sq);
}

AdaptiveFieldInterpolator::InterpolationOrder 
AdaptiveFieldInterpolator::select_order(double grad_magnitude) const {
    if (grad_magnitude < config_.gradient_threshold_low) {
        return InterpolationOrder::Linear;
    } else if (grad_magnitude > config_.gradient_threshold_high) {
        return config_.max_order;
    } else {
        // Intermediate gradient → use Cubic (good balance)
        return InterpolationOrder::Cubic;
    }
}

Vec3 AdaptiveFieldInterpolator::interpolate_linear(const Vec3& pos) const {
    int i, j, k;
    double fx, fy, fz;
    if (!pos_to_grid(pos, i, j, k, fx, fy, fz)) {
        return Vec3{0.0, 0.0, 0.0};  // Outside grid
    }
    
    // Trilinear weights
    double w000 = (1 - fx) * (1 - fy) * (1 - fz);
    double w100 = fx * (1 - fy) * (1 - fz);
    double w010 = (1 - fx) * fy * (1 - fz);
    double w110 = fx * fy * (1 - fz);
    double w001 = (1 - fx) * (1 - fy) * fz;
    double w101 = fx * (1 - fy) * fz;
    double w011 = (1 - fx) * fy * fz;
    double w111 = fx * fy * fz;
    
    // Interpolate each component
    double Ex = w000 * get_field_at(Ex_, i, j, k) + w100 * get_field_at(Ex_, i+1, j, k) +
                w010 * get_field_at(Ex_, i, j+1, k) + w110 * get_field_at(Ex_, i+1, j+1, k) +
                w001 * get_field_at(Ex_, i, j, k+1) + w101 * get_field_at(Ex_, i+1, j, k+1) +
                w011 * get_field_at(Ex_, i, j+1, k+1) + w111 * get_field_at(Ex_, i+1, j+1, k+1);
    
    double Ey = w000 * get_field_at(Ey_, i, j, k) + w100 * get_field_at(Ey_, i+1, j, k) +
                w010 * get_field_at(Ey_, i, j+1, k) + w110 * get_field_at(Ey_, i+1, j+1, k) +
                w001 * get_field_at(Ey_, i, j, k+1) + w101 * get_field_at(Ey_, i+1, j, k+1) +
                w011 * get_field_at(Ey_, i, j+1, k+1) + w111 * get_field_at(Ey_, i+1, j+1, k+1);
    
    double Ez = w000 * get_field_at(Ez_, i, j, k) + w100 * get_field_at(Ez_, i+1, j, k) +
                w010 * get_field_at(Ez_, i, j+1, k) + w110 * get_field_at(Ez_, i+1, j+1, k) +
                w001 * get_field_at(Ez_, i, j, k+1) + w101 * get_field_at(Ez_, i+1, j, k+1) +
                w011 * get_field_at(Ez_, i, j+1, k+1) + w111 * get_field_at(Ez_, i+1, j+1, k+1);
    
    return Vec3{Ex, Ey, Ez};
}

Vec3 AdaptiveFieldInterpolator::interpolate_cubic(const Vec3& pos) const {
    // Tricubic interpolation using Catmull-Rom splines
    // TODO(v1.1): Full implementation (complex - use library or future phase)
    // For now, fall back to linear
    return interpolate_linear(pos);
}

Vec3 AdaptiveFieldInterpolator::interpolate_quintic(const Vec3& pos) const {
    // Triquintic interpolation
    // TODO(v1.1): Full implementation (very complex - future phase)
    // For now, fall back to linear
    return interpolate_linear(pos);
}

bool AdaptiveFieldInterpolator::pos_to_grid(
    const Vec3& pos, int& i, int& j, int& k,
    double& fx, double& fy, double& fz
) const {
    // Convert position to grid coordinates
    Vec3 rel_pos = pos - origin_;
    double gx = rel_pos.x * inv_spacing_.x;
    double gy = rel_pos.y * inv_spacing_.y;
    double gz = rel_pos.z * inv_spacing_.z;
    
    // Check bounds
    if (gx < 0 || gx >= nx_-1 || gy < 0 || gy >= ny_-1 || gz < 0 || gz >= nz_-1) {
        return false;
    }
    
    // Integer and fractional parts
    i = static_cast<int>(gx);
    j = static_cast<int>(gy);
    k = static_cast<int>(gz);
    
    fx = gx - i;
    fy = gy - j;
    fz = gz - k;
    
    return true;
}

double AdaptiveFieldInterpolator::get_field_at(
    const std::vector<double>& field, int i, int j, int k
) const {
    // Bounds check
    if (i < 0 || i >= nx_ || j < 0 || j >= ny_ || k < 0 || k >= nz_) {
        return 0.0;
    }
    
    int idx = i * (ny_ * nz_) + j * nz_ + k;
    return field[idx];
}

}  // namespace gpu
}  // namespace icarion
