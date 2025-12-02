// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "RK45Strategy.h"
#include "core/physics/forces/ForceContext.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace ICARION {
namespace integrator {

namespace {
    // Dormand-Prince 5(4) Butcher tableau coefficients
    // Reference: Hairer, Norsett, Wanner (1993) "Solving Ordinary Differential Equations I"
    
    // Time nodes for stages
    constexpr double c2 = 1.0/5.0;
    constexpr double c3 = 3.0/10.0;
    constexpr double c4 = 4.0/5.0;
    constexpr double c5 = 8.0/9.0;
    constexpr double c6 = 1.0;
    constexpr double c7 = 1.0;
    
    // Runge-Kutta matrix (a_ij coefficients)
    // Stage 2
    constexpr double a21 = 1.0/5.0;
    
    // Stage 3
    constexpr double a31 = 3.0/40.0;
    constexpr double a32 = 9.0/40.0;
    
    // Stage 4
    constexpr double a41 = 44.0/45.0;
    constexpr double a42 = -56.0/15.0;
    constexpr double a43 = 32.0/9.0;
    
    // Stage 5
    constexpr double a51 = 19372.0/6561.0;
    constexpr double a52 = -25360.0/2187.0;
    constexpr double a53 = 64448.0/6561.0;
    constexpr double a54 = -212.0/729.0;
    
    // Stage 6
    constexpr double a61 = 9017.0/3168.0;
    constexpr double a62 = -355.0/33.0;
    constexpr double a63 = 46732.0/5247.0;
    constexpr double a64 = 49.0/176.0;
    constexpr double a65 = -5103.0/18656.0;
    
    // Stage 7 (FSAL: used for 4th order solution and next step's k1)
    constexpr double a71 = 35.0/384.0;
    constexpr double a72 = 0.0;
    constexpr double a73 = 500.0/1113.0;
    constexpr double a74 = 125.0/192.0;
    constexpr double a75 = -2187.0/6784.0;
    constexpr double a76 = 11.0/84.0;
    
    // 4th order solution weights (b_i coefficients) - SAME as stage 7
    constexpr double b41 = 35.0/384.0;
    constexpr double b42 = 0.0;
    constexpr double b43 = 500.0/1113.0;
    constexpr double b44 = 125.0/192.0;
    constexpr double b45 = -2187.0/6784.0;
    constexpr double b46 = 11.0/84.0;
    constexpr double b47 = 0.0;  // FSAL property
    
    // 5th order solution weights (b*_i coefficients) - for error estimation
    constexpr double b51 = 5179.0/57600.0;
    constexpr double b52 = 0.0;
    constexpr double b53 = 7571.0/16695.0;
    constexpr double b54 = 393.0/640.0;
    constexpr double b55 = -92097.0/339200.0;
    constexpr double b56 = 187.0/2100.0;
    constexpr double b57 = 1.0/40.0;
    
    // PI controller coefficients (optimized for order 5)
    constexpr double PI_BETA = 0.04;  // Stabilization parameter
    constexpr double PI_ALPHA = 0.2 - PI_BETA * 0.75;  // Proportional gain
    constexpr double REJECTION_EXPONENT = 0.2;  // 1/order for rejected steps (1/5 for order 5)
    
    // Step control parameters
    constexpr int MAX_REJECT_ATTEMPTS = 10;
    constexpr double SAFETY_MARGIN = 0.9;
    constexpr double MIN_ERROR_THRESHOLD = 1e-10;  // Below this, error estimate unreliable
    constexpr double ERROR_ACCEPTANCE_THRESHOLD = 1.0;  // Error <= 1.0 means step accepted
    constexpr double DT_MIN_TOLERANCE = 1.001;  // Accept step if dt near dt_min (1.001x factor)
    
    // Spatial dimensions
    constexpr int NUM_SPATIAL_DIMS = 3;  // x, y, z components
    constexpr int X_INDEX = 0;
    constexpr int Y_INDEX = 1;
    constexpr int Z_INDEX = 2;
}

RK45Strategy::RK45Strategy()
    : config_{}, stats_{}, last_error_(1.0), fsal_available_(false)
{
    // Use default configuration
}

RK45Strategy::RK45Strategy(const AdaptiveConfig& config)
    : config_(config), stats_{}, last_error_(1.0), fsal_available_(false)
{
    // Validate configuration
    if (config_.atol <= 0.0 || config_.rtol <= 0.0) {
        throw std::invalid_argument("RK45Strategy: tolerances must be positive");
    }
    if (config_.safety_factor <= 0.0 || config_.safety_factor >= 1.0) {
        throw std::invalid_argument("RK45Strategy: safety factor must be in (0, 1)");
    }
}

void RK45Strategy::compute_acceleration(
    double& ax, double& ay, double& az,
    const IonState& ion,
    double t,
    const physics::ForceRegistry& force_registry,
    const std::vector<IonState>& all_ions
) {
    physics::ForceContext ctx;
    ctx.domain = force_registry.domain();  // Get domain from registry
    ctx.all_ions = &all_ions;
    ctx.field_provider = nullptr;
    
    Vec3 F = force_registry.compute_total_force(ion, t, ctx);
    Vec3 a = F / ion.mass_kg;
    
    ax = a.x;
    ay = a.y;
    az = a.z;
}

double RK45Strategy::estimate_error(
    const IonState& y4,
    const IonState& y5,
    const IonState& y_current
) const {
    // Compute scaled error for position and velocity components
    // error = max_i |e_i| / (atol + rtol * |y_i|)
    
    double max_error = 0.0;
    
    // Position error (x, y, z)
    for (int i = 0; i < NUM_SPATIAL_DIMS; ++i) {
        double y5_val = (i == X_INDEX) ? y5.pos.x : (i == Y_INDEX) ? y5.pos.y : y5.pos.z;
        double y4_val = (i == X_INDEX) ? y4.pos.x : (i == Y_INDEX) ? y4.pos.y : y4.pos.z;
        double y_val = (i == X_INDEX) ? y_current.pos.x : (i == Y_INDEX) ? y_current.pos.y : y_current.pos.z;
        
        double err_abs = std::fabs(y5_val - y4_val);
        double scale = config_.atol + config_.rtol * std::fabs(y_val);
        double err_scaled = err_abs / scale;
        
        max_error = std::max(max_error, err_scaled);
    }
    
    // Velocity error (vx, vy, vz)
    for (int i = 0; i < NUM_SPATIAL_DIMS; ++i) {
        double y5_val = (i == X_INDEX) ? y5.vel.x : (i == Y_INDEX) ? y5.vel.y : y5.vel.z;
        double y4_val = (i == X_INDEX) ? y4.vel.x : (i == Y_INDEX) ? y4.vel.y : y4.vel.z;
        double y_val = (i == X_INDEX) ? y_current.vel.x : (i == Y_INDEX) ? y_current.vel.y : y_current.vel.z;
        
        double err_abs = std::fabs(y5_val - y4_val);
        double scale = config_.atol + config_.rtol * std::fabs(y_val);
        double err_scaled = err_abs / scale;
        
        max_error = std::max(max_error, err_scaled);
    }
    
    // Avoid division by zero in step control
    return std::max(max_error, MIN_ERROR_THRESHOLD);
}

double RK45Strategy::compute_new_step(
    double current_dt,
    double error,
    double dt_min,
    double dt_max
) {
    // PI controller for step size adjustment
    // Gustafsson, K. (1991): "Control theoretic techniques for stepsize selection"
    // dt_new = dt * (error_prev / error)^(beta/order) * (1 / error)^(alpha/order)
    
    double factor;
    
    if (error > 1.0) {
        // Step rejected: use conservative factor
        factor = config_.safety_factor * std::pow(1.0 / error, REJECTION_EXPONENT);
    } else {
        // Step accepted: use PI controller
        factor = config_.safety_factor 
                * std::pow(last_error_ / error, PI_BETA)
                * std::pow(1.0 / error, PI_ALPHA);
    }
    
    // Clamp growth/shrinkage per step
    factor = std::max(factor, config_.max_step_decrease);
    factor = std::min(factor, config_.max_step_increase);
    
    double dt_new = current_dt * factor;
    
    // Clamp to absolute limits
    dt_new = std::max(dt_new, dt_min);
    dt_new = std::min(dt_new, dt_max);
    
    return dt_new;
}

void RK45Strategy::step(
    IonState& ion,
    double t,
    double dt,
    const physics::ForceRegistry& force_registry,
    const std::vector<IonState>& all_ions
) {
    // Non-adaptive interface: use fixed dt
    double dt_variable = dt;
    step_adaptive(ion, t, dt_variable, force_registry, all_ions);
    // Note: dt_variable may change, but we ignore it for fixed-step interface
}

void RK45Strategy::step_adaptive(
    IonState& ion,
    double t,
    double& dt_inout,
    const physics::ForceRegistry& force_registry,
    const std::vector<IonState>& all_ions
) {
    const double dt_initial = dt_inout;
    
    // Domain pointer is retrieved for consistency; current implementation
    // does not use domain fields directly.
    const config::DomainConfig* domain = force_registry.domain();
    if (!domain) {
        throw std::runtime_error("RK45Strategy: ForceRegistry has no domain configured");
    }
    const double dt_min = dt_initial * config_.min_step_factor;
    const double dt_max = dt_initial * config_.max_step_factor;
    
    double dt = dt_inout;
    bool step_accepted = false;
    int attempts = 0;
    
    IonState y0 = ion;  // Save initial state
    
    while (!step_accepted && attempts < MAX_REJECT_ATTEMPTS) {
        // =====================================================================
        // Compute 7 RK stages (k1...k7)
        // For 2nd order ODE: x' = v, v' = a = f(x, v, t) / m
        // Each stage k_i stores (velocity, acceleration) pair
        // =====================================================================
        
        Vec3 k1_v, k1_a;
        Vec3 k2_v, k2_a;
        Vec3 k3_v, k3_a;
        Vec3 k4_v, k4_a;
        Vec3 k5_v, k5_a;
        Vec3 k6_v, k6_a;
        Vec3 k7_v, k7_a;
        
        // Stage 1: k1 = f(t, y)
        k1_v = y0.vel;  // dx/dt = v
        if (fsal_available_) {
            // Reuse k7 from previous step (FSAL property)
            k1_a.x = k1_stored_.ax;
            k1_a.y = k1_stored_.ay;
            k1_a.z = k1_stored_.az;
        } else {
            double ax, ay, az;
            compute_acceleration(ax, ay, az, y0, t, force_registry, all_ions);
            k1_a = Vec3{ax, ay, az};  // dv/dt = a
        }
        
        // Stage 2: k2 = f(t + c2*dt, y + dt*(a21*k1))
        IonState y2 = y0;
        y2.pos += k1_v * (dt * a21);  // x + v*dt*a21
        y2.vel += k1_a * (dt * a21);  // v + a*dt*a21
        k2_v = y2.vel;
        double ax2, ay2, az2;
        compute_acceleration(ax2, ay2, az2, y2, t + c2*dt, force_registry, all_ions);
        k2_a = Vec3{ax2, ay2, az2};
        
        // Stage 3: k3 = f(t + c3*dt, y + dt*(a31*k1 + a32*k2))
        IonState y3 = y0;
        y3.pos += (k1_v * a31 + k2_v * a32) * dt;
        y3.vel += (k1_a * a31 + k2_a * a32) * dt;
        k3_v = y3.vel;
        double ax3, ay3, az3;
        compute_acceleration(ax3, ay3, az3, y3, t + c3*dt, force_registry, all_ions);
        k3_a = Vec3{ax3, ay3, az3};
        
        // Stage 4: k4 = f(t + c4*dt, y + dt*(a41*k1 + a42*k2 + a43*k3))
        IonState y4_temp = y0;
        y4_temp.pos += (k1_v * a41 + k2_v * a42 + k3_v * a43) * dt;
        y4_temp.vel += (k1_a * a41 + k2_a * a42 + k3_a * a43) * dt;
        k4_v = y4_temp.vel;
        double ax4, ay4, az4;
        compute_acceleration(ax4, ay4, az4, y4_temp, t + c4*dt, force_registry, all_ions);
        k4_a = Vec3{ax4, ay4, az4};
        
        // Stage 5: k5 = f(t + c5*dt, y + dt*(a51*k1 + ... + a54*k4))
        IonState y5_temp = y0;
        y5_temp.pos += (k1_v * a51 + k2_v * a52 + k3_v * a53 + k4_v * a54) * dt;
        y5_temp.vel += (k1_a * a51 + k2_a * a52 + k3_a * a53 + k4_a * a54) * dt;
        k5_v = y5_temp.vel;
        double ax5, ay5, az5;
        compute_acceleration(ax5, ay5, az5, y5_temp, t + c5*dt, force_registry, all_ions);
        k5_a = Vec3{ax5, ay5, az5};
        
        // Stage 6: k6 = f(t + c6*dt, y + dt*(a61*k1 + ... + a65*k5))
        IonState y6 = y0;
        y6.pos += (k1_v * a61 + k2_v * a62 + k3_v * a63 + k4_v * a64 + k5_v * a65) * dt;
        y6.vel += (k1_a * a61 + k2_a * a62 + k3_a * a63 + k4_a * a64 + k5_a * a65) * dt;
        k6_v = y6.vel;
        double ax6, ay6, az6;
        compute_acceleration(ax6, ay6, az6, y6, t + c6*dt, force_registry, all_ions);
        k6_a = Vec3{ax6, ay6, az6};
        
        // Stage 7: k7 = f(t + dt, y + dt*(a71*k1 + ... + a76*k6))
        // This becomes the 4th-order solution AND next step's k1 (FSAL)
        IonState y7 = y0;
        y7.pos += (k1_v * a71 + k2_v * a72 + k3_v * a73 + k4_v * a74 + k5_v * a75 + k6_v * a76) * dt;
        y7.vel += (k1_a * a71 + k2_a * a72 + k3_a * a73 + k4_a * a74 + k5_a * a75 + k6_a * a76) * dt;
        k7_v = y7.vel;
        double ax7, ay7, az7;
        compute_acceleration(ax7, ay7, az7, y7, t + c7*dt, force_registry, all_ions);
        k7_a = Vec3{ax7, ay7, az7};
        
        // =====================================================================
        // Construct 4th and 5th order solutions
        // =====================================================================
        
        // 4th order solution (b4 coefficients) - SAME as y7 due to FSAL
        IonState y4 = y7;
        
        // 5th order solution (b5 coefficients)
        IonState y5 = y0;
        y5.pos += (k1_v * b51 + k2_v * b52 + k3_v * b53 + k4_v * b54 + k5_v * b55 + k6_v * b56 + k7_v * b57) * dt;
        y5.vel += (k1_a * b51 + k2_a * b52 + k3_a * b53 + k4_a * b54 + k5_a * b55 + k6_a * b56 + k7_a * b57) * dt;
        
        // =====================================================================
        // Error estimation and step acceptance
        // =====================================================================
        
        double error = estimate_error(y4, y5, y0);
        
        if (error <= ERROR_ACCEPTANCE_THRESHOLD || dt <= dt_min * DT_MIN_TOLERANCE) {
            // Accept step
            ion = y4;  // Use 4th-order solution
            step_accepted = true;
            
            // Store k7 for next step's k1 (FSAL)
            k1_stored_.ax = k7_a.x;
            k1_stored_.ay = k7_a.y;
            k1_stored_.az = k7_a.z;
            fsal_available_ = true;
            
            // Update statistics
            stats_.accepted_steps++;
            stats_.min_step_used = std::min(stats_.min_step_used, dt);
            stats_.max_step_used = std::max(stats_.max_step_used, dt);
            stats_.avg_error = (stats_.avg_error * (stats_.accepted_steps - 1) + error) 
                             / stats_.accepted_steps;
            
            // Compute new step size for next step
            double dt_next = compute_new_step(dt, error, dt_min, dt_max);
            dt_inout = dt_next;
            
            last_error_ = error;
            
        } else {
            // Reject step, reduce dt, retry
            stats_.rejected_steps++;
            
            double dt_next = compute_new_step(dt, error, dt_min, dt_max);
            dt = dt_next;
            
            attempts++;
            
            // Reset FSAL flag (must recompute k1)
            fsal_available_ = false;
        }
    }
    
    if (!step_accepted) {
        throw std::runtime_error(
            "RK45Strategy: Failed to converge after " 
            + std::to_string(MAX_REJECT_ATTEMPTS) + " attempts"
        );
    }
}

} // namespace integrator
} // namespace ICARION
