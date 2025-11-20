#pragma once
#include <cmath>
#include "core/physics/fields/physics_math_shared.h"

// Helper used by both host and device code to compute instantaneous AC omega.
// start_Hz and slope_Hz_per_s describe a linear sweep in Hz. If enable_frequency_sweep
// is zero, the caller should handle the non-sweep case (e.g., use a precomputed omega).
static inline __host__ __device__ double compute_ac_omega_from_hz(double start_Hz,
                                                                 double slope_Hz_per_s,
                                                                 double t_rel_s,
                                                                 int enable_frequency_sweep)
{
    if (!enable_frequency_sweep) {
        // Caller may prefer to use an explicit angular frequency; returning startHz*2pi
        // keeps the behaviour predictable when start_Hz is set.
        return TWO_PI() * start_Hz;
    }
    double cur = start_Hz + slope_Hz_per_s * t_rel_s;
    if (cur < 0.0) cur = 0.0;
    return TWO_PI() * cur;
}
