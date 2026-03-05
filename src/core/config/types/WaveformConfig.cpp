// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "WaveformConfig.h"
#include <stdexcept>
#include <algorithm>

namespace ICARION::config {

double ArbitraryWaveform::evaluate(double t) const {
    if (times_s.empty() || values.empty()) {
        throw std::runtime_error("ArbitraryWaveform: empty time/value arrays");
    }
    
    if (times_s.size() != values.size()) {
        throw std::runtime_error("ArbitraryWaveform: times and values must have same length");
    }
    
    // Before first point
    if (t <= times_s.front()) {
        return values.front();
    }
    
    // After last point
    if (t >= times_s.back()) {
        return values.back();
    }
    
    // Find interval [t_i, t_{i+1}] containing t
    auto it = std::lower_bound(times_s.begin(), times_s.end(), t);
    size_t idx = std::distance(times_s.begin(), it);
    
    if (idx == 0) {
        return values[0];
    }
    
    size_t i0 = idx - 1;
    size_t i1 = idx;
    
    double t0 = times_s[i0];
    double t1 = times_s[i1];
    double v0 = values[i0];
    double v1 = values[i1];
    
    switch (interp) {
        case Interpolation::Step:
            return v0;
            
        case Interpolation::Linear: {
            double alpha = (t - t0) / (t1 - t0);
            return v0 + alpha * (v1 - v0);
        }
            
        case Interpolation::Cubic: {
            // Cubic Hermite interpolation
            // For simplicity, use Catmull-Rom spline (finite differences for derivatives)
            double alpha = (t - t0) / (t1 - t0);
            
            // Get adjacent points for derivative estimates
            double m0, m1;
            
            if (i0 == 0) {
                m0 = (v1 - v0) / (t1 - t0);
            } else {
                double t_prev = times_s[i0 - 1];
                double v_prev = values[i0 - 1];
                m0 = (v1 - v_prev) / (t1 - t_prev);
            }
            
            if (i1 == times_s.size() - 1) {
                m1 = (v1 - v0) / (t1 - t0);
            } else {
                double t_next = times_s[i1 + 1];
                double v_next = values[i1 + 1];
                m1 = (v_next - v0) / (t_next - t0);
            }
            
            // Hermite basis functions
            double h00 = 2*alpha*alpha*alpha - 3*alpha*alpha + 1;
            double h10 = alpha*alpha*alpha - 2*alpha*alpha + alpha;
            double h01 = -2*alpha*alpha*alpha + 3*alpha*alpha;
            double h11 = alpha*alpha*alpha - alpha*alpha;
            
            double dt = t1 - t0;
            return h00 * v0 + h10 * dt * m0 + h01 * v1 + h11 * dt * m1;
        }
            
        default:
            throw std::runtime_error("ArbitraryWaveform: unknown interpolation mode");
    }
}

double ValueOrWaveform::evaluate(double t, const std::map<std::string, Waveform>& library) const {
    if (constant_value.has_value()) {
        return *constant_value;
    } else if (waveform.has_value()) {
        return waveform->evaluate(t);
    } else if (waveform_ref.has_value()) {
        auto it = library.find(*waveform_ref);
        if (it == library.end()) {
            throw std::runtime_error("Waveform reference not found: " + *waveform_ref);
        }
        return it->second.evaluate(t);
    } else {
        throw std::runtime_error("ValueOrWaveform: no value set (invalid state)");
    }
}

} // namespace ICARION::config
