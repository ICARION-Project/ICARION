// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "core/config/types/OutputConfig.h"
#include "core/utils/mathUtils.h"
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

namespace ICARION::integrator {

/**
 * @brief Optional deep collision diagnostics for mechanistic transport studies.
 *
 * Tracks per-ion collision summary statistics and, optionally, a sampled/full
 * collision event stream with axial momentum transfer diagnostics.
 */
class DeepCollisionDiagnosticsTracker {
public:
    explicit DeepCollisionDiagnosticsTracker(const config::OutputConfig& output);

    void configure(const config::OutputConfig& output);
    void reset(size_t n_ions);

    bool enabled() const noexcept { return enabled_; }

    void note_collision(size_t ion_idx,
                        int domain_idx,
                        bool collided,
                        double collision_time_s,
                        double ion_mass_kg,
                        const Vec3& pos_m,
                        double v_rel_before_m_s,
                        double sigma_mt_m2,
                        const Vec3& vel_before_m_s,
                        const Vec3& vel_after_m_s,
                        bool ion_ejected);

    void write_hdf5(const std::string& hdf5_filename) const;

private:
    bool enabled_ = false;
    bool store_events_ = false;
    bool full_events_ = false;
    std::string mode_ = "off";
    int domain_filter_idx_ = -1;
    size_t sample_every_n_ = 10;
    size_t max_events_per_ion_ = 0;

    std::vector<int32_t> collisions_total_;
    std::vector<int32_t> cooling_axial_count_;
    std::vector<int32_t> heating_axial_count_;
    std::vector<double> sum_delta_px_kgms_;
    std::vector<double> sum_delta_px2_kg2m2s2_;
    std::vector<double> sum_delta_py_kgms_;
    std::vector<double> sum_delta_py2_kg2m2s2_;
    std::vector<double> sum_delta_pz_kgms_;
    std::vector<double> sum_delta_pz2_kg2m2s2_;
    std::vector<double> sum_delta_ke_eV_;

    std::vector<uint32_t> sampled_events_per_ion_;

    std::vector<double> event_time_s_;
    std::vector<uint32_t> event_ion_index_;
    std::vector<int32_t> event_domain_index_;
    std::vector<double> event_delta_px_kgms_;
    std::vector<double> event_delta_py_kgms_;
    std::vector<double> event_delta_pz_kgms_;
    std::vector<double> event_delta_ke_eV_;
    std::vector<double> event_vx_before_ms_;
    std::vector<double> event_vx_after_ms_;
    std::vector<double> event_vy_before_ms_;
    std::vector<double> event_vy_after_ms_;
    std::vector<double> event_vz_before_ms_;
    std::vector<double> event_vz_after_ms_;
    std::vector<double> event_v_rel_before_ms_;
    std::vector<double> event_sigma_mt_m2_;
    std::vector<double> event_radius_m_;
    std::vector<double> event_vr_before_ms_;
    std::vector<uint8_t> event_ejected_flag_;
    mutable std::mutex events_mutex_;
};

}  // namespace ICARION::integrator
