// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "DeepCollisionDiagnosticsTracker.h"

#include "core/config/conversion/EnumMapper.h"
#include "core/io/hdf5Writer.h"
#include "utils/constants.h"
#include <algorithm>
#include <cmath>

namespace ICARION::integrator {

DeepCollisionDiagnosticsTracker::DeepCollisionDiagnosticsTracker(const config::OutputConfig& output) {
    configure(output);
}

void DeepCollisionDiagnosticsTracker::configure(const config::OutputConfig& output) {
    mode_ = config::EnumMapper::deep_analysis_mode_to_string(output.deep_analysis.mode_type);
    domain_filter_idx_ = output.deep_analysis.domain_filter_index;
    sample_every_n_ = std::max<size_t>(1, output.deep_analysis.sample_every_n);
    max_events_per_ion_ = output.deep_analysis.max_events_per_ion;

    enabled_ = (mode_ != "off");
    full_events_ = (mode_ == "full_events");
    store_events_ = (mode_ == "sampled_events" || mode_ == "full_events");
}

void DeepCollisionDiagnosticsTracker::reset(size_t n_ions) {
    if (!enabled_) {
        collisions_total_.clear();
        cooling_axial_count_.clear();
        heating_axial_count_.clear();
        sum_delta_px_kgms_.clear();
        sum_delta_px2_kg2m2s2_.clear();
        sum_delta_py_kgms_.clear();
        sum_delta_py2_kg2m2s2_.clear();
        sum_delta_pz_kgms_.clear();
        sum_delta_pz2_kg2m2s2_.clear();
        sum_delta_ke_eV_.clear();
        sampled_events_per_ion_.clear();
        event_time_s_.clear();
        event_ion_index_.clear();
        event_domain_index_.clear();
        event_delta_px_kgms_.clear();
        event_delta_py_kgms_.clear();
        event_delta_pz_kgms_.clear();
        event_delta_ke_eV_.clear();
        event_vx_before_ms_.clear();
        event_vx_after_ms_.clear();
        event_vy_before_ms_.clear();
        event_vy_after_ms_.clear();
        event_vz_before_ms_.clear();
        event_vz_after_ms_.clear();
        event_v_rel_before_ms_.clear();
        event_sigma_mt_m2_.clear();
        event_radius_m_.clear();
        event_vr_before_ms_.clear();
        event_ejected_flag_.clear();
        return;
    }

    collisions_total_.assign(n_ions, 0);
    cooling_axial_count_.assign(n_ions, 0);
    heating_axial_count_.assign(n_ions, 0);
    sum_delta_px_kgms_.assign(n_ions, 0.0);
    sum_delta_px2_kg2m2s2_.assign(n_ions, 0.0);
    sum_delta_py_kgms_.assign(n_ions, 0.0);
    sum_delta_py2_kg2m2s2_.assign(n_ions, 0.0);
    sum_delta_pz_kgms_.assign(n_ions, 0.0);
    sum_delta_pz2_kg2m2s2_.assign(n_ions, 0.0);
    sum_delta_ke_eV_.assign(n_ions, 0.0);

    sampled_events_per_ion_.assign(n_ions, 0);

    event_time_s_.clear();
    event_ion_index_.clear();
    event_domain_index_.clear();
    event_delta_px_kgms_.clear();
    event_delta_py_kgms_.clear();
    event_delta_pz_kgms_.clear();
    event_delta_ke_eV_.clear();
    event_vx_before_ms_.clear();
    event_vx_after_ms_.clear();
    event_vy_before_ms_.clear();
    event_vy_after_ms_.clear();
    event_vz_before_ms_.clear();
    event_vz_after_ms_.clear();
    event_v_rel_before_ms_.clear();
    event_sigma_mt_m2_.clear();
    event_radius_m_.clear();
    event_vr_before_ms_.clear();
    event_ejected_flag_.clear();
}

void DeepCollisionDiagnosticsTracker::note_collision(size_t ion_idx,
                                                     int domain_idx,
                                                     bool collided,
                                                     double collision_time_s,
                                                     double ion_mass_kg,
                                                     const Vec3& pos_m,
                                                     double v_rel_before_m_s,
                                                     double sigma_mt_m2,
                                                     const Vec3& vel_before_m_s,
                                                     const Vec3& vel_after_m_s,
                                                     bool ion_ejected) {
    if (!enabled_ || !collided) {
        return;
    }
    if (domain_filter_idx_ >= 0 && domain_idx != domain_filter_idx_) {
        return;
    }

    // Ion momentum change from this collision. Positive delta_pz accelerates the ion
    // along +z; axial cooling/heating below is classified from that signed component.
    const double delta_px = ion_mass_kg * (vel_after_m_s.x - vel_before_m_s.x);
    const double delta_py = ion_mass_kg * (vel_after_m_s.y - vel_before_m_s.y);
    const double delta_pz = ion_mass_kg * (vel_after_m_s.z - vel_before_m_s.z);
    const double v2_before = vel_before_m_s.x * vel_before_m_s.x +
                             vel_before_m_s.y * vel_before_m_s.y +
                             vel_before_m_s.z * vel_before_m_s.z;
    const double v2_after = vel_after_m_s.x * vel_after_m_s.x +
                            vel_after_m_s.y * vel_after_m_s.y +
                            vel_after_m_s.z * vel_after_m_s.z;
    const double delta_ke_eV = 0.5 * ion_mass_kg * (v2_after - v2_before) / ELEM_CHARGE_C;
    // Coordinate convention: z is axial, x/y span the transverse radial plane.
    const double radius_m = std::sqrt(pos_m.x * pos_m.x + pos_m.y * pos_m.y);
    double radial_velocity_before_ms = 0.0;
    if (radius_m > 1e-15) {
        radial_velocity_before_ms = (pos_m.x * vel_before_m_s.x + pos_m.y * vel_before_m_s.y) / radius_m;
    }

    collisions_total_[ion_idx] += 1;
    if (delta_pz < 0.0) {
        cooling_axial_count_[ion_idx] += 1;
    } else if (delta_pz > 0.0) {
        heating_axial_count_[ion_idx] += 1;
    }
    sum_delta_px_kgms_[ion_idx] += delta_px;
    sum_delta_px2_kg2m2s2_[ion_idx] += delta_px * delta_px;
    sum_delta_py_kgms_[ion_idx] += delta_py;
    sum_delta_py2_kg2m2s2_[ion_idx] += delta_py * delta_py;
    sum_delta_pz_kgms_[ion_idx] += delta_pz;
    sum_delta_pz2_kg2m2s2_[ion_idx] += delta_pz * delta_pz;
    sum_delta_ke_eV_[ion_idx] += delta_ke_eV;

    if (!store_events_) {
        return;
    }

    bool keep_event = false;
    if (full_events_) {
        keep_event = true;
    } else {
        const int32_t c = collisions_total_[ion_idx];
        keep_event = (c > 0) && (static_cast<size_t>(c - 1) % sample_every_n_ == 0);
    }

    if (keep_event && max_events_per_ion_ > 0 && sampled_events_per_ion_[ion_idx] >= max_events_per_ion_) {
        keep_event = false;
    }

    if (!keep_event) {
        return;
    }

    std::lock_guard<std::mutex> lock(events_mutex_);
    sampled_events_per_ion_[ion_idx] += 1;
    event_time_s_.push_back(collision_time_s);
    event_ion_index_.push_back(static_cast<uint32_t>(ion_idx));
    event_domain_index_.push_back(domain_idx);
    event_delta_px_kgms_.push_back(delta_px);
    event_delta_py_kgms_.push_back(delta_py);
    event_delta_pz_kgms_.push_back(delta_pz);
    event_delta_ke_eV_.push_back(delta_ke_eV);
    event_vx_before_ms_.push_back(vel_before_m_s.x);
    event_vx_after_ms_.push_back(vel_after_m_s.x);
    event_vy_before_ms_.push_back(vel_before_m_s.y);
    event_vy_after_ms_.push_back(vel_after_m_s.y);
    event_vz_before_ms_.push_back(vel_before_m_s.z);
    event_vz_after_ms_.push_back(vel_after_m_s.z);
    event_v_rel_before_ms_.push_back(v_rel_before_m_s);
    event_sigma_mt_m2_.push_back(sigma_mt_m2);
    event_radius_m_.push_back(radius_m);
    event_vr_before_ms_.push_back(radial_velocity_before_ms);
    event_ejected_flag_.push_back(ion_ejected ? 1u : 0u);
}

void DeepCollisionDiagnosticsTracker::write_hdf5(const std::string& hdf5_filename) const {
    if (!enabled_) {
        return;
    }

    std::vector<double> mean_delta_px(collisions_total_.size(), 0.0);
    std::vector<double> mean_delta_py(collisions_total_.size(), 0.0);
    std::vector<double> mean_delta_pz(collisions_total_.size(), 0.0);
    std::vector<double> mean_delta_ke(collisions_total_.size(), 0.0);

    for (size_t i = 0; i < collisions_total_.size(); ++i) {
        if (collisions_total_[i] <= 0) {
            continue;
        }
        const double inv_n = 1.0 / static_cast<double>(collisions_total_[i]);
        mean_delta_px[i] = sum_delta_px_kgms_[i] * inv_n;
        mean_delta_py[i] = sum_delta_py_kgms_[i] * inv_n;
        mean_delta_pz[i] = sum_delta_pz_kgms_[i] * inv_n;
        mean_delta_ke[i] = sum_delta_ke_eV_[i] * inv_n;
    }

    io::HDF5Writer::write_deep_collision_diagnostics(
        hdf5_filename,
        mode_,
        domain_filter_idx_,
        static_cast<int>(sample_every_n_),
        static_cast<int>(max_events_per_ion_),
        collisions_total_,
        cooling_axial_count_,
        heating_axial_count_,
        sum_delta_px_kgms_,
        sum_delta_px2_kg2m2s2_,
        sum_delta_py_kgms_,
        sum_delta_py2_kg2m2s2_,
        sum_delta_pz_kgms_,
        sum_delta_pz2_kg2m2s2_,
        sum_delta_ke_eV_,
        mean_delta_px,
        mean_delta_py,
        mean_delta_pz,
        mean_delta_ke,
        event_time_s_,
        event_ion_index_,
        event_domain_index_,
        event_delta_px_kgms_,
        event_delta_py_kgms_,
        event_delta_pz_kgms_,
        event_delta_ke_eV_,
        event_vx_before_ms_,
        event_vx_after_ms_,
        event_vy_before_ms_,
        event_vy_after_ms_,
        event_vz_before_ms_,
        event_vz_after_ms_,
        event_v_rel_before_ms_,
        event_sigma_mt_m2_,
        event_radius_m_,
        event_vr_before_ms_,
        event_ejected_flag_);
}

}  // namespace ICARION::integrator
