// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "InteractionPotentialCollisionHandler.h"
#include "core/physics/collisions/core/VelocitySampling.h"
#include "core/physics/collisions/core/CollisionGeometry.h"
#include "core/utils/mathUtils.h"
#include "core/log/Logger.h"
#include "utils/constants.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>

namespace ICARION::physics {

namespace {
    constexpr double MIN_RELATIVE_VELOCITY = 1e-10;

    struct alignas(64) InteractionPotentialStatsSlot {
        std::uint64_t generation = 0;
        CollisionStats stats{};
        size_t rate_samples = 0;
        double rate_sum = 0.0;
    };

    thread_local const InteractionPotentialCollisionStatsState* interaction_potential_cached_state = nullptr;
    thread_local std::uint64_t interaction_potential_cached_state_id = 0;
    thread_local InteractionPotentialStatsSlot* interaction_potential_cached_slot = nullptr;

    struct InteractionPotentialSlotCacheEntry {
        std::uint64_t state_id = 0;
        InteractionPotentialStatsSlot* slot = nullptr;
    };

    thread_local std::unordered_map<
        const InteractionPotentialCollisionStatsState*,
        InteractionPotentialSlotCacheEntry
    > interaction_potential_slot_cache;
}

struct InteractionPotentialCollisionStatsState {
    InteractionPotentialCollisionStatsState();
    mutable std::mutex mutex;
    const std::uint64_t id;
    std::atomic<std::uint64_t> generation{1};
    std::vector<std::unique_ptr<InteractionPotentialStatsSlot>> slots;
};

InteractionPotentialCollisionStatsState::InteractionPotentialCollisionStatsState()
    : id([]() {
        static std::atomic<std::uint64_t> next_id{1};
        return next_id.fetch_add(1, std::memory_order_relaxed);
    }())
{}

InteractionPotentialCollisionHandler::InteractionPotentialCollisionHandler(
    bool enable_logging,
    const config::SpeciesDatabase* species_db,
    const std::string& orientation_mode,
    int fixed_orientation_index,
    const std::string& vrel_log_prefix,
    const std::string& momentum_log_prefix
) : enable_logging_(enable_logging)
  , species_db_(species_db)
  , stats_state_(std::make_shared<InteractionPotentialCollisionStatsState>())
{
    if (!vrel_log_prefix.empty()) {
        vrel_logging_enabled_ = true;
        vrel_log_prefix_ = vrel_log_prefix;
    }

    std::string mode = orientation_mode;
    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
    if (mode == "fixed") {
        orientation_mode_ = OrientationMode::Fixed;
    } else if (mode == "random") {
        orientation_mode_ = OrientationMode::Random;
    } else {
        throw std::invalid_argument(
            "[InteractionPotential] Unknown orientation mode '" + orientation_mode +
            "' (expected 'random' or 'fixed')"
        );
    }
    fixed_orientation_index_ = (fixed_orientation_index > 0)
        ? static_cast<size_t>(fixed_orientation_index)
        : 0;

    if (!momentum_log_prefix.empty()) {
        momentum_logging_enabled_ = true;
        momentum_log_prefix_ = momentum_log_prefix;
    }

    if (!species_db_) {
        return;
    }

    for (const auto& [id, props] : species_db_->species) {
        if (!props.ipm_samples_file) {
            continue;
        }
        InteractionPotentialOfflineSampleSet samples;
        std::string err;
        if (!load_interaction_potential_offline_sample_set_file(*props.ipm_samples_file, samples, &err)) {
            ICARION::log::Logger::get("collision")->warn(
                "[InteractionPotential] Failed to load offline samples for '{}': {}",
                id, err
            );
            continue;
        }
        if (!samples.valid()) {
            ICARION::log::Logger::get("collision")->warn(
                "[InteractionPotential] Offline samples for '{}' are invalid.",
                id
            );
            continue;
        }
        const bool has_full_cdf =
            !samples.cdf_offsets.empty() &&
            !samples.cdf_counts.empty() &&
            !samples.cdf_values.empty() &&
            !samples.dp_samples.empty();
        if (!has_full_cdf && !samples.dp_stats.empty()) {
            ICARION::log::Logger::get("collision")->warn(
                "[InteractionPotential] Offline samples for '{}' do not contain full CDF momentum kicks; "
                "runtime will use the lower-fidelity dp_stats fallback. Regenerate with --store-full-cdf "
                "for scientific production runs.",
                id
            );
        }
        logged_use_once_.emplace(id, std::make_unique<std::once_flag>());
        samples_.emplace(id, std::move(samples));
    }
}

void InteractionPotentialCollisionHandler::record_attempt(
    double rate,
    bool rejected,
    bool collision
) const {
    auto* state = stats_state_.get();
    InteractionPotentialStatsSlot* slot = interaction_potential_cached_slot;
    if (interaction_potential_cached_state != state || interaction_potential_cached_state_id != state->id) {
        auto cache_it = interaction_potential_slot_cache.find(state);
        if (cache_it != interaction_potential_slot_cache.end() && cache_it->second.state_id == state->id) {
            slot = cache_it->second.slot;
        } else {
            std::lock_guard<std::mutex> lock(state->mutex);
            auto owned_slot = std::make_unique<InteractionPotentialStatsSlot>();
            slot = owned_slot.get();
            state->slots.push_back(std::move(owned_slot));
            interaction_potential_slot_cache[state] = {state->id, slot};
        }
        interaction_potential_cached_state = state;
        interaction_potential_cached_state_id = state->id;
        interaction_potential_cached_slot = slot;
    }

    const std::uint64_t generation = state->generation.load(std::memory_order_acquire);
    if (slot->generation != generation) {
        slot->generation = generation;
        slot->stats = {};
        slot->rate_samples = 0;
        slot->rate_sum = 0.0;
    }
    slot->rate_sum += rate;
    slot->rate_samples++;
    slot->stats.rejected_collisions += rejected ? 1 : 0;
    slot->stats.total_collisions += collision ? 1 : 0;
}

InteractionPotentialCollisionHandler::~InteractionPotentialCollisionHandler() {
    try {
        flush_logs();
    } catch (const std::exception& e) {
        ICARION::log::Logger::get("collision")->warn(
            "[InteractionPotential] Destructor log flush failed: {}", e.what()
        );
    } catch (...) {
        ICARION::log::Logger::get("collision")->warn(
            "[InteractionPotential] Destructor log flush failed with unknown error"
        );
    }
}

void InteractionPotentialCollisionHandler::flush_logs() {
    std::unordered_map<std::string, VrelLog> vrel_logs_snapshot;
    std::unordered_map<std::string, MomentumLog> momentum_logs_snapshot;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        vrel_logs_snapshot = vrel_logs_;
        momentum_logs_snapshot = momentum_logs_;
    }

    if (vrel_logging_enabled_ && !vrel_logs_snapshot.empty()) {
        std::filesystem::path prefix_path(vrel_log_prefix_.empty() ? "analysis/output/ipm_vrel" : vrel_log_prefix_);
        const bool base_is_csv = prefix_path.extension() == ".csv";
        if (base_is_csv) {
            prefix_path = prefix_path.parent_path() / prefix_path.stem();
        }

        for (const auto& [species, log] : vrel_logs_snapshot) {
            if (log.samples == 0 || log.logv_bins.empty()) {
                continue;
            }

            std::filesystem::path out_path;
            if (base_is_csv && vrel_logs_snapshot.size() == 1) {
                out_path = std::filesystem::path(vrel_log_prefix_);
            } else {
                out_path = prefix_path;
                out_path += "_" + species + ".csv";
            }

            if (!out_path.parent_path().empty()) {
                std::filesystem::create_directories(out_path.parent_path());
            }

            std::ofstream ofs(out_path);
            if (!ofs) {
                ICARION::log::Logger::get("collision")->warn(
                    "[InteractionPotential] Failed to write v_rel log to '{}'", out_path.string()
                );
                continue;
            }

            const double v_mean = log.v_sum / static_cast<double>(log.samples);
            const double v_rms = std::sqrt(log.v2_sum / static_cast<double>(log.samples));

            ofs << "# species: " << species << "\n";
            ofs << "# samples: " << log.samples << "\n";
            ofs << "# v_mean_m_s: " << v_mean << "\n";
            ofs << "# v_rms_m_s: " << v_rms << "\n";
            ofs << "logv,v_rel_m_s,count\n";

            for (size_t i = 0; i < log.logv_bins.size(); ++i) {
                const double logv = log.logv_bins[i];
                const double v = std::exp(logv);
                const auto count = (i < log.bin_counts.size()) ? log.bin_counts[i] : 0ULL;
                ofs << logv << "," << v << "," << count << "\n";
            }
        }
    }

    if (!momentum_logging_enabled_ || momentum_logs_snapshot.empty()) {
        return;
    }

    std::filesystem::path prefix_path(momentum_log_prefix_.empty() ? "analysis/output/ipm_momentum" : momentum_log_prefix_);
    const bool base_is_csv = prefix_path.extension() == ".csv";
    if (base_is_csv) {
        prefix_path = prefix_path.parent_path() / prefix_path.stem();
    }

    for (const auto& [species, log] : momentum_logs_snapshot) {
        if (log.attempts == 0 || log.collisions == 0) {
            continue;
        }

        std::filesystem::path out_path;
        if (base_is_csv && momentum_logs_snapshot.size() == 1) {
            out_path = std::filesystem::path(momentum_log_prefix_);
        } else {
            out_path = prefix_path;
            out_path += "_" + species + ".csv";
        }

        if (!out_path.parent_path().empty()) {
            std::filesystem::create_directories(out_path.parent_path());
        }

        std::ofstream ofs(out_path);
        if (!ofs) {
            ICARION::log::Logger::get("collision")->warn(
                "[InteractionPotential] Failed to write momentum log to '{}'", out_path.string()
            );
            continue;
        }

        const double mean_k = log.k_sum / static_cast<double>(log.attempts);
        const double mean_v = log.v_sum / static_cast<double>(log.attempts);
        const double v_rms = std::sqrt(log.v2_sum / static_cast<double>(log.attempts));
        const double mean_dp_par = log.dp_par_sum / static_cast<double>(log.collisions);
        const double force_est = mean_k * mean_dp_par;
        const double e_est = (log.q_c != 0.0) ? (force_est / log.q_c) : 0.0;

        ofs << "species,attempts,collisions,mean_k_s-1,mean_vrel_m_s,vrel_rms_m_s,mean_dp_par_SI,force_est_N,E_est_V_per_m,ion_charge_C\n";
        ofs << species << "," << log.attempts << "," << log.collisions << ","
            << mean_k << "," << mean_v << "," << v_rms << ","
            << mean_dp_par << "," << force_est << "," << e_est << "," << log.q_c << "\n";
    }
}

bool InteractionPotentialCollisionHandler::handle_collision(
    core::IonCollisionData& view,
    double dt,
    PhysicsRng& rng,
    const config::EnvironmentConfig& env,
    CollisionEventDiagnostics* diagnostics
) {
    const std::string& species_id = view.species_id();
    auto it = samples_.find(species_id);
    if (it == samples_.end()) {
        throw std::runtime_error("[InteractionPotential] Missing offline samples for species '" + species_id + "'");
    }
    const InteractionPotentialOfflineSampleSet& samples = it->second;
    // Constructor validation is sufficient because loaded sample tables are immutable.
    // Revalidating full CDF vectors here would scan every sample on every ion step.

    if (enable_logging_) {
        auto once_it = logged_use_once_.find(species_id);
        if (once_it != logged_use_once_.end()) {
            std::call_once(*once_it->second, [&]() {
                ICARION::log::Logger::get("collision")->info(
                    "[InteractionPotential] Using offline samples for '{}' (gas '{}', N_orient={}, bins={})",
                    species_id,
                    samples.gas.empty() ? env.gas_species : samples.gas,
                    samples.n_orient,
                    samples.n_bins
                );
            });
        }
    }

    const config::GasMixtureComponent* active_component = nullptr;
    if (!env.gas_mixture.empty()) {
        for (const auto& comp : env.gas_mixture) {
            if (!comp.participates_in_collisions) {
                continue;
            }
            if (active_component != nullptr) {
                throw std::runtime_error("[InteractionPotential] Gas mixtures with multiple active components are not supported yet");
            }
            active_component = &comp;
        }
        if (active_component == nullptr) {
            static std::once_flag no_active_mixture_warn_once;
            std::call_once(no_active_mixture_warn_once, []() {
                ICARION::log::Logger::main()->warn(
                    "[InteractionPotential] gas_mixture is configured but has no active collision components; "
                    "collisions are disabled for this domain until mixture flags/densities are fixed.");
            });
            return false;
        }
    }

    const std::string gas_id = (active_component == nullptr) ? env.gas_species : active_component->species;
    if (!samples.gas.empty() && samples.gas != gas_id) {
        throw std::runtime_error(
            "[InteractionPotential] Gas mismatch for species '" + species_id +
            "': samples tagged gas '" + samples.gas +
            "' but environment gas is '" + gas_id + "'"
        );
    }

    const double m_ion = view.kin.get_mass();
    const double m_gas = (active_component == nullptr) ? env.gas_mass_kg : active_component->mass_kg;
    const double n_gas = (active_component == nullptr) ? env.particle_density_m_3 : active_component->density_m3;
    const Vec3 v_gas = env.gas_velocity_m_s;
    const double T = env.temperature_K;

    const Vec3 v_neutral = collision_core::VelocitySampling::sample_neutral_velocity(
        T, m_gas, v_gas, rng
    );
    const Vec3 v_rel = view.kin.vel() - v_neutral;
    const double v_rel_mag = norm(v_rel);
    if (v_rel_mag < MIN_RELATIVE_VELOCITY) {
        return false;
    }

    const double logv = std::log(v_rel_mag);
    const auto& bins = samples.logv_bins;
    if (bins.empty()) {
        return false;
    }
    const double v_min_sample = std::exp(bins.front());
    const double v_max_sample = std::exp(bins.back());
    const double low_bin_step = (bins.size() > 1) ? (bins[1] - bins[0]) : 0.0;
    const double high_bin_step = (bins.size() > 1) ? (bins.back() - bins[bins.size() - 2]) : 0.0;
    const double low_warn_threshold = std::exp(bins.front() - low_bin_step);
    const double high_warn_threshold = std::exp(bins.back() + high_bin_step);

    size_t bin_idx = 0;
    if (logv <= bins.front()) {
        bin_idx = 0;
        bool do_warn = false;
        if (v_rel_mag < low_warn_threshold) {
            std::lock_guard<std::mutex> lock(state_mutex_);
            do_warn = warned_low_coverage_.insert(species_id).second;
        }
        if (do_warn) {
            ICARION::log::Logger::get("collision")->warn(
                "[InteractionPotential] v_rel below sample coverage for species='{}' gas='{}': "
                "v_rel={:.3e} m/s, sampled_range=[{:.3e}, {:.3e}] m/s. "
                "Using lowest sample bin; regenerate offline samples for this T/E_N range or use HSS fallback.",
                species_id,
                gas_id,
                v_rel_mag,
                v_min_sample,
                v_max_sample
            );
        }
    } else if (logv >= bins.back()) {
        bin_idx = bins.size() - 1;
        bool do_warn = false;
        if (v_rel_mag > high_warn_threshold) {
            std::lock_guard<std::mutex> lock(state_mutex_);
            do_warn = warned_high_coverage_.insert(species_id).second;
        }
        if (do_warn) {
            ICARION::log::Logger::get("collision")->warn(
                "[InteractionPotential] v_rel above sample coverage for species='{}' gas='{}': "
                "v_rel={:.3e} m/s, sampled_range=[{:.3e}, {:.3e}] m/s. "
                "Using highest sample bin; regenerate offline samples for this T/E_N range or use HSS fallback.",
                species_id,
                gas_id,
                v_rel_mag,
                v_min_sample,
                v_max_sample
            );
        }
    } else {
        auto itb = std::upper_bound(bins.begin(), bins.end(), logv);
        size_t hi = static_cast<size_t>(std::distance(bins.begin(), itb));
        size_t lo = hi - 1;
        const double t = (logv - bins[lo]) / (bins[hi] - bins[lo]);
        bin_idx = (rng.uniform01() < t) ? hi : lo;
    }

    const size_t n_orient = samples.n_orient;
    size_t orient_idx = 0;
    if (orientation_mode_ == OrientationMode::Fixed) {
        orient_idx = fixed_orientation_index_ % std::max<size_t>(1, n_orient);
    } else {
        orient_idx = static_cast<size_t>(rng.uniform01() * static_cast<double>(n_orient));
    }
    const size_t idx = samples.index(std::min(orient_idx, n_orient - 1), bin_idx);

    const double sigma_mt = samples.sigma_mt_m2[idx];
    const double sigma_event = samples.sigma_event_m2[idx];
    if (sigma_event <= 0.0 || n_gas <= 0.0) {
        return false;
    }

    const double k = n_gas * v_rel_mag * sigma_event;
    if (momentum_logging_enabled_) {
        std::lock_guard<std::mutex> lock(state_mutex_);
        auto& log = momentum_logs_[species_id];
        log.attempts += 1;
        log.k_sum += k;
        log.v_sum += v_rel_mag;
        log.v2_sum += v_rel_mag * v_rel_mag;
        log.q_c = view.kin.get_charge();
    }
    const double P = 1.0 - std::exp(-k * dt);
    if (rng.uniform01() >= P) {
        record_attempt(k, true, false);
        return false;
    }

    if (vrel_logging_enabled_) {
        std::lock_guard<std::mutex> lock(state_mutex_);
        auto& log = vrel_logs_[species_id];
        if (log.logv_bins.empty()) {
            log.logv_bins = samples.logv_bins;
            log.bin_counts.assign(samples.logv_bins.size(), 0);
        }
        log.samples += 1;
        log.v_sum += v_rel_mag;
        log.v2_sum += v_rel_mag * v_rel_mag;
        if (bin_idx < log.bin_counts.size()) {
            log.bin_counts[bin_idx] += 1;
        }
    }

    Vec3 dp_collision{0.0, 0.0, 0.0};
    bool sampled = false;
    if (!samples.cdf_offsets.empty() && !samples.cdf_counts.empty() && !samples.cdf_values.empty()) {
        const long long offset = samples.cdf_offsets[idx];
        const long long count = samples.cdf_counts[idx];
        if (offset >= 0 && count > 0) {
            const size_t offset_u = static_cast<size_t>(offset);
            const size_t count_u = static_cast<size_t>(count);
            if (offset_u > samples.cdf_values.size() ||
                count_u > (samples.cdf_values.size() - offset_u)) {
                throw std::runtime_error("[InteractionPotential] Invalid CDF offset/count table for species '" + species_id + "'");
            }
            const double u = rng.uniform01();
            auto begin = samples.cdf_values.begin() + static_cast<std::ptrdiff_t>(offset_u);
            auto end = begin + static_cast<std::ptrdiff_t>(count_u);
            auto itc = std::lower_bound(begin, end, u);
            const size_t local = (itc == end)
                ? (count_u - 1)
                : static_cast<size_t>(std::distance(begin, itc));
            const size_t sample_idx = offset_u + local;
            const size_t base = sample_idx * 3;
            if (base + 2 < samples.dp_samples.size()) {
                dp_collision = Vec3(
                    samples.dp_samples[base],
                    samples.dp_samples[base + 1],
                    samples.dp_samples[base + 2]
                );
                sampled = true;
            }
        }
    }

    if (!sampled && !samples.dp_stats.empty()) {
        const size_t stats_idx = idx * 4;
        const double mean_par = samples.dp_stats[stats_idx + 0];
        const double mean_perp = samples.dp_stats[stats_idx + 1];
        const double var_par = samples.dp_stats[stats_idx + 2];
        const double var_perp = samples.dp_stats[stats_idx + 3];
        const double std_par = std::sqrt(std::max(0.0, var_par));
        const double std_perp = std::sqrt(std::max(0.0, var_perp));
        const double dp_par = mean_par + std_par * rng.normal();
        const double dp_perp = std::abs(mean_perp + std_perp * rng.normal());
        const double phi = 2.0 * M_PI * rng.uniform01();
        dp_collision = Vec3(dp_perp * std::cos(phi), dp_perp * std::sin(phi), dp_par);
        sampled = true;
    }

    if (!sampled) {
        record_attempt(k, true, false);
        return false;
    }

    if (momentum_logging_enabled_) {
        std::lock_guard<std::mutex> lock(state_mutex_);
        auto& log = momentum_logs_[species_id];
        log.collisions += 1;
        log.dp_par_sum += dp_collision.z;
    }

    Vec3 collision_axis = v_rel / v_rel_mag;
    Vec3 t1, t2;
    collision_core::CollisionGeometry::construct_orthonormal_basis(collision_axis, t1, t2);
    const Vec3 dp_lab = t1 * dp_collision.x + t2 * dp_collision.y + collision_axis * dp_collision.z;

    Vec3 v_post = view.kin.vel() + dp_lab / m_ion;
    view.kin.set_vel(v_post);
    if (diagnostics) {
        diagnostics->v_rel_before_m_s = v_rel_mag;
        diagnostics->sigma_mt_m2 = sigma_mt;
    }
    record_attempt(k, false, true);
    return true;
}

CollisionStats InteractionPotentialCollisionHandler::get_stats() const {
    auto* state = stats_state_.get();
    std::lock_guard<std::mutex> lock(state->mutex);
    const std::uint64_t generation = state->generation.load(std::memory_order_acquire);
    CollisionStats reduced{};
    size_t rate_samples = 0;
    double rate_sum = 0.0;
    for (const auto& slot : state->slots) {
        if (!slot || slot->generation != generation) {
            continue;
        }
        reduced.total_collisions += slot->stats.total_collisions;
        reduced.rejected_collisions += slot->stats.rejected_collisions;
        rate_samples += slot->rate_samples;
        rate_sum += slot->rate_sum;
    }
    reduced.average_collision_rate = rate_samples > 0
        ? rate_sum / static_cast<double>(rate_samples)
        : 0.0;
    return reduced;
}

void InteractionPotentialCollisionHandler::reset_stats() {
    stats_state_->generation.fetch_add(1, std::memory_order_acq_rel);
}

} // namespace ICARION::physics
