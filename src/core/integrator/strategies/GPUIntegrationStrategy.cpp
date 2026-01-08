// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "GPUIntegrationStrategy.h"

namespace ICARION::integrator {

GPUIntegrationStrategy::GPUIntegrationStrategy(
    Kind kind,
    std::unique_ptr<IIntegrationStrategy> cpu_fallback,
    size_t gpu_threshold)
    : kind_(kind),
      threshold_(gpu_threshold),
      cpu_fallback_(std::move(cpu_fallback)) {
    adaptive_ = cpu_fallback_ ? cpu_fallback_->is_adaptive() : false;
#ifdef ICARION_USE_GPU
    try {
        gpu_context_ = icarion::gpu::GPUContext::create(0);
        if (gpu_context_) {
            gpu_helper_ = icarion::gpu::GPUIntegrationHelper::create(*gpu_context_, gpu_threshold);
            if (!gpu_helper_) {
                gpu_context_.reset();
            }
        }
    } catch (const std::exception& e) {
        log::Logger::main()->warn("GPUIntegrationStrategy: Failed to create GPU context/helper: {}", e.what());
        gpu_context_.reset();
        gpu_helper_.reset();
    }
#else
    (void)gpu_threshold;
#endif
}

void GPUIntegrationStrategy::set_gpu_damping_constant(double nu_const) {
#ifdef ICARION_USE_GPU
    if (gpu_helper_) {
        gpu_helper_->set_damping_constant(nu_const);
        gpu_damping_enabled_ = true;
    }
#else
    (void)nu_const;
#endif
}

void GPUIntegrationStrategy::set_gpu_damping_per_ion(const std::vector<float>& nu_per_ion) {
#ifdef ICARION_USE_GPU
    if (gpu_helper_) {
        gpu_helper_->set_damping_per_ion(nu_per_ion);
        gpu_damping_enabled_ = true;
    }
#else
    (void)nu_per_ion;
#endif
}

void GPUIntegrationStrategy::step(core::IonEnsemble& ensemble,
                                  size_t ion_idx,
                                  double t,
                                  double dt,
                                  const physics::ForceRegistry& force_registry) {
    if (cpu_fallback_) {
        cpu_fallback_->step(ensemble, ion_idx, t, dt, force_registry);
    }
}

std::string GPUIntegrationStrategy::name() const {
    std::string suffix = "";
#ifdef ICARION_USE_GPU
    if (gpu_helper_) {
        suffix = "+GPU";
    }
#endif
    return cpu_fallback_ ? cpu_fallback_->name() + suffix : "GPUIntegrator";
}

bool GPUIntegrationStrategy::is_adaptive() const {
    return adaptive_;
}

void GPUIntegrationStrategy::set_parallel_enabled(bool enabled) {
    if (cpu_fallback_) {
        cpu_fallback_->set_parallel_enabled(enabled);
    }
}

bool GPUIntegrationStrategy::supports_batch() const {
#ifdef ICARION_USE_GPU
    return gpu_helper_ && gpu_helper_->is_enabled();
#else
    return false;
#endif
}

bool GPUIntegrationStrategy::step_batch(
    core::IonEnsemble& ensemble,
    double t,
    double dt,
    const std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
    const std::vector<int>& domain_indices) {
#ifndef ICARION_USE_GPU
    (void)ensemble;
    (void)t;
    (void)dt;
    (void)registries;
    (void)domain_indices;
    return false;
#else
    if (!warned_experimental_) {
        log::Logger::main()->warn("GPU integration is experimental; unsupported forces will force CPU fallback.");
        warned_experimental_ = true;
    }

    if (!supports_batch()) {
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }

    const size_t n = ensemble.size();
    std::vector<size_t> mapping;
    mapping.reserve(n);
    int selected_domain = -1;
    bool multi_domain = false;

    for (size_t i = 0; i < n; ++i) {
        if (domain_indices[i] < 0) {
            continue;
        }
        if (!ensemble.is_active(i)) {
            continue;
        }
        mapping.push_back(i);
        if (selected_domain < 0) {
            selected_domain = domain_indices[i];
        } else if (domain_indices[i] != selected_domain) {
            multi_domain = true;
        }
    }

    if (mapping.size() < threshold_ || selected_domain < 0) {
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }

    if (multi_domain) {
        if (!warned_multiple_domains_) {
            log::Logger::main()->warn(
                "GPUIntegrationStrategy: multiple domains present in batch – falling back to CPU");
            warned_multiple_domains_ = true;
        }
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }

    if (selected_domain >= static_cast<int>(registries.size()) || !registries[selected_domain]) {
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }

    const auto& forces = registries[selected_domain]->forces();
    size_t electric_count = 0;
    bool has_space_charge = registries[selected_domain]->space_charge_model() != nullptr;
    bool has_damping = false;
    bool has_magnetic = false;
    bool has_other = false;
    for (const auto& f : forces) {
        if (dynamic_cast<const physics::ElectricFieldForce*>(f.get())) {
            ++electric_count;
            continue;
        }
        if (dynamic_cast<const physics::DampingForce*>(f.get())) {
            has_damping = true;
            continue;
        }
        if (dynamic_cast<const physics::MagneticFieldForce*>(f.get())) {
            has_magnetic = true;
            continue;
        }
        has_other = true;
    }

    auto log_gpu_reject = [&](const std::string& reason) {
        if (!warned_force_mix_) {
            log::Logger::main()->warn("GPUIntegrationStrategy: {} – using CPU for domain '{}'",
                                      reason,
                                      registries[selected_domain]->domain() ? registries[selected_domain]->domain()->name : "<unknown>");
            warned_force_mix_ = true;
        }
    };

    if (has_space_charge) {
        log_gpu_reject("space-charge present");
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }
    if (has_damping && !gpu_damping_enabled_) {
        log_gpu_reject("damping force present");
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }
    if (has_magnetic) {
        log_gpu_reject("magnetic force present");
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }
    if (has_other) {
        log_gpu_reject("unsupported force type present");
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }
    if (electric_count != 1) {
        log_gpu_reject("requires exactly one ElectricFieldForce");
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }

    const auto* e_force = extract_electric_force(registries[selected_domain]);
    if (!e_force || !e_force->get_field_provider()) {
        if (!warned_force_mix_) {
            log::Logger::main()->warn(
                "GPUIntegrationStrategy: domain '{}' lacks GPU-supportable ElectricFieldForce; using CPU",
                registries[selected_domain]->domain() ? registries[selected_domain]->domain()->name : "<unknown>");
            warned_force_mix_ = true;
        }
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }

    std::vector<IonState> gpu_ions;
    gpu_ions.reserve(mapping.size());
    for (size_t idx : mapping) {
        gpu_ions.push_back(ensemble.ion_state(idx));
    }

    bool gpu_success = false;
    switch (kind_) {
        case Kind::RK4:
            gpu_success = gpu_helper_->integrate_batch_rk4(gpu_ions, dt, t, e_force->get_field_provider());
            break;
        case Kind::RK45:
            gpu_success = gpu_helper_->integrate_batch_rk45(gpu_ions, dt, t, e_force->get_field_provider());
            break;
        case Kind::Boris:
            gpu_success = gpu_helper_->integrate_batch_boris(gpu_ions, dt, t, e_force->get_field_provider());
            break;
    }

    if (!gpu_success) {
        return run_cpu_fallback(ensemble, t, dt, registries, domain_indices);
    }

    for (size_t i = 0; i < mapping.size(); ++i) {
        ensemble.apply_ion_state(mapping[i], gpu_ions[i]);
    }

    return true;
#endif
}

bool GPUIntegrationStrategy::run_cpu_fallback(
    core::IonEnsemble& ensemble,
    double t,
    double dt,
    const std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
    const std::vector<int>& domain_indices) {
    if (!cpu_fallback_) {
        return false;
    }

    const size_t n = ensemble.size();
    for (size_t i = 0; i < n; ++i) {
        if (domain_indices[i] < 0 || !ensemble.is_active(i)) {
            continue;
        }
        const int dom = domain_indices[i];
        if (dom < 0 || dom >= static_cast<int>(registries.size())) {
            continue;
        }
        const auto& registry = registries[dom];
        if (!registry) {
            continue;
        }
        cpu_fallback_->step(ensemble, i, t, dt, *registry);
    }
    return false;
}

const physics::ElectricFieldForce* GPUIntegrationStrategy::extract_electric_force(
    const std::shared_ptr<physics::ForceRegistry>& registry) const {
    if (!registry || registry->forces().size() != 1) {
        return nullptr;
    }
    return dynamic_cast<const physics::ElectricFieldForce*>(registry->forces().front().get());
}

} // namespace ICARION::integrator
