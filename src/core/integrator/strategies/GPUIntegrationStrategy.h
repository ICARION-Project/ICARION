// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IIntegrationStrategy.h"
#include "core/log/Logger.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/physics/forces/MagneticFieldForce.h"
#include <vector>
#include <memory>

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUIntegrationHelper.h"
#endif

namespace ICARION::integrator {

class GPUIntegrationStrategy : public IIntegrationStrategy {
public:
    enum class Kind { RK4, RK45, Boris };

    GPUIntegrationStrategy(Kind kind,
                           std::unique_ptr<IIntegrationStrategy> cpu_fallback,
                           size_t gpu_threshold = 5000);

    // Optional damping configuration (constant or per-ion nu) for GPU path.
    // No effect if GPU helper is unavailable.
    void set_gpu_damping_constant(double nu_const);
    void set_gpu_damping_per_ion(const std::vector<float>& nu_per_ion);

    void step(core::IonEnsemble& ensemble,
              size_t ion_idx,
              double t,
              double dt,
              const physics::ForceRegistry& force_registry) override;

    std::string name() const override;
    bool is_adaptive() const override;

    bool supports_batch() const override;
    bool step_batch(core::IonEnsemble& ensemble,
                    double t,
                    double dt,
                    const std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
                    const std::vector<int>& domain_indices) override;

private:
    bool run_cpu_fallback(core::IonEnsemble& ensemble,
                          double t,
                          double dt,
                          const std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
                          const std::vector<int>& domain_indices);

    const physics::ElectricFieldForce* extract_electric_force(const std::shared_ptr<physics::ForceRegistry>& registry) const;

    Kind kind_;
    size_t threshold_;
    std::unique_ptr<IIntegrationStrategy> cpu_fallback_;
    bool adaptive_ = false;
    bool warned_multiple_domains_ = false;
    bool warned_force_mix_ = false;
    bool warned_experimental_ = false;
#ifdef ICARION_USE_GPU
    std::shared_ptr<icarion::gpu::GPUContext> gpu_context_;
    std::shared_ptr<icarion::gpu::GPUIntegrationHelper> gpu_helper_;
#endif
};

}  // namespace ICARION::integrator
