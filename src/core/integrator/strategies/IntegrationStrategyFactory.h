// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

#include "IIntegrationStrategy.h"
#include "RK4Strategy.h"
#include "RK45Strategy.h"
#include "BorisStrategy.h"
#include "core/config/types/FullConfig.h"
#include "GPUIntegrationStrategy.h"

namespace ICARION {
namespace integrator {

class IntegrationStrategyFactory {
public:
    static std::unique_ptr<IIntegrationStrategy> create(const std::string& strategy_name) {
        return create(strategy_name, nullptr);
    }

    static std::unique_ptr<IIntegrationStrategy> create(
        const std::string& strategy_name,
        const config::FullConfig* config) {
        auto cpu_strategy = create_cpu(strategy_name);

        if (!config || !config->simulation.enable_gpu) {
            return cpu_strategy;
        }

#ifdef ICARION_USE_GPU
        GPUIntegrationStrategy::Kind kind = GPUIntegrationStrategy::Kind::RK4;
        if (strategy_name == "RK45") {
            kind = GPUIntegrationStrategy::Kind::RK45;
        } else if (strategy_name == "Boris") {
            kind = GPUIntegrationStrategy::Kind::Boris;
        }
        return std::make_unique<GPUIntegrationStrategy>(kind, std::move(cpu_strategy));
#else
        return cpu_strategy;
#endif
    }

    static std::vector<std::string> supported_strategies() {
        return {"RK4", "RK45", "Boris"};
    }

private:
    static std::unique_ptr<IIntegrationStrategy> create_cpu(const std::string& strategy_name) {
        if (strategy_name == "RK4") {
            return std::make_unique<RK4Strategy>();
        } else if (strategy_name == "RK45") {
            return std::make_unique<RK45Strategy>();
        } else if (strategy_name == "Boris") {
            return std::make_unique<BorisStrategy>();
        }
        throw std::invalid_argument(
            "IntegrationStrategyFactory: Unknown strategy '" + strategy_name +
            "'. Supported: 'RK4', 'RK45', 'Boris'.");
    }
};

} // namespace integrator
} // namespace ICARION
