// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IFieldModel.h"
#include <memory>

// Forward declaration of grid/map field provider (lives in global namespace)
class IFieldProvider;

namespace ICARION::config {

/**
 * @brief Field model wrapper around an IFieldProvider (grid/map fields).
 *
 * Provides E-field sampling through the same interface as analytical models,
 * enabling SSOT for field access in forces/integrators.
 */
class FieldProviderModel : public IFieldModel {
public:
    explicit FieldProviderModel(std::shared_ptr<::IFieldProvider> provider)
        : provider_(std::move(provider)) {}

    Vec3 E(const Vec3& global_pos, double t) const override;

private:
    std::shared_ptr<::IFieldProvider> provider_;
};

} // namespace ICARION::config
