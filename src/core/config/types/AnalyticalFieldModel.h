// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IFieldModel.h"
#include "DomainConfig.h"
#include <map>

namespace ICARION::config {

/**
 * @brief Analytical electric field model mirroring ElectricFieldForce formulas.
 *
 * Evaluates instrument-specific analytical fields based on DomainConfig
 * (IMS, LQIT, TOF, FTICR, Orbitrap, QuadrupoleRF).
 */
class AnalyticalFieldModel : public IFieldModel {
public:
    explicit AnalyticalFieldModel(const DomainConfig& domain) : domain_(&domain) {}

    Vec3 E(const Vec3& global_pos, double t) const override;

private:
    const DomainConfig* domain_;

    Vec3 compute_lqit_field(const Vec3& pos, double t) const;
    Vec3 compute_ims_field(const Vec3& pos, double t) const;
    Vec3 compute_tof_field(const Vec3& pos) const;
    Vec3 compute_fticr_field(const Vec3& pos) const;
    Vec3 compute_orbitrap_field(const Vec3& pos, double t) const;
    Vec3 compute_quadrupole_rf_field(const Vec3& pos, double t) const;
};

} // namespace ICARION::config
