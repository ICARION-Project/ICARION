// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_LEGACY_ADAPTER_H
#define ICARION_CONFIG_LEGACY_ADAPTER_H

#include "../types/FullConfig.h"
#include "core/param/paramUtils.h"
#include <vector>

namespace ICARION::config {

/**
 * @brief Adapter to convert new FullConfig to legacy GlobalParams + InstrumentDomain
 * 
 * ⚠️ TEMPORARY ADAPTER - TO BE REMOVED IN FUTURE REFACTORING ⚠️
 * 
 * This adapter exists solely to bridge the new modular config system with the
 * legacy parameter structures still used by the integrator, physics, and IO code.
 * 
 * Future work:
 * - Refactor integrator to use SimulationConfig + DomainConfig directly
 * - Refactor physics modules to use PhysicsConfig + EnvironmentConfig
 * - Refactor HDF5 writer to serialize FullConfig instead of GlobalParams
 * - Remove this adapter entirely once migration is complete
 * 
 * @see FullConfig - Modern modular configuration
 * @see GlobalParams - Legacy global parameters (to be deprecated)
 * @see InstrumentDomain - Legacy domain parameters (to be deprecated)
 */
class LegacyAdapter {
public:
    /**
     * @brief Convert FullConfig to legacy GlobalParams
     * 
     * Maps modern config sections to legacy flat structure:
     * - SimulationConfig → time parameters, execution flags
     * - PhysicsConfig → collision model, feature flags
     * - OutputConfig → output file paths
     * 
     * @param config Modern configuration
     * @return GlobalParams Legacy parameter structure
     */
    static core::GlobalParams to_global_params(const FullConfig& config);
    
    /**
     * @brief Convert FullConfig.domains to legacy InstrumentDomain vector
     * 
     * Maps each DomainConfig to InstrumentDomain:
     * - GeometryConfig → Geometry
     * - EnvironmentConfig → Environment
     * - FieldsConfig → DC/RF/AC/Magnetic voltages
     * 
     * @param config Modern configuration
     * @return Vector of legacy InstrumentDomain structures
     */
    static std::vector<core::InstrumentDomain> to_instrument_domains(const FullConfig& config);

private:
    // Helper: Convert single DomainConfig to InstrumentDomain
    static core::InstrumentDomain convert_domain(const DomainConfig& domain);
    
    // Helper: Convert GeometryConfig to legacy Geometry
    static core::Geometry convert_geometry(const GeometryConfig& geom);
    
    // Helper: Convert EnvironmentConfig to legacy Environment
    static core::Environment convert_environment(const EnvironmentConfig& env);
    
    // Helper: Convert FieldsConfig to legacy voltage structures
    static void convert_fields(const FieldsConfig& fields, core::InstrumentDomain& dom);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_LEGACY_ADAPTER_H
