// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "core/config/types/IonConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/validation/ValidationResult.h"
#include "core/types/IonState.h"
#include <vector>
#include <random>
#include <filesystem>

namespace ICARION::config {

/**
 * @brief Result of ion generation
 */
struct IonGenerationResult {
    std::vector<core::IonState> ions;
    ValidationResult validation;
    
    bool success() const { return validation.valid; }
};

/**
 * @brief Ion initialization from configuration
 * 
 * Generates ion ensemble from IonConfig using:
 * - Species database for physical properties
 * - Position/velocity distributions (per species)
 * - Random number generator for sampling
 */
class IonLoader {
public:
    /**
     * @brief Generate ions from configuration
     * 
     * @param config Ion configuration
     * @param species_db Species database for mass/charge lookup
     * @param rng Random number generator
     * @return IonGenerationResult with ions and validation status
     * 
     * Behavior:
     * - If config.from_file: Load from JSON file
     * - Else: Generate from species list (each species has own boundaries)
     * - Returns validation errors/warnings instead of throwing
     */
    static IonGenerationResult generate_ions(
        const IonConfig& config,
        const SpeciesDatabase& species_db,
        std::mt19937& rng
    );
    
    /**
     * @brief Load ions from JSON file
     * 
     * @param filepath Path to ion cloud JSON
     * @param species_db Species database for property lookup
     * @return IonGenerationResult with ions and validation status
     * 
     * Expected format:
     * ```json
     * {
     *   "ions": [
     *     {
     *       "species": "H3O+",
     *       "pos": [0, 0, 0],
     *       "vel": [0, 0, 0],
     *       "birth_time": 0.0
     *     }
     *   ]
     * }
     * ```
     */
    static IonGenerationResult load_from_file(
        const std::filesystem::path& filepath,
        const SpeciesDatabase& species_db
    );

private:
    /**
     * @brief Generate ions for single species (with its own boundaries)
     */
    static std::vector<core::IonState> generate_species(
        const IonSpeciesConfig& spec_config,
        const SpeciesDatabase& species_db,
        std::mt19937& rng
    );
    
    /**
     * @brief Sample position from distribution
     */
    static Vec3 sample_position(
        const PositionConfig& config,
        std::mt19937& rng
    );
    
    /**
     * @brief Sample velocity from distribution
     * 
     * @param config Velocity configuration
     * @param ion_mass_kg Ion mass [kg]
     * @param rng Random number generator
     * 
     * Note: For Thermal distribution, directions are isotropic (no directed drift)
     */
    static Vec3 sample_velocity(
        const VelocityConfig& config,
        double ion_mass_kg,
        std::mt19937& rng
    );
};

} // namespace ICARION::config
