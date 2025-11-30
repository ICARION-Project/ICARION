// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>
#include <vector>

#include "core/utils/mathUtils.h"
#include "utils/constants.h"

namespace ICARION {
namespace io {

/**
 * @brief Single atom in a molecule
 * 
 * Contains position, mass, charge, and Lennard-Jones parameters.
 * Typically derived from DFT calculations.
 * 
 * **Unit Conventions:**
 * - JSON input: positions and LJ_sigma in Angström [Å]
 * - Struct storage: positions and LJ_sigma in SI units [m]
 * - Conversion: 1 Å = 1e-10 m (applied automatically by loader)
 */
struct Atom {
    std::string element;          ///< Element symbol (e.g., "C", "H", "O", "N")
    Vec3 pos_m;                   ///< Position in molecule coordinate system [m] (JSON: Å)
    double mass_u;                ///< Atomic mass [u]
    double partial_charge_e;      ///< Partial charge [elementary charges]
    double LJ_sigma_m;            ///< Lennard-Jones sigma parameter [m] (JSON: Å)
    double LJ_epsilon_eV;         ///< Lennard-Jones epsilon parameter [eV]
    
    /**
     * @brief Validate atom parameters
     * 
     * @throws std::runtime_error if parameters are non-physical
     */
    void validate() const;
};

/**
 * @brief Complete molecular structure
 * 
 * Contains all atoms with their positions, charges, and LJ parameters.
 * Also includes molecular-level properties like total CCS.
 */
struct Molecule {
    std::string name;             ///< Molecule identifier
    double diameter_m;            ///< Approximate molecular diameter [m]
    double CCS_m2;                ///< Collision cross-section [m²]
    std::vector<Atom> atoms;      ///< List of all atoms
    
    // Computed properties
    double total_mass_u;          ///< Sum of all atomic masses [u]
    double total_charge_e;        ///< Sum of all partial charges [e]
    Vec3 center_of_mass_m;        ///< Center of mass position [m]
    
    /**
     * @brief Calculate derived molecular properties
     * 
     * Computes total mass, total charge, and center of mass.
     * Should be called after loading all atoms.
     */
    void calculate_properties();
    
    /**
     * @brief Validate molecule structure
     * 
     * @param charge_tolerance Maximum allowed deviation from integer charge
     * @throws std::runtime_error if validation fails
     * 
     * Checks:
     * - At least one atom present
     * - All atoms have valid parameters
     * - Total charge is close to an integer value
     * - CCS is positive
     */
    void validate(double charge_tolerance = 0.1) const;
    
    /**
     * @brief Get number of atoms
     * 
     * @return Number of atoms in molecule
     */
    size_t num_atoms() const { return atoms.size(); }
};

/**
 * @brief Load molecular geometry from JSON file
 * 
 * @param filepath Path to molecule JSON file
 * @return Molecule structure with all atoms and properties
 * @throws std::runtime_error on file access errors or validation failures
 * 
 * Expected JSON format:
 * ```json
 * {
 *   "molecule": {
 *     "name": "H3O+",
 *     "diameter_m": 3.2e-10,
 *     "CCS_m2": 7.8e-19,
 *     "atoms": [
 *       {
 *         "element": "O",
 *         "pos": [0.0, 0.0, 0.0],          // Positions in Angstrom [Å]
 *         "mass_u": 15.999,
 *         "partial_charge_e": -0.5,
 *         "LJ_sigma_m": 3.07e-10,
 *         "LJ_epsilon_eV": 0.0957
 *       },
 *       {
 *         "element": "H",
 *         "pos": [0.958, 0.0, 0.0],        // Positions in Angstrom [Å]
 *         "mass_u": 1.008,
 *         "partial_charge_e": 0.5,
 *         "LJ_sigma_m": 2.42e-10,
 *         "LJ_epsilon_eV": 0.0157
 *       }
 *     ]
 *   }
 * }
 * ```
 * 
 * **Units:**
 * - `pos`: Angstrom [Å] (automatically converted to meters)
 * - `mass_u`: Atomic mass units [u]
 * - `partial_charge_e`: Elementary charge units [e]
 * - `LJ_sigma_m`: Meters [m]
 * - `LJ_epsilon_eV`: Electron volts [eV]
 * - `diameter_m`, `CCS_m2`: SI units [m], [m²]
 * 
 * Required fields:
 * - molecule.name
 * - molecule.atoms array
 * - For each atom: element, pos, mass_u, partial_charge_e
 * 
 * Optional fields:
 * - molecule.diameter_m (calculated if missing)
 * - molecule.CCS_m2 (must be provided or calculated externally)
 * - atom.LJ_sigma_m, atom.LJ_epsilon_eV (default to 0)
 */
Molecule load_molecule(const std::string& filepath);

/**
 * @brief Load molecule with automatic property calculation
 * 
 * @param filepath Path to molecule JSON file
 * @param calculate_diameter If true, calculates diameter from atom positions
 * @return Molecule structure with calculated properties
 * @throws std::runtime_error on file access errors or validation failures
 * 
 * Convenience function that loads molecule and automatically calculates
 * molecular diameter if not provided in file.
 */
Molecule load_molecule(const std::string& filepath, bool calculate_diameter);

}  // namespace io
}  // namespace ICARION
