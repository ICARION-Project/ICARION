// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        geometryReader.h
 *   @brief       Loads geometry file for EHSS collision model.
 *
 *   @details
 *   Parses a JSON file containing molecular geometry definitions and returns one or more
 *   `MoleculeRecord` entries. Each molecule includes its name, diameter, and a list of
 *   atom records with mass, charge, position, and Lennard-Jones parameters. If a specific
 *   `targetName` is provided, only the matching molecule is returned.
 *
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once

#include <cmath>
#include <fstream>
#include <json/json.h> 
#include <stdexcept>
#include <string>
#include <vector>
#include "../utils/lj_params.h"

// Simple atom record
/**
 * @struct AtomRecord
 * @brief Single atom in molecular geometry for EHSS collisions
 * 
 * Contains atomic properties needed for collision calculations:
 * - Position in molecule frame
 * - Mass and charge
 * - Lennard-Jones parameters for interaction potentials
 */
struct AtomRecord {
    std::string type;          ///< Element symbol (e.g., "C", "H", "N")
    double      mass_amu;      ///< Atomic mass [amu]
    double      charge_e;      ///< Total atomic charge [elementary charges]
    double      partCharge_e;  ///< Partial charge [e] for electrostatics
    double      posx_m;        ///< X position [m] in molecule frame
    double      posy_m;        ///< Y position [m] in molecule frame
    double      posz_m;        ///< Z position [m] in molecule frame
    double      LJ_sigma_m;    ///< Lennard-Jones sigma parameter [m]
    double      LJ_eps_J;      ///< Lennard-Jones epsilon parameter [J]
};

// Simple molecule record
/**
 * @struct MoleculeRecord
 * @brief Complete molecular geometry for EHSS collision model
 * 
 * Contains molecular structure and reference data:
 * - Atom positions and properties
 * - Effective diameter for collision estimates
 * - Optional reference CCS from MOBCAL calculations
 */
struct MoleculeRecord {
    std::string             name;        ///< Molecule identifier (e.g., "Pentanal", "H2O")
    double                  diameter_m;  ///< Effective molecular diameter [m]
    double                  CCS_m2;      ///< Reference collision cross-section [m²] (if available)
    bool                    has_CCS;     ///< true if reference CCS was provided in geometry file
    std::vector<AtomRecord> atoms;       ///< List of atoms with positions and properties
};

/**
 * @brief Read molecular geometry from JSON file for EHSS collisions
 * 
 * @param filename Path to JSON geometry file
 * @param targetName Optional: specific molecule to extract (if empty, returns all)
 * 
 * @return Vector of molecule records with atomic coordinates and properties
 * 
 * JSON file format:
 * [
 *   {
 *     "name": "Pentanal",
 *     "diameter": 6.5e-10,
 *     "CCS_m2": 9.8e-19,  // optional reference CCS
 *     "atoms": [
 *       {
 *         "type": "C",
 *         "mass": 12.011,
 *         "charge": 0.0,
 *         "partCharge": -0.1,
 *         "posx": 0.0,  // Angstrom
 *         "posy": 0.0,
 *         "posz": 0.0,
 *         "LJsigma": 3.4,  // Angstrom (optional, defaults from element table)
 *         "LJeps": 0.47    // kJ/mol
 *       },
 *       ...
 *     ]
 *   }
 * ]
 * 
 * Coordinates are converted from Angstrom to meters internally.
 * LJ parameters converted from kJ/mol to Joules per particle.
 * 
 * @throws std::runtime_error if file not found or targetName not found
 */
inline std::vector<MoleculeRecord> read_geometry_file(const std::string& filename,
                                                      const std::string& targetName = "") {
    std::ifstream in(filename);
    if (!in.good()) {
        throw std::runtime_error("Geometry file not found: " + filename);
    }

    Json::Value root;
    in >> root;
    in.close();

    std::vector<MoleculeRecord> molecules;

    for (Json::ArrayIndex i = 0; i < root.size(); ++i) {
        std::string name = root[i]["name"].asString();

        if (!targetName.empty() && name != targetName) {
            continue;
        }

        MoleculeRecord mol;
        mol.name = name;

        mol.diameter_m = root[i]["diameter"].asDouble();
        
        // Read optional CCS field (MobCal reference value)
        if (root[i].isMember("CCS_m2")) {
            mol.CCS_m2 = root[i]["CCS_m2"].asDouble();
            mol.has_CCS = true;
        } else {
            mol.CCS_m2 = 0.0;
            mol.has_CCS = false;
        }

        for (Json::ArrayIndex j = 0; j < root[i]["atoms"].size(); ++j) {
            const auto& atom = root[i]["atoms"][j];
            AtomRecord  a;

            a.type         = atom["type"].asString();
            a.mass_amu     = atom["mass"].asDouble();
            a.charge_e     = atom["charge"].asDouble();
            a.partCharge_e = atom["partCharge"].asDouble();

            a.posx_m     = atom["posx"].asDouble() * 1e-10;
            a.posy_m     = atom["posy"].asDouble() * 1e-10;
            a.posz_m     = atom["posz"].asDouble() * 1e-10;
            // Read optional LJ sigma (Angstrom in JSON). If missing, fall back to
            // the element lookup table provided by `lj_params.h`.
            if (atom.isMember("LJsigma")) {
                a.LJ_sigma_m = atom["LJsigma"].asDouble() * 1e-10;
            } else {
                // Use lookup defaults for common elements (meters)
                a.LJ_sigma_m = get_lj_sigma_m(a.type);
            }
            a.LJ_eps_J   = atom["LJeps"].asDouble() * 1e3 / 6.02214076e23;  // kJ/mol → J

            mol.atoms.push_back(a);
        }

        molecules.push_back(mol);
    }

    if (!targetName.empty() && molecules.empty()) {
        throw std::runtime_error("Molecule '" + targetName + "' not found in " + filename);
    }

    return molecules;
}
