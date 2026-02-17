// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include "core/io/moleculeLoader.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <json/json.h>

namespace ICARION {
namespace io {

// -----------------------------
// Atom member functions
// -----------------------------

void Atom::validate() const {
    if (element.empty()) {
        throw std::runtime_error("Atom element symbol cannot be empty");
    }
    if (mass_u <= 0.0) {
        throw std::runtime_error("Atom '" + element + "': mass must be positive, got " + std::to_string(mass_u));
    }
    if (std::abs(partial_charge_e) > 5.0) {
        throw std::runtime_error("Atom '" + element + "': partial charge seems unreasonable (|q| > 5), got " + std::to_string(partial_charge_e));
    }
    if (LJ_sigma_m < 0.0) {
        throw std::runtime_error("Atom '" + element + "': LJ sigma cannot be negative, got " + std::to_string(LJ_sigma_m));
    }
    if (LJ_epsilon_eV < 0.0) {
        throw std::runtime_error("Atom '" + element + "': LJ epsilon cannot be negative, got " + std::to_string(LJ_epsilon_eV));
    }
}

// -----------------------------
// Molecule member functions
// -----------------------------

void Molecule::calculate_properties() {
    if (atoms.empty()) {
        throw std::runtime_error("Cannot calculate properties: molecule '" + name + "' has no atoms");
    }
    
    // Calculate total mass and charge
    total_mass_u = 0.0;
    total_charge_e = 0.0;
    for (const auto& atom : atoms) {
        total_mass_u += atom.mass_u;
        total_charge_e += atom.partial_charge_e;
    }
    
    // Calculate center of mass
    Vec3 com_weighted = Vec3{0.0, 0.0, 0.0};
    for (const auto& atom : atoms) {
        com_weighted = com_weighted + atom.pos_m * atom.mass_u;
    }
    center_of_mass_m = com_weighted / total_mass_u;
    
    // Calculate diameter if not set (max pairwise distance)
    if (diameter_m <= 0.0 && atoms.size() > 1) {
        double max_dist = 0.0;
        for (size_t i = 0; i < atoms.size(); ++i) {
            for (size_t j = i + 1; j < atoms.size(); ++j) {
                Vec3 diff = atoms[i].pos_m - atoms[j].pos_m;
                double dist = std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
                max_dist = std::max(max_dist, dist);
            }
        }
        diameter_m = max_dist;
    }
}

void Molecule::validate(double charge_tolerance) const {
    if (name.empty()) {
        throw std::runtime_error("Molecule name cannot be empty");
    }
    if (atoms.empty()) {
        throw std::runtime_error("Molecule '" + name + "' has no atoms");
    }
    
    // Validate all atoms
    for (const auto& atom : atoms) {
        atom.validate();
    }
    
    // Check charge neutrality (or integer charge for ions)
    double charge_remainder = std::abs(total_charge_e - std::round(total_charge_e));
    if (charge_remainder > charge_tolerance) {
        throw std::runtime_error("Molecule '" + name + "': total charge (" + 
                                std::to_string(total_charge_e) + 
                                " e) is not close to an integer value");
    }
    
    if (CCS_m2 < 0.0) {
        throw std::runtime_error("Molecule '" + name + "': CCS cannot be negative");
    }
    
    if (diameter_m < 0.0) {
        throw std::runtime_error("Molecule '" + name + "': diameter cannot be negative");
    }
}

// -----------------------------
// Loader functions
// -----------------------------

Molecule load_molecule(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open molecule file: " + filepath);
    }
    
    Json::Value root;
    try {
        file >> root;
    } catch (const Json::Exception& e) {
        throw std::runtime_error("JSON parse error in " + filepath + ": " + e.what());
    }
    
    if (!root.isMember("molecule")) {
        throw std::runtime_error("Molecule file missing 'molecule' object: " + filepath);
    }
    
    const Json::Value& mol_json = root["molecule"];
    
    // Required fields
    if (!mol_json.isMember("name")) {
        throw std::runtime_error("Molecule missing required field 'name' in: " + filepath);
    }
    if (!mol_json.isMember("atoms")) {
        throw std::runtime_error("Molecule missing required field 'atoms' in: " + filepath);
    }
    if (!mol_json["atoms"].isArray()) {
        throw std::runtime_error("'atoms' must be an array in: " + filepath);
    }
    
    Molecule mol;
    mol.name = mol_json["name"].asString();
    mol.diameter_m = mol_json.get("diameter_m", 0.0).asDouble();
    mol.CCS_m2 = mol_json.get("CCS_m2", 0.0).asDouble();
    
    // Load atoms
    const Json::Value& atoms_json = mol_json["atoms"];
    mol.atoms.reserve(atoms_json.size());
    
    for (Json::ArrayIndex i = 0; i < atoms_json.size(); ++i) {
        const Json::Value& atom_json = atoms_json[i];
        
        // Required fields
        if (!atom_json.isMember("element")) {
            throw std::runtime_error("Atom " + std::to_string(i) + " missing 'element' in: " + filepath);
        }
        if (!atom_json.isMember("pos")) {
            throw std::runtime_error("Atom " + std::to_string(i) + " missing 'pos' in: " + filepath);
        }
        if (!atom_json.isMember("mass_u")) {
            throw std::runtime_error("Atom " + std::to_string(i) + " missing 'mass_u' in: " + filepath);
        }
        if (!atom_json.isMember("partial_charge_e")) {
            throw std::runtime_error("Atom " + std::to_string(i) + " missing 'partial_charge_e' in: " + filepath);
        }
        
        Atom atom;
        atom.element = atom_json["element"].asString();
        atom.mass_u = atom_json["mass_u"].asDouble();
        atom.partial_charge_e = atom_json["partial_charge_e"].asDouble();
        
        // Parse position array (JSON positions are in Angstrom, convert to meters)
        const Json::Value& pos = atom_json["pos"];
        if (!pos.isArray() || pos.size() != 3) {
            throw std::runtime_error("Atom " + std::to_string(i) + " 'pos' must be array of 3 numbers in: " + filepath);
        }
        atom.pos_m = Vec3{
            pos[0].asDouble() * ANGSTROM_TO_M,
            pos[1].asDouble() * ANGSTROM_TO_M,
            pos[2].asDouble() * ANGSTROM_TO_M
        };
        
        // Optional LJ parameters (JSON in Ångström, convert to meters)
        atom.LJ_sigma_m = atom_json.get("LJ_sigma_angstrom", 0.0).asDouble() * ANGSTROM_TO_M;
        atom.LJ_epsilon_eV = atom_json.get("LJ_epsilon_eV", 0.0).asDouble();
        
        mol.atoms.push_back(atom);
    }
    
    // Calculate derived properties
    mol.calculate_properties();
    
    // Validate
    mol.validate();
    
    std::cout << "Loaded molecule '" << mol.name << "' with " << mol.atoms.size() 
              << " atoms (mass: " << mol.total_mass_u << " u, charge: " 
              << mol.total_charge_e << " e)\n";
    
    return mol;
}

Molecule load_molecule(const std::string& filepath, bool calculate_diameter) {
    auto mol = load_molecule(filepath);
    
    if (calculate_diameter && mol.diameter_m <= 0.0) {
        mol.calculate_properties();  // Recalculate to update diameter
    }
    
    return mol;
}

}  // namespace io
}  // namespace ICARION
