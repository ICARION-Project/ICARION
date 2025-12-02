#!/usr/bin/env python3
"""
Generate thermalization test configs for validation suite.

Test Matrix:
- 3 temperatures: 150K, 300K, 1000K
- 5 pressures: 0.2, 2, 20, 200, 2000 Pa
- 3 ion species: H3O+, PentanalH+, 26DTBPH+
- 2 collision models: HSS, EHSS
Total: 90 configs

Author: ICARION Validation Suite
"""

import json
import os
from pathlib import Path
import math

# Physical constants
AMU_TO_KG = 1.66053906660e-27
ELEM_CHARGE = 1.602176634e-19
K_BOLTZMANN = 1.380649e-23
HE_MASS_AMU = 4.002602

def calc_collision_time(mass_amu, ccs_m2, pressure_pa, temperature_k):
    """
    Calculate expected collision time using kinetic theory.
    
    τ_collision = 1 / (n * σ * v_rel)
    where:
    - n = number density [m^-3]
    - σ = collision cross section [m^2]
    - v_rel = relative velocity from Maxwell-Boltzmann distribution
    """
    # Number density from ideal gas law
    n_he = pressure_pa / (K_BOLTZMANN * temperature_k)
    
    # Reduced mass (ion + He)
    mu_kg = (mass_amu * HE_MASS_AMU) / (mass_amu + HE_MASS_AMU) * AMU_TO_KG
    
    # Relative velocity (Maxwell-Boltzmann)
    v_rel = math.sqrt(8 * K_BOLTZMANN * temperature_k / (math.pi * mu_kg))
    
    # Collision frequency
    nu_collision = n_he * ccs_m2 * v_rel
    
    # Collision time
    tau_collision = 1.0 / nu_collision
    
    return tau_collision, v_rel, n_he

def generate_config(ion_species, collision_model, temperature_k, pressure_pa):
    """Generate single thermalization test config."""
    
    # Ion species data
    species_data = {
        "H3O+": {"mass_amu": 19.0, "ccs_m2": 24.9e-20, "duration_factor": 20},
        "PentanalH+": {"mass_amu": 87.0, "ccs_m2": 53.7e-20, "duration_factor": 100},
        "26DTBPH+": {"mass_amu": 192.0, "ccs_m2": 87.02e-20, "duration_factor": 150}
    }
    
    ion_data = species_data[ion_species]
    
    # Calculate collision time and timestep
    tau_coll, v_rel, n_he = calc_collision_time(
        ion_data["mass_amu"], 
        ion_data["ccs_m2"], 
        pressure_pa, 
        temperature_k
    )
    
    # dt = τ_collision / 50
    dt_s = tau_coll / 50.0
    
    # Total time = duration_factor * τ_collision
    total_time_s = ion_data["duration_factor"] * tau_coll
    
    # Write interval = 1 collision time
    write_interval = max(1, int(tau_coll / dt_s))
    
    # Generate config
    config = {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,  # Use 4 threads for faster execution
            "rng_seed": 42
        },
        "physics": {
            "collision_model": collision_model,
            "enable_space_charge": False,
            "enable_reactions": False
        },
        "species_database": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": f"/home/chsch95/ICARION/results/v1.0_test/physics/thermalization",
            "trajectory_file": f"therm_{collision_model.lower()}_{ion_species.replace('+', 'p').replace(',', '').replace('-', '')}_{temperature_k}K_{pressure_pa}Pa.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": ion_species,
                    "count": 10000,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.01],
                        "std": [0.001, 0.001, 0.0005]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": 0.1
                    }
                }
            ]
        },
        "domains": [
            {
                "name": "thermalization_chamber",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, -5.0],
                    "radius_m": 10,
                    "length_m": 10
                },
                "environment": {
                    "pressure_Pa": pressure_pa,
                    "temperature_K": temperature_k,
                    "gas_species": ["He"]
                },
                "fields": {
                    "electric": {"type": "none"},
                    "magnetic": {"type": "none"}
                }
            }
        ]
    }
    
    # Add metadata for analysis
    config["_metadata"] = {
        "test_type": "thermalization",
        "collision_model": collision_model,
        "ion_species": ion_species,
        "temperature_K": temperature_k,
        "pressure_Pa": pressure_pa,
        "expected_collision_time_s": tau_coll,
        "timestep_factor": 50,
        "duration_collision_times": ion_data["duration_factor"],
        "relative_velocity_ms": v_rel,
        "he_number_density_m3": n_he
    }
    
    return config

def main():
    """Generate all thermalization test configs."""
    
    # Test parameters
    temperatures = [150, 300, 1000]  # K
    pressures = [0.2, 2.0, 20.0, 200.0, 2000.0]  # Pa
    ion_species = ["H3O+", "PentanalH+", "26DTBPH+"]
    collision_models = ["HSS", "EHSS"]
    
    # Output directory
    output_dir = Path(__file__).parent.parent / "configs" / "physics" / "thermalization"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate configs
    total_configs = 0
    
    for temp in temperatures:
        for pressure in pressures:
            for ion in ion_species:
                for model in collision_models:
                    # Generate config
                    config = generate_config(ion, model, temp, pressure)
                    
                    # Generate filename
                    ion_clean = ion.replace('+', 'p').replace(',', '').replace('-', '')
                    filename = f"{model.lower()}_{ion_clean}_{temp}K_{pressure}Pa.json"
                    filepath = output_dir / filename
                    
                    # Write config
                    with open(filepath, 'w') as f:
                        json.dump(config, f, indent=2)
                    
                    total_configs += 1
                    print(f"Generated: {filename}")
                    
                    # Print collision time info for first few
                    if total_configs <= 5:
                        meta = config["_metadata"]
                        print(f"  τ_coll = {meta['expected_collision_time_s']:.2e} s")
                        print(f"  dt = {config['simulation']['dt_s']:.2e} s")
                        print(f"  total_time = {config['simulation']['total_time_s']:.2e} s")
                        print()
    
    print(f"\nGenerated {total_configs} thermalization configs")
    print(f"Matrix: {len(temperatures)} temps × {len(pressures)} pressures × {len(ion_species)} ions × {len(collision_models)} models")
    
    # Generate summary
    summary_file = output_dir / "README.md"
    with open(summary_file, 'w') as f:
        f.write("# Thermalization Test Configs\n\n")
        f.write(f"**Generated:** {total_configs} configs\n")
        f.write(f"**Test Matrix:**\n")
        f.write(f"- Temperatures: {temperatures} K\n")
        f.write(f"- Pressures: {pressures} Pa\n")
        f.write(f"- Ion species: {ion_species}\n")
        f.write(f"- Collision models: {collision_models}\n\n")
        f.write("**Expected Results:**\n")
        f.write("- Final temperature should match environment temperature\n")
        f.write("- Thermalization time ~ few collision times\n")
        f.write("- HSS and EHSS should give similar results\n")

if __name__ == "__main__":
    main()