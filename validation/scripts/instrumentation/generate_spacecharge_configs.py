#!/usr/bin/env python3
"""
Generate space charge validation suite

Validates:
1. Coulomb expansion from point source (1D, 3D)
2. Direct solver (N<1000) vs Grid solver (N≥1000) comparison
3. Charge conservation and force accuracy

Tests automatic solver selection threshold at N=1000.
"""

import json
import numpy as np
from pathlib import Path

# Output directory
OUT_DIR = Path("../configs/physics/spacecharge")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Physical constants
Q_E = 1.602176634e-19  # Elementary charge
M_PROTON = 1.67262192e-27  # Proton mass
EPSILON_0 = 8.8541878128e-12  # Vacuum permittivity

# Test configurations
print("\n" + "="*80)
print("Space Charge Validation Suite Generator")
print("="*80)
print("Validates Coulomb expansion and solver accuracy")
print("Automatic solver selection: N<1000 → Direct, N≥1000 → Grid")
print("="*80 + "\n")

configs_generated = 0

# ============================================================================
# TEST 1: Coulomb Expansion (Small N, Direct Solver)
# ============================================================================

print("TEST 1: Coulomb Expansion (Direct Solver, N<1000)")
print("-" * 60)

# Small ion counts for direct solver validation
EXPANSION_ION_COUNTS = [100, 500]

for n_ions in EXPANSION_ION_COUNTS:
    # Initial cloud size (tight Gaussian)
    sigma_init_m = 1e-4  # 0.1 mm initial spread
    
    # Domain size (large enough to avoid boundary effects)
    domain_size_m = 0.05  # 5 cm cube
    
    # Simulation time (watch expansion)
    # Time scale: τ ~ (m*σ²)/(q²/4πε₀) for Coulomb expansion
    # For H3O+: τ ~ 1 µs for σ=0.1mm
    total_time_s = 50e-6  # 50 µs
    
    # Timestep (conservative for direct Coulomb)
    dt_s = 1e-9  # 1 ns
    
    # Write frequently to track expansion
    write_interval = 100  # Every 100 ns
    
    config = {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "None",  # Collisionless (pure Coulomb)
            "enable_space_charge": True,
            "enable_reactions": False
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "validation/results/v1.0_test/physics/spacecharge",
            "trajectory_file": f"coulomb_expansion_N{n_ions}_direct.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": "H3O+",
                    "count": n_ions,
                    "position": {
                        "type": "gaussian",
                        "center": [domain_size_m/2, domain_size_m/2, domain_size_m/2],
                        "std": [sigma_init_m, sigma_init_m, sigma_init_m]
                    },
                    "velocity": {
                        "type": "maxwell",
                        "temperature_K": 1.0  # Very cold (minimal thermal motion)
                    }
                }
            ]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"Coulomb expansion (N={n_ions}, Direct solver)",
            "instrument": "Generic",
            "geometry": {
                "type": "box",
                "dimensions": [domain_size_m, domain_size_m, domain_size_m],
                "origin_m": [0.0, 0.0, 0.0]
            },
            "env": {
                "gas_species": "He",
                "temperature_K": 1.0,
                "pressure_Pa": 0.0,  # Perfect vacuum (no collisions)
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {}  # No external fields
        }]
    }
    
    # Write config
    out_file = OUT_DIR / f"coulomb_expansion_N{n_ions}_direct.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1
    
    print(f"✅ N = {n_ions} ions (Direct solver)")
    print(f"   Initial σ: {sigma_init_m*1e3:.2f} mm")
    print(f"   Domain: {domain_size_m*100:.1f} cm cube")
    print(f"   Time: {total_time_s*1e6:.1f} µs")
    print(f"   Expected: Radial expansion ∝ √t")
    print(f"   Output: {out_file.name}")
    print()

# ============================================================================
# TEST 2: Grid Solver (N≥1000)
# ============================================================================

print("TEST 2: Grid-Based Solver (N≥1000)")
print("-" * 60)

GRID_ION_COUNTS = [1000, 5000, 10000]

for n_ions in GRID_ION_COUNTS:
    sigma_init_m = 1e-4  # 0.1 mm
    domain_size_m = 0.05  # 5 cm
    total_time_s = 50e-6  # 50 µs
    dt_s = 1e-9  # 1 ns
    write_interval = 100
    
    config = {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "None",
            "enable_space_charge": True,
            "enable_reactions": False
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "validation/results/v1.0_test/physics/spacecharge",
            "trajectory_file": f"coulomb_expansion_N{n_ions}_grid.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": "H3O+",
                    "count": n_ions,
                    "position": {
                        "type": "gaussian",
                        "center": [domain_size_m/2, domain_size_m/2, domain_size_m/2],
                        "std": [sigma_init_m, sigma_init_m, sigma_init_m]
                    },
                    "velocity": {
                        "type": "maxwell",
                        "temperature_K": 1.0
                    }
                }
            ]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"Coulomb expansion (N={n_ions}, Grid solver)",
            "instrument": "Generic",
            "geometry": {
                "type": "box",
                "dimensions": [domain_size_m, domain_size_m, domain_size_m],
                "origin_m": [0.0, 0.0, 0.0]
            },
            "env": {
                "gas_species": "He",
                "temperature_K": 1.0,
                "pressure_Pa": 0.0,
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {}
        }]
    }
    
    out_file = OUT_DIR / f"coulomb_expansion_N{n_ions}_grid.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1
    
    print(f"✅ N = {n_ions} ions (Grid solver, auto-selected)")
    print(f"   Grid method: CIC charge deposition")
    print(f"   Poisson solver: Auto-select (likely Multigrid)")
    print(f"   Expected: σ(t) ∝ √t (same as Direct)")
    print(f"   Output: {out_file.name}")
    print()

# ============================================================================
# TEST 3: Direct vs Grid Comparison (N=1000, threshold test)
# ============================================================================

print("TEST 3: Direct vs Grid Comparison (N=1000, threshold)")
print("-" * 60)

n_ions = 1000  # Exact threshold
sigma_init_m = 1e-4
domain_size_m = 0.05
total_time_s = 20e-6  # Shorter for comparison
dt_s = 1e-9
write_interval = 50

# Both methods should give similar results at threshold
# This validates the crossover point

config_comparison = {
    "simulation": {
        "total_time_s": total_time_s,
        "dt_s": dt_s,
        "write_interval": write_interval,
        "integrator": "RK4",
        "enable_gpu": False,
        "enable_openmp": True,
        "rng_seed": 42
    },
    "physics": {
        "collision_model": "None",
        "enable_space_charge": True,
        "enable_reactions": False
    },
    "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
    "output": {
        "folder": "validation/results/v1.0_test/physics/spacecharge",
        "trajectory_file": f"threshold_comparison_N{n_ions}.h5",
        "print_progress": True
    },
    "ions": {
        "species": [
            {
                "species_id": "H3O+",
                "count": n_ions,
                "position": {
                    "type": "gaussian",
                    "center": [domain_size_m/2, domain_size_m/2, domain_size_m/2],
                    "std": [sigma_init_m, sigma_init_m, sigma_init_m]
                },
                "velocity": {
                    "type": "maxwell",
                    "temperature_K": 1.0
                }
            }
        ]
    },
    "domains": [{
        "domain_index": 0,
        "name": f"Threshold test (N={n_ions}, auto-select solver)",
        "instrument": "Generic",
        "geometry": {
            "type": "box",
            "dimensions": [domain_size_m, domain_size_m, domain_size_m],
            "origin_m": [0.0, 0.0, 0.0]
        },
        "env": {
            "gas_species": "He",
            "temperature_K": 1.0,
            "pressure_Pa": 0.0,
            "gas_velocity_m_s": [0.0, 0.0, 0.0]
        },
        "fields": {}
    }]
}

out_file = OUT_DIR / f"threshold_comparison_N{n_ions}.json"
with open(out_file, 'w') as f:
    json.dump(config_comparison, f, indent=2)

configs_generated += 1

print(f"✅ N = {n_ions} ions (Threshold test)")
print(f"   Auto-selection: Should choose Grid solver")
print(f"   Time: {total_time_s*1e6:.1f} µs")
print(f"   Purpose: Validate solver crossover")
print(f"   Output: {out_file.name}")
print()

# ============================================================================
# TEST 4: IMS with Space Charge (Realistic scenario)
# ============================================================================

print("TEST 4: IMS with Space Charge (Realistic scenario)")
print("-" * 60)

# IMS parameters
drift_length_m = 0.05  # 5 cm
E_field_V_m = 5000.0  # 5 kV/m
pressure_Pa = 1000.0  # 1 kPa He
temperature_K = 300.0

# Two cases: with and without space charge
for enable_sc in [False, True]:
    n_ions = 5000  # Grid solver
    
    config_ims = {
        "simulation": {
            "total_time_s": 5e-3,  # 5 ms (full drift time)
            "dt_s": 1e-7,  # 100 ns
            "write_interval": 100,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "HSS",
            "enable_space_charge": enable_sc,
            "enable_reactions": False
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "validation/results/v1.0_test/physics/spacecharge",
            "trajectory_file": f"ims_spacecharge_{'on' if enable_sc else 'off'}_N{n_ions}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": "H3O+",
                    "count": n_ions,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.001],  # Start at z=1mm
                        "std": [0.002, 0.002, 0.0005]  # Compact initial cloud
                    },
                    "velocity": {
                        "type": "maxwell",
                        "temperature_K": temperature_K
                    }
                }
            ]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"IMS with space charge {'ON' if enable_sc else 'OFF'}",
            "instrument": "IMS",
            "geometry": {
                "type": "cylinder",
                "length_m": drift_length_m,
                "radius_m": 0.025,
                "origin_m": [0.0, 0.0, 0.0]
            },
            "env": {
                "gas_species": "He",
                "temperature_K": temperature_K,
                "pressure_Pa": pressure_Pa,
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {
                "DC": {
                    "type": "uniform",
                    "field_V_m": [0.0, 0.0, E_field_V_m]
                }
            }
        }]
    }
    
    out_file = OUT_DIR / f"ims_spacecharge_{'on' if enable_sc else 'off'}_N{n_ions}.json"
    with open(out_file, 'w') as f:
        json.dump(config_ims, f, indent=2)
    
    configs_generated += 1
    
    print(f"✅ IMS with space charge {'ON' if enable_sc else 'OFF'}")
    print(f"   N = {n_ions} ions (Grid solver)")
    print(f"   E-field: {E_field_V_m} V/m")
    print(f"   Expected: Space charge → slower drift + broadening")
    print(f"   Output: {out_file.name}")
    print()

print("="*80)
print(f"Total configurations: {configs_generated}")
print(f"Output: {OUT_DIR}")
print("\nTest categories:")
print("  1. Coulomb expansion (N=100, 500) - Direct solver")
print("  2. Coulomb expansion (N=1000, 5000, 10000) - Grid solver")
print("  3. Threshold test (N=1000) - Auto-selection validation")
print("  4. IMS comparison (space charge ON/OFF) - Realistic scenario")
print("\nExpected results:")
print("  - Radial spread: σ(t) ∝ √t (Coulomb expansion)")
print("  - Direct ≈ Grid at N=1000 (within 30% tolerance)")
print("  - IMS with SC: Longer drift time + broader peak")
print("\nNext steps:")
print("  1. Run: for cfg in configs/physics/spacecharge/*.json; do")
print("             icarion_main $cfg; done")
print("  2. Analyze: python3 analyze_spacecharge.py")
print("  3. Validate: σ(t) scaling, Direct vs Grid agreement")
print("="*80 + "\n")
