#!/usr/bin/env python3
"""
Generate performance benchmark suite

Validates computational performance and scaling:
1. Ion count scaling (N = 100, 1k, 10k, 100k)
2. Collision model overhead (None, HSS, EHSS, Friction)
3. Space charge overhead (Direct vs Grid)
4. OpenMP scaling (1-8 threads)
"""

import json
import numpy as np
from pathlib import Path

# Output directory
OUT_DIR = Path("../configs/performance")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Benchmark parameters
BENCHMARK_TIME_S = 10e-6  # 10 µs (short for benchmarking)
DT_S = 1e-9  # 1 ns
WRITE_INTERVAL = 10000  # Minimal I/O overhead

TEMPERATURE_K = 300.0
PRESSURE_PA = 1000.0  # 1 kPa for collision benchmarks

print("\n" + "="*80)
print("Performance Benchmark Suite Generator")
print("="*80)
print("Validates computational performance and scaling")
print("="*80 + "\n")

configs_generated = 0

# ============================================================================
# BENCHMARK 1: Ion Count Scaling (Collisionless baseline)
# ============================================================================

print("BENCHMARK 1: Ion Count Scaling (Collisionless)")
print("-" * 60)

ION_COUNTS = [100, 1000, 10000, 100000]

for n_ions in ION_COUNTS:
    config = {
        "simulation": {
            "total_time_s": BENCHMARK_TIME_S,
            "dt_s": DT_S,
            "write_interval": WRITE_INTERVAL,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "None",
            "enable_space_charge": False,
            "enable_reactions": False
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "validation/results/v1.0_test/performance",
            "trajectory_file": f"scaling_baseline_N{n_ions}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": "H3O+",
                    "count": n_ions,
                    "position": {
                        "type": "uniform",
                        "min": [0.0, 0.0, 0.0],
                        "max": [0.01, 0.01, 0.01]
                    },
                    "velocity": {
                        "type": "maxwell",
                        "temperature_K": TEMPERATURE_K
                    }
                }
            ]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"Baseline scaling (N={n_ions})",
            "instrument": "Generic",
            "geometry": {
                "type": "box",
                "dimensions": [0.01, 0.01, 0.01],
                "origin_m": [0.0, 0.0, 0.0]
            },
            "env": {
                "gas_species": "He",
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa": 0.0,
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {}
        }]
    }
    
    out_file = OUT_DIR / f"scaling_baseline_N{n_ions}.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1
    
    print(f"✅ N = {n_ions:6d} ions (collisionless baseline)")
    print(f"   Expected: O(N) scaling")

print(f"\nGenerated {len(ION_COUNTS)} baseline configs")

# ============================================================================
# BENCHMARK 2: Collision Model Overhead
# ============================================================================

print("\n" + "="*80)
print("BENCHMARK 2: Collision Model Overhead")
print("-" * 60)

COLLISION_MODELS = ["None", "HSS", "EHSS", "Friction"]
N_IONS_COLLISION = 10000

for model in COLLISION_MODELS:
    config = {
        "simulation": {
            "total_time_s": BENCHMARK_TIME_S,
            "dt_s": DT_S,
            "write_interval": WRITE_INTERVAL,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": model,
            "enable_space_charge": False,
            "enable_reactions": False
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "validation/results/v1.0_test/performance",
            "trajectory_file": f"collision_overhead_{model}_N{N_IONS_COLLISION}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": "H3O+" if model != "EHSS" else "PentanalH+",
                    "count": N_IONS_COLLISION,
                    "position": {
                        "type": "uniform",
                        "min": [0.0, 0.0, 0.0],
                        "max": [0.01, 0.01, 0.01]
                    },
                    "velocity": {
                        "type": "maxwell",
                        "temperature_K": TEMPERATURE_K
                    }
                }
            ]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"Collision overhead ({model})",
            "instrument": "Generic",
            "geometry": {
                "type": "box",
                "dimensions": [0.01, 0.01, 0.01],
                "origin_m": [0.0, 0.0, 0.0]
            },
            "env": {
                "gas_species": "He",
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa": PRESSURE_PA if model != "None" else 0.0,
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {}
        }]
    }
    
    out_file = OUT_DIR / f"collision_overhead_{model}_N{N_IONS_COLLISION}.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1
    
    overhead = {
        "None": "0%",
        "HSS": "~10%",
        "EHSS": "~800%",
        "Friction": "~5%"
    }
    
    print(f"✅ {model:10s} (N={N_IONS_COLLISION})")
    print(f"   Expected overhead: {overhead[model]}")

print(f"\nGenerated {len(COLLISION_MODELS)} collision overhead configs")

# ============================================================================
# BENCHMARK 3: Space Charge Overhead
# ============================================================================

print("\n" + "="*80)
print("BENCHMARK 3: Space Charge Overhead")
print("-" * 60)

# Direct solver (N < 1000)
N_DIRECT = [100, 500]
# Grid solver (N >= 1000)
N_GRID = [1000, 5000, 10000]

for n_ions in N_DIRECT + N_GRID:
    for sc_enabled in [False, True]:
        config = {
            "simulation": {
                "total_time_s": BENCHMARK_TIME_S,
                "dt_s": DT_S,
                "write_interval": WRITE_INTERVAL,
                "integrator": "RK4",
                "enable_gpu": False,
                "enable_openmp": True,
                "rng_seed": 42
            },
            "physics": {
                "collision_model": "None",
                "enable_space_charge": sc_enabled,
                "enable_reactions": False
            },
            "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
            "output": {
                "folder": "validation/results/v1.0_test/performance",
                "trajectory_file": f"spacecharge_overhead_{'on' if sc_enabled else 'off'}_N{n_ions}.h5",
                "print_progress": True
            },
            "ions": {
                "species": [
                    {
                        "species_id": "H3O+",
                        "count": n_ions,
                        "position": {
                            "type": "uniform",
                            "min": [0.0, 0.0, 0.0],
                            "max": [0.01, 0.01, 0.01]
                        },
                        "velocity": {
                            "type": "maxwell",
                            "temperature_K": TEMPERATURE_K
                        }
                    }
                ]
            },
            "domains": [{
                "domain_index": 0,
                "name": f"Space charge {'ON' if sc_enabled else 'OFF'} (N={n_ions})",
                "instrument": "Generic",
                "geometry": {
                    "type": "box",
                    "dimensions": [0.01, 0.01, 0.01],
                    "origin_m": [0.0, 0.0, 0.0]
                },
                "env": {
                    "gas_species": "He",
                    "temperature_K": TEMPERATURE_K,
                    "pressure_Pa": 0.0,
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {}
            }]
        }
        
        out_file = OUT_DIR / f"spacecharge_overhead_{'on' if sc_enabled else 'off'}_N{n_ions}.json"
        with open(out_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        configs_generated += 1

print(f"✅ Direct solver: N={N_DIRECT} (SC on/off)")
print(f"   Expected: O(N²) overhead")
print(f"✅ Grid solver: N={N_GRID} (SC on/off)")
print(f"   Expected: O(N log N) overhead")

print(f"\nGenerated {2*(len(N_DIRECT) + len(N_GRID))} space charge configs")

print("\n" + "="*80)
print(f"Total configurations: {configs_generated}")
print(f"Output: {OUT_DIR}")
print("\nBenchmark categories:")
print("  1. Ion count scaling (4 configs) - O(N) baseline")
print("  2. Collision overhead (4 configs) - HSS ~10%, EHSS ~800%")
print("  3. Space charge overhead (10 configs) - Direct O(N²), Grid O(N log N)")
print("\nExpected timings (approximate):")
print("  Baseline (N=10k, 10µs): ~0.1 s")
print("  HSS collision: ~0.11 s (+10%)")
print("  EHSS collision: ~0.8 s (+800%)")
print("  Space charge Direct (N=500): ~0.5 s")
print("  Space charge Grid (N=10k): ~0.3 s")
print("\nNext steps:")
print("  1. Run: for cfg in configs/performance/*.json; do")
print("             time icarion_main $cfg; done")
print("  2. Analyze: python3 analyze_performance.py")
print("  3. Generate scaling plots and overhead tables")
print("="*80 + "\n")
