#!/usr/bin/env python3
"""
Generate GPU performance benchmark suite

Validates GPU acceleration performance:
1. GPU vs CPU comparison (N = 1k, 5k, 10k, 50k, 100k)
2. Integrator comparison (RK4, RK45, Boris) on GPU
3. GPU threshold validation (N near threshold)
4. Multi-timestep GPU efficiency
"""

import json
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent


def find_validation_dir() -> Path:
    for parent in SCRIPT_DIR.parents:
        if parent.name == "validation":
            return parent
    raise RuntimeError("Unable to locate 'validation' directory relative to script")


VALIDATION_DIR = find_validation_dir()
OUT_DIR = VALIDATION_DIR / "configs" / "performance" / "gpu"
RESULTS_BASE = VALIDATION_DIR / "results" / "performance" / "gpu"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Benchmark parameters
BENCHMARK_TIME_S = 10e-6  # 10 µs (short for benchmarking)
DT_S = 1e-9  # 1 ns
WRITE_INTERVAL = 10000  # Minimal I/O overhead

TEMPERATURE_K = 300.0
PRESSURE_PA = 1e-6  # Minimal pressure (collisionless regime, MFP >> domain)

print("\n" + "="*80)
print("GPU Performance Benchmark Suite Generator")
print("="*80)
print("Validates GPU acceleration and scaling")
print("="*80 + "\n")

configs_generated = 0

# ============================================================================
# BENCHMARK 1: GPU vs CPU Scaling Comparison
# ============================================================================

print("BENCHMARK 1: GPU vs CPU Scaling Comparison")
print("-" * 60)

ION_COUNTS = [1000, 5000, 10000, 50000, 100000]
INTEGRATORS = ["RK4"]  # Start with RK4, add others later

for integrator in INTEGRATORS:
    for n_ions in ION_COUNTS:
        # CPU version
        config_cpu = {
            "simulation": {
                "total_time_s": BENCHMARK_TIME_S,
                "dt_s": DT_S,
                "write_interval": WRITE_INTERVAL,
                "integrator": integrator,
                "enable_gpu": False,
                "enable_openmp": True,
                "rng_seed": 42
            },
            "physics": {
                "collision_model": "NoCollisions",
                "enable_space_charge": False
            },
            "output": {
                "folder": str(RESULTS_BASE / integrator),
                "trajectory_file": f"cpu_N{n_ions}.h5",
                "print_progress": True
            },
            "ions": {
                "species": [
                    {
                        "id": "H3O+",
                        "count": n_ions,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.005],
                            "std": [0.001, 0.001, 0.001]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": TEMPERATURE_K
                        }
                    }
                ]
            },
            "domains": [{
                "name": f"CPU Benchmark (N={n_ions})",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, 0.0],
                    "length_m": 0.05,
                    "radius_m": 0.01
                },
                "env": {
                    "temperature_K": TEMPERATURE_K,
                    "pressure_Pa": PRESSURE_PA,
                    "gas_species": "He",
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "DC": {
                        "EN_Td": 0.0,
                        "axial_V": 0.0
                    }
                }
            }]
        }
        
        # GPU version
        config_gpu = json.loads(json.dumps(config_cpu))
        config_gpu["simulation"]["enable_gpu"] = True
        config_gpu["simulation"]["gpu_threshold"] = 1000  # Force GPU usage
        config_gpu["output"]["trajectory_file"] = f"gpu_N{n_ions}.h5"
        config_gpu["domains"][0]["name"] = f"GPU Benchmark (N={n_ions})"
        
        # Save configs
        out_cpu = OUT_DIR / f"{integrator}_cpu_N{n_ions}.json"
        out_gpu = OUT_DIR / f"{integrator}_gpu_N{n_ions}.json"
        
        with open(out_cpu, 'w') as f:
            json.dump(config_cpu, f, indent=2)
        with open(out_gpu, 'w') as f:
            json.dump(config_gpu, f, indent=2)
        
        configs_generated += 2
        
        print(f"✅ {integrator:6s}: N = {n_ions:6d} ions (CPU + GPU)")

print(f"\nGenerated {configs_generated} configs (CPU+GPU pairs)")

# ============================================================================
# BENCHMARK 2: GPU Integrator Comparison
# ============================================================================

print("\n" + "="*80)
print("BENCHMARK 2: GPU Integrator Comparison (RK4 vs RK45 vs Boris)")
print("-" * 60)

INTEGRATORS_ALL = ["RK4", "RK45", "Boris"]
N_IONS_INTEGRATOR = 10000

for integrator in INTEGRATORS_ALL:
    config = {
        "simulation": {
            "total_time_s": BENCHMARK_TIME_S,
            "dt_s": DT_S,
            "write_interval": WRITE_INTERVAL,
            "integrator": integrator,
            "enable_gpu": True,
            "gpu_threshold": 1000,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "NoCollisions",
            "enable_space_charge": False
        },
        "output": {
            "folder": str(RESULTS_BASE / "integrators"),
            "trajectory_file": f"{integrator.lower()}_gpu_N{N_IONS_INTEGRATOR}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": N_IONS_INTEGRATOR,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.005],
                        "std": [0.001, 0.001, 0.001]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": TEMPERATURE_K
                    }
                }
            ]
        },
        "domains": [{
            "name": f"{integrator} GPU Benchmark",
            "instrument": "IMS",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": 0.05,
                "radius_m": 0.01
            },
            "env": {
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa": PRESSURE_PA,
                "gas_species": "He",
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {
                "DC": {
                    "EN_Td": 0.0,
                    "axial_V": 0.0
                }
            }
        }]
    }
    
    # Add magnetic field for Boris test
    if integrator == "Boris":
        config["domains"][0]["fields"]["magnetic"] = {
            "enabled": True,
            "field_strength_T": [0.0, 0.0, 1.0]
        }
    
    out_file = OUT_DIR / f"integrator_{integrator.lower()}_gpu_N{N_IONS_INTEGRATOR}.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1
    
    print(f"✅ {integrator:6s}: GPU benchmark @ N={N_IONS_INTEGRATOR}")

print(f"\nGenerated {len(INTEGRATORS_ALL)} integrator comparison configs")

# ============================================================================
# BENCHMARK 3: GPU Threshold Validation
# ============================================================================

print("\n" + "="*80)
print("BENCHMARK 3: GPU Threshold Validation")
print("-" * 60)

# Test around default thresholds:
# - Boris: 2500 (threshold/2)
# - RK4/RK45: 5000 (threshold)

THRESHOLD_TESTS = [
    ("RK4", [4000, 4500, 5000, 5500, 6000]),
    ("RK45", [4000, 4500, 5000, 5500, 6000]),
    ("Boris", [2000, 2250, 2500, 2750, 3000])
]

for integrator, ion_counts in THRESHOLD_TESTS:
    for n_ions in ion_counts:
        config = {
            "simulation": {
                "total_time_s": BENCHMARK_TIME_S,
                "dt_s": DT_S,
                "write_interval": WRITE_INTERVAL,
                "integrator": integrator,
                "enable_gpu": True,
                "gpu_threshold": 5000,  # Default
                "enable_openmp": True,
                "rng_seed": 42
            },
            "physics": {
                "collision_model": "NoCollisions",
                "enable_space_charge": False
            },
            "output": {
                "folder": str(RESULTS_BASE / "threshold"),
                "trajectory_file": f"{integrator.lower()}_threshold_N{n_ions}.h5",
                "print_progress": True
            },
            "ions": {
                "species": [
                    {
                        "id": "H3O+",
                        "count": n_ions,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.005],
                            "std": [0.001, 0.001, 0.001]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": TEMPERATURE_K
                        }
                    }
                ]
            },
            "domains": [{
                "name": f"{integrator} Threshold Test (N={n_ions})",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, 0.0],
                    "length_m": 0.05,
                    "radius_m": 0.01
                },
                "env": {
                    "temperature_K": TEMPERATURE_K,
                    "pressure_Pa": PRESSURE_PA,
                    "gas_species": "He",
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "DC": {
                        "EN_Td": 0.0,
                        "axial_V": 0.0
                    }
                }
            }]
        }
        
        if integrator == "Boris":
            config["domains"][0]["fields"]["magnetic"] = {
                "enabled": True,
                "field_strength_T": [0.0, 0.0, 1.0]
            }
        
        out_file = OUT_DIR / f"threshold_{integrator.lower()}_N{n_ions}.json"
        with open(out_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        configs_generated += 1

print(f"✅ Generated threshold validation configs for RK4, RK45, Boris")

# ============================================================================
# BENCHMARK 4: Long Simulation GPU Efficiency
# ============================================================================

print("\n" + "="*80)
print("BENCHMARK 4: Long Simulation GPU Efficiency")
print("-" * 60)

LONG_BENCHMARK_TIME_S = 100e-6  # 100 µs (10× longer)
N_IONS_LONG = 10000

for integrator in INTEGRATORS_ALL:
    config = {
        "simulation": {
            "total_time_s": LONG_BENCHMARK_TIME_S,
            "dt_s": DT_S,
            "write_interval": 10000,
            "integrator": integrator,
            "enable_gpu": True,
            "gpu_threshold": 1000,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "NoCollisions",
            "enable_space_charge": False
        },
        "output": {
            "folder": str(RESULTS_BASE / "long"),
            "trajectory_file": f"{integrator.lower()}_long_gpu_N{N_IONS_LONG}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": N_IONS_LONG,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.005],
                        "std": [0.001, 0.001, 0.001]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": TEMPERATURE_K
                    }
                }
            ]
        },
        "domains": [{
            "name": f"{integrator} Long GPU Test",
            "instrument": "IMS",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": 0.05,
                "radius_m": 0.01
            },
            "env": {
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa": PRESSURE_PA,
                "gas_species": "He",
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {
                "DC": {
                    "EN_Td": 0.0,
                    "axial_V": 0.0
                }
            }
        }]
    }
    
    if integrator == "Boris":
        config["domains"][0]["fields"]["magnetic"] = {
            "enabled": True,
            "field_strength_T": [0.0, 0.0, 1.0]
        }
    
    out_file = OUT_DIR / f"long_{integrator.lower()}_gpu_N{N_IONS_LONG}.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1

print(f"✅ Generated long simulation configs for {len(INTEGRATORS_ALL)} integrators")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "="*80)
print("SUMMARY: GPU Performance Benchmark Suite")
print("="*80)
print(f"Total configs generated: {configs_generated}")
print(f"Output directory: {OUT_DIR}")
print("\nBenchmark Categories:")
print("  1. GPU vs CPU Scaling:     10 configs (5 ion counts × 2 modes)")
print("  2. Integrator Comparison:   3 configs (RK4, RK45, Boris)")
print("  3. Threshold Validation:   15 configs (3 integrators × 5 counts)")
print("  4. Long Simulation:         3 configs (RK4, RK45, Boris)")
print(f"\nTotal: {configs_generated} configs")
print("="*80)
