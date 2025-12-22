#!/usr/bin/env python3
"""
GPU Performance Benchmark Analysis

Parses HDF5 results from GPU performance tests and generates:
1. CPU vs GPU speedup plots
2. Integrator comparison charts
3. Threshold validation analysis
4. Long simulation efficiency metrics

Requires:
- h5py
- numpy
- matplotlib
"""

import sys
import json
import argparse
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

try:
    import h5py
except ImportError:
    print("ERROR: h5py not installed. Install with: pip install h5py")
    sys.exit(1)

# Configuration
SCRIPT_DIR = Path(__file__).parent
VALIDATION_ROOT = SCRIPT_DIR.parent.parent
PROJECT_ROOT = VALIDATION_ROOT.parent
RESULTS_ROOT = VALIDATION_ROOT / "results" / "v1.0_test" / "performance"
CPU_RESULTS_DIR = RESULTS_ROOT / "logs"
GPU_RESULTS_DIR = RESULTS_ROOT / "gpu_logs"
OUTPUT_DIR = RESULTS_ROOT / "gpu_analysis"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("\n" + "="*80)
print("GPU Performance Benchmark Analysis")
print("="*80)
print(f"CPU logs: {CPU_RESULTS_DIR}")
print(f"GPU logs: {GPU_RESULTS_DIR}")
print(f"Output directory:  {OUTPUT_DIR}")
print("="*80 + "\n")

# ============================================================================
# Helper Functions
# ============================================================================

def parse_hdf5_timing(h5_file):
    """Extract timing information from HDF5 file"""
    try:
        with h5py.File(h5_file, 'r') as f:
            # Try to extract wall time from metadata
            if 'metadata' in f.attrs:
                metadata = json.loads(f.attrs['metadata'])
                if 'wall_time_s' in metadata:
                    return metadata['wall_time_s']
            
            # Fallback: estimate from trajectories
            if 'trajectories' in f and 'time' in f['trajectories']:
                times = f['trajectories']['time'][:]
                n_steps = len(times)
                return None  # Can't determine wall time from trajectories alone
    except Exception as e:
        print(f"⚠️  Error parsing {h5_file}: {e}")
        return None
    
    return None

TIME_PATTERNS = [
    re.compile(r"CPU time:\s*([0-9.+\-Ee]+)"),
    re.compile(r"Wall time:\s*([0-9.+\-Ee]+)"),
    re.compile(r"Total time:\s*([0-9.+\-Ee]+)")
]

def parse_log_timing(log_file):
    """Extract timing information from log file"""
    try:
        text = Path(log_file).read_text()
    except Exception as e:
        print(f"⚠️  Error reading log {log_file}: {e}")
        return None

    for pattern in TIME_PATTERNS:
        match = pattern.search(text)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                continue
    return None

def parse_time_file(time_file):
    """Read wall time (seconds) from the companion .time file."""
    try:
        text = Path(time_file).read_text().strip()
        if text:
            return float(text)
    except Exception as e:
        print(f"⚠️  Error parsing time file {time_file}: {e}")
    return None

def resolve_artifact(base_name, extension, prefer_gpu=False):
    search_dirs = [GPU_RESULTS_DIR, CPU_RESULTS_DIR] if prefer_gpu else [CPU_RESULTS_DIR, GPU_RESULTS_DIR]
    for directory in search_dirs:
        candidate = directory / f"{base_name}.{extension}"
        if candidate.exists():
            return candidate
    return None

def get_wall_time(base_name, prefer_gpu=False):
    log_file = resolve_artifact(base_name, "log", prefer_gpu)
    time_file = resolve_artifact(base_name, "time", prefer_gpu)
    wall_time = None
    if log_file:
        wall_time = parse_log_timing(log_file)
    if wall_time is None and time_file:
        wall_time = parse_time_file(time_file)
    return log_file, wall_time

def extract_timing_data(pattern, results_dir=GPU_RESULTS_DIR):
    """Extract timing data for configs matching pattern"""
    timings = {}
    
    for log_file in results_dir.glob(f"{pattern}.log"):
        config_name = log_file.stem
        wall_time = parse_log_timing(log_file)
        
        if wall_time is not None:
            timings[config_name] = wall_time
        else:
            print(f"⚠️  No timing data for {config_name}")
    
    return timings

# ============================================================================
# ANALYSIS 1: GPU vs CPU Speedup
# ============================================================================

print("ANALYSIS 1: GPU vs CPU Speedup")
print("-" * 60)

ION_COUNTS = [1000, 5000, 10000, 50000, 100000]
cpu_times = []
gpu_times = []
speedups = []

for n_ions in ION_COUNTS:
    cpu_log, cpu_time = get_wall_time(f"RK4_cpu_N{n_ions}")
    gpu_log, gpu_time = get_wall_time(f"RK4_gpu_N{n_ions}", prefer_gpu=True)
    
    if cpu_time is not None and gpu_time is not None:
        speedup = cpu_time / gpu_time
        cpu_times.append(cpu_time)
        gpu_times.append(gpu_time)
        speedups.append(speedup)
        
        print(f"N={n_ions:6d}: CPU={cpu_time:8.3f}s, GPU={gpu_time:8.3f}s, Speedup={speedup:5.2f}×")
    else:
        cpu_times.append(np.nan)
        gpu_times.append(np.nan)
        speedups.append(np.nan)
        
        status = []
        if not cpu_log:
            status.append("CPU log missing")
        elif cpu_time is None:
            status.append("CPU timing missing")
        if not gpu_log:
            status.append("GPU log missing")
        elif gpu_time is None:
            status.append("GPU timing missing")
        
        print(f"N={n_ions:6d}: ⚠️  {', '.join(status)}")

# Plot speedup
if any(not np.isnan(s) for s in speedups):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Speedup plot
    valid_idx = [i for i, s in enumerate(speedups) if not np.isnan(s)]
    valid_ions = [ION_COUNTS[i] for i in valid_idx]
    valid_speedups = [speedups[i] for i in valid_idx]
    
    ax1.plot(valid_ions, valid_speedups, 'o-', linewidth=2, markersize=8, label='RK4 GPU Speedup')
    ax1.axhline(y=1.0, color='gray', linestyle='--', label='No Speedup')
    ax1.set_xscale('log')
    ax1.set_xlabel('Number of Ions', fontsize=12)
    ax1.set_ylabel('Speedup (CPU time / GPU time)', fontsize=12)
    ax1.set_title('GPU vs CPU Speedup (RK4)', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Absolute timing plot
    valid_cpu = [cpu_times[i] for i in valid_idx]
    valid_gpu = [gpu_times[i] for i in valid_idx]
    
    ax2.plot(valid_ions, valid_cpu, 'o-', linewidth=2, markersize=8, label='CPU')
    ax2.plot(valid_ions, valid_gpu, 's-', linewidth=2, markersize=8, label='GPU')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Number of Ions', fontsize=12)
    ax2.set_ylabel('Wall Time (s)', fontsize=12)
    ax2.set_title('Absolute Performance (RK4)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'gpu_vs_cpu_speedup.png', dpi=150)
    print(f"✅ Saved speedup plot: {OUTPUT_DIR / 'gpu_vs_cpu_speedup.png'}")
else:
    print("⚠️  No valid speedup data, skipping plot")

# ============================================================================
# ANALYSIS 2: Integrator Comparison
# ============================================================================

print("\n" + "="*80)
print("ANALYSIS 2: GPU Integrator Comparison")
print("-" * 60)

INTEGRATORS = ["rk4", "rk45", "boris"]
integrator_times = {}
N_IONS = 10000

for integrator in INTEGRATORS:
    log_file, wall_time = get_wall_time(f"integrator_{integrator}_gpu_N{N_IONS}", prefer_gpu=True)
    
    if log_file and wall_time is not None:
        if wall_time:
            integrator_times[integrator.upper()] = wall_time
            print(f"{integrator.upper():6s}: {wall_time:8.3f} s")
        else:
            print(f"{integrator.upper():6s}: ⚠️  No timing data")
    else:
        print(f"{integrator.upper():6s}: ⚠️  Log file missing")

# Plot integrator comparison
if integrator_times:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    integrators_list = list(integrator_times.keys())
    times_list = list(integrator_times.values())
    
    bars = ax.bar(integrators_list, times_list, color=['#1f77b4', '#ff7f0e', '#2ca02c'], alpha=0.8)
    
    # Add value labels on bars
    for bar, time in zip(bars, times_list):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{time:.3f}s',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.set_ylabel('Wall Time (s)', fontsize=12)
    ax.set_title(f'GPU Integrator Comparison (N={N_IONS})', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'integrator_comparison.png', dpi=150)
    print(f"✅ Saved integrator comparison: {OUTPUT_DIR / 'integrator_comparison.png'}")
else:
    print("⚠️  No integrator timing data, skipping plot")

# ============================================================================
# ANALYSIS 3: Threshold Validation
# ============================================================================

print("\n" + "="*80)
print("ANALYSIS 3: GPU Threshold Validation")
print("-" * 60)

THRESHOLD_DATA = {
    "RK4": [4000, 4500, 5000, 5500, 6000],
    "RK45": [4000, 4500, 5000, 5500, 6000],
    "Boris": [2000, 2250, 2500, 2750, 3000]
}

threshold_results = {}

for integrator, ion_counts in THRESHOLD_DATA.items():
    integrator_lower = integrator.lower()
    times = []
    
    print(f"\n{integrator}:")
    for n_ions in ion_counts:
        log_file, wall_time = get_wall_time(f"threshold_{integrator_lower}_N{n_ions}", prefer_gpu=True)
    
        if log_file and wall_time is not None:
            if wall_time:
                times.append(wall_time)
                print(f"  N={n_ions:5d}: {wall_time:8.3f} s")
            else:
                times.append(np.nan)
                print(f"  N={n_ions:5d}: ⚠️  No timing data")
        else:
            times.append(np.nan)
            print(f"  N={n_ions:5d}: ⚠️  Log missing")
    
    threshold_results[integrator] = (ion_counts, times)

# Plot threshold validation
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, (integrator, (ion_counts, times)) in zip(axes, threshold_results.items()):
    valid_idx = [i for i, t in enumerate(times) if not np.isnan(t)]
    valid_ions = [ion_counts[i] for i in valid_idx]
    valid_times = [times[i] for i in valid_idx]
    
    if valid_ions:
        ax.plot(valid_ions, valid_times, 'o-', linewidth=2, markersize=8)
        
        # Mark default threshold
        if integrator == "Boris":
            ax.axvline(x=2500, color='red', linestyle='--', label='Default Threshold (2500)')
        else:
            ax.axvline(x=5000, color='red', linestyle='--', label='Default Threshold (5000)')
        
        ax.set_xlabel('Number of Ions', fontsize=12)
        ax.set_ylabel('Wall Time (s)', fontsize=12)
        ax.set_title(f'{integrator} Threshold Test', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', fontsize=14)
        ax.set_title(f'{integrator} Threshold Test', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'threshold_validation.png', dpi=150)
print(f"\n✅ Saved threshold validation: {OUTPUT_DIR / 'threshold_validation.png'}")

# ============================================================================
# ANALYSIS 4: Long Simulation Efficiency
# ============================================================================

print("\n" + "="*80)
print("ANALYSIS 4: Long Simulation Efficiency")
print("-" * 60)

long_times = {}

for integrator in ["rk4", "rk45", "boris"]:
    log_file, wall_time = get_wall_time(f"long_{integrator}_gpu_N10000", prefer_gpu=True)
    
    if log_file and wall_time is not None:
        if wall_time:
            long_times[integrator.upper()] = wall_time
            print(f"{integrator.upper():6s}: {wall_time:8.3f} s (100 µs simulation)")
        else:
            print(f"{integrator.upper():6s}: ⚠️  No timing data")
    else:
        print(f"{integrator.upper():6s}: ⚠️  Log missing")

# Plot long simulation efficiency
if long_times:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    integrators_list = list(long_times.keys())
    times_list = list(long_times.values())
    
    bars = ax.bar(integrators_list, times_list, color=['#1f77b4', '#ff7f0e', '#2ca02c'], alpha=0.8)
    
    for bar, time in zip(bars, times_list):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{time:.3f}s',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.set_ylabel('Wall Time (s)', fontsize=12)
    ax.set_title('Long Simulation Efficiency (100 µs, N=10000)', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'long_simulation_efficiency.png', dpi=150)
    print(f"✅ Saved long simulation plot: {OUTPUT_DIR / 'long_simulation_efficiency.png'}")
else:
    print("⚠️  No long simulation data, skipping plot")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Output directory: {OUTPUT_DIR}")
print("\nGenerated plots:")
if (OUTPUT_DIR / 'gpu_vs_cpu_speedup.png').exists():
    print("  ✅ gpu_vs_cpu_speedup.png")
if (OUTPUT_DIR / 'integrator_comparison.png').exists():
    print("  ✅ integrator_comparison.png")
if (OUTPUT_DIR / 'threshold_validation.png').exists():
    print("  ✅ threshold_validation.png")
if (OUTPUT_DIR / 'long_simulation_efficiency.png').exists():
    print("  ✅ long_simulation_efficiency.png")
print("="*80 + "\n")
