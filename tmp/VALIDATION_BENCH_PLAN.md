# 🎯 ICARION Validation Bench - Implementation Plan

**Goal:** Publication-quality physics validation suite  
**Strategy:** Small increments, test after each config, document progress  
**Location:** `validation/` directory  
**Documentation:** `validation/README.md` (updated after each session)

---

## 📊 TEST OVERVIEW (20 Validation Tests)

| Category | Tests | Runtime | Priority |
|----------|-------|---------|----------|
| **Instrument Physics** | 5 (IMS, TOF, Orbitrap, LQIT, Energy) | ~15 min | ⭐⭐⭐ Critical |
| **Collision Physics** | 6 (HSS/EHSS thermalization) | ~30 min | ⭐⭐⭐ Critical |
| **Transport** | 2 (gas flow, combined) | ~5 min | ⭐⭐ Important |
| **Space Charge** | 2 (1D, radial) | ~10 min | ⭐⭐ Important |
| **Reactions** | 1 (A+ → B+ decay) | ~2 min | ⭐⭐ Important |
| **Performance** | 4 (scaling benchmarks) | ~45 min | ⭐ Useful |
| **TOTAL** | **20 tests** | **~2 hours** | - |

**Note:** This is a publication-quality validation suite combining physics accuracy + performance benchmarking

---

## 📁 DIRECTORY STRUCTURE

```
validation/
├── README.md                           # Progress tracker (updated after each session)
├── run_validation_suite.sh             # Master script (Session 4)
│
├── configs/
│   ├── instruments/                    # Session 1: Instrument Physics
│   │   ├── ims_mobility.json           # IMS drift velocity
│   │   ├── tof_flight_time.json        # TOF mass separation
│   │   ├── orbitrap_frequency.json     # Orbitrap axial frequency
│   │   ├── lqit_stability.json         # LQIT Mathieu stability
│   │   └── energy_conservation.json    # Vacuum acceleration
│   │
│   ├── physics/
│   │   ├── thermalization/             # Session 2: Collision Physics
│   │   │   ├── hss_300K_1atm.json      # HSS baseline
│   │   │   ├── hss_450K_1atm.json      # HSS high T
│   │   │   ├── hss_300K_lowP.json      # HSS low P
│   │   │   ├── ehss_300K_1atm.json     # EHSS baseline
│   │   │   ├── ehss_450K_1atm.json     # EHSS high T
│   │   │   └── ehss_300K_lowP.json     # EHSS low P
│   │   │
│   │   ├── transport/                  # Session 3: Transport
│   │   │   ├── gas_flow_only.json      # Pure gas drag
│   │   │   └── gas_flow_drift.json     # Gas + E-field
│   │   │
│   │   ├── spacecharge/                # Session 3: Space Charge
│   │   │   ├── 1d_expansion.json       # Axial expansion
│   │   │   └── radial_defocusing.json  # Radial spread
│   │   │
│   │   └── reactions/                  # Session 3: Reactions
│   │       ├── simple_reaction.json    # A+ → B+ kinetics
│   │       └── reaction_validation.json # Reaction database
│   │
│   └── performance/                    # Session 4: Performance
│       ├── scaling_N_100.json          # N=100 ions
│       ├── scaling_N_1000.json         # N=1k ions
│       ├── scaling_N_10000.json        # N=10k ions
│       └── scaling_N_100000.json       # N=100k ions
│
├── scripts/
│   ├── test_single_config.sh           # Quick test script
│   ├── run_instrument_tests.sh         # Session 1: Instrument physics
│   ├── run_thermalization_tests.sh     # Session 2: Collision physics
│   ├── run_transport_tests.sh          # Session 3: Transport
│   ├── run_spacecharge_tests.sh        # Session 3: Space charge
│   ├── run_reactions_tests.sh          # Session 3: Reactions
│   ├── run_performance_tests.sh        # Session 4: Performance benchmarks
│   ├── run_all_validation.sh           # Master script (all sessions)
│   │
│   └── utils/
│       ├── extract_instruments.py      # IMS mobility, TOF timing, orbit frequency
│       ├── extract_thermalization.py   # Extract T_final
│       ├── extract_transport.py        # Extract v_drift
│       ├── extract_spacecharge.py      # Extract cloud size
│       ├── extract_reactions.py        # Extract species counts
│       ├── extract_performance.py      # CPU time scaling, integrator comparison
│       └── generate_report.py          # Generate Markdown report
│
├── reference_data/
│   ├── instruments_theory.csv          # IMS K₀, TOF t_flight, Orbitrap ω_z
│   ├── thermalization_theory.csv       # Expected T_final
│   ├── transport_theory.csv            # Expected v_drift
│   ├── spacecharge_theory.csv          # Expected expansion
│   ├── reactions_theory.csv            # Expected N(t)
│   └── performance_scaling.csv         # Expected CPU time scaling
│
├── results/
│   └── v1.0_YYYYMMDD_HHMMSS/          # Timestamped runs
│       ├── instruments/
│       ├── thermalization/
│       ├── transport/
│       ├── spacecharge/
│       ├── reactions/
│       ├── performance/
│       └── summary.json
│
└── templates/
    └── VALIDATION_REPORT_TEMPLATE.md
```

---

## 🎯 SESSION 1: INSTRUMENT PHYSICS (5 configs)

**Duration:** 15 minutes  
**Goal:** Validate analytical ion optics models (no collisions/space charge)  
**Status:** ⏳ TODO

### **Test 1.1: IMS Drift Velocity** ⏳

**File:** `validation/configs/instruments/ims_mobility.json`

**Physics:** Mason-Schamp mobility equation: `v_drift = K₀ * (E/N) * (P/P₀) * (T₀/T)`

**Setup:**
- N=100 H3O+ ions, thermal 300K
- Uniform E = 200 V/cm over 10cm
- N2 gas, 1 atm, 300K, HSS collisions
- Duration: 10 ms
- Expected: v_drift ≈ 220 m/s (K₀≈2.7 cm²/(V·s))

**Test:**
```bash
./build/icarion_cli --config validation/configs/instruments/ims_mobility.json \
  --output validation/results/v1.0_test/instruments/ims_mobility.h5
```

**Success Criteria:**
- ✅ Terminal velocity within 1 ms
- ✅ v_drift within 5% of theory
- ✅ Maxwellian velocity distribution

---

### **Test 1.2: TOF Flight Time** ⏳

**File:** `validation/configs/instruments/tof_flight_time.json`

**Physics:** Ballistic flight time: `t = L * sqrt(m / (2*q*U))`

**Setup:**
- N=100 ions (m/z = 100, 200, 300)
- Extraction: 1000V over 5mm
- Free flight: 1m (no fields/collisions)
- Duration: 100 µs
- Expected: t₁₀₀≈32µs, t₂₀₀≈45µs, t₃₀₀≈55µs

**Test:**
```bash
./build/icarion_cli --config validation/configs/instruments/tof_flight_time.json \
  --output validation/results/v1.0_test/instruments/tof_flight_time.h5
```

**Success Criteria:**
- ✅ Flight times within 1% of theory
- ✅ Mass separation Δt/t > 0.05
- ✅ No scattering

---

### **Test 1.3: Orbitrap Axial Frequency** ⏳

**File:** `validation/configs/instruments/orbitrap_frequency.json`

**Physics:** Harmonic oscillation: `ω_z = sqrt(q*k/m)`, k=field curvature

**Setup:**
- N=10 ions (m/z = 200)
- k = 2000 V/mm²
- Initial z = 1mm offset
- Duration: 10 ms
- Expected: f_z ≈ 4.5 kHz

**Test:**
```bash
./build/icarion_cli --config validation/configs/instruments/orbitrap_frequency.json \
  --output validation/results/v1.0_test/instruments/orbitrap_frequency.h5
```

**Success Criteria:**
- ✅ Frequency within 0.5% of theory
- ✅ Radial orbit stable
- ✅ Frequency amplitude-independent

---

### **Test 1.4: LQIT Stability Boundary** ⏳

**File:** `validation/configs/instruments/lqit_stability.json`

**Physics:** Mathieu stability: q=0.908 marginally stable

**Setup:**
- N=10 ions (m/z = 200)
- r₀=5mm, Ω=1MHz, q=0.908
- Initial r = 0.1mm
- Duration: 100 RF cycles (100µs)
- Expected: Bounded oscillation

**Test:**
```bash
./build/icarion_cli --config validation/configs/instruments/lqit_stability.json \
  --output validation/results/v1.0_test/instruments/lqit_stability.h5
```

**Success Criteria:**
- ✅ Ion confined r < r₀ for 100 cycles
- ✅ Amplitude growth matches theory
- ✅ Loss if q > 0.908

---

### **Test 1.5: Energy Conservation** ⏳

**File:** `validation/configs/instruments/energy_conservation.json`

**Physics:** Symplectic integrator: ΔE/E < 10⁻¹²

**Setup:**
- N=100 ions, uniform E=100V/cm
- No collisions, no space charge
- Duration: 1 ms
- Expected: Energy drift < 10⁻¹⁰ (Boris)

**Test:**
```bash
./build/icarion_cli --config validation/configs/instruments/energy_conservation.json \
  --output validation/results/v1.0_test/instruments/energy_conservation.h5
```

**Success Criteria:**
- ✅ ΔE/E < 10⁻¹⁰ (Boris) or 10⁻⁸ (RK4)
- ✅ Linear momentum conserved
- ✅ Reproducible trajectory

---

**Session 1 Deliverables:**
- 5 instrument configs
- 5 HDF5 outputs
- `extract_instruments.py` analysis script
- Pass/fail table

---

## 🎯 SESSION 2: THERMALIZATION (6 configs)

**Duration:** 30 minutes  
**Goal:** Validate HSS/EHSS collision models thermalize correctly  
**Status:** ⏳ TODO

### **Test 1.1: HSS @ 300K, 1 atm (Baseline)** ⏳

**File:** `validation/configs/physics/thermalization/hss_300K_1atm.json`

**Physics:**
- N=1000 H3O+ ions, T_initial=1000 K (hot)
- Buffer gas: N2, T=300 K, P=101325 Pa
- Collision model: HSS
- Runtime: 50 µs (≈10 mean collision times)
- **Expected:** T_final → 300 K (±5%)

**Config Summary:**
```json
{
  "title": "HSS Thermalization: 300K, 1 atm (Baseline)",
  "simulation": {"total_time_s": 5e-5, "dt_s": 1e-9, "rng_seed": 42},
  "physics": {"collision_model": "HSS"},
  "ions": {
    "species": [{
      "id": "H3O+", "count": 1000,
      "position": {"type": "gaussian", "center": [0,0,0.05], "std": [0.001,0.001,0.001]},
      "velocity": {"type": "thermal", "temperature_K": 1000.0}
    }]
  },
  "domains": [{
    "instrument": "IMS",
    "geometry": {"length_m": 0.1, "radius_m": 0.02},
    "environment": {"pressure_Pa": 101325, "temperature_K": 300, "gas_species": "N2"},
    "fields": {"DC": {"axial_V": 0.0}}
  }]
}
```

**Test Command:**
```bash
cd validation
bash scripts/test_single_config.sh configs/physics/thermalization/hss_300K_1atm.json
```

**Success Criteria:**
- ✅ Simulation completes without crash
- ✅ HDF5 file created (size > 1 MB)
- ✅ Ion count = 1000 at all timesteps
- ✅ Final temperature ≈ 300 K (visual check of velocities)

**Common Issues:**
- Species database path wrong → Fix: Use absolute path or `../../../data/species_database_v1.json`
- H3O+ not in database → Add to species_database_v1.json
- Output folder permissions → Create `validation/results/` directory

---

### **Test 1.2: HSS @ 450K, 1 atm** ⏳

**File:** `validation/configs/physics/thermalization/hss_450K_1atm.json`

**Changes from 1.1:**
- `temperature_K`: 300 → 450

**Expected:** T_final → 450 K (±5%)

**Test:** Copy baseline, modify, run

---

### **Test 1.3: HSS @ 300K, Low P** ⏳

**File:** `validation/configs/physics/thermalization/hss_300K_lowP.json`

**Changes from 1.1:**
- `pressure_Pa`: 101325 → 1000 (1% of STP)
- `total_time_s`: 5e-5 → 5e-4 (10x longer, fewer collisions)

**Expected:** T_final → 300 K (±10%, slower convergence)

**Test:** Copy baseline, modify, run

---

### **Test 1.4: EHSS @ 300K, 1 atm** ⏳

**File:** `validation/configs/physics/thermalization/ehss_300K_1atm.json`

**Changes from 1.1:**
- `collision_model`: "HSS" → "EHSS"
- Add: `"molecule_database": "../../../data/molecules"`

**Expected:** T_final → 300 K (±5%, similar to HSS)

**Test:** Copy HSS baseline, change collision model, run

---

### **Test 1.5: EHSS @ 450K, 1 atm** ⏳

**File:** `validation/configs/physics/thermalization/ehss_450K_1atm.json`

**Changes from 1.4:**
- `temperature_K`: 300 → 450

**Expected:** T_final → 450 K (±5%)

---

### **Test 1.6: EHSS @ 300K, Low P** ⏳

**File:** `validation/configs/physics/thermalization/ehss_300K_lowP.json`

**Changes from 1.4:**
- `pressure_Pa`: 101325 → 1000
- `total_time_s`: 5e-5 → 5e-4

**Expected:** T_final → 300 K (±10%)

---

### **Session 1 Actions:**

1. ✅ Create `validation/` directory structure
2. ✅ Create `scripts/test_single_config.sh` (quick test helper)
3. ⏳ Create 6 thermalization configs
4. ⏳ Test each config individually
5. ⏳ Debug any config errors (species, paths, fields)
6. ⏳ Document results in `validation/README.md`

**Deliverables:**
- 6 thermalization JSON configs
- All 6 configs tested and working
- `validation/README.md` updated with Session 1 results

---

## 🎯 SESSION 2: TRANSPORT + SPACE CHARGE (4 configs)

**Duration:** 2 hours  
**Goal:** Validate gas flow transport + Coulomb expansion  
**Status:** ⏳ TODO (after Session 1)

### **Test 2.1: Gas Flow Only (No E-field)** ⏳

**File:** `validation/configs/physics/transport/gas_flow_only.json`

**Physics:**
- N=500 H3O+ ions, thermal equilibrium (T=300K)
- Buffer gas: N2, v_gas = [0, 0, 100 m/s] (axial flow)
- No E-field (DC.axial_V = 0)
- Collision model: HSS
- Runtime: 100 µs
- **Expected:** v_drift → 100 m/s (±10%)

**Key Config:**
```json
{
  "environment": {
    "gas_velocity_m_s": [0.0, 0.0, 100.0]
  },
  "fields": {"DC": {"axial_V": 0.0}}
}
```

**Success Criteria:**
- ✅ Ions reach z ≈ 0.01 m (100 m/s × 100 µs)
- ✅ Mean velocity ≈ 100 m/s

---

### **Test 2.2: Gas Flow + E-field** ⏳

**File:** `validation/configs/physics/transport/gas_flow_drift.json`

**Changes from 2.1:**
- Add: `"axial_V": 100.0` (100V over 0.15m = 667 V/m)

**Expected:** v_total = v_gas + K₀×E ≈ 100 + 187 = 287 m/s

---

### **Test 2.3: 1D Space Charge Expansion** ⏳

**File:** `validation/configs/physics/spacecharge/1d_expansion.json`

**Physics:**
- N=500 H3O+ ions, tight Gaussian (σ=0.5 mm)
- Vacuum (no collisions)
- No E-field (Coulomb only)
- Runtime: 10 µs
- **Expected:** σ_z(t=10µs) > σ_z(t=0)

**Key Config:**
```json
{
  "physics": {
    "collision_model": "NoCollisions",
    "enable_space_charge": true
  },
  "ions": {
    "species": [{
      "count": 500,
      "position": {"type": "gaussian", "std": [0.0005, 0.0005, 0.0005]},
      "velocity": {"type": "fixed", "value": [0, 0, 0]}
    }]
  }
}
```

**Success Criteria:**
- ✅ Cloud expands (σ increases)
- ✅ Expansion rate ~ 1/√t (Coulomb scaling)

---

### **Test 2.4: Radial Space Charge Defocusing** ⏳

**File:** `validation/configs/physics/spacecharge/radial_defocusing.json`

**Changes from 2.3:**
- Initial cloud: cylindrical (σ_r=0.1mm, σ_z=5mm)
- Initial velocity: v_z=1000 m/s (moving beam)

**Expected:** σ_r(t) increases due to Coulomb repulsion

---

### **Session 2 Actions:**

1. ⏳ Create 2 transport configs
2. ⏳ Create 2 space charge configs
3. ⏳ Test each individually
4. ⏳ Verify HDF5 output (velocity/position trends)
5. ⏳ Document in `validation/README.md`

**Potential Issues:**
- Gas flow may not work if not implemented in force registry
- Space charge slow for N=500 (may need smaller N or GPU)
- Need very small dt (1e-10) for space charge accuracy

---

## 🎯 SESSION 3: REACTIONS (1 config)

**Duration:** 1 hour  
**Goal:** Validate simple reaction kinetics  
**Status:** ⏳ TODO (after Session 2)

### **Test 3.1: Simple A+ → B+ Decay** ⏳

**File:** `validation/configs/physics/reactions/simple_reaction.json`

**Physics:**
- N=1000 A+ ions
- Reaction: A+ → B+ with k=1e6 s⁻¹
- No collisions, no E-field (pure kinetics)
- Runtime: 10 µs
- **Expected:** N_A(t) = N₀ × exp(-kt)

**Key Config:**
```json
{
  "physics": {
    "collision_model": "NoCollisions",
    "enable_reactions": true
  },
  "reaction_database": "./reaction_validation.json"
}
```

**Auxiliary File:** `validation/configs/physics/reactions/reaction_validation.json`

```json
{
  "reactions": [
    {
      "id": "decay_A_to_B",
      "reactant": "A+",
      "product": "B+",
      "rate_model": "FirstOrder",
      "k_first_order_s": 1e6
    }
  ]
}
```

**Prerequisites:**
- Add A+ and B+ to `data/species_database_v1.json`
  - Mass: 19 amu (same as H3O+)
  - CCS: 25 Ų (dummy value)

**Success Criteria:**
- ✅ N_A decays exponentially
- ✅ N_A + N_B = N₀ (conservation)
- ✅ Half-life ≈ ln(2)/k = 0.693 µs

---

### **Session 3 Actions:**

1. ⏳ Create reaction config
2. ⏳ Create reaction database file
3. ⏳ Add A+, B+ to species database
4. ⏳ Test config
5. ⏳ Verify exponential decay in HDF5
6. ⏳ Document in `validation/README.md`

---

## 🎯 SESSION 4: PERFORMANCE BENCHMARKS (4 tests)

**Duration:** 45 minutes  
**Goal:** Validate CPU scaling, integrator performance, collision overhead  
**Status:** ⏳ TODO

### **Test 4.1: CPU Scaling with N** ⏳

**Files:**
- `validation/configs/performance/scaling_N_100.json`
- `validation/configs/performance/scaling_N_1000.json`
- `validation/configs/performance/scaling_N_10000.json`
- `validation/configs/performance/scaling_N_100000.json`

**Physics:** Measure CPU time vs N ions (no collisions, no space charge)

**Setup:**
- Fixed timestep dt=1e-8s, duration=1ms
- Uniform E-field (100 V/cm)
- Integrator: RK4
- Expected: CPU time ∝ N (linear scaling)

**Test:**
```bash
for N in 100 1000 10000 100000; do
  ./build/icarion_cli --config validation/configs/performance/scaling_N_${N}.json \
    --output validation/results/v1.0_test/performance/scaling_N_${N}.h5
done
```

**Success Criteria:**
- ✅ Linear scaling: t_CPU(10N) ≈ 10 × t_CPU(N)
- ✅ No memory leak: RAM usage stable
- ✅ Throughput > 1M ion-steps/sec on single core

**Expected Timing:**
- N=100: <1 sec
- N=1k: ~5 sec
- N=10k: ~1 min
- N=100k: ~10 min

---

### **Test 4.2: Integrator Comparison** ⏳

**Physics:** Compare accuracy vs speed for E+B field

**Setup:**
- N=1000 ions in uniform E + B field
- Timesteps: dt=1e-8, 1e-9, 1e-10 s
- Duration: 100 µs
- Measure: Energy drift, CPU time

**Success Criteria:**
- ✅ Boris fastest for E+B (symplectic)
- ✅ Energy conservation: Boris < 10⁻¹⁰, RK4 < 10⁻⁸

---

### **Test 4.3: Collision Overhead** ⏳

**Physics:** Measure collision model overhead

**Setup:**
- N=10k ions, 1 atm N2, 300K
- Duration: 1 ms
- Models: NoCollisions, HSS, EHSS, OU

**Success Criteria:**
- ✅ HSS overhead < 10%
- ✅ EHSS overhead: 50-100%
- ✅ OU overhead: <5%

---

### **Test 4.4: Space Charge Scaling** ⏳

**Physics:** Compare direct N² vs grid N log N

**Setup:**
- Direct method: N=100, 1000
- Grid method: N=10k
- Duration: 100 µs

**Success Criteria:**
- ✅ Direct: t_CPU(1000) ≈ 100 × t_CPU(100) (N² scaling)
- ✅ Grid: sub-quadratic scaling
- ✅ Grid accuracy: <5% error vs direct

---

**Session 4 Deliverables:**
- 11 performance configs total
- Analysis script: `extract_performance.py`
- Performance table in report

---

## 🎯 SESSION 5: SCRIPTS + ORCHESTRATION (30 minutes)

**Duration:** 2 hours  
**Goal:** Automate all tests + generate report  
**Status:** ⏳ TODO (after Session 3)

### **Script 4.1: Quick Test Helper** ✅

**File:** `validation/scripts/test_single_config.sh`

```bash
#!/bin/bash
# Quick test for a single config

CONFIG=$1
OUTPUT_DIR="results/test_$(date +%Y%m%d_%H%M%S)"

echo "Testing: $CONFIG"
mkdir -p "$OUTPUT_DIR"

../../build/icarion_core "$CONFIG" \
    --output-folder "$OUTPUT_DIR" \
    > "$OUTPUT_DIR/simulation.log" 2>&1

if [ $? -eq 0 ]; then
    echo "✅ SUCCESS"
    ls -lh "$OUTPUT_DIR/trajectories.h5"
else
    echo "❌ FAILED - Check log: $OUTPUT_DIR/simulation.log"
fi
```

---

### **Script 4.2: Run All Thermalization** ⏳

**File:** `validation/scripts/run_thermalization_tests.sh`

```bash
#!/bin/bash
OUTPUT_DIR=${1:-"results/validation_$(date +%Y%m%d_%H%M%S)"}
mkdir -p "$OUTPUT_DIR/thermalization"

CONFIGS=(
    "hss_300K_1atm" "hss_450K_1atm" "hss_300K_lowP"
    "ehss_300K_1atm" "ehss_450K_1atm" "ehss_300K_lowP"
)

for cfg in "${CONFIGS[@]}"; do
    echo "Running: $cfg..."
    mkdir -p "$OUTPUT_DIR/thermalization/$cfg"
    
    ../../build/icarion_core \
        "configs/physics/thermalization/${cfg}.json" \
        --output-folder "$OUTPUT_DIR/thermalization/$cfg" \
        > "$OUTPUT_DIR/thermalization/$cfg/simulation.log" 2>&1
    
    echo "  $([ $? -eq 0 ] && echo '✅ PASS' || echo '❌ FAIL')"
done
```

---

### **Script 4.3: Extract Thermalization Metrics** ⏳

**File:** `validation/scripts/utils/extract_thermalization.py`

```python
#!/usr/bin/env python3
import h5py, numpy as np, json, sys
from pathlib import Path

def compute_temperature(velocities, mass_kg=19.0*1.66054e-27):
    """T from velocity distribution via equipartition"""
    v2 = np.sum(velocities**2, axis=1)
    KE_avg = 0.5 * mass_kg * np.mean(v2)
    return (2.0/3.0) * KE_avg / 1.380649e-23

def extract(h5_file, T_expected):
    with h5py.File(h5_file, 'r') as f:
        species = list(f['trajectories'].keys())[0]
        vel = f[f'trajectories/{species}/velocity'][-1, :, :]
        T_final = compute_temperature(vel)
        error = abs(T_final - T_expected) / T_expected
        
        return {
            "T_final_K": float(T_final),
            "T_expected_K": float(T_expected),
            "relative_error": float(error),
            "status": "PASS" if error < 0.10 else "FAIL"
        }

if __name__ == "__main__":
    test_dir = Path(sys.argv[1])
    T_expected = float(sys.argv[2])
    
    metrics = extract(test_dir / "trajectories.h5", T_expected)
    
    with open(test_dir / "metrics.json", 'w') as f:
        json.dump(metrics, f, indent=2)
    
    print(json.dumps(metrics))
```

---

### **Script 4.4: Generate Report** ⏳

**File:** `validation/scripts/utils/generate_report.py`

```python
#!/usr/bin/env python3
import json
from pathlib import Path
import sys
from datetime import datetime

def load_metrics(results_dir):
    metrics = {}
    for mf in Path(results_dir).rglob("metrics.json"):
        path = str(mf.parent.relative_to(results_dir))
        with open(mf) as f:
            metrics[path] = json.load(f)
    return metrics

def generate_report(results_dir):
    metrics = load_metrics(results_dir)
    n_pass = sum(1 for m in metrics.values() if m.get("status") == "PASS")
    n_total = len(metrics)
    
    report = f"""# ICARION v1.0 Physics Validation Report

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Summary

- **Total Tests:** {n_total}
- **Passed:** {n_pass} ({100*n_pass//n_total if n_total else 0}%)
- **Failed:** {n_total - n_pass}

## Results

| Test | Status | Metric |
|------|--------|--------|
"""
    
    for test, data in sorted(metrics.items()):
        emoji = "✅" if data["status"] == "PASS" else "❌"
        metric = f"{data.get('T_final_K', 'N/A'):.1f} K" if 'T_final_K' in data else "N/A"
        report += f"| `{test}` | {emoji} | {metric} |\n"
    
    return report

if __name__ == "__main__":
    print(generate_report(sys.argv[1]))
```

---

### **Script 4.5: Master Suite** ⏳

**File:** `validation/run_validation_suite.sh`

```bash
#!/bin/bash
set -e

OUTPUT_DIR="results/validation_$(date +%Y%m%d_%H%M%S)"

echo "ICARION Physics Validation Suite"
echo "Output: $OUTPUT_DIR"

# Build
echo "==> Building..."
cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j$(nproc) && cd ../validation

# Run tests
echo "==> Running tests..."
bash scripts/run_thermalization_tests.sh "$OUTPUT_DIR"
bash scripts/run_transport_tests.sh "$OUTPUT_DIR"
bash scripts/run_spacecharge_tests.sh "$OUTPUT_DIR"
bash scripts/run_reactions_tests.sh "$OUTPUT_DIR"

# Extract metrics
echo "==> Extracting metrics..."
for test in "$OUTPUT_DIR"/thermalization/*/; do
    name=$(basename "$test")
    T=$(echo "$name" | grep -oP '\d+K' | grep -oP '\d+')
    python3 scripts/utils/extract_thermalization.py "$test" "$T"
done

# Generate report
echo "==> Generating report..."
python3 scripts/utils/generate_report.py "$OUTPUT_DIR" > VALIDATION_REPORT_v1.0.md

echo "✅ Done! Report: VALIDATION_REPORT_v1.0.md"
```

---

### **Session 4 Actions:**

1. ⏳ Create all bash scripts
2. ⏳ Create Python utilities
3. ⏳ Test full validation suite end-to-end
4. ⏳ Fix any bugs in orchestration
5. ⏳ Generate first VALIDATION_REPORT_v1.0.md
6. ⏳ Document in `validation/README.md`

---

## 📋 PROGRESS TRACKING

### **Completed:**
- ✅ Plan created (`tmp/VALIDATION_BENCH_PLAN.md`)

### **TODO:**
- ⏳ Session 1: Instrument Physics (5 configs, ~15 min)
- ⏳ Session 2: Thermalization (6 configs, ~30 min)
- ⏳ Session 3: Transport + Space Charge + Reactions (5 configs, ~17 min)
- ⏳ Session 4: Performance Benchmarks (11 configs, ~45 min)
- ⏳ Session 5: Scripts + Orchestration (~30 min)

### **Expected Issues & Solutions:**

| Issue | Solution |
|-------|----------|
| Species database path | Use absolute: `/home/chsch95/ICARION/data/species_database_v1.json` |
| H3O+ not in database | Add with mass=19, CCS=104 Ų (N2) |
| Output folder permissions | `mkdir -p validation/results` |
| Gas flow not implemented | Check `DomainConfig::gas_velocity_m_s` integration |
| Space charge too slow | Reduce N to 200 ions or add timeout |
| Reaction format wrong | Check schema: `k_first_order_s` vs `rate_constant_cm3s` |
| EHSS needs molecule files | Copy H3O.xyz to `data/molecules/` |

---

## 🎯 NEXT STEPS

**Ready to start Session 1?**

```bash
# Create directory structure
cd /home/chsch95/ICARION
mkdir -p validation/{configs/{instruments,physics/{thermalization,transport,spacecharge,reactions},performance},scripts/utils,results,reference_data,templates}

# Create test script
cat > validation/scripts/test_single_config.sh << 'EOF'
#!/bin/bash
CONFIG=$1
OUTPUT_DIR="results/test_$(date +%Y%m%d_%H%M%S)"
echo "Testing: $CONFIG"
mkdir -p "$OUTPUT_DIR"
../../build/icarion_cli "$CONFIG" --output "$OUTPUT_DIR/trajectories.h5" > "$OUTPUT_DIR/simulation.log" 2>&1
[ $? -eq 0 ] && echo "✅ SUCCESS" || echo "❌ FAILED"
ls -lh "$OUTPUT_DIR/"
EOF
chmod +x validation/scripts/test_single_config.sh

# Start with first instrument config (IMS mobility)
nano validation/configs/instruments/ims_mobility.json
```

**Let's go! 🚀**

