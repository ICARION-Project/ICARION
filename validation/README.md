# ICARION Validation Suite

**Version:** 1.0  
**Last Updated:** 2025-11-28  
**Branch:** `feature/validation-suite`  
**Goal:** Systematic validation of physics and performance after SoA Foundation v1.0

---

## 📊 Overall Progress

| Session | Category | Configs | Status | Notes |
|---------|----------|---------|--------|-------|
| **1** | **Thermalization** | 90 | ✅ **COMPLETE** | Cold start (0.1K) → target temp, HSS/EHSS validated |
| **2** | Instrument Physics | TBD | 🔄 IN PROGRESS | IMS, Orbitrap, LQIT validation |
| **3** | Transport/SpaceCharge | TBD | ⏳ PLANNED | Drift velocity, Coulomb expansion |
| **4** | Reactions | TBD | ⏳ PLANNED | Kinetics validation |
| **5** | Performance | TBD | ⏳ PLANNED | Scaling, integrators, overhead |

**Completed:** 90 thermalization configs + analysis tools  
**Current Focus:** Instrument validation (IMS, Orbitrap, LQIT)

---

## 🔄 SESSION 2: INSTRUMENT PHYSICS (Status: 🔄 IN PROGRESS)

### **Planned Tests:**

#### **IMS (Ion Mobility Spectrometry)**
- Drift velocity vs E/N (reduced field)
- Mobility verification against Mason-Schamp equation
- Resolution vs drift length and voltage

#### **Orbitrap**
- Axial frequency vs m/z (hyperlogarithmic geometry)
- Mass resolution validation
- Stability of ion trajectories

#### **LQIT (Linear Quadrupole Ion Trap)**
- Ion confinement at various q/a parameters
- RF trapping secular motion frequency
- Ion motion described by Pseudopotential approximation

#### **Quadrupole Mass Filter**
- Stability map validation (a-q diagram)
- Transmission efficiency vs m/z
- RF/DC voltage ratio effects

#### **FT-ICR (Fourier Transform Ion Cyclotron Resonance)**
- Cyclotron frequency vs m/z
- Mass resolution validation
- Energy conservation in collisionless environment

#### **Cross-Instrument**
- Energy conservation (collisionless)
- Boundary detection accuracy
- Field solver validation

### **Expected Results:**

| Instrument | Test | Measurement | Expected | Tolerance |
|------------|------|-------------|----------|-----------|
| IMS | Drift velocity | v_drift | K₀ × E/N | ±5% |
| Orbitrap | Axial frequency | f_z(m/z) | k·√(m/z) | ±1% |
| LQIT | Stability | secular frequency | Pseudopotential approximation | ±3% |
| Quadrupole | Stability-map | a-q parameter | Mathieu region | Stable |
| FT-ICR | Stability | cyclotron frequency | Theoretical model | Stable |
| All | Energy | ΔE/E | <10⁻¹⁰ | Collisionless |

### **Status:**

```
🔄 Directories created, awaiting test config generation
```

---

## SESSION 1: THERMALIZATION (Status: COMPLETE)

### **Design:**

**Systematic parameter sweep:** 3 temperatures × 5 pressures × 3 ion species × 2 collision models = **90 configurations**

- **Temperatures:** 150K, 300K, 1000K
- **Pressures:** 0.2, 2.0, 20.0, 200.0, 2000.0 Pa
- **Ion species:** H3O+, PentanalH+, 2,6-DTBPH+
- **Collision models:** HSS (Hard-Sphere Scattering), EHSS (Elastic Hard-Sphere Scattering)
- **Initial conditions:** 0.1K ion temperature (cold start for clear thermalization signal)
- **Duration:** 20-150 collision times (species-dependent)
- **Timestep:** dt = τ_collision/50 (kinetic theory calculation)

### **Key Results:**

| Model | T_initial | T_final | T_target | Error | CPU Time | Status |
|-------|-----------|---------|----------|-------|----------|--------|
| **HSS**  | 4.1 K | 288.2 K | 300.0 K | 3.9% | 4.3 s | PASS |
| **EHSS** | 1.3 K | 319.6 K | 300.0 K | 6.5% | 34.6 s | PASS |

**Key Findings:**
- Both collision models thermalize correctly from cold initial conditions
- Final temperatures within 10% of target across all conditions
- EHSS ~8x slower than HSS due to molecular structure calculations
- Cold start (0.1K → 300K) provides clear thermalization validation
- Collision time calculations accurate (dt = τ_coll/50 resolution sufficient)

### **Tools:**

- `scripts/generate_thermalization_configs.py` - Generate 90 systematic test configs
- `scripts/run_thermalization_tests.sh` - Test orchestration (quick/subset/full modes)
- `scripts/final_therm_check.py` - HDF5 velocity analysis for thermalization validation
- `scripts/test_single_config.sh` - Single test runner with logging

### **Issues Resolved:**

1. **EHSS GeometryMap loading** - Fixed with absolute path to `/home/chsch95/ICARION/data/molecules/`
2. **ConfigLoader species_id field** - Added backward compatibility for old "id" format
3. **Initial temperature selection** - Changed from T_target to 0.1K for clear thermalization signal

---

## ⏳ SESSION 3: TRANSPORT PHYSICS (Status: ⏳ PLANNED)

### **Planned Tests:**

#### **Drift Velocity**
- Mason-Schamp equation validation
- Mobility vs E/N verification
- Temperature dependence

#### **Gas Flow**
- Pure gas drag (no E-field)
- Combined gas flow + E-field
- Flow profile effects

#### **Diffusion**
- Einstein relation (D = μ·k_B·T/q)
- Cloud spreading vs time
- Temperature effects

### **Expected Results:**

| Test | Measurement | Expected | Notes |
|------|-------------|----------|-------|
| Drift velocity | v_drift | K₀ × E | Mason-Schamp |
| Gas drag | <v> | v_gas | Pure flow |
| Diffusion | σ(t) | √(2Dt) | Einstein |

---

## ⏳ SESSION 4: SPACE CHARGE (Status: ⏳ PLANNED)

### **Planned Tests:**

#### **Coulomb Expansion**
- 1D expansion from point source
- Radial defocusing in drift tube
- Self-bunching effects

#### **Space Charge Algorithms**
- Direct summation (N²) accuracy
- Grid-based methods validation
- Performance comparison

### **Expected Results:**

| Test | Measurement | Expected | Notes |
|------|-------------|----------|-------|
| 1D expansion | σ(t) | Coulomb growth | Linear |
| Radial spread | σ_r(t) | Defocusing rate | IMS |

---

## ⏳ SESSION 5: REACTION KINETICS (Status: ⏳ PLANNED)

### **Planned Tests:**

#### **First-Order Reactions**
- A+ → B+ decay validation
- Rate constant verification
- Temperature dependence

#### **Bimolecular Reactions**
- A+ + B → C+ kinetics
- Collision frequency accuracy

### **Expected Results:**

| Test | Measurement | Expected | Notes |
|------|-------------|----------|-------|
| First-order | N_A(t) | exp(-kt) | Decay |
| Bimolecular | N_C(t) | Rate equation | Collision |

---

## ⏳ SESSION 6: PERFORMANCE BENCHMARKS (Status: ⏳ PLANNED)

### **Planned Tests:**

#### **Scaling**
- Ion count scaling: N = 100, 1k, 10k, 100k
- OpenMP thread scaling (1-32 threads)
- GPU vs CPU performance

#### **Integrators**
- RK4 vs RK45 vs Boris accuracy and speed
- Adaptive timestep efficiency
- Energy conservation

#### **Physics Overhead**
- Collisionless baseline
- HSS collision overhead
- EHSS collision overhead
- Space charge algorithms (direct vs grid)

### **Expected Results:**

| Test | Expected | Notes |
|------|----------|-------|
| Ion scaling | O(N) | Linear |
| OpenMP | Speedup for N>1000 | Overhead for small N |
| Boris | Fastest for E+B | Energy conserving |
| HSS overhead | <10% | Simple model |
| EHSS overhead | ~8x slower | Molecular structure |
| Grid space charge | O(N log N) | vs O(N²) direct |

### **Known Results:**

From thermalization validation:
- **HSS**: 4.3s for 10,000 ions, 12μs simulation
- **EHSS**: 34.6s for 10,000 ions, 12μs simulation (8x slower)
- **Performance ratio**: EHSS/HSS ≈ 8.0

---

## 📦 VALIDATION INFRASTRUCTURE

### **Available Scripts:**

✅ **Thermalization (Session 1):**
- `scripts/generate_thermalization_configs.py` - Generate 90 systematic configs
- `scripts/run_thermalization_tests.sh` - Test orchestration (quick/subset/full)
- `scripts/test_single_config.sh` - Single test runner with logging
- `scripts/final_therm_check.py` - HDF5 velocity analysis

⏳ **Instruments (Session 2):**
- `scripts/run_instrument_tests.sh` - IMS, Orbitrap, LQIT validation
- `scripts/analyze_ims.py` - IMS drift velocity extraction
- `scripts/analyze_orbitrap.py` - Orbitrap frequency analysis

⏳ **Transport/Space Charge/Reactions (Sessions 3-5):**
- `scripts/run_transport_tests.sh` - Drift and gas flow validation
- `scripts/run_spacecharge_tests.sh` - Coulomb expansion tests
- `scripts/run_reactions_tests.sh` - Kinetics validation

⏳ **Performance (Session 6):**
- `scripts/run_performance_tests.sh` - Scaling and overhead benchmarks

⏳ **Orchestration:**
- `scripts/run_all_validation.sh` - Master script for full suite

### **Directory Structure:**

```
validation/
├── configs/
│   ├── physics/
│   │   └── thermalization/     # ✅ 90 configs (HSS/EHSS)
│   ├── instruments/            # 🔄 IMS, Orbitrap, LQIT
│   │   ├── ims/
│   │   ├── orbitrap/
│   │   └── lqit/
│   ├── transport/              # ⏳ Drift, gas flow
│   ├── spacecharge/            # ⏳ Coulomb expansion
│   ├── reactions/              # ⏳ Kinetics
│   └── performance/            # ⏳ Scaling benchmarks
├── scripts/                    # ✅ Test orchestration and analysis
├── results/                    # Test outputs (HDF5, logs)
└── README.md                   # This file
```

---

## 📝 USAGE

### **Quick Start:**

```bash
# Run single thermalization test
cd /home/chsch95/ICARION/validation
../build/src/icarion_main configs/physics/thermalization/hss_H3Op_300K_20.0Pa.json

# Analyze results
source /home/chsch95/dev_venv/bin/activate
python3 scripts/final_therm_check.py

# Run full thermalization suite (90 configs)
./scripts/run_thermalization_tests.sh full
```

### **Test Modes:**

- **Quick mode** (6 configs): Representative subset for rapid validation
- **Subset mode** (18 configs): One ion species, all conditions
- **Full mode** (90 configs): Complete systematic sweep

---

## 🎯 NEXT STEPS

### **Current Focus: Session 2 - Instrument Validation**

**Planned work:**
1. Create IMS validation configs (drift velocity, mobility)
2. Create Orbitrap validation configs (axial frequency, m/z)
3. Create LQIT validation configs (Mathieu stability)
4. Implement analysis scripts for instrument metrics
5. Validate against theoretical predictions

**Commands:**

```bash
cd /home/chsch95/ICARION/validation

# Create IMS test configs
scripts/generate_ims_configs.py

# Run IMS validation
./scripts/run_instrument_tests.sh ims

# Analyze results
python3 scripts/analyze_ims.py results/instruments/ims/
```

---

## 📚 REFERENCES

- **Thermalization validation**: Maxwell-Boltzmann distribution, kinetic theory
- **IMS theory**: Mason-Schamp equation, reduced mobility
- **Orbitrap**: Makarov (2000), hyperlogarithmic geometry
- **LQIT**: Mathieu stability diagram, Paul trap theory
- **Space charge**: Coulomb expansion, Poisson equation
- **Collision models**: Hard-sphere scattering (HSS), Elastic HSS (EHSS)

---

**Last updated:** 2025-11-28  
**Status:** Session 1 (Thermalization) complete ✅, Session 2 (Instruments) in progress 🔄

