# ICARION Validation Bench - Progress Tracker

**Last Updated:** 2025-11-26  
**Plan:** See `../tmp/VALIDATION_BENCH_PLAN.md`  
**Goal:** Publication-quality validation suite (20 tests)

---

## 📊 Overall Progress

| Session | Category | Tests | Time | Status | Notes |
|---------|----------|-------|------|--------|-------|
| **1** | Instrument Physics | 5 | ~15 min | ⏳ TODO | IMS, TOF, Orbitrap, LQIT, Energy |
| **2** | Thermalization | 6 | ~30 min | ⏳ TODO | HSS/EHSS at 3 conditions each |
| **3** | Transport/SpaceCharge/Reactions | 5 | ~17 min | ⏳ TODO | Gas flow + Coulomb + kinetics |
| **4** | Performance Benchmarks | 11 | ~45 min | ⏳ TODO | Scaling + integrators + overhead |
| **5** | Scripts | - | ~30 min | ⏳ TODO | Orchestration + report generator |

**Total:** 20 tests (5 instrument + 11 physics + 4 performance)  
**Runtime:** ~2 hours (full suite)

---

## ✅ SESSION 1: INSTRUMENT PHYSICS (Status: ⏳ TODO)

### **Tests:**

1. ✅ `ims_mobility.json` - IMS drift velocity (Mason-Schamp) - **PASS** (1m50s)
2. ⏳ `tof_flight_time.json` - TOF ballistic flight
3. ⏳ `orbitrap_frequency.json` - Orbitrap axial oscillation
4. ⏳ `lqit_stability.json` - LQIT Mathieu stability
5. ⏳ `energy_conservation.json` - Energy conservation test

### **Expected Results:**

| Test | Measurement | Expected | Tolerance | Status |
|------|-------------|----------|-----------|--------|
| IMS mobility | v_drift | ~220 m/s | ±5% | ⏳ |
| TOF flight | t(m/z=200) | ~45 µs | ±1% | ⏳ |
| Orbitrap | f_z | ~4.5 kHz | ±0.5% | ⏳ |
| LQIT stability | confinement | bounded | 100 cycles | ⏳ |
| Energy conservation | ΔE/E | <10⁻¹⁰ | Boris | ⏳ |

### **Test Results:**

```
⏳ Not yet run
```

### **Issues Encountered:**

```
None yet
```

---

## ✅ SESSION 2: THERMALIZATION (Status: ⏳ TODO)

### **Tests:**

1. ⏳ `hss_300K_1atm.json` - HSS @ 300K, 1 atm (baseline)
2. ⏳ `hss_450K_1atm.json` - HSS @ 450K, 1 atm (high T)
3. ⏳ `hss_300K_lowP.json` - HSS @ 300K, 1000 Pa (low P)
4. ⏳ `ehss_300K_1atm.json` - EHSS @ 300K, 1 atm (baseline)
5. ⏳ `ehss_450K_1atm.json` - EHSS @ 450K, 1 atm (high T)
6. ⏳ `ehss_300K_lowP.json` - EHSS @ 300K, 1000 Pa (low P)

### **Expected Results:**

| Test | T_initial | T_buffer | T_final (expected) | Tolerance |
|------|-----------|----------|-------------------|-----------|
| hss_300K_1atm | 1000 K | 300 K | 300 K | ±5% |
| hss_450K_1atm | 1000 K | 450 K | 450 K | ±5% |
| hss_300K_lowP | 1000 K | 300 K | 300 K | ±10% |
| ehss_300K_1atm | 1000 K | 300 K | 300 K | ±5% |
| ehss_450K_1atm | 1000 K | 450 K | 450 K | ±5% |
| ehss_300K_lowP | 1000 K | 300 K | 300 K | ±10% |

### **Test Results:**

```
⏳ Not yet run
```

### **Issues Encountered:**

```
None yet
```

---

## ✅ SESSION 3: TRANSPORT / SPACE CHARGE / REACTIONS (Status: ⏳ TODO)

### **Tests:**

1. ⏳ `gas_flow_only.json` - Pure gas drag (no E-field)
2. ⏳ `gas_flow_drift.json` - Gas + E-field combined
3. ⏳ `1d_expansion.json` - 1D Coulomb expansion
4. ⏳ `radial_defocusing.json` - Radial space charge spread
5. ⏳ `simple_reaction.json` - A+ → B+ first-order decay

### **Expected Results:**

| Test | Measurement | Expected | Status |
|------|-------------|----------|--------|
| Gas flow only | <v> | ~flow velocity | ⏳ |
| Gas+E drift | v_drift | mobility × E | ⏳ |
| 1D expansion | σ(t) | Coulomb growth | ⏳ |
| Radial defocus | σ_r(t) | Radial spread | ⏳ |
| Reaction | N_A(t) | exp(-kt) | ⏳ |

### **Test Results:**

```
⏳ Not yet run
```

---

## ✅ SESSION 4: PERFORMANCE BENCHMARKS (Status: 🔄 IN PROGRESS)

### **Tests:**

1. ✅ **OpenMP scaling test** - N=100 ions - **COMPLETE**
   - Finding: N=100 too small for OpenMP benefit (overhead dominates)
   - 1 thread: 17s, 2-32 threads: ~100s (5.8x slower!)
   - Need N≥1000 for parallel speedup
2. ⏳ CPU scaling (N=100, 1k, 10k, 100k) - 4 configs
3. ⏳ Integrator comparison (RK4, RK45, Boris) - 3 configs
4. ⏳ Collision overhead (none, HSS, EHSS, OU) - 4 configs
5. ⏳ Space charge scaling (direct N², grid N log N) - 3 configs

**Total:** 11 performance configs

### **Expected Results:**

| Test | Expected | Status |
|------|----------|--------|
| CPU scaling | Linear O(N) | ⏳ |
| Boris integrator | Fastest for E+B | ⏳ |
| HSS overhead | <10% | ⏳ |
| EHSS overhead | 50-100% | ⏳ |
| Grid space charge | Sub-quadratic | ⏳ |

### **Test Results:**

```
⏳ Not yet run
```

---

## ✅ SESSION 5: ORCHESTRATION (Status: ⏳ TODO)

### **Scripts:**

- ⏳ `scripts/test_single_config.sh` - Quick test helper
- ⏳ `scripts/run_instrument_tests.sh` - Session 1: Instruments
- ⏳ `scripts/run_thermalization_tests.sh` - Session 2: Thermalization
- ⏳ `scripts/run_transport_tests.sh` - Session 3: Transport
- ⏳ `scripts/run_spacecharge_tests.sh` - Session 3: Space charge
- ⏳ `scripts/run_reactions_tests.sh` - Session 3: Reactions
- ⏳ `scripts/run_performance_tests.sh` - Session 4: Performance
- ⏳ `scripts/run_all_validation.sh` - Master script (all sessions)
- ⏳ `scripts/utils/extract_instruments.py` - Extract IMS/TOF/Orbitrap metrics
- ⏳ `scripts/utils/extract_thermalization.py` - Extract T_final
- ⏳ `scripts/utils/extract_transport.py` - Extract drift velocity
- ⏳ `scripts/utils/extract_spacecharge.py` - Extract cloud size
- ⏳ `scripts/utils/extract_reactions.py` - Extract species counts
- ⏳ `scripts/utils/extract_performance.py` - Extract CPU time, scaling
- ⏳ `scripts/utils/generate_report.py` - Generate Markdown report

### **Deliverables:**

- ⏳ `VALIDATION_REPORT_v1.0.md` (auto-generated)
- ⏳ `results/validation_YYYYMMDD_HHMMSS/` (timestamped runs)

---

## 📝 NEXT STEPS

**Current Task:** Session 1 - Create first thermalization config

**Commands:**

```bash
# 1. Create directory structure
cd /home/chsch95/ICARION
mkdir -p validation/{configs/physics/{thermalization,transport,spacecharge,reactions},scripts/utils,results,reference_data,templates}

# 2. Create test script
cat > validation/scripts/test_single_config.sh << 'EOF'
#!/bin/bash
CONFIG=$1
OUTPUT_DIR="results/test_$(date +%Y%m%d_%H%M%S)"
echo "Testing: $CONFIG"
mkdir -p "$OUTPUT_DIR"
../../build/icarion_core "$CONFIG" --output-folder "$OUTPUT_DIR" > "$OUTPUT_DIR/simulation.log" 2>&1
[ $? -eq 0 ] && echo "✅ SUCCESS" || echo "❌ FAILED"
ls -lh "$OUTPUT_DIR/"
EOF
chmod +x validation/scripts/test_single_config.sh

# 3. Create first config
cd validation
# Now ready to create configs/physics/thermalization/hss_300K_1atm.json
```

---

## 📚 DOCUMENTATION

After each session, update this file with:
- ✅/❌ status for each test
- Actual vs. expected results
- Issues encountered + solutions
- Runtime measurements

**Format:**

```markdown
### **Test Results:**

| Test | Status | T_final | Expected | Error | Runtime |
|------|--------|---------|----------|-------|---------|
| hss_300K_1atm | ✅ PASS | 298.5 K | 300 K | 0.5% | 45s |
| ... | ... | ... | ... | ... | ... |
```

---

**Ready to start! 🚀**

