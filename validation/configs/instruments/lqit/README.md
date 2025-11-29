# LQIT Validation Test Suite

## Overview
Tests Linear Quadrupole Ion Trap (LQIT) stability, secular motion, and resonant excitation.

## Test Parameters
- **Ion Species:** H3O+ (19.02 u, +1 e)
- **Ion Count:** 1000
- **Trap Geometry:** r₀ = 5.0 mm, L = 50.0 mm
- **RF Frequency:** 1.0 MHz
- **Environment:** He @ 0.1 Pa (UHV)

## Test Suite

### Part 1: Mathieu Stability (6 tests)
1. **Stable Region (q=0.4):** All collision models
2. **Stability Boundary (q=0.7):** Near-critical stability
3. **Unstable Region (q=0.95):** Expected ion ejection
4. **DC Offset:** Stability with axial confinement

### Part 2: AC Resonant Excitation (4 tests)
1. **On-Resonance (f_AC = f_sec):** Strong heating/ejection
2. **Off-Resonance (f_AC = 0.5×f_sec):** Minimal effect
3. **Off-Resonance (f_AC = 2.0×f_sec):** Minimal effect
4. **Collision Model Comparison:** HSS vs Friction damping

## Expected Results

### Stability Tests
- **Stable (q=0.4):** All ions confined, oscillate at f_sec ≈ 141 kHz
- **Stable (q=0.7):** Ions confined, f_sec ≈ 247 kHz  
- **Unstable (q=0.95):** Ions ejected within 5-10 RF periods

### AC Excitation Tests
- **On-Resonance:** Ion amplitude grows exponentially → ejection
- **Off-Resonance:** Ion amplitude stable (no parametric resonance)
- **Collision Damping:** HSS shows stronger damping than Friction in UHV

## Analysis
```bash
# Stability analysis
python3 ../scripts/analyze_lqit_stability.py results/v1.0_test/instruments/lqit/lqit_*_q*.h5

# AC excitation analysis
python3 ../scripts/analyze_lqit_excitation.py results/v1.0_test/instruments/lqit/lqit_*_ac*.h5
```

## Physics Background

### Mathieu Stability
For 2D quadrupole (LQIT):
```
q = 4·Q·V_rf / (m·Ω²·r₀²)
a = -8·Q·V_dc / (m·Ω²·r₀²)
```
Stability: q < 0.908, |a| < 0.2

### Secular Frequency
Pseudopotential approximation (q << 1):
```
ω_sec ≈ q·Ω / (2√2)
```

### Resonant Excitation
AC field at f_sec causes parametric resonance → amplitude growth → ejection

## References
- March & Todd, "Quadrupole Ion Trap Mass Spectrometry" (2005)
- Dehmelt pseudopotential approximation
- Mathieu equation solutions
