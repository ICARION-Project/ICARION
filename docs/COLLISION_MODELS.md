# ICARION Collision Models - Status and Recommendations

## Overview

ICARION implements multiple collision models for ion-neutral interactions. This document provides guidance on which models to use for production simulations vs experimental testing.

## Model Status Summary

| Model | Status | Use Case | Validation |
|-------|--------|----------|------------|
| **NoCollisions** | **PRODUCTION READY** | Vacuum simulations | N/A |
| **Friction** | **PRODUCTION READY** | All simulations | Validated |
| **HSS** | **PRODUCTION READY** | Stochastic collisions | Validated |
| **EHSS** | **PRODUCTION (mobility bias noted)** | Molecular structure resolved scattering | Thermalization validated; mobility overestimates |
| **HSD** | **EXPERIMENTAL** | Deterministic hard sphere | Partially validated |
| **Langevin** | **EXPERIMENTAL** | Polarizable gases only | Not validated for N2 |

---

## Production-Ready Models

### 1. Friction (CollisionModel::Friction) **Recommended for high pressure simulations**

**Formula:** γ = q/(K·m) where K = K₀·(n₀/n)

**Advantages:**
- Based on experimentally measured or literature mobility (K₀)
- Mason-Schamp equation foundation
- **Exact agreement** with IMS drift measurements (356 m/s measured vs 356 m/s theory in N2)
- Works for all gas types (polar and non-polar)
- Density correction included automatically

**Validation:**
- H3O+ in N2 at 1000 Pa: Within ±40% tolerance (test_ims_drift.cpp)
- Consistent with literature mobility values in N2
- Model provides good approximation, but does not account for diffusion (no randomization)

**Usage:**
```cpp
config::CollisionModel::Friction
// Requires: ion.reduced_mobility_cm2_Vs (from experiment/simulation/literature)
```

**When to use:**
- At high pressures (> 100 mbar) and when only mobility is of interest, not the actual ion cloud shape
- IMS, DTIMS
- Any simulation where experimental mobility is known and diffusion is not critical

---

### 2. HSS (Hard-Sphere Stochastic)

**Formula:** Stochastic, isotropic collision -> velocity randomization

**Advantages:**
- Full stochastic treatment (diffusion included)
- Validated thermalization behavior (see [VALIDATION_REPORT_v1.0.md](../validation/VALIDATION_REPORT_v1.0.md))
- Good agreement with mobility (±1%)

**Validation:**
- H3O+ in N2: 360 m/s measured vs 357 m/s theory (0.8% error)
- Diffusion coefficient correct
- Full stochastic treatment provides realistic thermal behavior

**Note**
- Simplified elastic collision model
- Does not account for molecular structure (in contras to EHSS model)
- Does not account for inelastic collisions

**Usage:**
```cpp
config::CollisionModel::HSS
// Requires: ion.CCS_m2 (collision cross-section)
```

**When to use:**
- Need accurate thermal diffusion
- Trajectory-level accuracy required
- Modeling real experimental conditions
- No DFT-optimized structures available 

---

### 3. EHSS (Exact hard-sphere scattering)

**Formula:** Exact scattering angles + stochastic rteatment

**Validation:**
- H3O+ in N2: ~23% higher drift velocity than theory in simple CTest
- Thermalization behavior correct; mobility bias remains
- More accurate angular scattering than HSS when molecular structure is important

**Usage:**
```cpp
config::CollisionModel::EHSS
// Requires: ion.CCS_m2
```

**Note**
- Does not account for long-range interactions
- Does not account for inelastic collisions
- GPU helper supports HSS/EHSS batch collisions; dispatch uses active-ion threshold (default 5000).

**When to use:**
- Studying structure-dependent effects
- When detailed scattering physics is required
- DFT or ab initio simulations available 

---

## Experimental Models

### 4. HSD (CollisionModel::HSD) - Hard-Sphere Deterministic - EXPERIMENTAL

**Formula:** γ = n·σ·v_th (kinetic theory)

**Status:** Under validation, use with caution

**Known Issues:**
- Gives γ = 1.20e8 Hz (16% lower than Friction: 1.43e8 Hz)
- Requires **accurate CCS for target gas** !
- Example: H3O+ CCS = 104 Å in N2 
- CCS must be calculated from Mason-Schamp for each gas
- Missing thermal diffusion (needs OU thermalization)

**When NOT yet to use:**
- Production simulations
- When mobility validation is critical
- With literature CCS values from different gases

**When to use (research only):**
- Fundamental collision physics research
- Comparing kinetic theory to experiment
- When you have gas-specific CCS measurements

**CCS Calculation:**
```python
# Calculate CCS from measured mobility using Mason-Schamp
CCS = (3*q)/(16*N₀*K₀) × √(2π/(μ*kB*T))
# For H3O+ in N2: CCS = 104 Å (not 24.9 Å from He!)
```

---

### 5. Langevin (CollisionModel::Langevin) - EXPERIMENTAL

**Formula:** γ = n·σ_Langevin(v)·v_th·m_red/m_ion

**Status:** **Not validated**

**Critical Limitations:**
- Predicts K = 0.042 cm²/(V·s) for H3O+ in N2 vs measured 3.2 cm²/(V·s) at 100 Td
- Model should better describe mobility at low E/N!
- Overpredicts damping by factor of 6

**When NOT to use:**
- Gases with low polarizability (N2, O2, Ar, He)
- Production simulations
- IMS mobility measurements

**When to use (if at all):**
- Ion-polar molecule interactions (H2O, NH3, etc.)
- Research on polarization effects
- **Validate against experiment first!**

---

## Recommendations by Simulation Type

### IMS/DTIMS Drift Time Measurements
**Use:** `Friction` model
- Exact mobility match
- Validated extensively
- Most reliable for production

### Trajectory Simulations with Diffusion
**Use:** `HSS` model
- Full stochastic treatment
- Accurate diffusion
- Good mobility agreement

### Structure-Dependent Scattering Studies
**Use:** `EHSS` model
- Detailed scattering physics
- Requires CCS from DFT/ab initio

### Research / Model Development
**Use with caution:** `HSD`, `Langevin`
- Validate against known systems first
- Document all assumptions
- Compare to Friction model baseline

---

## For Developers

### Adding New Collision Models

1. Implement in `DampingForce.cpp` or as stochastic handler
2. Validate against at least 3 experimental systems
3. Document limitations and valid parameter ranges
4. Mark as EXPERIMENTAL until extensively validated
5. Compare to Friction model baseline

### Validation Checklist

- [ ] Mobility matches experiment (±5%)
- [ ] Works across pressure range (100-100000 Pa), or specialized for low/high P
- [ ] Works across temperature range (200-400 K)
- [ ] Tested with multiple ion-gas pairs
- [ ] Diffusion coefficient validated (if stochastic)
- [ ] Documented failure modes

---

## References

- Mason, E.A. & McDaniel, E.W. "Transport Properties of Ions in Gases" (1988)
- Revercomb, H.E. & Mason, E.A. "Theory of plasma chromatography/gaseous electrophoresis" Anal. Chem. (1975)
- Langevin, P. "Une formule fondamentale de théorie cinétique" Ann. Chim. Phys. (1905)

---

**Last Updated:** December 2025  
**Version:** v1.0
**Maintainer:** ICARION Development Team
