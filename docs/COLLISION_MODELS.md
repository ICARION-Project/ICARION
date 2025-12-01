# ICARION Collision Models - Status and Recommendations

## Overview

ICARION implements multiple collision models for ion-neutral interactions. This document provides guidance on which models to use for production simulations vs experimental testing.

## Model Status Summary

| Model | Status | Use Case | Validation |
|-------|--------|----------|------------|
| **NoCollisions** | **PRODUCTION READY** | Vacuum simulations | N/A |
| **Friction** | **PRODUCTION READY** | All simulations | Validated (±40% tolerance) |
| **HSS** | **PRODUCTION READY** | Stochastic collisions | Validated (±25% tolerance) |
| **EHSS** | **PRODUCTION READY** | Molecular structure resolved scattering | Validated (±50% tolerance) |
| **HSD** | **EXPERIMENTAL** | Deterministic hard sphere | Partially validated |
| **Langevin** | **EXPERIMENTAL** | Polarizable gases only | Not validated for N2 |

---

## Production-Ready Models

### 1. Friction (CollisionModel::Friction) **RECOMMENDED**

**Formula:** γ = q/(K·m) where K = K₀·(n₀/n)

**Advantages:**
- Based on experimentally measured mobility (K₀)
- Mason-Schamp equation foundation
- **Exact agreement** with IMS drift measurements (356 m/s measured vs 356 m/s theory in N2)
- Works for all gas types (polar and non-polar)
- Density correction included automatically

**Validation:**
- H3O+ in N2 at 1000 Pa: Within ±40% tolerance (test_ims_drift.cpp)
- Consistent with literature mobility values in N2
- Theory provides good approximation, but does not account for diffusion (no randomization)

**Usage:**
```cpp
config::CollisionModel::Friction
// Requires: ion.reduced_mobility_cm2_Vs (from experiment/literature)
```

**When to use:**
- At high pressures (> 50 mbar) and when only mobility is of interest, not ion cloud shape
- IMS, DTIMS
- Any simulation where experimental mobility is known and diffusion is not critical

---

### 2. HSS (Hard-Sphere Stochastic)

**Formula:** Stochastic velocity randomization after collision

**Advantages:**
- Full stochastic treatment (diffusion included)
- Validated thermalization behavior
- Good agreement with mobility (±1%)

**Validation:**
- H3O+ in N2: 360 m/s measured vs 357 m/s theory (0.8% error)
- Diffusion coefficient correct
- Full stochastic treatment provides realistic thermal behavior

**Usage:**
```cpp
config::CollisionModel::HSS
// Requires: ion.CCS_m2 (collision cross-section)
```

**When to use:**
- Need accurate thermal diffusion
- Trajectory-level accuracy required
- Modeling real experimental conditions
- No DFT simmulations available 

---

### 3. EHSS (Exact hard sphere scattering)

**Formula:** Exact scattering angles + stochastic

**Validation:**
- H3O+ in N2: 440 m/s vs 357 m/s theory (23% deviation)
- Note: Higher velocity deviation expected due to molecular geometry effects
- Should be the most accurate physical model available when molecular structure is important, as it includes detailed scattering physics

**Usage:**
```cpp
config::CollisionModel::EHSS
// Requires: ion.CCS_m2
```

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
- Example: H3O+ CCS = 104 Ų in N2 (not 24.9 Ų from He!)
- CCS must be calculated from Mason-Schamp for each gas
- Missing thermal diffusion (needs OU thermalization)

**When NOT to use:**
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
# For H3O+ in N2: CCS = 104 Ų (not 24.9 Ų from He!)
```

---

### 5. Langevin (CollisionModel::Langevin) - EXPERIMENTAL

**Formula:** γ = n·σ_Langevin(v)·v_th·m_red/m_ion

**Status:** **Not validated for non-polar gases**

**Critical Limitations:**
- **Only valid for high polarizability molecules!**
- N2 has low polarisability → model gives 76x wrong result!
- Predicts K = 0.042 cm²/(V·s) vs measured 3.2 cm²/(V·s)
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

## Testing and Validation

Current test results (H3O+ in N2 at 1000 Pa, 300K):

```
Model          | Expected    | Tolerance   | Status
---------------|-------------|-------------|--------
NoCollisions   | N/A         | N/A         | PRODUCTION
Friction       | 356 m/s     | ±40%        | PRODUCTION
HSS            | 357 m/s     | ±25%        | PRODUCTION
EHSS           | Variable    | ±50%        | PRODUCTION (geometry-dependent)
HSD            | (untested)  | Unknown     | EXPERIMENTAL
Langevin       | (fails N2)  | N/A         | EXPERIMENTAL (polar gases only)
```

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

**Last Updated:** 2025-11-26  
**Version:** v1.0
**Maintainer:** ICARION Development Team
