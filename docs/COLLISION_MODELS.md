# ICARION Collision Models - Status and Recommendations

## Overview

ICARION implements multiple collision models for ion-neutral interactions. This document provides guidance on which models to use for production simulations vs experimental testing.

## Model Status Summary

| Model | Status | Use Case | Validation |
|-------|--------|----------|------------|
| **NoCollisions** | **STABLE** | Vacuum simulations | N/A |
| **Friction** | **USABLE (mobility required; deterministic)** | Drift-only when K0 is known; no diffusion | IMS drift sanity test; validation report includes low E/N comparisons |
| **HSS** | **VALIDATED (thermalization + drift)** | Stochastic collisions with diffusion | Validation suite + IMS drift tests (see validation report); constant-CCS limits at high E/N |
| **EHSS** | **VALIDATED (thermalization), drift bias** | Geometry-resolved stochastic collisions | Thermalization validated; drift shows systematic bias (see validation report) |
| **HSD** | **EXPERIMENTAL** | Deterministic hard sphere (no diffusion) | Basic IMS test only |
| **Langevin** | **EXPERIMENTAL** | Polarizable gases only | N2 drift test shows large error (manual) |

---

## Recommended Models

### 1. Friction (CollisionModel::Friction) **Recommended when mobility is known and deterministic drift is acceptable**

**Formula:** γ = q/(K·m) where K = K₀·(n₀/n)

**Advantages:**
- Based on experimentally measured, modelled (e.g., IMoS [https://www.imospedia.com/imos/], MobCal-MPI 2.0 [https://github.com/HopkinsLaboratory/MobCal-MPI]) or literature mobility (K₀)
- Mason-Schamp equation foundation
- Deterministic drift speed when K0 is known
- Number density correction included automatically

**Validation:**
- IMS drift sanity test exists in `tests/instruments/test_ims_drift.cpp` (single config, broad tolerance)
- Validation report compares Friction against Mason-Schamp across low-field sweeps; treat it as a mobility surrogate
- Model provides deterministic drift speed; diffusion is not modeled

**Usage:**
```cpp
config::CollisionModel::Friction
// Requires: ion.reduced_mobility_cm2_Vs (from experiment/simulation/literature)
```

**When to use:**
- When reduced mobility is known and you want deterministic drift speed (i.e., correct arrival times)
- As a baseline in IMS/DTIMS comparisons
- When diffusion is not critical

**Important limitations:**
- Requires `ion.reduced_mobility_cm2_Vs`; if missing, gamma is zero (no damping).
- For gas mixtures, CCS_HSS (or `ion.CCS_m2` fallback) is only used to scale per-gas mobility, not to replace K0.
- Not a stochastic model; diffusion is not modeled. Can therefore model peak position, but no reasonable peak shapes/distributions.

---

### 2. HSS (Hard-Sphere Stochastic)

**Formula:** Stochastic, isotropic collision -> velocity randomization

**Advantages:**
- Full stochastic treatment (diffusion included)
- Validated thermalization behavior (see [VALIDATION_REPORT_v1.0.0.md](../validation/VALIDATION_REPORT_v1.0.0.md))
- IMS drift tests use broad tolerances; validation report shows ~10% accuracy once collisions are frequent

**Validation:**
- IMS drift test for H3O+ is manual (`[.]` tag) and reports close agreement at one point
- Validation report shows a low-field accuracy corridor once enough collisions occur
- Cross-tool comparisons to IDSimF/SIMION are ongoing; treat any external agreement as qualitative until a public report is available
- Thermalization and diffusion behavior are validated in the collision test suite

**Note**
- Simplified elastic collision model
- Does not account for molecular structure (contrast to EHSS)
- Does not account for inelastic collisions

**Usage:**
```cpp
config::CollisionModel::HSS
// Requires: ion.CCS_m2 (collision cross-section)
```

**When to use:**
- Need accurate thermal diffusion/peak shapes
- Trajectory-level accuracy required
- Modeling real experimental conditions
- No DFT-optimized structures available 

---

### 3. EHSS (Explicit hard-sphere scattering)

**Formula:** Explicit scattering with stochastic treatment

**Validation:**
- Thermalization behavior is validated in the collision test suite and validation report
- IMS drift test for H3O+ is manual (`[.]` tag) and shows large deviation; the validation report notes a systematic drift bias; however, note that this is still a simplified model that neglects long-range interactions
- More accurate angular scattering than HSS when molecular structure is important

**Usage:**
```cpp
config::CollisionModel::EHSS
// Prefer: CCS_EHSS map or geometry; falls back to derived/reference CCS with warnings
```

**Note**
- Does not account for long-range interactions
- Does not account for inelastic collisions
- GPU helper supports HSS/EHSS batch collisions; dispatch uses active-ion threshold (default 5000).

**Optional orientation-sampled mode**
- Provide `EHSS_samples_file` in the species database to sample orientation-dependent
  projection areas per collision. The sampled orientation is reused in the EHSS
  kernel (force-hit path) to keep rate and scattering consistent.

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
- Requires accurate CCS for the target gas (not a reference gas)
- Gas-specific CCS should be derived per gas (Mason-Schamp or `ccs_precompute`)
- Missing thermal diffusion

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

**Note**
- might only provide good agreement with experimental or literature mobility values under certain conditions (e.g., when long-range interactions can be neglected)

---

### 5. Langevin (CollisionModel::Langevin) - EXPERIMENTAL

**Formula:** γ = n·σ_Langevin(v)·v_th·m_red/m_ion

**Status:** **Not validated**

**Critical Limitations:**
- IMS drift fixture for N2 shows order-of-magnitude mobility error (manual test, `[.]` tag)
- Model should better describe mobility at low E/N in polarizable gases, but N2 is non-polar
- Overpredicts damping in the N2 fixture compared to Friction

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
**Use:** `HSS` for stochastic drift + diffusion; `Friction` when K0 is known and you want deterministic drift
- HSS is the primary validated stochastic model for drift in the validation report
- Friction is a mobility surrogate without diffusion

### Trajectory Simulations with Diffusion
**Use:** `HSS` model
- Full stochastic treatment
- Accurate diffusion
- Drift accuracy in the low-field validation envelope

### Structure-Dependent Scattering Studies
**Use:** `EHSS` model
- Detailed scattering physics
- Requires geometry or CCS_EHSS map (DFT/ab initio recommended)

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
**Version:** 1.0.0
**Maintainer:** Christoph Schaefer
