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

### High-Pressure Stochastic Runs

At high pressure or large `dt_s`, one macro-step can contain more than one physical collision. ICARION provides two approximation knobs:

- `physics.collision_subcycles_per_step`: splits the collision update into equal micro-steps and recomputes rates in each micro-step.
- `physics.collision_multi_event_mode`: enables a practical multi-collision approximation and uses at least `collision_max_events_per_step` micro subcycles.
- `physics.collision_max_events_per_step`: legacy field name; in multi-event mode this is a minimum micro-subcycle count, not a guaranteed upper bound on physical continuous-time collision events.

These options reduce bias from large steps, but they are still approximations. The most accurate setup is to choose `dt_s` such that less than one collision per ion per step is likely.

There are two separate accuracy criteria:

- Collision-count/statistics accuracy: each collision micro step should keep
  `mu = lambda * dt_collision` below the range validated for the selected
  model and pressure regime, where `lambda` is the instantaneous collision
  rate. If a sub-step allows at most one accepted event, the missed multiple
  event error is second order in `mu`, while the relative bias in the expected
  event count is roughly proportional to `mu / 2` for small `mu`. As a generic
  orientation, `mu <= 0.02` is around percent level, `mu <= 0.05` is a few
  percent, and `mu >= 0.1` is no longer a strict setting. Prefer measured
  validation thresholds over these generic numbers.
- Trajectory accuracy: `dt_s` must still resolve deterministic motion. Collision
  subcycling improves collision counts and velocity relaxation inside a
  macro-step, but it does not make an arbitrarily large RK/global step valid.
  The global step must resolve RF or field timescales, momentum relaxation,
  relevant field gradients, time to walls or domain boundaries, and intervals
  with strong velocity changes.

In practice, use both checks: keep `max(lambda * dt_collision)` below the
validated collision threshold, and keep the global `dt_s` small enough for the
fields, geometry, boundary crossings, and deterministic acceleration.

Current public validation coverage for this release:

| Model | Tested Regime | Maximum `mu` | Deviation | Runtime Gain |
| --- | --- | --- | --- | --- |
| HSS multi-event dispatch | CTest synthetic handlers with fixed `dt_s`, no physical rate model | n/a | Verifies handler call count, collision microstep size, event accounting, and that multi-event mode is not a global event cap | Not measured |
| HSS drift convergence gate | `validation/scripts/physics/validate_multi_event_mode.py`, short IMS drift envelope against a smaller-`dt` reference | Not reported by the gate | Checks median arrival and active ion fraction envelope; not a Poisson count or relaxation limit proof | Reported in generated gate output, not fixed in docs |
| Full high pressure physical envelope | Poisson collision counts, relaxation against a small-`dt` reference, mobility/temperature convergence, and error vs. `mu` sweeps | Not established as release evidence | Not established as release evidence | Not established as release evidence |

Do not interpret collision subcycling as permission to take a large global
integration step. The subcycles improve stochastic collision statistics inside
one macro-step; they do not automatically feed collision impulses back into the
deterministic trajectory at exact event times.

### InteractionPotentialModel Offline Tables

`InteractionPotentialModel` uses precomputed offline sample files referenced from the species database via `ipm_samples_file`. The runtime interpolates the table for the current relative speed and samples from stored momentum-transfer possibilities. For one relative speed and orientation, multiple momentum-transfer outcomes remain possible because the impact parameter is still sampled; the table stores that distribution rather than a single deterministic transfer.

#### Methodological context

Classical trajectory calculations for CCS, mobility, and collision integrals are established methods represented by MOBCAL/MobCal-MPI, IMoS, Collidoscope, and CoSIMS. ICARION does not claim those methods or offline scattering calculations as novel. Its use is to retain stochastic momentum-transfer tables for sampling within time-resolved instrument simulations. See [Related and Complementary Software](RELATED_SOFTWARE.md) and the [IPM precomputation guide](../rtd/ipm-precomputation.md).

Files produced by `interaction_potential_precompute` also contain structured
reproducibility information under `/metadata`: software/build/system facts, RNG
and cell-seed provenance, resolved species and neutral parameters, resolved
precompute settings, input filenames and SHA-256 hashes, embedded text inputs, and
checkpoint/completion state. This hierarchy is additive; the version-1 root
attributes and numerical lookup datasets used by the runtime are unchanged.
See `tools/README.md` for `h5dump` and h5py inspection examples.

Format version 1 separates the event cross section from the diagnostic momentum-transfer cross section. Runtime collision probabilities use `sigma_event_m2`, while `sigma_mt_m2` records the corresponding momentum-transfer integral for auditing. Stored `dp_samples` and `dp_stats` are sampled over impact-parameter area, so the momentum-transfer factor is not applied a second time during runtime sampling.

Use full-CDF sample tables for scientific production runs. The compact
`dp_stats` fallback reconstructs momentum kicks from per-cell moments only; it
does not preserve parallel/perpendicular correlations, non-Gaussian tails, or
per-event kinematic bounds from the offline trajectories. The precompute tool
therefore writes full-CDF tables by default; `--compact-dp-stats` is reserved
for lightweight debugging or legacy compatibility.

The precompute tool uses `qmc` orientation sampling by default. `lebedev` remains
available as an explicit option, but v1.1 runtime sampling draws stored
orientations uniformly by index and does not consume Lebedev quadrature weights.
Tables store `attempted_trajectories`, `accepted_trajectories`, and
`rejected_non_asymptotic` per cell, plus `max_energy_step_error` and
`max_energy_cumulative_error` diagnostics for the offline trajectory integrator.
Precompute aborts if the non-asymptotic fraction exceeds
`--max-non-asymptotic-frac` (default: `0.005`). The stored CDF and `dp_stats`
payloads contain the accepted asymptotic trajectories only; for strict
production tables, require zero rejects. A nonzero reject fraction introduces
an unresolved conditioning bias; the threshold bounds only the fraction of
unresolved trajectories, not the resulting error in transport moments.

For relative speeds between stored bins, the runtime randomly chooses one of
the neighboring log-speed bins. The sampled absolute momentum kick belongs to
the selected bin speed, so production tables should be checked for convergence
with respect to velocity-bin density as well as orientations, trials, `b_max_m`,
and CDF sample count.

Supported runtime controls include:

- `physics.ipm_orientation_mode`: `random` samples a new stored orientation per collision; `fixed` uses `ipm_fixed_orientation_index`.
- `physics.ipm_vrel_log_prefix`: optional relative-speed histogram CSV output.
- `physics.ipm_momentum_log_prefix`: optional momentum-transfer diagnostic CSV output.

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
