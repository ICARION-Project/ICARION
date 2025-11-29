# ICARION Troubleshooting Guide

**Version:** 1.0  
**Last Updated:** November 26, 2025

This guide helps users diagnose and resolve common issues when running ICARION simulations.

---

## Table of Contents

1. [Ion Loss and Deactivation Issues](#ion-loss-and-deactivation-issues)
2. [Performance Issues](#performance-issues)
3. [Numerical Stability](#numerical-stability)
4. [Collision Physics](#collision-physics)

---

## Ion Loss and Deactivation Issues

### Ions Immediately Deactivated at Entrance

**Symptom:**
- Ions starting at entrance boundary (e.g., z=0) are marked inactive after first collision
- Ion final position shows small negative z-coordinate (e.g., z = -0.0001 mm)
- Problem worse at low pressures (< 1000 Pa)
- More common in IMS, TOF, or other drift instruments

**Example:**
```
Final ion position: z = -0.081 mm
Final ion velocity: vz = -202 m/s
Ion active: NO
Simulation time: < 1 μs
```

**Cause:**

When ions start exactly at a domain boundary (e.g., ion at z=0, domain origin at z=0):

1. Ion is technically at the entrance boundary
2. First collision randomizes velocity with thermal component
3. Random velocity may have backward component (vz < 0)
4. Ion moves slightly backward (z < 0)
5. Boundary check fails → ion deactivated

This is especially problematic at **low pressure** where:
- Thermal velocities from collisions are large compared to drift velocity
- Mean free path is long → few collisions → large random kicks
- Diffusion dominates over directed drift

**Solution: Add Domain Buffer Zones**

**QUICK FIX:** Shift the domain origin backward to create an entrance buffer.

**Example for 10 cm IMS drift tube:**

```json
{
  "domains": [{
    "name": "ims_drift",
    "instrument": "IMS",
    "geometry": {
      "origin_m": [0.0, 0.0, -0.02],  ← Shift 2cm backward
      "length_m": 0.14,                ← Original 10cm + 4cm buffer (2cm each side)
      "radius_m": 0.5                  ← Wide radius to prevent radial losses
    },
    "fields": {
      "dc": {
        "axial_V": {
          "constant_value": 1000.0      ← Voltage across 10cm effective length
        }
      }
    }
  }]
}
```

**How it works:**
- Ion starting position: z = 0 (lab frame)
- Domain extends from z = -0.02 m to z = 0.12 m
- Ion is now **2 cm inside** the domain entrance
- Backward collisions at z=0 → z=-0.001 m are still within bounds
- Provides safety margin for thermal diffusion

**Best Practices:**

| Pressure Range | Recommended Buffer | Reason |
|---------------|-------------------|---------|
| < 100 Pa | ≥ 2-5 cm | Long mean free path, large thermal kicks |
| 100-1000 Pa | ≥ 1-2 cm | Moderate collisions |
| > 10 kPa | ≥ 0.5-1 cm | Short mean free path, thermal diffusion small |

**Additional Tips:**
- Use **wide radius** (≥ 5× length) to prevent radial losses from diffusion
- For very low pressure (< 10 Pa): consider even larger buffers or higher pressure
- Check `tests/instruments/test_ims_drift.cpp` for working example

**Physical Intuition:**

Think of it like a river entering a lake:
- **Without buffer**: Ions start right at the shore → any backward motion = out of bounds
- **With buffer**: Ions start 2m into the lake → backward fluctuations stay within bounds
- Buffer size should match the "diffusion length" of your ions

**Future Development (v1.1+):**

Alternative boundary conditions planned:
- **Reflecting boundaries**: Elastic collisions at walls (no ion loss)
- **Emitting boundaries**: Continuous ion injection at entrance
- **Periodic boundaries**: Wraparound for bulk gas simulations

Currently (v1.0), only **absorbing boundaries** are implemented.

---

### High Radial Losses in Drift Tubes

**Symptom:**
- Large fraction of ions hit radial boundaries and deactivate
- Final ion positions at r ≈ radius_m
- Problem worse at low pressures or long drift times

**Cause:**

Thermal diffusion causes radial spreading:
- Diffusion coefficient: D ∝ 1/pressure
- Radial displacement grows as √(Dt)
- Low pressure + long drift time → ions reach walls

**Solution:**

1. **Increase domain radius:**
   ```json
   "geometry": {
     "radius_m": 0.5  // 50 cm radius for 10 cm length (5:1 ratio)
   }
   ```

2. **Reduce drift time:**
   - Increase pressure (more collisions = slower drift BUT less diffusion per collision)
   - Shorten drift length
   - Increase electric field (faster transit)

3. **Add radial confinement:**
   - Use RF radial confinement (ion guides)
   - Add magnetic field (cyclotron confinement)

---

## Performance Issues

### Simulation Running Extremely Slow

**Symptom:**
- Simulation appears to hang
- Very slow progress (< 1% per minute)
- High CPU usage but no output

**Common Causes:**

#### 1. Too Small Timestep with Long Simulation Time

**Check:**
```json
"simulation": {
  "dt_s": 1e-12,           ← 1 ps timestep
  "total_time_s": 0.001    ← 1 ms simulation
}
// → 1,000,000,000 steps!
```

**Solution:**
- Increase timestep: `dt_s = 1e-9` (1 ns) for typical ion dynamics
- Reduce simulation time to minimum needed
- Use RK45 adaptive integrator for efficiency

#### 2. Very Low Pressure with Small Timestep

**Issue:**
- Low pressure (< 1 Pa) → very rare collisions
- Simulation still runs at full timestep resolution
- Most timesteps do nothing except advance time

**Solution:**
- Increase pressure if physically reasonable
- Use larger timestep for low-pressure simulations
- Consider collisionless simulation if pressure < 1e-6 Pa

#### 3. Too Many Ions with Space Charge

**Check:**
```json
"physics": {
  "space_charge": {
    "enabled": true,
    "method": "direct"      ← O(N²) scaling!
  }
}
"ions": {
  "species": [{
    "count": 10000          ← 10k ions → 100M interactions per step!
  }]
}
```

**Solution:**
- Use grid-based space charge: `"method": "grid"`
- Reduce ion count for testing
- Disable space charge if not critical: `"enabled": false`

---

## Numerical Stability

### Ions Gain Excessive Energy

**Symptom:**
- Ion kinetic energy grows exponentially
- Final velocities >> expected thermal or drift velocities
- Simulation eventually crashes with NaN positions

**Cause:**

Timestep too large for:
- RF frequencies (Ω × dt >> 1)
- Strong magnetic fields (ω_c × dt >> 1)
- Steep electric field gradients

**Solution:**

1. **Reduce timestep:**
   - RF: `dt ≤ 0.01 / frequency_Hz`
   - Magnetic: `dt ≤ 0.01 / (q B / m)`

2. **Use appropriate integrator:**
   - **RK4**: Good general-purpose, fixed timestep
   - **RK45**: Adaptive timestep, automatically adjusts
   - **Boris**: Best for strong magnetic fields (symplectic)

3. **Check field sanity:**
   ```json
   "fields": {
     "rf": {
       "frequency_Hz": 1e6,     ← 1 MHz
       "amplitude_V": 1000      ← Check not too strong
     }
   }
   ```

---

## Collision Physics

### Ions Don't Thermalize Correctly

**Symptom:**
- Ion kinetic energy equilibrates to ~70-80% of expected thermal energy
- Problem persists across different seeds, timesteps, gas molecules
- Energy distribution is not Maxwell-Boltzmann

**Cause:**

**Fixed in v1.0** (commit 92d29c1):

Older versions calculated collision probability with bulk gas velocity, but applied collision with different sampled neutral velocity. This velocity inconsistency created systematic bias.

**Solution:**

Update to v1.0 or later. The fix:
1. Sample neutral velocity FIRST from Maxwell-Boltzmann distribution
2. Calculate collision probability using SAME sampled velocity
3. Apply collision with that consistent velocity

**Verification:**

Run thermalization test:
```bash
./build/tests/test_physics_validation "OU thermalization"
```

Expected: Mean kinetic energy = (3/2) k_B T within 5% after ~100 collisions.

---

### Collision Rates Too High or Too Low

**Symptom:**
- Mean free time << or >> expected from kinetic theory
- Mobility measurements differ significantly from literature

**Check:**

1. **Collision cross-section (CCS):**
   ```json
   "species": [{
     "id": "H3O+",
     "CCS_A2": 24.9        ← Check literature value (Å²)
   }]
   ```

2. **Gas pressure and temperature:**
   ```json
   "environment": {
     "gas_molecule": "N2",
     "pressure_Pa": 101325,   ← 1 atm
     "temperature_K": 300
   }
   ```

3. **Number density calculation:**
   ```
   n = P / (k_B T)
   n = 101325 / (1.38e-23 × 300) = 2.45e25 m⁻³
   ```

**Expected mean free path:**
```
λ = 1 / (n σ)
λ = 1 / (2.45e25 × 24.9e-20) = 16.4 μm  (for H3O+ in N2 at 1 atm)
```

**Expected collision rate:**
```
ν = v_mean / λ
v_mean ≈ √(8 k_B T / π m) ≈ 580 m/s  (for H3O+)
ν ≈ 580 / 16.4e-6 ≈ 35 MHz  (1 collision every 28 ns)
```

---

## Getting More Help

**If problem persists:**

1. **Check examples:**
   - Browse `examples/` for similar configurations
   - Run example configs to verify installation

2. **Enable debug output:**
   ```json
   "simulation": {
     "write_interval": 10,    ← Write every 10 steps for diagnosis
     "verbose": true
   }
   ```

3. **Minimal reproducer:**
   - Simplify config to smallest failing case
   - Single ion, short time, simple geometry
   - Check if problem is in physics or config

4. **Validation tests:**
   ```bash
   cd build
   ctest -R validation  # Run validation suite
   ```

5. **GitHub Issues:**
   - Provide minimal config JSON
   - Include error messages and output
   - Specify ICARION version and platform

**Reference Documentation:**
- Configuration: `docs/INPUT_FORMAT_SPECIFICATION.md`
- Physics models: `docs/VALIDATION_METHODS.md`
- Developer info: `docs/DEVELOPERS_GUIDE.md`

---

**Document Status:** Living document, updated as new issues are identified and resolved.
