# CRITICAL BUG: HSS Collision Model Does Not Limit Drift Velocity

**Status:** RELEASE BLOCKER  
**Date:** 2025-11-25  
**Discovered in:** feature/instrument-physics-tests branch

## Problem

When simulating ion drift in IMS with HSS collisions enabled, ions gain excessive velocities (1000-7000 m/s) instead of reaching the expected drift velocity limited by collisions.

## Expected Behavior

For IMS drift with collisions, the drift velocity should be:

```
v_d = K0 × (E/N) × N0
```

Where:
- K0 = reduced mobility (e.g., 2.8 cm²/Vs for H3O+ in N2)
- E/N = reduced field strength in Townsend (Td)
- N0 = Loschmidt constant = 2.6868×10²⁵ m⁻³

**Example:** For E/N = 10 Td, K0 = 2.8 cm²/Vs:
```
v_d = 2.8×10⁻⁴ m²/Vs × 10 × 2.6868×10²⁵ m⁻³ × 10⁻²¹ = 75.2 m/s
```

## Observed Behavior

Test case: `test_ims_drift.cpp` - IMS with HSS collisions at E/N=10 Td
- **Expected:** v_d ≈ 75 m/s, ion reaches 50mm in ~0.67 ms
- **Observed:** Ion velocity ~1256 m/s, gains ~6984 m/s radial velocity
- **Result:** Ion hits radial wall after only 14 μs due to excessive thermal/drift velocities

## Evidence

From test output:
```
Final vel: (-1176.91, 6984.34, 1623.49) m/s
Final time: 0.014195 ms
Active: 0 (hit wall)
```

Velocities are 10-100x too high!

## Root Cause (Hypothesis)

The HSS collision handler (`HSSCollisionHandler.cpp`) likely has one of these issues:

1. **Collision frequency too low:** Not enough collisions happening to limit velocity
2. **Post-collision velocity wrong:** Isotropic scattering not properly thermalizing
3. **Integration issue:** SimulationEngine not applying collision handler correctly to single ions
4. **Missing drift velocity limit:** No mechanism to enforce v_d = K0 × E/N × N0

Note: Standalone `test_hss_collision_handler` tests pass, suggesting issue is in SimulationEngine integration.

## Impact

- **IMS simulations invalid** with collisions enabled
- **Drift time measurements wrong** by orders of magnitude
- **Ion losses unrealistic** due to excessive diffusion
- **Cannot validate mobility measurements**

## Workaround

For now, tests use `NoCollisions` model for validation. This works for:
- Electric field acceleration tests
- Field scaling tests
- TOF flight time (vacuum)

But blocks validation of:
- IMS mobility measurements
- Collision-limited drift
- Diffusion/thermalization

## Action Items

- [ ] Debug HSS collision handler in SimulationEngine context
- [ ] Check collision probability calculation (P = 1 - exp(-n·σ·v·dt))
- [ ] Verify ion domain properties are set correctly (gas mass, temperature, density)
- [ ] Add detailed logging to collision handler
- [ ] Create minimal reproduction case (single ion, HSS, monitor velocity vs time)
- [ ] Fix before v1.0 release

## References

- Test file: `tests/instruments/test_ims_drift.cpp`
- Collision handler: `src/core/physics/collisions/HSSCollisionHandler.cpp`
- Integration: `src/core/integrator/SimulationEngine.cpp`
- Working unit tests: `tests/physics/collisions/test_hss_collision_handler.cpp`
