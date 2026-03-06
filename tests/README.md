# ICARION Test Suite

**IMPORTANT: CTests vs Validation Suite**

This directory contains **CTests** - regression tests for CI/CD (generally short; some benchmarks/GPU tests can run longer).

For scientific validation (long, high-fidelity tests for publications), see `validation/` directory.

## CTests vs Validation Suite

| Aspect | CTests (`tests/`) | Validation Suite (`validation/`) |
|--------|-------------------|----------------------------------|
| **Purpose** | Regression testing, CI/CD | Scientific validation for papers |
| **Runtime** | Seconds to minutes (varies; benchmarks/GPU tests can be longer) | Up to ~30 minutes per test |
| **Ensemble Size** | Varies (10-50k ions depending on test) | Large (1k-10k ions typical) |
| **Output** | Pass/Fail | Plots, tables, quantitative metrics |
| **Frequency** | In CI when configured | Before releases/publications |

---

Test organization and quick pointers for v1.0.0.

## Layout

```
tests/
├── instruments/    # Instrument physics (IMS, Orbitrap, TOF, LQIT, Quadrupole, FTICR, domain transitions)
├── physics/        # Forces, collisions, reactions, space charge
├── integrator/     # RK4/RK45/Boris, SimulationEngine, GPU parity
├── config/         # Loaders/validators for config, species/reactions, waveforms, field arrays
├── io/             # HDF5 writer tests (v1/v2)
├── integration/    # End-to-end main.cpp pipeline
├── gpu/            # GPU-specific kernels (boundaries, space charge, interpolation)
├── unit/           # SoA benchmarks/tests (IonEnsemble, SimulationEngine SoA, OpenMP benchmarks)
├── utils/          # Hashing utilities
├── helpers/        # Test utilities (headers only)
```

## Notable Test Files (by area)

- **Instruments:** `test_instrument_basic.cpp`, `test_ims_drift.cpp`, `test_lqit_stability.cpp`, `test_orbitrap_frequency.cpp`, `test_quadrupole_filter.cpp`, `test_tof_flight_time.cpp`, `test_fticr_cyclotron.cpp`, `test_domain_transition.cpp`
- **Physics – Collisions:** `test_collision_energy_conservation.cpp`, `test_hss_collision_handler.cpp`, `test_ehss_collision_handler.cpp`, `test_ou_collision_handler.cpp`, `test_temperature_scaling.cpp`, `test_multi_gas_collision.cpp`, CPU SoA parity (`test_collision_soa_parity.cpp`), GPU parity/thermalization (`test_gpu_collision_parity.cpp`, `test_gpu_thermalization.cpp`, `test_gpu_ehss_thermalization.cpp`)
- **Physics – Collisions (GPU EHSS):** `test_gpu_ehss_mapping.cpp` (geometry→species map), `test_gpu_ehss_parity.cpp` (CPU vs GPU thermalization sanity)
- **Physics – Forces:** `test_electric_field_force.cpp`, `test_field_model_parity.cpp`, `test_field_model_provider.cpp`, `test_magnetic_damping_forces.cpp`, `test_force_registry.cpp`
- **Physics – Reactions:** `test_reaction_factory.cpp`, `test_stochastic_reaction_handler.cpp`, `test_multi_gas_reaction.cpp`, SoA parity (`test_reaction_soa_parity.cpp`)
- **Physics – Reactions (GPU):** `test_gpu_reaction_parity.cpp` (GPU vs CPU stochastic parity for constant/Arrhenius/modified rates incl. multi-gas; requires `USE_GPU_ACCEL=ON`)
- **Physics – Space Charge:** `test_poisson_solver.cpp`, `test_charge_deposition.cpp`, `test_space_charge_integration.cpp`, `test_space_charge_model_direct.cpp`, `test_space_charge_model_parity.cpp`, `test_space_charge_gpu_model.cpp`
- **Physics – Gas Flow Transport:** `test_gas_flow_transport.cpp`
- **Integrator:** `test_rk4_strategy.cpp`, `test_rk45_strategy.cpp`, `test_boris_strategy.cpp`, `test_domain_manager.cpp`, `test_domain_geometry.cpp`, `test_output_manager.cpp`, `test_simulation_engine.cpp` (SoA parity, birth/transition), GPU integration/parity tests (`test_gpu_integration.cpp`, `test_gpu_rk45.cpp`, `test_gpu_boris.cpp`, `test_rk45_boris_parity.cpp`, `test_gpu_field_interpolation.cpp`), `test_simulation_engine_soa.cpp` (SoA unit/parity)
  - **RK45 per-ion dt/OpenMP:** `test_rk45_per_ion_dt.cpp`
  - **Integration batch dt fallback:** `test_integration_batch_dt.cpp`
  - **SimulationEngine per-ion dt determinism:** `test_simulation_engine_per_ion_dt.cpp`
  - **SimulationEngine RNG/compaction determinism:** `test_engine_rng_compaction.cpp`
  - **SimulationEngine delayed-birth regression:** covered in `test_simulation_engine.cpp` (`All ions delayed birth must still advance time`) to prevent zero-time stalls when all ions are born in the future
- **Config:** `test_config_loader.cpp`, `test_field_array_terms_loader.cpp`, `test_field_array_e2e.cpp`, `test_ion_loader.cpp`, `test_species_loader_unit.cpp`, `test_reaction_loader_unit.cpp`, `test_reaction_validation.cpp`, `test_waveform_loader.cpp`, `test_waveform_types.cpp`, `test_database_integration.cpp`
- **I/O:** `test_hdf5_writer.cpp`, `test_hdf5_writer_v2.cpp`
- **GPU (misc):** `test_gpu_boundaries.cpp`, `test_gpu_space_charge.cpp`, `test_adaptive_interpolation.cpp`
- **Unit/SoA:** `test_ion_ensemble.cpp`, `test_simulation_engine_soa.cpp`, `benchmark_openmp_soa.cpp`, `benchmark_soa_performance.cpp`
- **Utils:** `test_hash.cpp`

## Coverage Snapshot (v1.0.0)

- **Instruments:** IMS, LQIT, Orbitrap, TOF, Quadrupole, FTICR, multi-domain transitions
- **Integrators:** RK4, RK45 (adaptive), Boris; SimulationEngine CPU/GPU parity and field interpolation
- **Forces:** Electric, magnetic, damping, space charge; registry superposition
- **Collisions:** HSS, EHSS, OU, deterministic damping; geometry/mixture cases; GPU parity/thermalization
- **Reactions:** Stochastic handler, factory, multi-gas mixtures, validation
- **Space Charge:** Direct Poisson solver, charge deposition, integration tests (CPU); GPU helper tests
- **Gas Flow Transport:** SIFT-MS like ion transport without field  
- **Config:** Loaders/validators for configs, species, reactions, waveforms, field arrays
- **I/O:** HDF5 writer v2 (primary) and legacy writer regression
- **SoA/Unit:** IonEnsemble, SimulationEngine SoA, OpenMP/SoA benchmarks

## Running Tests

From `build/`:
```bash
ctest --output-on-failure
```

Run a specific binary:
```bash
./tests/instruments/test_ims_drift
./tests/physics/forces/test_electric_field_force
```

Catch2 tags (where present):
```bash
./tests/instruments/test_instrument_basic "[fast]"
ctest -E "slow|performance"   # skip marked slow/perf tests if any
```

If CTest times out in constrained environments, run the desired binary directly (all CTests are standalone executables under `build/tests/...`). Some SimulationEngine tests write HDF5 files to `/tmp`.

## Notes and Known Status

- GPU tests are only built with `-DUSE_GPU_ACCEL=ON` (defines `ICARION_USE_GPU`); some assume a CUDA device and do not skip if none is available.
- EHSS/HSS mobility/thermalization tests have stochastic tolerance; occasional seed sensitivity may appear—rerun if needed.
- Space-charge GPU model is opt-in (`physics.enable_space_charge_gpu`) and falls back automatically; GPU tests validate helper behavior while CPU-only builds exercise the stub.
- HDF5 writer v2 is the primary path; legacy writer tests remain for regression coverage.

For full test intent/details, see the individual test files.
