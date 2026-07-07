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

Test organization and quick pointers for the current v1.1.0 branch.

## Layout

```
tests/
├── instruments/    # Instrument physics (IMS, Orbitrap, TOF, LQIT, Quadrupole, FTICR, domain transitions)
├── physics/        # Forces, collisions, reactions, space charge
├── integrator/     # RK4/RK45/Boris, SimulationEngine, domains, boundaries, RNG
├── config/         # Loaders/validators for config, species/reactions, waveforms, field arrays
├── io/             # HDF5 writer tests
├── integration/    # End-to-end main.cpp pipeline
├── gpu/            # GPU-specific kernels (boundaries, space charge, interpolation)
├── unit/           # SoA benchmarks/tests (IonEnsemble, SimulationEngine SoA, OpenMP benchmarks)
├── utils/          # Hashing utilities
├── helpers/        # Test utilities (headers only)
```

## CTests by Area

Current CPU-only build (`ctest -N`) exposes 70 tests. Each entry lists the CTest name followed by the source file.

- **Unit/SoA/benchmarks:** `test_ion_ensemble` (`test_ion_ensemble.cpp`), `test_simulation_engine_soa` (`test_simulation_engine_soa.cpp`), `benchmark_soa_performance` (`benchmark_soa_performance.cpp`), `benchmark_openmp_soa` (`benchmark_openmp_soa.cpp`), `test_domain_manager_lookup` (`test_domain_manager_lookup.cpp`), `benchmark_domain_lookup` (`benchmark_domain_lookup.cpp`)
- **Config:** `ConfigLoader` (`test_config_loader.cpp`), `ConfigOverride` (`test_config_override.cpp`), `SpeciesLoaderSmoke` (`test_species_loader.cpp`), `SpeciesLoaderUnit` (`test_species_loader_unit.cpp`), `ReactionLoaderUnit` (`test_reaction_loader_unit.cpp`), `ReactionValidation` (`test_reaction_validation.cpp`), `DatabaseIntegration` (`test_database_integration.cpp`), `IonLoader` (`test_ion_loader.cpp`), `FieldArrayTermsLoader` (`test_field_array_terms_loader.cpp`)
- **I/O and utils:** `HDF5Writer` (`test_hdf5_writer.cpp`), `test_hash` (`test_hash.cpp`)
- **Forces and fields:** `ForceRegistry` (`test_force_registry.cpp`), `ElectricFieldForce` (`test_electric_field_force.cpp`), `TIMSAxialFieldModel` (`test_tims_axial_field_model.cpp`), `FieldModelParity` (`test_field_model_parity.cpp`), `FieldModelProvider` (`test_field_model_provider.cpp`), `MagneticAndDampingForces` (`test_magnetic_damping_forces.cpp`), `ForceIntegration` (`test_force_integration.cpp`), `GasFlowTransport` (`test_gas_flow_transport.cpp`)
- **Collisions:** `CollisionGeometry` (`core/test_collision_geometry.cpp`), `VelocitySampling` (`core/test_velocity_sampling.cpp`), `CollisionKernels` (`core/test_collision_kernels.cpp`), `CollisionFactory` (`test_collision_factory.cpp`), `EHSSCollisionHandler` (`test_ehss_collision_handler.cpp`), `HSSCollisionHandler` (`test_hss_collision_handler.cpp`), `OUCollisionHandler` (`test_ou_collision_handler.cpp`), `GeometryUtils` (`test_geometry_utils.cpp`), `CollisionEnergyConservation` (`test_collision_energy_conservation.cpp`), `EHSSSamplesLoader` (`test_ehss_samples_loader.cpp`), `InteractionPotentialCollisionHandler` (`test_interaction_potential_collision_handler.cpp`), `TemperatureScaling` (`test_temperature_scaling.cpp`), `MultiGasCollision` (`test_multi_gas_collision.cpp`), `CollisionSoAParity` (`test_collision_soa_parity.cpp`)
- **Reactions:** `StochasticReactionHandler` (`test_stochastic_reaction_handler.cpp`), `ReactionFactory` (`test_reaction_factory.cpp`), `MultiGasReaction` (`test_multi_gas_reaction.cpp`), `ReactionSoAParity` (`test_reaction_soa_parity.cpp`)
- **Space charge:** `PoissonSolver` (`test_poisson_solver.cpp`), `ChargeDeposition` (`test_charge_deposition.cpp`), `SpaceChargeIntegration` (`test_space_charge_integration.cpp`), `SpaceChargeDirectModel` (`test_space_charge_model_direct.cpp`), `SpaceChargeFactory` (`test_space_charge_factory.cpp`), `SpaceChargeGridModelParity` (`test_space_charge_model_parity.cpp`), `SpaceChargeGPUModelTest` (`test_space_charge_gpu_model.cpp`)
- **Integrator/engine:** `test_rk4_strategy` (`test_rk4_strategy.cpp`), `test_rk45_strategy` (`test_rk45_strategy.cpp`), `test_boris_strategy` (`test_boris_strategy.cpp`), `test_domain_manager` (`test_domain_manager.cpp`), `test_output_manager` (`test_output_manager.cpp`), `test_simulation_engine` (`test_simulation_engine.cpp`), `test_simulation_engine_per_ion_dt` (`test_simulation_engine_per_ion_dt.cpp`), `test_simulation_engine_adaptive_sc` (`test_simulation_engine_adaptive_sc.cpp`), `test_domain_geometry` (`test_domain_geometry.cpp`), `test_boundary_actions` (`test_boundary_actions.cpp`), `test_engine_rng_compaction` (`test_engine_rng_compaction.cpp`)
- **Integration:** `test_main_integration` (`test_main_integration.cpp`)
- **Instruments:** `InstrumentBasic` (`test_instrument_basic.cpp`), `IMSDrift` (`test_ims_drift.cpp`), `TOFFlightTime` (`test_tof_flight_time.cpp`), `LQITStability` (`test_lqit_stability.cpp`), `OrbitrapFrequency` (`test_orbitrap_frequency.cpp`), `DomainTransition` (`test_domain_transition.cpp`), `FTICRCyclotron` (`test_fticr_cyclotron.cpp`), `QuadrupoleFilter` (`test_quadrupole_filter.cpp`)
- **GPU-only tests:** additional CUDA tests are built only with `-DUSE_GPU_ACCEL=ON`; examples include `test_gpu_collision_parity.cpp`, `test_gpu_reaction_parity.cpp`, `test_gpu_boundaries.cpp`, `test_gpu_space_charge.cpp`, and `test_gpu_field_interpolation.cpp`.

## Coverage Snapshot (v1.1.0 branch)

- **Instruments:** IMS, TIMS axial ramp fields, LQIT, Orbitrap, TOF, Quadrupole, FTICR, multi-domain transitions
- **Integrators:** RK4, RK45 adaptive/fixed behavior, Boris, per-ion dt handling, adaptive space-charge stage refresh, boundary actions, RNG/compaction determinism
- **Forces:** Electric and field-model providers, magnetic, damping, gas-flow transport, space charge, registry superposition
- **Collisions:** HSS, EHSS, OU, interaction-potential collisions, EHSS offline sample loading, geometry kernels, velocity sampling, energy diagnostics, temperature scaling, multi-gas cases, CPU SoA parity
- **Reactions:** Stochastic handler, factory setup, multi-gas mixtures, loader/validation rules, CPU SoA parity
- **Space Charge:** Direct/grid model selection, Poisson solver, charge deposition, integration tests, CPU/GPU-model availability behavior
- **Config:** Loaders/validators for configs, species, reactions, ions, field-array terms, database integration
- **I/O:** HDF5 writer regression coverage
- **SoA/Unit:** IonEnsemble, SimulationEngine SoA, OpenMP/SoA benchmarks, domain lookup benchmarks

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

- GPU-specific tests are only built with `-DUSE_GPU_ACCEL=ON` (defines `ICARION_USE_GPU`); some assume a CUDA device and do not skip if none is available. `SpaceChargeGPUModelTest` is present in CPU-only builds because it validates CPU-side availability/fallback behavior.
- EHSS/HSS mobility/thermalization tests have stochastic tolerance; occasional seed sensitivity may appear—rerun if needed.
- Space-charge GPU model is opt-in (`physics.enable_space_charge_gpu`) and falls back automatically; GPU tests validate helper behavior while CPU-only builds exercise the stub.

For full test intent/details, see the individual test files.
