# ICARION Release Notes

## v1.0.0 (December 2025)

**Highlights**
- Core SimulationEngine (AoS/SoA) with multi-domain support and OutputManager
- GPU acceleration: RK4/RK45/Boris batch integrators; HSS/EHSS collision helper (active-ion threshold default 5000)
- SoA (`IonEnsemble`) path and integrator `process_timestep` hooks
- HDF5 output with selected config/system metadata, species/reaction subsets, waveform library (v1.1 addition)
- Comprehensive CLI (`icarion_main`): validation, schema dumps, profiling/benchmark flags
- Validated collision models: Friction (deterministic friction force), HSS (isotropic collisions), EHSS (structure-aware, mobility bias noted)
- Numerical safety logging/reporting, OpenMP-aware threading defaults

**Known Limitations**
- GPU space-charge P³M helper exists but is not yet dispatched from the main loop (CPU direct summation only)
- GPU boundary helper (absorption/cylindrical) not wired into timestep loop
- EHSS mobility overestimates drift speed (~23% in simple H3O+/N2 test); thermalization is correct
- Full JSON config is written as a snapshot file next to output; HDF5 still stores selected fields only

**Breaking/Behavior Notes**
- Domain configs use `environment`, `geometry.origin_m`, and lowercase field keys (`dc`, `rf`, `ac`)
- GPU dispatch thresholds: integrators 5000 ions (Boris ~2500), collisions 5000 active ions by default

**Upgrade Tips**
- Validate configs with `--validate-config` before runs
- For EHSS, ensure species include geometry or gas-specific CCS; otherwise fallback/throws
- For GPU runs, build with `-DUSE_GPU_ACCEL=ON` and set `enable_gpu: true` in config

**Future (planned)**
- Wire GPU space charge and boundary helpers into the main loop
- Broader reaction type persistence in HDF5
- FieldSolver/Optimizer modules (scaffolded)
