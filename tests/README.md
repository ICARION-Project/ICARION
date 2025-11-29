# ICARION Test Suite

Comprehensive test suite for ICARION v1.0. Tests are organized by component type and functionality.

## Test Structure

```
tests/
├── instruments/     # Instrument-specific physics validation
├── physics/         # Core physics components
│   ├── collisions/  # Collision models and thermalization
│   ├── forces/      # Force calculations
│   ├── reactions/   # Chemical reactions
│   └── spacecharge/ # Space charge effects
├── integrator/      # Numerical integration strategies
├── config/          # Configuration loading and validation
├── io/              # HDF5 file I/O
├── integration/     # End-to-end integration tests
├── utils/           # Utility functions
└── debug/           # Debug/troubleshooting tests
```

---

## Test Classes by Category

### 1. Instruments (`instruments/`)

**Purpose:** Validate instrument-specific physics and mathematical models

#### `test_instrument_basic.cpp` - Fast Sanity Checks
- IMS domain instantiation and field setup
- Orbitrap domain instantiation and k-factor
- LQIT domain instantiation and Mathieu parameters
- Quadrupole domain instantiation and stability
- TOF domain instantiation and drift region
- FTICR domain instantiation and B-field
- All instruments have valid collision-free environments

#### `test_ims_drift.cpp` - Ion Mobility Spectrometry Physics
- Electric field acceleration (no collisions): Basic E-field acceleration in vacuum
- HSS collisions do not crash: Stability check with stochastic collisions
- Field scaling (no collisions): Higher E-field produces faster drift
- Mobility measurement with HSS collisions: Drift velocity vs Mason-Schamp equation (±25% accuracy)
- Mobility measurement with EHSS collisions: Geometry-dependent collision model (±50% accuracy)
- Friction mobility (DampingForce only): Deterministic damping model (0% error, PRODUCTION READY)
- Langevin mobility (DampingForce only): Polarization-based damping (EXPERIMENTAL, polar molecules only)
- HSD mobility (DampingForce only): Hard-sphere deterministic model (EXPERIMENTAL)

#### `test_lqit_stability.cpp` - Linear Quadrupole Ion Trap
- Stable trapping (low q parameter): Ion confined in stable Mathieu parameter region
- Unstable region (high q parameter): Ion escapes at unstable parameters
- Mass-selective stability: Different masses have different stability regions

#### `test_tof_flight_time.cpp` - Time-of-Flight
- Ion accelerates through drift tube: Basic acceleration and drift
- Heavier ions arrive later: Mass-dependent flight time
- Flight time scales with sqrt(m): Kinematic relationship validation

#### `test_orbitrap_frequency.cpp` - Orbitrap
- Axial oscillation frequency: Frequency measurement from trajectory
- Frequency scales with sqrt(1/m): Mass-dependent frequency scaling
- Mass-dependent detection: Frequency resolution

#### `test_fticr_cyclotron.cpp` - Fourier Transform Ion Cyclotron Resonance
- Cyclotron frequency: ω_c = qB/m validation
- Magnetic field strength: B-field effects on cyclotron motion

#### `test_quadrupole_filter.cpp` - Quadrupole Mass Filter
- Mass-selective filtering: Stable vs unstable trajectories
- Mathieu stability diagram: a-q parameter space validation

#### `test_domain_transition.cpp` - Multi-Domain Transitions
- Ion crosses from IMS to TOF: Seamless domain boundary crossing
- Aperture blocks off-axis ions: Geometric acceptance
- Properties update between domains: Field and environment changes

---

### 2. Physics (`physics/`)

#### 2.1 Collisions (`physics/collisions/`)

**Purpose:** Validate collision models, energy conservation, thermalization

##### `test_geometry_utils.cpp` - Molecular Geometry Tools
- convert_molecule_to_geometry: Basic conversion, single atom, empty molecule
- load_geometry_map: Single file mode, directory mode, missing species handling
- Physical constants verification

##### `test_collision_energy_conservation.cpp` - Energy Conservation
- Hard-sphere collision conserves total energy
- EHSS collision with geometry conserves energy
- Center-of-mass frame energy conservation
- Thermalization to correct equilibrium validates collision physics

##### `test_temperature_scaling.cpp` - Temperature Dependence
- Temperature scaling: Does equilibrium scale with T_gas?
- Low energy start: Heating to equilibrium

##### `test_hss_collision_handler.cpp` - Hard Sphere Stochastic (HSS)
- Thermalization of H3O+: Room temperature equilibrium
- Thermalization from high energy: Cooling to equilibrium
- Isotropic velocity distribution: Angular isotropy check

##### `test_ehss_collision_handler.cpp` - Enhanced HSS (Geometry-Dependent)
- Thermalization of H3O+: With molecular geometry
- Thermalization from high energy: Geometry effects on cooling

##### `test_ou_collision_handler.cpp` - Ornstein-Uhlenbeck Thermalization
- Thermalization of H3O+: OU process equilibrium
- Thermalization from high energy: OU process cooling
- Isotropic velocity distribution: OU isotropy
- Fluctuation-Dissipation balance: Statistical mechanics validation

##### `test_collision_factory.cpp` - Collision Handler Factory
- EHSS handler creation, requires geometry
- HSS handler creation
- OU handler with Friction, requires positive gamma
- Deterministic without OU returns null
- HSD without OU returns null
- NoCollisions returns null
- Unknown model throws
- OU with HSD is valid
- Logging enabled works

##### `test_multi_gas_collision.cpp` - Gas Mixture Handling
- EHSS throws on missing geometry
- EHSS uses CCS_EHSS map in mixture, missing CCS falls back to geometry
- HSS uses gas-specific CCS map in mixture
- HSS throws when mixture has no sigma
- HSS mixture thermalization proxy via collision counts

##### `test_collision_energy_diagnostic.cpp` - Energy Diagnostics
- (Detailed energy tracking for debugging)

##### `test_realistic_transfer.cpp` - Realistic Transfer Conditions
- (Complex multi-stage simulations)

##### `test_thermalization_convergence.cpp` - Thermalization Convergence
- (Convergence rate analysis)

#### 2.2 Forces (`physics/forces/`)

**Purpose:** Validate force calculation accuracy

##### `test_electric_field_force.cpp` - Electric Field Force
- Constructor validation, field provider mode
- IMS (Ion Mobility Spectrometry): Uniform axial field
- TOF (Time-of-Flight): Field-free drift
- LQIT (Linear Quadrupole Ion Trap): Quadrupole potential
- FTICR: Trapping potential
- Orbitrap: Logarithmic radial potential with axial harmonic
- QuadrupoleRF: RF quadrupole fields
- Edge cases

##### `test_magnetic_damping_forces.cpp` - Magnetic and Damping Forces
- MagneticFieldForce: Constructor validation, Lorentz force F=q(v×B), cyclotron motion, linear gradient field, disabled force
- DampingForce: Friction model (mobility-based), HardSphere model (deterministic), no damping
- Combined Magnetic + Damping

##### `test_force_registry.cpp` - Force Registry (Superposition)
- Empty registry, single force, multiple forces
- Conditional force
- Realistic physics (gravity)
- Clear functionality
- Null force handling
- Force context (SSOT)

##### `test_space_charge_force.cpp` - Space Charge (Coulomb Interactions)
- SpaceChargeDirect: Name, applies to all ions, negative softening throws
- No ions returns zero force, null ion ensemble, single ion
- Self-interaction excluded
- Two positive charges repel along X, opposite charges attract
- Newton's third law
- Three-dimensional repulsion, three ions linear
- Softening prevents infinite force, converges to Coulomb at large distances
- Overlapping ions return zero force
- Neutral ion produces no force
- Time independent

##### `test_force_integration.cpp` - Force Integration Tests
- Electric + Magnetic forces combine
- All forces (Electric + Magnetic + Damping + SpaceCharge)
- Force superposition principle
- Performance with 100 ions

#### 2.3 Reactions (`physics/reactions/`)

**Purpose:** Validate stochastic chemical reaction handling

##### `test_reaction_factory.cpp` - Reaction Handler Factory
- Create StochasticReactionHandler
- Reactions disabled
- Enable logging
- SSOT compliance

##### `test_multi_gas_reaction.cpp` - Multi-Gas Reactions
- StochasticReactionHandler uses mixture partial densities

##### `test_stochastic_reaction_handler.cpp` - Stochastic Reaction Handler
- SSOT compliance
- Second-order reaction, third-order reaction
- Buffer gas fallback
- Species database lookup
- No applicable reactions
- Reaction statistics
- Handles competing channels correctly
- Handles zero reactions gracefully
- Handles very large k_eff, very small k_eff
- Early exit for k_total < 1e-60
- Numerical safety for k*dt > 50
- Arrhenius temperature dependence
- Modified Arrhenius (capture)

#### 2.4 Space Charge (`physics/spacecharge/`)

**Purpose:** Validate Poisson solver and charge deposition

##### `test_poisson_solver.cpp` - Poisson Equation Solver
- Point charge in vacuum: Analytical comparison
- Uniform charged sphere: Analytical comparison
- Charge conservation
- Solver method comparison: Performance analysis
- Field computation from potential: Gradient calculation

##### `test_charge_deposition.cpp` - Charge Deposition
- Charge conservation
- Single ion at grid point: NGP (Nearest Grid Point) method
- CIC method validation: Cloud-In-Cell interpolation
- OpenMP thread safety: Parallel deposition
- Inactive ions excluded
- Boundary handling
- CIC smoothness validation
- CIC trilinear weights
- CIC vs analytical for uniform distribution

##### `test_space_charge_integration.cpp` - Space Charge Integration
- Two-ion Coulomb repulsion (Direct vs Grid)
- E-field smoothness near ion
- Force symmetry

---

### 3. Integrator (`integrator/`)

**Purpose:** Validate numerical integration accuracy and stability

#### `test_rk4_strategy.cpp` - 4th Order Runge-Kutta
- Basic properties
- Free fall (constant acceleration): Kinematic accuracy
- Harmonic oscillator (periodic motion): Energy conservation
- Exponential decay (damping): Exponential convergence
- Factory creation
- SSOT compliance

#### `test_rk45_strategy.cpp` - Adaptive RK45 (Dormand-Prince)
- Constructor and configuration
- Free fall with adaptive timestep
- Harmonic oscillator with adaptive timestep
- Exponential decay with error control
- Convergence order verification: 5th order accuracy
- Step rejection with tight tolerance
- FSAL optimization: First Same As Last
- Statistics collection

#### `test_boris_strategy.cpp` - Boris Integrator (Magnetic Fields)
- Constructor and properties
- Factory creation
- Constant electric field acceleration
- Cyclotron motion in uniform B-field
- E×B drift motion
- Long-term energy conservation: Symplectic property
- Agreement with RK4 for E-field only

#### `test_domain_manager.cpp` - Domain Management
- (Domain boundary handling, transitions)

#### `test_output_manager.cpp` - Output Management
- (Trajectory output scheduling)

#### `test_simulation_engine.cpp` - Simulation Engine
- (High-level simulation orchestration)

---

### 4. Configuration (`config/`)

**Purpose:** Validate JSON config loading, parsing, and validation

#### `test_config_loader.cpp` - Main Config Loader
- Loads minimal valid config
- Loads config with two domains
- Handles file not found
- Requires dt_s and total_time_s

#### `test_field_array_terms_loader.cpp` - Field Array Terms
- Load single constant term, DC_Axial term, RF term with phase
- Load multiple terms
- Default values
- All scale types
- Error: missing file, unknown scale_type
- Empty array is valid, no field_array_terms key is valid

#### `test_field_array_e2e.cpp` - Field Array End-to-End
- Load dc_axial_unit.h5 and verify structure
- Interpolate dc_axial_unit at center
- Verify DC_Axial scaling logic, Constant scaling, RF scaling (time-dependent)

#### `test_reaction_loader_unit.cpp` - Reaction Loader
- Load simple reaction without order terms, with order terms
- Multiple reactions
- Species validation with database
- Invalid species reference throws
- Missing required fields throws, file not found throws
- Effective rate calculation, multiple order terms
- ReactionDatabase: Get reactions for species

#### `test_reaction_validation.cpp` - Reaction Validation
- Reaction validation rules
- Dimensional consistency warnings

#### `test_ion_loader.cpp` - Ion Distribution Loader
- Point position distribution
- Gaussian position distribution
- Thermal velocity always random
- Multiple species with different boundaries

#### `test_species_loader_unit.cpp` - Species Database Loader
- Load valid ion species, neutral species
- Load multiple species
- Missing required field throws, invalid mass throws
- File not found throws, invalid JSON throws
- Unit conversions are correct
- SpeciesDatabase: Dictionary lookup

#### `test_waveform_loader.cpp` - Waveform Loader
- Load constant, linear (with defaults), quadratic, sinusoidal waveforms
- Load pulsed waveform, arbitrary waveform (with step interpolation)
- Load waveform library
- ValueOrWaveform: static value, inline waveform, waveform reference
- Missing required field throws, unknown type throws
- Invalid time range throws
- Arbitrary with mismatched arrays throws, unsorted times throws
- Invalid waveform reference throws

#### `test_waveform_types.cpp` - Waveform Types
- ConstantWaveform, LinearWaveform (before start, during ramp, after end clamped/unclamped, negative ramp)
- QuadraticWaveform (before start, during interval, after end)
- SinusoidalWaveform (at key phases, with phase offset)
- PulsedWaveform (before, during, after pulse)
- ArbitraryWaveform (linear interpolation, step interpolation, out of bounds)
- Waveform variant: evaluate via std::visit
- ValueOrWaveform: constant value, inline waveform, waveform reference
- Invalid cases (no option, multiple options, missing reference)
- Empty arrays throw, mismatched array sizes throw

---

### 5. I/O (`io/`)

**Purpose:** Validate HDF5 trajectory file writing

#### `test_hdf5_writer_v2.cpp` - HDF5 Writer v2
- HDF5 library is available and functional
- Creates correct file structure
- Writes simulation metadata correctly
- Writes reproducibility metadata (SHA256 hash, git commit)
- Writes species metadata in tabular format
- Writes domain configuration
- Writes ion initial conditions
- Appends trajectory data correctly
- Finalization writes completion metadata
- Handles empty reaction database
- Writes reaction metadata correctly
- Writes multiple domains correctly
- Writes system metadata correctly
- Handles large ion ensembles (performance)
- SHA256 hashing with ConfigLoader integration

#### `test_hdf5_writer.cpp` - HDF5 Writer (Legacy)
- (Same tests as v2, older version)

---

### 6. Integration (`integration/`)

**Purpose:** End-to-end simulation pipeline tests

#### `test_main_integration.cpp` - Main Integration
- Complete simulation pipeline: Config loading → initialization → simulation → output
- Error handling: Invalid configs, missing files

---

### 7. Utils (`utils/`)

**Purpose:** Utility function tests

#### `test_hash.cpp` - SHA256 Hashing
- SHA256 hashing works correctly
- SHA256 error handling
- SHA256 integration with HDF5Writer

---

### 8. Debug (`debug/`)

**Purpose:** Bug reproduction and troubleshooting

#### `test_orbitrap_engine_bug.cpp` - Orbitrap Bug Reproducer
- SimulationEngine: Orbitrap bug reproducer

#### `debug_orbitrap_engine.cpp` - Debug Executable
- (Standalone debugging tool)

#### `debug_hss_first_step.cpp` - HSS Debug
- (HSS collision first-step debugging)

---

## Running Tests

### Run all tests:
```bash
cd build
ctest --output-on-failure
```

### Run specific test class:
```bash
./build/tests/instruments/test_ims_drift
./build/tests/physics/forces/test_electric_field_force
```

### Run with tags (Catch2):
```bash
# Fast tests only
./build/tests/instruments/test_instrument_basic "[fast]"

# Physics validation
./build/tests/instruments/test_ims_drift "[physics]"

# Skip slow tests
ctest -E "slow|performance"
```

### Known Test Behavior

**Expected Failures (marked `[!mayfail]`):**
- `test_ims_drift.cpp`: HSS/EHSS mobility tests (RNG seed dependency in CTest)
- `test_ims_drift.cpp`: Langevin/HSD mobility (EXPERIMENTAL models)
- `test_ou_collision_handler.cpp`: OU thermalization (statistical convergence)
- `test_hss_collision_handler.cpp`: HSS thermalization (statistical convergence)

**Run manually:**
```bash
./build/tests/instruments/test_ims_drift "[!mayfail]"
./build/tests/physics/collisions/test_ou_collision_handler "[!mayfail]"
```

---

## Test Coverage

- **Instruments:** 7 instrument types, 30+ physics validation tests
- **Collisions:** 4 collision models (HSS, EHSS, OU, deterministic damping)
- **Forces:** Electric, magnetic, damping, space charge
- **Integrators:** RK4, RK45 (adaptive), Boris (symplectic)
- **Configuration:** Complete JSON schema validation
- **I/O:** HDF5 output schema v2.0
- **Energy Conservation:** All collision models validated
- **Thermalization:** Temperature scaling, equilibrium convergence

---

## References

- **Architecture:** `docs/ARCHITECTURE.md`
- **Collision Models:** `docs/COLLISION_MODELS.md`
- **Troubleshooting:** `docs/TROUBLESHOOTING.md`
- **Developer Guide:** `docs/DEVELOPERS_GUIDE.md`

