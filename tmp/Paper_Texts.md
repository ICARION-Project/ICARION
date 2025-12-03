# Design Overview

ICARION is a modular C++17 framework for simulating ion trajectories, ion–neutral collisions, and ion–molecule reactions under user-defined electric field configurations and background gas environments. The software separates physical models, numerical integration, domain descriptions, and optional GPU backends into independent components with clear interfaces. This modular structure enables the flexible combination of interaction models, domain types, field descriptions, and integrators within a single simulation workflow.
The framework follows a multi-domain design, representing an instrument as a sequence of regions with individual geometry, gas properties, and electric field descriptions. Native geometries currently supported by the simulation engine include cylindrical regions and a specialized Orbitrap-type radial geometry; Cartesian box grids exist only in auxiliary space-charge prototypes. Analytical electric field models are fully integrated, while three-dimensional field maps imported from external solvers are supported on the CPU backend.
Stochastic processes rely on separate pseudorandom generators on CPU and GPU. The CPU backend uses std::mt19937_64-based generators for stochastic collisions and chemical reactions, whereas the GPU backend employs XORWOW generators from cuRAND. Random-number streams are not persistent on a per-ion basis, and bitwise deterministic agreement between CPU and GPU backends is not guaranteed.
GPU acceleration is available for a subset of collision and force computations; however, several GPU code paths, including space charge and generalized field interpolation, remain experimental and are disabled by default. The default execution relies on the CPU backend, with optional offloading when a compatible GPU implementation is available for the selected physical models.
Simulation parameters are defined through structured JSON files that describe species, reactions, initial ion ensembles, domains, and global settings. HDF5 output files contain trajectories and selected metadata, but the full input specification is not embedded within the output. To support reproducibility, ICARION writes a configuration snapshot alongside each output file, capturing the exact JSON inputs used for the run.


# Domain and Environment Layer

A domain defines a continuous region with its own origin, axial length, and radial extent. The current engine natively supports cylindrical regions and an Orbitrap-type radial potential; there is no general boundary mask beyond radius/length limits. Domains may represent drift regions, RF quadrupoles, ion traps, time-of-flight sections, or Orbitraps and can be chained to model multi-stage instruments, and the DomainManager handles coordinate transforms, domain membership, and axial transitions.
Each domain specifies its gas environment (species, temperature, pressure, optional flow), which sets neutral density, collision rates, and reaction rates and may change discontinuously between domains. Environment data is cached in a DomainContext and provided to the physical modules during integration.
Electric fields are supplied either analytically (uniform DC, RF quadrupolar fields, Orbitrap potential) or via tabulated 3D potential arrays from external solvers. Tabulated maps are stored in HDF5 and sampled with trilinear gradients on the CPU; GPU use of imported maps is only partially wired. Both field types share a common access interface.

# Species and Ion Clouds 

Ion species are defined in structured JSON files that specify charge state, mass, mobility/CCS, and optional polarizability. For structure-resolved collisions (EHSS), per-species geometry can be provided in separate JSON files containing atom coordinates and radii; if absent, EHSS falls back to the stored CCS. 
Initial ion ensembles are defined inline in the configuration (point sources, uniform or Gaussian clouds, with fixed velocities, thermal velocities, or specified kinetic energies), and multi-species packets are supported. This enables realistic source conditions, gated injection into drift regions, and packet formation typical of IMS–MS workflows.

# Physical Interaction Modules

Ion motion in ICARION is governed by the combined effect of electric forces, deterministic collision forces, stochastic collisions, and chemical reactions. Physical modules operate on a local DomainContext object to ensure domain-specific fields and gas properties are used consistently. 

## Electric forces

Electric forces follow directly from the field modules associated with the current domain. Time-dependent fields are applied by scaling analytical or tabulated potentials with user-defined waveforms.

## Deterministic collision forces

Deterministic collision models approximate the average drag from frequent collisions. Three formulations are available: hard-sphere damping (geometric cross sections), a Langevin-type term (polarization), and mobility-based friction from literature mobilities. Hard-sphere and Langevin are experimental; the friction model is the recommended default. Thermal kicks are not applied here (they live in the stochastic pipeline).

## Stochastic collisions

Stochastic collisions represent discrete momentum-transfer events. ICARION provides an isotropic hard-sphere Monte Carlo model (HSS) and a structure-resolved EHSS model that samples geometry-dependent scattering based on molecular atom coordinates/radii. Both conserve momentum and kinetic energy.

## Chemical reactions

ICARION supports unimolecular and bimolecular reactions using Arrhenius and modified Arrhenius expressions, with optional higher-order terms. Reactions are evaluated alongside stochastic collisions; product species properties are read from the species database and velocities are assigned via momentum conservation.

# Numerical Integration

Time propagation is performed on the CPU using fixed-step fourth-order Runge–Kutta (RK4) or adaptive Runge–Kutta–Fehlberg (RK45) schemes; a Boris pusher is available for magnetic-field-dominated cases. Optional modules include OU thermalization (disabled by default) and space-charge forces: a direct O(N²) Coulomb term and a box-grid Poisson solver on the CPU, with an experimental P³M helper on the GPU. At each step the engine applies deterministic forces from the current domain, invokes stochastic collision/reaction handlers, applies optional OU kicks, and then performs boundary/domain checks. GPU integrators and field interpolation exist in experimental form but are not wired by default.

# Simulation Engine

The simulation engine coordinates all physical modules during time integration. At each timestep it determines the active domain, gathers the DomainContext, applies collision/reaction modules using the backend RNGs (std::mt19937_64 on CPU; XORWOW on GPU), evaluates forces, advances the state with the chosen integrator, checks boundaries/transitions, and writes outputs when scheduled. GPU offload is optional and limited to selected experimental paths; CPU and GPU results are not guaranteed bitwise-identical.

# Parallelisation

ICARION parallelises ion trajectories across CPU threads using OpenMP and offers CUDA acceleration for selected modules. A structure-of-arrays (SoA) layout improves CPU cache behaviour and GPU coalescing. RNG streams are not persisted per ion, and deterministic agreement across backends is not guaranteed.

# Input/Output and Reproducibility
ICARION uses a structured JSON input system (config, species, reactions, ions, optional geometry). Outputs are written to HDF5 with trajectories and selected metadata; the full inputs are not embedded. A separate configuration snapshot is written alongside each run to aid reproducibility. Persistent RNG seeds and Git metadata are recorded where available; exact bitwise replay across CPU/GPU is not claimed.
