# Related and Complementary Software

This list is not exhaustive. The packages solve overlapping but often distinct problems; inclusion does not imply numerical equivalence or direct feature parity. Specialized tools should be used as validation references where appropriate.

ICARION is part of a broader ecosystem of software for ion–neutral collision calculations, ion mobility prediction, ion optics, and numerical simulation of ion trajectories. The projects below have overlapping but generally distinct scopes. They are listed to document the methodological context of ICARION and to help users select the appropriate tool for a given problem.

## Collision cross sections and ion mobility

### MOBCAL and MobCal-MPI

MOBCAL established the trajectory method as a widely used approach for calculating ion mobility collision cross sections from molecular structures and explicit ion–neutral interaction potentials. The original MobCal-MPI publication describes the parallelized package; MobCal-MPI 2.0 is the later field-dependent extension.

These programs primarily calculate collision integrals, collision cross sections, and mobilities for molecular structures. ICARION uses related trajectory calculations to generate stochastic collision event tables that can subsequently be applied during time-resolved simulations of complete instruments.

### Collidoscope

Collidoscope is a parallel trajectory method collision cross section calculator with explicit He and N₂ collision gas models and optimized trajectory sampling.

Its trajectory sampling and interaction potential calculations are methodologically related to ICARION’s offline interaction potential precomputation. Collidoscope primarily reports structural collision cross sections, whereas ICARION preserves event cross sections and momentum-transfer distributions for later stochastic runtime sampling.

### Ion Mobility Suite

The Ion Mobility Suite (IMoS) provides methods for calculating ion mobilities and collision cross sections using detailed scattering models. IMoS 2 extends this framework to trajectory-resolved solutions of generalized Boltzmann transport equations for polyatomic ions in arbitrary electric fields.

IMoS therefore provides an important reference for ion transport theory and field-dependent mobility calculations. ICARION instead focuses on coupling collision models to explicit ion trajectories within user-defined IMS and MS instrument geometries.

### CoSIMS

CoSIMS is an optimized trajectory-based CCS simulator for ion mobility applications. It is a methodological reference for molecular ion-neutral scattering rather than a complete instrument trajectory framework.

## Ion trajectory and instrument simulation

### IDSimF

IDSimF is an open-source framework for simulations of nonrelativistic molecular ion dynamics in mass spectrometry and ion mobility spectrometry. It supports modular trajectory solvers, electric fields, ion–neutral collision models, space charge interactions, and chemical reaction kinetics.

IDSimF is among the closest existing projects to ICARION at the overall framework level. Both projects target modular ion dynamics simulations rather than only a single instrument or collision cross section calculation. Their architectures, configuration systems, numerical models, and intended workflows nevertheless differ, and results should not be assumed to be interchangeable without model-specific validation.

### ITSIM

ITSIM is a historically important ion trajectory simulation package developed primarily for quadrupole ion trap mass spectrometry. Later versions supported multiparticle simulations and arbitrary electrode geometries using externally calculated electric field arrays.

ITSIM represents an important predecessor for instrument-specific ion trap trajectory simulation. ICARION generalizes this type of workflow to multiple IMS and MS device classes, modular collision and reaction models, and reproducible HDF5-based simulation outputs.

### SIMION

SIMION is a widely used general-purpose platform for calculating electrostatic fields and simulating charged particle trajectories in user-defined electrode geometries. It has been extensively used in ion optics and mass spectrometer development.

See the [official SIMION user manual](https://simion.com/info/manual.html) for its documented scope and workflows.

## Scope of ICARION

ICARION is not intended to replace every specialized collision cross section, mobility, field solver, or ion optics program. Its primary goal is to provide a common, reproducible framework in which:

* electric and magnetic fields,
* molecular and effective ion–neutral collision models,
* stochastic reactions,
* space charge,
* multiple connected instrument domains, and
* structured simulation and provenance data

can be combined within one time-resolved ion trajectory simulation workflow.

Where applicable, ICARION should be validated against specialized software and analytical reference solutions rather than being treated as an independent source of ground truth.

CCS and mobility calculators evaluate molecular transport quantities; instrument frameworks propagate ions through time-dependent geometries and fields. ICARION belongs primarily to the latter category and consumes offline IPM scattering tables during instrument simulations. See [Collision models](COLLISION_MODELS.md), [IPM tool documentation](../tools/README.md#interaction_potential_precompute), and the [public IPM guide](../rtd/ipm-precomputation.md).

## Primary references

- M. F. Mesleh, J. M. Hunter, A. A. Shvartsburg, G. C. Schatz, and M. F. Jarrold, “Structural Information from Ion Mobility Measurements: Effects of the Long-Range Potential,” *J. Phys. Chem.* **100** (1996), 16082–16086. DOI: [10.1021/jp961623v](https://doi.org/10.1021/jp961623v).
- C. Ieritano, J. Crouse, J. L. Campbell, and W. S. Hopkins, “A Parallelized Molecular Collision Cross Section Package with Optimized Accuracy and Efficiency,” *Analyst* **144** (2019), 1660–1670. DOI: [10.1039/C8AN02150C](https://doi.org/10.1039/C8AN02150C).
- A. Haack, C. Ieritano, and W. S. Hopkins, “MobCal-MPI 2.0: an accurate and parallelized package for calculating field-dependent collision cross sections and ion mobilities,” *Analyst* **148** (2023), 3257–3273. DOI: [10.1039/D3AN00545C](https://doi.org/10.1039/D3AN00545C).
- S. A. Ewing, M. T. Donor, J. W. Wilson, and J. S. Prell, “Collidoscope: An Improved Tool for Computing Collisional Cross-Sections with the Trajectory Method,” *J. Am. Soc. Mass Spectrom.* **28** (2017), 587–596. DOI: [10.1007/s13361-017-1594-2](https://doi.org/10.1007/s13361-017-1594-2).
- V. Shrivastav, M. Nahin, C. J. Hogan Jr., and C. Larriba-Andaluz, “Benchmark Comparison for a Multi-Processing Ion Mobility Calculator in the Free Molecular Regime,” *J. Am. Soc. Mass Spectrom.* **28** (2017), 1540–1551. DOI: [10.1007/s13361-017-1661-8](https://doi.org/10.1007/s13361-017-1661-8).
- C. Larriba-Andaluz and C. J. Hogan Jr., “Collision cross section calculations for polyatomic ions considering rotating diatomic/linear gas molecules,” *J. Chem. Phys.* **141** (2014), 194107. DOI: [10.1063/1.4901890](https://doi.org/10.1063/1.4901890).
- C. Larriba-Andaluz, “Ion Mobility Suite (IMoS) 2: A General Computational Framework for Solving the Generalized Boltzmann Equation with Rotational Transport for Polyatomic Ions in Arbitrary Fields,” *Anal. Chem.* **98** (2026), 15730–15740. DOI: [10.1021/acs.analchem.6c01278](https://doi.org/10.1021/acs.analchem.6c01278).
- C. A. Myers, R. J. D'Esposito, D. Fabris, S. V. Ranganathan, and A. A. Chen, “CoSIMS: An Optimized Trajectory-Based Collision Simulator for Ion Mobility Spectrometry,” *J. Phys. Chem. B* **123** (2019), 4347–4357. DOI: [10.1021/acs.jpcb.9b01018](https://doi.org/10.1021/acs.jpcb.9b01018).
- M. Rajkovic, S. Benter, M. Hammelrath, M. Thinius, T. Benter, and W. Wißdorf, “IDSimF: An Open-Source Framework for the Simulation of Molecular Ion Dynamics in Mass Spectrometry and Ion Mobility Spectrometry,” *J. Am. Soc. Mass Spectrom.* **35** (2024), 1451-1460. DOI: [10.1021/jasms.4c00054](https://doi.org/10.1021/jasms.4c00054).
- G. Wu, R. G. Cooks, Z. Ouyang, M. Yu, W. J. Chappell, and W. R. Plass, “Ion trajectory simulation for electrode configurations with arbitrary geometries,” *J. Am. Soc. Mass Spectrom.* **17** (2006), 1216–1228. DOI: [10.1016/j.jasms.2006.05.004](https://doi.org/10.1016/j.jasms.2006.05.004).
