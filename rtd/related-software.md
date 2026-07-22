# Related and complementary software

This list is not exhaustive. The packages solve overlapping but often distinct problems; inclusion does not imply numerical equivalence or direct feature parity. Specialized tools should be used as validation references where appropriate.

## Collision cross sections and ion mobility

### MOBCAL and MobCal-MPI

MOBCAL established widely used trajectory method workflows for structure-derived collision cross sections and ion mobilities. The original MobCal-MPI publication describes the parallelized package; MobCal-MPI 2.0 is the later field-dependent extension. These are CCS/mobility tools rather than complete time-resolved instrument simulators.

### IMoS and IMoS 2

The Ion Mobility Suite (IMoS) covers momentum-transfer and trajectory approaches for CCS and mobility, including rotating diatomic gases. IMoS 2 formulates trajectory-resolved generalized Boltzmann transport with rotational degrees of freedom and arbitrary fields.

### Collidoscope

Collidoscope is a trajectory method CCS calculator with parallel and optimized sampling. Its structural CCS output has a different purpose from ICARION's persistent event-level momentum-transfer tables.

### CoSIMS

CoSIMS is an optimized trajectory-based CCS simulator for ion mobility applications. It is another methodological reference for molecular ion-neutral scattering, not an instrument trajectory framework.

## Ion trajectory and instrument simulation

### IDSimF

IDSimF is an open-source framework for molecular ion dynamics in mass spectrometry and ion mobility spectrometry, including fields, collision models, space charge, and reactions.

### ITSIM

ITSIM is an ion trajectory simulation package developed around quadrupole ion traps. ITSIM 6.0 added multiparticle simulations with arbitrary electrode geometries and externally calculated fields.

### SIMION

SIMION is a general-purpose commercial environment for electrostatic field and charged particle trajectory calculations in user-defined electrode geometries. See the [official SIMION user manual](https://simion.com/info/manual.html).

## Scope of ICARION

CCS and mobility calculators evaluate molecular transport quantities. Instrument frameworks propagate ions through time-dependent geometries and fields, often with gas collisions, reactions, and ensemble effects. ICARION belongs primarily to the latter category; its IPM workflow uses offline molecular scattering to create persistent stochastic momentum-transfer tables consumed inside time-resolved instrument simulations. This does not make the classical trajectory method or offline scattering calculation an ICARION novelty.

See [Collision models](collision-models.md) and [Interaction-potential precomputation](ipm-precomputation.md).

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
