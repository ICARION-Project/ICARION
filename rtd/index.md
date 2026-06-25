# ICARION

**ICARION** (*Ion Collision And Reaction IntegratiON*) is a modular framework for particle-based ion trajectory simulations in electric fields and neutral gas environments.

ICARION is designed for reproducible simulations of ion transport, stochastic ion-neutral collisions, ion-molecule reactions, and coupled multi-domain ion mobility–mass spectrometry (IMS-MS) instruments.

## What ICARION is for

ICARION can be used to model ion motion in instrument regions such as drift cells, RF quadrupoles, time-of-flight sections, linear ion traps, Orbitraps, and FT-ICR. Its main purpose is to combine deterministic electric or magnetic fields with collisional transport and optional chemistry in one reproducible simulation workflow.

## Core ideas

- **Configuration-driven simulations:** simulations are defined in structured JSON input files.
- **Multi-domain instruments:** different instrument regions can be combined seamlessly in one run.
- **Explicit physics modules:** collision models, reaction handling, force models, and numerical integrators are modular.
- **Reproducible output:** simulation results and metadata are written to structured HDF5 files with relevant metadata to ensure reproducibility.
- **Validation-oriented design:** examples and validation workflows are intended to make physical behavior traceable.

## Start here

New users should begin with:

1. [Installation](installation.md)
2. [Release packages and launcher](release-packages.md)
3. [First simulation](first-simulation.md)
4. [CLI reference](cli-reference.md)
5. [Configuration](configuration.md)
6. [Output files](output-files.md)

For common first-run problems, see [Troubleshooting](troubleshooting.md).

## Current status

ICARION v1.0.1 focuses on the validated CPU execution path. GPU-related helper code may be present in the repository, but GPU execution should be considered experimental unless stated otherwise for a specific release.

For reproducibility checks, see [Validation](validation.md).

!!! note
    This documentation is an initial user-facing manual. It is intended to complement the shorter repository README and the more detailed technical Markdown files already present in the repository.
