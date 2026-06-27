# Architecture

ICARION is organized as a modular simulation framework with a strict separation between [configuration](configuration.md), domain handling, physical models, numerical integration, and [output](output-files.md).

## High-level components

The central simulation loop is coordinated by the simulation engine. The engine delegates most model-specific work to specialized components:

- domain manager,
- field models,
- force registry,
- collision handler,
- reaction handler,
- numerical integrator,
- space charge model,
- and output manager.

## Separation of concerns

The simulation engine should not encode instrument-specific physics directly. Instead, instrument behavior is defined by configuration files and implemented through interchangeable strategy-like modules.

This makes it possible to extend ICARION with new collision models, integrators, force terms, or domain types without rewriting the main simulation loop.

## Configuration-driven execution

At runtime, ICARION reads the JSON configuration, constructs the required domains and physical model handlers, initializes the ion ensemble, and then advances active ions through the simulation loop. For more information, see [configuration](configuration.md).

## Output and metadata

The output manager writes trajectory data and metadata to HDF5 [outputs](output-files.md). The goal is to make simulation results reproducible and inspectable without relying on hidden runtime state.

## Developer notes

For implementation details, refer to the architecture documentation and source files in the repository.
