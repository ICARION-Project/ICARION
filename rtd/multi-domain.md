# Multi-domain simulations

A central feature of ICARION is the ability to combine multiple instrument domains in one simulation. 

## What is a domain?

A domain is a spatial region with its own geometry, fields, gas environment, pressure, temperature, and boundary behavior.

Examples of domains include:

- an IMS drift cell,
- an RF quadrupole mass filter or ion guide,
- a time-of-flight acceleration region,
- a linear ion trap,
- an Orbitrap,
- or an FT-ICR.

## Why multi-domain simulations?

Many IMS-MS instruments are composed of several regions with very different physical conditions. A simulation may need to describe a collisional drift region followed by a lower-pressure RF ion guide or mass analyzer.

In ICARION, ions can move from one domain into the next without manually exporting and reinitializing the ion ensemble.

## Domain transitions

Domain transitions are handled by the domain manager. The simulation resolves the active domain of each ion and applies the corresponding field model, environment, and boundary conditions.

!!! note
    In v1.0.x, `physics.collision_model` is global for the full run, not domain-local. Multi-domain simulations can vary pressure, temperature, gas species, fields, and geometry per domain. They cannot use `HSS` in one domain and `NoCollisions` in another domain within the same run.

## Example concept

Multi-domain instruments are defined through the `domains` array in the [configuration](configuration.md).

A simple coupled IMS-quadrupole setup may use:

1. an IMS domain with helium gas and a uniform axial electric field,
2. followed by a quadrupole mass filter domain with RF/DC fields and much lower pressure.

This allows one simulation to describe both ion mobility separation and mass filtering.
