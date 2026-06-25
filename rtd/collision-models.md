# Collision models

Collision models define how ICARION calculates the effect of ion-neutral collisions. They are selected in the `physics` section of the main [configuration](configuration.md), but their behavior also depends on the gas conditions in each domain and the ion properties in the [species database](species-database.md).

In ICARION, the electric and magnetic fields determine the deterministic acceleration of an ion, while the collision model determines how the interaction with the gas modifies that motion. This includes mobility-dependent drift, collision-induced damping or thermalization, diffusion, and Monte Carlo ion-neutral scattering.

```text
ion trajectory = field-driven motion + gas interaction model + optional reactions
```

This page explains the physical meaning of the available models, when to use them, and which input data they require.

---

## Mental model

A simulation step in ICARION can be thought of as three coupled operations:

```text
1. Evaluate fields
   E(r,t), B(r,t), RF/DC/AC waveforms, imported field arrays (and optional collisional damping)

2. Advance the ion trajectory
   RK4, RK45, or Boris

3. Apply gas physics
   stochastic thermal kicks or explicit ion-neutral collisions
```

Stochastic collision models control step 3 by applying explicit gas-event updates after trajectory propagation. Deterministic collision models instead contribute effective drag forces during field/force evaluation. This is why the same instrument geometry can behave very differently depending on whether it is simulated as vacuum, continuum mobility drift, or explicit stochastic collisions.

For example, a standalone time-of-flight simulation may often use `NoCollisions`, while an IMS drift simulation normally requires a gas model. In a coupled IMS-MS simulation, domains can use different gas pressures, temperatures, species, fields, and geometries; however, in v1.0.x the selected `physics.collision_model` is global for the whole run. See [Multi-domain simulations](multi-domain.md).

---

## Continuum vs. discrete-collision descriptions

ICARION contains two broad classes of gas interaction models.

### Continuum models

Continuum models do not simulate individual ion-neutral collisions. Instead, the gas is represented by an effective damping relation. These models are computationally efficient and useful when many collisions occur over the timescale of interest.

Typical use cases:

- mobility-dependent drift in high-pressure IMS regions,
- fast baseline calculations,
- deterministic drift time checks.

The main limitation is that the individual collision history is not available. A pure friction model gives a mean drift behavior but not a physically meaningful event-by-event collision sequence and thus no realistic arrival time distributions (only average values).

### Discrete stochastic models

Discrete stochastic models treat ion-neutral collisions as random binary collision events. The probability of a collision depends on gas density, relative velocity, and collision cross section. When a collision occurs, the post-collision velocity is sampled from the chosen scattering model.

Typical use cases:

- low- and intermediate-pressure regimes,
- diffusion and peak broadening,
- trajectory-level collision statistics.

The main cost is runtime: explicit collision models are more expensive than continuum damping.

---

## Model overview

| Model | Physical picture | Typical use | Main limitation |
|---|---|---|---|
| `NoCollisions` | Vacuum trajectory | TOF, Orbitrap, FT-ICR, pure ion optics | No gas damping or thermalization |
| `Friction` | Mobility-based continuum damping | IMS drift when `K0` is known | Deterministic; no explicit diffusion by itself |
| `Friction` + OU | Mobility damping plus thermal kicks | Mean drift with thermal velocity distribution | Still not event-resolved |
| `HSS` | Stochastic hard-sphere collisions | Explicit gas collisions, diffusion | Simplified spherical scattering |
| `EHSS` | Geometry-/orientation-aware hard-sphere scattering | Structure-sensitive stochastic scattering | More expensive; needs geometry, samples, or EHSS CCS data |
| `HSD` | Deterministic hard-sphere damping | Development/research comparisons | Experimental; no proper stochastic diffusion |
| `Langevin` | Polarization-based damping | Research on long-range ion-neutral effects | Experimental; validate before use |

!!! note
    For production simulations, the safest default choices are usually `NoCollisions` for vacuum regions, `Friction` for mobility-based drift when a reliable mobility is known, and `HSS` for stochastic transport with diffusion. `EHSS` is valuable when molecular shape matters, but it requires more careful input data and validation.

---

## Choosing a model

Use `NoCollisions` when the simulated region is intended to be operated in vacuum and gas scattering should not affect the trajectory.

Use `Friction` when the reduced mobility is known and you mainly need the correct mean drift behavior or arrival time. This is often a good first model for IMS drift sanity checks.

Use `Friction` plus an OU thermalization term when deterministic damping should be balanced by thermal velocity fluctuations. The OU term prevents all random motion from simply damping away.

Use `HSS` when individual collision timing, stochastic scattering, diffusion, or pressure-dependent transport matters.

Use `EHSS` when molecular shape, orientation, or geometry-derived scattering is important.

Use `HSD` or `Langevin` only for development, comparison, or research workflows where the model assumptions are explicitly checked.

---

## Required inputs

Most collision models need three classes of input:

| Input | Where it is configured | Used for |
|---|---|---|
| Gas pressure | Domain `environment.pressure_Pa` | Gas number density and collision frequency |
| Gas temperature | Domain `environment.temperature_K` | Neutral velocity distribution and thermalization |
| Gas species | Domain `environment.gas_species` | Gas-dependent CCS, mass, and collision parameters |
| Ion mass and charge | [Species database](species-database.md) | Dynamics and collision kinematics |
| Mobility `mobility_cm2Vs` | [Species database](species-database.md) | Friction model |
| CCS values | [Species database](species-database.md) | HSS collision frequency |
| Geometry or EHSS samples | [Species database](species-database.md) | EHSS scattering |

A minimal gas-containing domain therefore usually needs at least:

```json
{
  "physics": {
    "collision_model": "HSS"
  },
  "domains": [
    {
      "name": "drift_region",
      "instrument": "IMS",
      "environment": {
        "pressure_Pa": 2000.0,
        "temperature_K": 300.0,
        "gas_species": "He"
      }
    }
  ],
  "species_database": "data/species_database_v1.json"
}
```

See [Configuration](configuration.md) for the full JSON structure.

---

## `NoCollisions`

`NoCollisions` disables gas interactions. The ion motion is then determined only by fields, initial conditions, boundaries, and the integrator.

This is appropriate for vacuum ion optics or reduced model systems where gas interactions should intentionally be absent.

Typical examples:

- time-of-flight regions,
- idealized Orbitrap or FT-ICR trajectories,
- RF/DC ion optics in vacuum,
- debugging field maps or integrators.

A simple vacuum configuration can use:

```json
{
  "physics": {
    "collision_model": "NoCollisions"
  }
}
```

Do not use this model for gas-filled IMS regions unless the absence of gas is intentional.

---

## `Friction`: mobility-based continuum damping

The `Friction` model represents the gas by an effective damping force. Instead of sampling individual collisions, ICARION uses a mobility relationship between electric field, ion mobility, and drift velocity.

Conceptually:

```text
known mobility -> damping coefficient -> mobility-based drift
```

A common relation is:

```text
gamma = q / (K m)
```

where `q` is the ion charge, `m` is the ion mass, `K` is the mobility under the local gas conditions, and `gamma` is the damping coefficient.

The reduced mobility is usually specified in the species database:

```json
{
  "species": {
    "H3O+": {
      "mass_amu": 19.02,
      "charge": 1,
      "mobility_cm2Vs": 24.1,
      "CCS_A2": 24.9
    }
  }
}
```

### What the model gives you

The friction model is useful when you want a smooth, mobility-dependent trajectory. It can reproduce mean drift behavior when the mobility input is reliable.

It is especially useful for:

- drift time baselines,
- sanity checks against known mobilities,
- high-pressure regions where many collisions average out and explicit collision models are expensive,
- fast parameter scans where explicit collisions would be too expensive.

### What it does not give you

A pure friction model is not an event-resolved collision model. By itself, it does not generate realistic stochastic diffusion or collisional scattering.

This means that `Friction` can be a good model for peak positions or mean arrival times, but not necessarily for peak widths, arrival time distributions, or microscopic collision statistics.

---

## OU thermalization

A deterministic damping term alone removes random kinetic energy from the ions. Without a balancing stochastic term, random thermal motion would be damped away.

The Ornstein-Uhlenbeck (OU) thermalization term adds stochastic velocity kicks that balance damping and maintain a thermal velocity distribution at the domain temperature.

Conceptually:

```text
friction removes kinetic energy
OU kicks add thermal fluctuations
balance -> Maxwell-Boltzmann-like equilibrium
```

This is useful for deterministic damping models such as Friction when thermal spreading should be represented without explicit binary collisions.

!!! warning
    OU thermalization must be consistent with the damping coefficient used by the deterministic model. If the damping coefficient and OU term do not match, the equilibrium temperature can be wrong.

OU thermalization is not the same as HSS or EHSS. HSS and EHSS already sample thermal neutral velocities and explicit collision events, so an additional OU thermostat would generally double-count thermalization.

---

## `HSS`: stochastic hard-sphere collisions

`HSS` treats ion-neutral interactions as stochastic binary collisions. During each time step, ICARION estimates whether a collision occurs from Poisson sampling using the gas density, the relative ion-neutral velocity, and an effective collision cross section.

A simplified event loop is:

```text
for each ion and time step:
    sample neutral velocity from gas temperature
    compute relative velocity
    compute collision probability
    if a collision occurs:
        sample an isotropic scattering direction
        transform to the center-of-mass frame
        apply elastic hard-sphere scattering
        transform back to the laboratory frame
```

This makes HSS qualitatively different from Friction: the gas influence is not only a smooth damping force, but a sequence of random collision events.

### Required species data

HSS requires a collision cross section. In practice, ICARION can use gas-specific CCS information from the species database:

```json
{
  "species": {
    "H3O+": {
      "mass_amu": 19.02,
      "charge": 1,
      "CCS_A2": 24.9,
      "CCS_HSS": {
        "He": 24.9,
        "N2": 104.02
      }
    }
  }
}
```

Gas-specific values are important because a CCS measured or computed for one gas should not automatically be assumed to apply to another gas. Note, however, that the CCS should also depend on field and temperature. This is currently neglected in ICARION.

### When HSS is a good choice

Use HSS when:

- diffusion or peak broadening matters,
- the timing of individual collisions matters,
- you want to simulate stochastic collision events,
- no reliable molecular geometry or EHSS orientation data are available.

### Limitations

HSS is still a simplified model. It reduces the collision physics to effective hard-sphere scattering and does not explicitly include long-range ion-neutral interactions, inelastic collisions, reactive scattering, or detailed molecular anisotropy.

For many simulations this is a useful trade-off: HSS is more physical than a purely deterministic damping model, but less expensive and less input-dependent than a full geometry-aware scattering model.

---

## `EHSS`: explicit hard-sphere scattering

`EHSS` extends the hard-sphere collision picture by using molecular geometry, orientation information, or precomputed EHSS data. The goal is to represent shape-dependent scattering more explicitly than in the isotropic HSS model.

In simple terms:

```text
HSS:  ion represented by an effective spherical cross section
EHSS: ion shape/orientation can affect collision outcome
```

### Required species data

EHSS needs enough information to determine an orientation- or geometry-dependent collision target. Depending on the available species entry, ICARION may use:

- gas-specific `CCS_EHSS` values,
- precomputed `EHSS_samples_file` data,
- a molecular `geometry_file`,
- or other geometry-derived collision information available in the species database.

Example species entry:

```json
{
  "species": {
    "PentanalH+": {
      "mass_amu": 87.0,
      "charge": 1,
      "mobility_cm2Vs": 10.25,
      "CCS_A2": 53.7,
      "geometry_file": "molecules/PentanalH+.json",
      "EHSS_samples_file": "molecules/precomputed_ccs/PentanalH+_ehss_samples.json"
    }
  }
}
```

### When EHSS is useful

Use EHSS when:

- molecular shape is expected to affect scattering,
- you compare conformers or structure-dependent transport,
- you have DFT/ab initio geometries or precomputed orientation samples,
- HSS is too isotropic for the scientific question.

### Important limitations

EHSS is still a hard-sphere model. It does not automatically include long-range ion-induced dipole interactions, charge-transfer physics, inelasticity, or reactive chemistry. Reactions are handled separately; see [Reactions](reactions.md).

EHSS is also more sensitive to input quality. Poor geometries, inconsistent radii, inconsistent CCS values, or missing gas-specific data can produce misleading results. For production work, compare EHSS results against HSS and/or Friction baselines and run the relevant checks in [Validation](validation.md).

---

## Collision event sampling in HSS and EHSS

HSS and EHSS are collision models based on Monte Carlo event sampling. Instead of applying a continuous damping force at every step, ICARION checks, for each ion and time step, whether a discrete ion-neutral collision event occurs. This event decision is based on the collision rate

```text
k = N * sigma_eff * |v_ion - v_neutral|
```

where `N` is the neutral number density, `sigma_eff` is the effective collision cross section, and `|v_ion - v_neutral|` is the actual relative speed between the ion and a sampled neutral particle.

For a time step `dt`, the expected number of collisions is

```text
lambda = k * dt
```

ICARION treats collision events as a Poisson process and samples the probability that at least one collision occurs during the current time step:

```text
P(n >= 1) = 1 - exp(-lambda)
```

A uniform random number `u` in `[0, 1)` is drawn. If `u < P(n >= 1)`, a collision event is applied. Otherwise, the ion continues without a collision during that step.

!!! note
    In the current event handler, at most one stochastic HSS/EHSS collision is applied per ion and integration step. The probability `1 - exp(-lambda)` is the exact Poisson probability for one or more events, but multiple events within the same step are not individually resolved. For event-resolved simulations, choose the time step such that `lambda << 1` for the relevant pressure, cross section, and relative velocity range.

### Neutral velocity sampling

Before evaluating the collision probability, ICARION samples a neutral velocity from the gas temperature and bulk gas flow. This is important because the collision rate depends on the actual relative velocity, not only on the ion velocity or a thermal average.

Conceptually, the single gas event loop is:

```text
for each ion and time step:
    sample neutral velocity from gas temperature and gas flow
    compute relative velocity
    compute collision rate k = N * sigma_eff * |v_rel|
    compute event probability P = 1 - exp(-k * dt)

    if random_uniform() < P:
        apply one elastic collision event
    else:
        continue without collision
```

This means that even a slow ion in a warm gas can collide because the neutral particle has thermal motion. Conversely, a fast drifting ion experiences an increased collision rate because the relative velocity increases.

### Gas mixtures

For gas mixtures, ICARION evaluates each gas component separately. For component `i`,

```text
k_i = N_i * sigma_i * |v_ion - v_neutral_i|
```

is computed from the component density `N_i`, the component-specific effective cross section `sigma_i`, and a separately sampled neutral velocity for that component.

The total collision rate is

```text
k_tot = sum_i(k_i)
```

The event probability is then

```text
P_tot = 1 - exp(-k_tot * dt)
```

If a collision occurs, the gas component is selected with probability

```text
P(component i | collision) = k_i / k_tot
```

The collision is then performed with the same neutral velocity that was used to compute `k_i`. This keeps the event-rate calculation and the actual collision kinematics consistent.

### HSS collision event

In HSS, the ion-neutral interaction is reduced to an isotropic hard-sphere collision. Once an event has been accepted:

1. The ion and sampled neutral velocities are transformed conceptually into the center-of-mass frame.
2. The relative speed is preserved.
3. A new scattering direction is sampled isotropically.
4. The post-collision ion velocity is transformed back to the laboratory frame.

The event is elastic: kinetic energy and momentum are conserved for the binary collision. HSS therefore gives stochastic velocity randomization, thermalization, diffusion broadening, and pressure-dependent collision statistics, but it does not represent molecular shape explicitly.

### EHSS collision event

EHSS uses the same Poisson event sampling logic as HSS, but the accepted collision is resolved with molecular geometry information.

The effective event rate still has the form

```text
k = N * sigma_eff * |v_rel|
```

but `sigma_eff` can come from geometry-derived or precomputed EHSS information. If orientation samples are available, ICARION can use one of these orientations by randomly drawing an orientation-dependent projected area for the event rate. The same sampled orientation is then used for calculating scattering in the collision step, so that the collision probability and the collision geometry remain consistent.

Once an EHSS event has been accepted:

1. The relative velocity defines the collision axis.
2. A molecular orientation of the ion geometry is selected or sampled.
3. An impact parameter is sampled in the plane perpendicular to the relative velocity.
4. The algorithm checks whether the incoming gas hard sphere intersects the molecular hard-sphere representation of the ion.
5. If a contact is found, the contact normal defines the scattering direction. If multiple contacts are found, the first local contact normal in the sampled ray-sphere test is used, as it represents the first contact with an atom of the ion.
6. The relative velocity is reflected specularly in the center-of-mass frame.
7. The ion velocity is transformed back to the laboratory frame.

Compared with HSS, EHSS therefore keeps the same stochastic event picture but replaces isotropic scattering by geometry-dependent scattering. It is useful when molecular shape, orientation, or structure-dependent collision behavior matters.

### Practical timestep criterion

For stochastic collision models, the time step is not only an integrator setting. It also controls how well individual collision events are resolved.

A useful diagnostic quantity is

```text
lambda = N * sigma_eff * |v_rel| * dt
```

For event-resolved simulations, `lambda` should usually be much smaller than 1 for typical ions and gas conditions. If `lambda` approaches or exceeds 1, a collision becomes likely in almost every step and unresolved multiple-collision events become possible within one integration step. In that regime, reducing `dt` or switching to a continuum model such as Friction may be more appropriate, depending on the scientific question.

A simple rule of thumb is:

```text
lambda << 1      event-level stochastic collisions are well resolved
lambda ~ 0.1     usually acceptable threshold criterion for many trajectory simulations
lambda ~ 1       collisions occur nearly every step; reduce dt if event timing matters
lambda >> 1      multiple collisions per step are likely; event-level interpretation becomes coarse
```

### Consequences for interpretation

The stochastic event models should be interpreted as Monte Carlo realizations of ion-neutral collisions. Individual trajectories are random, but ensemble properties such as thermalization, arrival time distributions, transmission, and collision-induced broadening should converge with sufficient ion count and appropriate time step.

This is different from the Friction model, where the gas acts continuously through an effective mobility-dependent damping term. HSS and EHSS can represent discrete collision timing and scattering, while Friction represents the averaged transport limit.


## CCS, mobility, and gas-specific data

A common source of confusion is the relation between mobility, CCS, and the collision model.

### Mobility

Mobility is mainly used by the Friction model. It tells ICARION how fast the ion should drift in response to an electric field under given gas conditions.

### CCS

CCS values are mainly used by stochastic collision models to determine collision frequencies and effective target sizes. In the species database, a generic `CCS_A2` may be available, but gas-specific maps are preferable when using HSS or EHSS.

### Practical hierarchy

For the HSS model, the preferred data hierarchy is:

```text
gas-specific CCS_HSS[gas]
    -> model/config override if provided
        -> generic CCS_A2 fallback
```

For the EHSS model, the preferred data hierarchy is:

```text
gas-specific CCS_EHSS[gas]
    -> precomputed EHSS orientation samples
        -> molecular geometry
            -> derived or generic CCS fallback with warning/error if data are insufficient
```

This is why the [species database](species-database.md) is not just a list of masses and charges. It is part of the physical model definition.

---

## Pressure regimes

The appropriate collision model depends strongly on pressure and on the physical observable of interest.

At high pressure, many collisions occur during a short time interval. Their detailed timing may average out, and a continuum mobility model can be efficient and sufficient for mean drift.

At lower pressure, collision events become sparse enough that their timing and angular scattering matter. In this regime, HSS or EHSS can capture effects that a smooth mobility model cannot.

In intermediate regimes, comparing multiple models is often useful:

```text
Friction  -> mean mobility reference
HSS       -> stochastic collision and diffusion reference
EHSS      -> shape-sensitive stochastic reference
```

This comparison is especially useful in coupled instruments where one domain behaves like an IMS region and another behaves like a transfer or m/z-analysis region.

---

## Collision models and reactions

Collision models and reaction models are separate layers in ICARION.

The collision model determines how ion velocities and trajectories are modified by the neutral gas. The reaction model determines whether an ion changes species according to configured kinetic channels.

For coupled simulations of transport and reactions, the collision model still matters because it determines how long ions remain in each domain, how trajectories spread, and how often ions sample different gas conditions. The actual reaction definitions are configured separately; see [Reactions](reactions.md).

---

## Output and diagnostics

All collision models affect the trajectory output. The resulting positions, velocities, species indices, domain indices, and metadata are stored in HDF5 output; see [Output files](output-files.md).

Continuum models mainly affect mean trajectories, drift times, and velocity damping. Stochastic models additionally affect velocity distributions, diffusion, peak shapes, and collision-related diagnostics when those are enabled.

For practical post-processing, see [Analysis](analysis.md). For standardized checks of thermalization, transport, reactions, space charge, and instruments, see [Validation](validation.md).

---

## Common mistakes

### Using `NoCollisions` in a gas-filled region

If an IMS or drift region should contain gas, `NoCollisions` is usually not appropriate. The simulated ions will not experience gas damping, diffusion, or scattering.

### Expecting per-domain collision model selection

In v1.0.x, `physics.collision_model` applies to the full simulation. You can vary pressure, temperature, gas species, fields, and geometry by domain, but you cannot set `HSS` for one domain and `NoCollisions` for another in the same config.

### Using `Friction` and expecting collision statistics

`Friction` is not an event-resolved model. It can provide a mobility-based drift baseline, but it does not produce a physical collision sequence.

### Using stochastic models without checking CCS units

CCS values in the species database are given in Å² and converted internally. Make sure the value corresponds to the gas and model you intend to use.

### Running EHSS without sufficient molecular data

EHSS needs geometry, orientation samples, or EHSS-specific CCS information. If these are missing or inconsistent, the model may fail or fall back to less specific data.

### Adding OU thermalization on top of stochastic collisions

HSS and EHSS already include thermal neutral velocities and stochastic collision outcomes. Adding an independent OU thermostat can double-count thermalization unless the code path is explicitly designed for it.

---

## Recommended workflow

For a new instrument or new species, a robust workflow is:

```text
1. Run NoCollisions if the fields or boundaries need debugging.
2. Run Friction if a reliable mobility is available.
3. Run HSS to include stochastic collisions and diffusion.
4. Run EHSS only when geometry- or orientation-dependent scattering is needed.
5. Compare arrival times, losses, velocity distributions, and final species/domain states.
6. Validate the chosen model against a known reference case before production use.
```

This staged approach helps separate field and geometry problems from ion transport problems.
