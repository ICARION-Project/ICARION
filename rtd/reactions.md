# Reactions

ICARION can model stochastic conversion of one ion species into another during a trajectory simulation. Reactions are configured through a reaction database and enabled in the simulation config.

Use reactions when you want to simulate processes such as:

- first-order ion decay or fragmentation,
- pseudo-first-order ion-neutral reactions,
- bimolecular ion-molecule reactions with explicit neutral concentration,
- competing reaction channels,
- simple reaction chains.

---

## Enable reactions in a config

A reaction run needs three things:

1. a species database containing all reactants and products, (see [species database](species-database.md))
2. a reaction database defining possible channels,
3. `physics.enable_reactions = true`.

```json
{
  "species_database": "data/species_database_v1.json",
  "reaction_database": "data/reactions_database_v1.json",
  "physics": {
    "collision_model": "HSS",
    "enable_reactions": true
  }
}
```

The species IDs used in the reaction database must exactly match the IDs in the species database.

---

## Minimal reaction database

A reaction database contains a `reactions` array:

```json
{
  "reactions": [
    {
      "id": "h3o_to_pentanal_first_order",
      "reactant": "H3O+",
      "product": "PentanalH+",
      "rate_constant": 1200.0
    }
  ]
}
```

With no `order` term, the rate constant is interpreted as a first-order rate in `s^-1`.

The probability for one ion to react during one time step is:

```text
P = 1 - exp(-k_eff * dt)
```

For small probabilities this is approximately `P ≈ k_eff * dt`.

The resulting species change can be analyzed from the HDF5 output files. See also [Output files](output-files.md) and [analysis](analysis.md).

---

## Worked example: first-order conversion

### `species_reaction_demo.json`

```json
{
  "species": {
    "H3O+": {
      "name": "Hydronium ion",
      "mass_amu": 19.02,
      "charge": 1,
      "CCS_A2": 24.9,
      "mobility_cm2Vs": 24.1
    },
    "PentanalH+": {
      "name": "Pentanal protonated",
      "mass_amu": 87.0,
      "charge": 1,
      "CCS_A2": 53.7,
      "mobility_cm2Vs": 10.25
    }
  }
}
```

### `reactions_first_order.json`

```json
{
  "reactions": [
    {
      "id": "h3o_to_pentanal_first_order",
      "reactant": "H3O+",
      "product": "PentanalH+",
      "rate_constant": 1200.0,
      "description": "Example first-order conversion"
    }
  ]
}
```

### Simulation config excerpt

The following is an excerpt and not a complete runnable configuration. A full
simulation file also needs at least an `output` section and one `domains` entry;
see [Configuration files](configuration.md).

```json
{
  "species_database": "species_reaction_demo.json",
  "reaction_database": "reactions_first_order.json",
  "simulation": {
    "total_time_s": 0.002,
    "dt_s": 1e-7,
    "integrator": "RK4",
    "write_interval": 100,
    "rng_seed": 7
  },
  "physics": {
    "collision_model": "NoCollisions",
    "enable_reactions": true
  },
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 10000,
        "position": {
          "type": "point",
          "center": [0.0, 0.0, 0.0]
        },
        "velocity": {
          "type": "fixed",
          "value": [0.0, 0.0, 0.0]
        }
      }
    ]
  }
}
```

For a pure first-order conversion, the expected remaining reactant fraction is:

```text
N_reactant(t) / N0 = exp(-k * t)
```

For `k = 1200 s^-1` and `t = 0.002 s`, the expected remaining fraction is about `exp(-2.4) ≈ 0.091`.

---

## Bimolecular / pseudo-first-order reactions

To model an ion reacting with a neutral concentration, add an `order` term.

```json
{
  "reactions": [
    {
      "id": "h3o_plus_pentanal",
      "reactant": "H3O+",
      "product": "PentanalH+",
      "rate_constant": 3.0e-15,
      "order": [
        {
          "species": "Pentanal",
          "exponent": 1,
          "concentration_m3": 1.0e18
        }
      ]
    }
  ]
}
```

Here:

```text
k_eff = k * [Pentanal]
```

If `k = 3.0e-15 m^3/s` and `[Pentanal] = 1.0e18 m^-3`, then:

```text
k_eff = 3000 s^-1
```

The ion then reacts stochastically with that effective first-order rate.

---

## Buffer-gas fallback concentration

For pseudo-first-order reactions using the configured buffer gas density, an order term can use the neutral/buffer-gas fallback convention:

```json
{
  "id": "buffer_gas_channel",
  "reactant": "H3O+",
  "product": "PentanalH+",
  "rate_constant": 1.0e-16,
  "order": [
    {
      "species": "neutral",
      "exponent": 1,
      "concentration_m3": -1
    }
  ]
}
```

In that case ICARION uses the gas density implied by the domain environment. This is useful for generic buffer-gas processes, but explicit concentrations are usually clearer.

---

## Temperature-dependent rates

The reaction database supports constant, Arrhenius, and modified Arrhenius forms.

### Constant

```json
{
  "id": "constant_example",
  "reactant": "H3O+",
  "product": "PentanalH+",
  "rate_constant": 1200.0,
  "rate_model": "Constant"
}
```

### Arrhenius

```json
{
  "id": "arrhenius_example",
  "reactant": "H3O+",
  "product": "PentanalH+",
  "rate_constant": 200000.0,
  "rate_model": "Arrhenius",
  "activation_energy_eV": 0.08
}
```

### Modified Arrhenius

```json
{
  "id": "modified_arrhenius_example",
  "reactant": "H3O+",
  "product": "H3O+_cluster",
  "rate_constant": 2.0e-9,
  "rate_model": "ModifiedArrhenius",
  "temperature_exponent": -0.5,
  "reference_temperature_K": 300.0,
  "activation_energy_eV": 0.0,
  "order": [
    {
      "species": "H2O",
      "exponent": 1,
      "concentration_m3": 1.0e20
    }
  ]
}
```

Make sure all product species, including cluster products such as `H3O+_cluster`, are also present in the species database.

For validation of these rate constant implementations, see [validation](validation.md).

---

## Competing channels

Multiple reactions can share the same reactant:

```json
{
  "reactions": [
    {
      "id": "channel_a",
      "reactant": "H3O+",
      "product": "PentanalH+",
      "rate_constant": 1000.0
    },
    {
      "id": "channel_b",
      "reactant": "H3O+",
      "product": "CaffeineH+",
      "rate_constant": 500.0
    }
  ]
}
```

The total reaction probability is determined from the sum of effective channel rates. If a reaction occurs, the channel is selected according to the relative rates.

For the example above, channel A should occur about twice as often as channel B, ignoring statistical noise.

---

## Validation checklist for reaction configs

Before running a large reaction simulation:

- Run `icarion --dry-run config.json`.
- Check that `physics.enable_reactions` is `true`.
- Check that `reaction_database` points to the intended file.
- Check that every `reactant` and `product` exists in the species database.
- Check rate-constant units: `s^-1` for first-order, `m^3/s` for one concentration term, `m^6/s` for two concentration factors.
- Start with a simple one-channel reaction and compare to `exp(-k_eff t)`.
- Use a fixed `rng_seed` for reproducible debugging.

---

## Common mistakes

| Symptom | Likely cause |
|---|---|
| No species conversion | `enable_reactions` is false or reaction database not loaded. |
| Dry-run fails | Reactant/product missing from species database. |
| Reaction too fast | `rate_constant` entered in cm³/s instead of m³/s, or concentration too high. |
| Reaction too slow | Missing concentration term or wrong order. |
| Product has wrong transport behavior | Product species lacks CCS/mobility values. |
