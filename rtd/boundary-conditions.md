# Boundary conditions

Boundary conditions define what happens when an ion reaches a domain wall. They
are configured per domain under `domains[].boundary`.

```json
"boundary": {
  "type": "Absorption"
}
```

If no boundary block is given, ICARION defaults to absorption.

## Supported types

| Type | Alias | Effect |
| --- | --- | --- |
| `Absorption` | `absorption` | Removes the ion from the simulation. |
| `SpecularReflection` | `specular` | Mirrors the velocity at the boundary normal. |
| `DiffuseReflection` | `diffuse` | Emits into a randomized direction with partial thermalization. |
| `ThermalReflection` | `thermal` | Re-emits from a Maxwell-Boltzmann distribution at wall temperature. |

Use the canonical names in new configurations. The lowercase aliases are
accepted for convenience.

## Reflection parameters

Diffuse and thermal reflections use:

```json
"boundary": {
  "type": "DiffuseReflection",
  "accommodation_coeff": 1.0,
  "temperature_K": 300.0
}
```

- `accommodation_coeff` is constrained to `0.0 ... 1.0`.
- `0.0` behaves more elastically.
- `1.0` means full thermal accommodation.
- `temperature_K` is the wall temperature.
- If `temperature_K` is omitted, the domain environment temperature is used.

`accommodation_coeff` is relevant for diffuse reflection. For pure absorption it
has no practical effect.

## Output interpretation

Absorbed ions are no longer active. Their removal is visible in trajectory and
ion status output, including elimination timing when that dataset is
written. Reflection keeps the ion alive, but changes its velocity at the wall.

For detector domains, absorption is usually the clearest model. For wall 
collision studies, reflection can be useful, but it changes residence 
times and energy distributions. Treat those results as part of the
physical model, not as a numerical detail.

## Common mistakes

- Starting ions outside the domain causes immediate boundary handling.
- A timestep that crosses through a narrow domain can make losses look sudden.
- Reflective boundaries can hide geometry mistakes because ions remain alive.
- Wall temperature and gas temperature are independent if `temperature_K` is set
  explicitly.

When debugging unexpected losses, first run with `--validate-config`, inspect
the initial ion cloud, and reduce `simulation.dt_s`.
