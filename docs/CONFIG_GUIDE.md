# ICARION Configuration Guide

This guide explains how to create and structure JSON configuration files for ICARION simulations.

## Table of Contents

- [Quick Start](#quick-start)
- [Configuration Structure](#configuration-structure)
- [Required Fields](#required-fields)
- [Schema Validation](#schema-validation)
- [Helper Tools](#helper-tools)
- [Examples](#examples)

---

## Quick Start

### Minimal Configuration

```json
{
  "simulation": {
    "total_time_s": 1e-3,
    "dt_s": 1e-9,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "NoCollisions"
  },
  "output": {
    "folder": "./results",
    "trajectory_file": "trajectories.h5"
  },
  "domains": [
    {
      "name": "my_domain",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.05,
        "radius_m": 0.01
      },
      "environment": {
        "pressure_Pa": 101325.0,
        "temperature_K": 300.0,
        "gas_species": "He"
      }
    }
  ]
}
```

### Creating a Config with the Helper Tool

```bash
# Interactive mode
python3 scripts/create_config.py

# Command-line mode
python3 scripts/create_config.py \
  --template ims \
  --output my_config.json \
  --time 1e-3 \
  --timestep 1e-9
```

---

## Configuration Structure

ICARION uses a hierarchical JSON configuration with four main sections:

```json
{
  "simulation": { ... },    // Time parameters, integrator, execution mode
  "physics": { ... },       // Collision model, reactions, space charge
  "output": { ... },        // Output paths and options
  "domains": [ ... ]        // Array of instrument domains
}
```

### Optional Top-Level Fields

```json
{
  "title": "My Simulation",                          // Human-readable description
  "species_database": "data/species_database_v1.json",  // Species properties database (optional)
  "reaction_database": "data/reactions_database_v1.json", // Reaction rates database (optional)
  "ion_cloud": "data/initial_cloud.h5"               // Initial ion distribution
}
```

**Note on Database Paths:**

- Both `species_database` and `reaction_database` are **optional**
- If not specified, ICARION will automatically search for global fallback databases:
  - `data/species_database_v1.json` (global species database)
  - `data/reactions_database_v1.json` (global reactions database)
- The fallback search starts from the config file's directory and searches up to 5 levels
- A log message `[INFO] No [species/reaction] database specified, using global fallback: ...` indicates fallback usage
- If no database is found (neither specified nor fallback), databases remain empty (which is valid for simulations without reactions)

### Gas-specific CCS maps (HSS/EHSS)

You can store gas-dependent CCS values (generated via `ccs_precompute`):
```json
\"species\": {
  \"H3O+\": {
    \"mass_amu\": 19.0,
    \"charge\": 1,
    \"CCS_reference_gas\": \"He\",
    \"CCS_model\": \"HSS\",
    \"CCS_HSS\": { \"He\": 110.0, \"N2\": 130.0, \"O2\": 140.0 },   // Å²
    \"CCS_EHSS\": { \"He\": 120.0 }  // optional, Å² (nutzt geometry_file wenn nicht vorhanden)
  }
}
```

- HSS: uses σ per gas from `CCS_HSS[gas]`, else mixture override `cross_section_m2`, else `ion.CCS_m2`.
- EHSS: uses `CCS_EHSS[gas]` if present, else orientation samples (if provided), else geometry, else (without geometry) throws.
- EHSS offline runs can also reference a precomputed single-gas sample table with `EHSS_offline_samples_file`.
- InteractionPotentialModel runs reference offline tables with the canonical `ipm_samples_file` key.
- Tool: `./ccs_precompute --input species.json --output out.json --species H3O+ --ref-gas He --ref-ccs-A2 110.0 [--model HSS|EHSS] [--override] [--n-orientations 300]`.
- Orientation samples tool: `./ehss_samples_precompute --input species.json --output h3o_samples.json --species H3O+ [--n-orientations 300] [--n-samples 8000]`.

---

## Species and Reaction Databases (v1.1.0)

ICARION supports external databases for species properties and reaction rates.

### Species Database

Species databases define physical properties of ions and neutrals in user-friendly units:

```json
{
  "species": {
    "H3O+": {
      "name": "Hydronium ion",
      "mass_amu": 19.02,
      "charge": 1,
      "mobility_cm2Vs": 2.8,
      "CCS_A2": 11.0,
      "reference_temperature_K": 300.0,
      "reference_pressure_Pa": 101325.0,
      "ccs_method": "trajectory"
    },
    "H2O": {
      "name": "Water",
      "mass_amu": 18.015,
      "charge": 0,
      "polarizability_A3": 1.47
    }
  }
}
```

**Required fields:**

- `mass_amu`: Molecular mass in atomic mass units
- `charge`: Charge state (integer: -2, -1, 0, +1, +2, ...)

**Optional fields (ions):**

- `mobility_cm2Vs`: Reduced ion mobility at STP [cm²/(V·s)]
- `CCS_A2`: Collision cross-section [Ų]

**Optional fields (neutrals):**

- `polarizability_A3`: Electric polarizability [ų]

**Optional metadata:**

- `name`: Human-readable name
- `geometry_file`: Path to molecular geometry (e.g., .xyz, .pdb, required for EHSS)
- `EHSS_samples_file`: Path to precomputed EHSS orientation samples (JSON)
- `reference_temperature_K`: Temperature for mobility/CCS measurements (just for reproducibility)
- `reference_pressure_Pa`: Pressure for mobility measurements (just for reproducibility)
- `ccs_method`: Method used to determine CCS (`"trajectory"`, `"projection"`, etc., for reproducibility)

**Unit conversions:** ICARION automatically converts to SI internally:

- `mass_amu` → kg (× 1.66054e-27)
- `charge` → C (× 1.60218e-19)
- `mobility_cm2Vs` → m²/(V·s) (× 1e-4)
- `CCS_A2` → m² (× 1e-20)
- `polarizability_A3` → m³ (× 1e-30)

### Reaction Database

Reaction databases define ion-molecule reactions with **temperature-dependent** and concentration-dependent rates:

```json
{
  "reactions": [
    {
      "id": "rxn_001_constant",
      "reactant": "H3O+",
      "product": "H5O2+",
      "rate_constant": 3.5e-15,
      "rate_model": "Constant",
      "order": [
        {
          "species": "H2O",
          "exponent": 1,
          "concentration_m3": 2.5e25
        }
      ],
      "description": "Constant rate (no T-dependence)"
    },
    {
      "id": "rxn_002_arrhenius",
      "reactant": "H3O+",
      "product": "NH4+",
      "rate_constant": 1.5e-15,
      "rate_model": "Arrhenius",
      "activation_energy_eV": 0.12,
      "order": [
        {
          "species": "NH3",
          "exponent": 1
        }
      ],
      "description": "Proton transfer with activation barrier",
      "reference": "Doe et al., Int. J. Mass Spectrom. (2021)"
    },
    {
      "id": "rxn_003_capture",
      "reactant": "H3O+",
      "product": "H3O+·H2O",
      "rate_constant": 2.0e-15,
      "rate_model": "ModifiedArrhenius",
      "temperature_exponent": -0.5,
      "reference_temperature_K": 300.0,
      "activation_energy_eV": 0.0,
      "order": [
        {
          "species": "H2O",
          "exponent": 1
        }
      ],
      "description": "Ion-dipole capture (T^-0.5, anti-Arrhenius)"
    }
  ]
}
```

**Required fields:**

- `id`: Unique reaction identifier
- `reactant`: Species ID (must exist in species database)
- `product`: Species ID (must exist in species database)
- `rate_constant`: Base rate constant (k₀ or A). **Units depend on reaction order** (see below)

**Optional fields (Temperature Dependence):**

- `rate_model`: Temperature dependence model (default: `"Constant"`)
  - **`"Constant"`**: k(T) = k₀ (no T-dependence)
  - **`"Arrhenius"`**: k(T) = A × exp(-Eₐ / (kB·T))
  - **`"ModifiedArrhenius"`**: k(T) = A × (T/T₀)ⁿ × exp(-Eₐ / (kB·T))
- `activation_energy_eV`: Activation energy Eₐ in eV (Arrhenius/ModifiedArrhenius)
- `temperature_exponent`: Temperature exponent n (ModifiedArrhenius only)
- `reference_temperature_K`: Reference temperature T₀ in Kelvin (ModifiedArrhenius only, default: 300 K)

**Optional fields (Concentration Dependence):**

- `order`: Array of concentration-dependent terms
  - `species`: Species ID for concentration dependence, `"neutral"` for non-equilibrium buffer gas fallback, or `"M"` for a third-body placeholder
  - `exponent`: Concentration exponent (allowed: **0, 1, or 2**)
  - `concentration_m3`: Fixed concentration [m⁻³] or **-1** for buffer gas fallback
- `description`: Human-readable description
- `reference`: Literature reference

---

### Supported Reaction Cases & Validation Rules

ICARION supports the following reaction types with **strict validation**:

| **Case** | **JSON Configuration** | **Meaning** | **Example** |
|----------|------------------------|-------------|-------------|
| **1. Pseudo-first-order** | `"species": "neutral"` <br> `"exponent": 1` <br> `"concentration_m3": -1` | Use buffer gas density `n_gas` from simulation config | H₃O⁺ + N₂ → products <br> (k in [m³/s], rate = k·n_gas) |
| **2. Explicit concentration** | `"species": "O2"` <br> `"exponent": 1` <br> `"concentration_m3": 2.5e25` | User-defined fixed concentration | H₃O⁺ + O₂ → H₃O⁺·O₂ <br> (k in [m³/s], n_O₂ = 2.5×10²⁵ m⁻³) |
| **3. Bimolecular (Ion + X)** | `"species": "O2"` <br> `"exponent": 1` | Concentration from species_db or runtime | H₃O⁺ + O₂ → products <br> (typical 2-body ion-molecule) |
| **4. Termolecular (3-body)** | `"species": "H2O"` <br> `"exponent": 2` | Quadratic concentration dependence | H₃O⁺ + 2 H₂O → H₃O⁺·(H₂O)₂ <br> (k in [m⁶/s]) |
| **5. Autocatalytic** | `"species": "H3O+"` <br> `"exponent": 1` | Product species appears in order term | Ion → Ion collision chains <br> (mathematical only, rare) |

**IMPORTANT: Order-Dependent Units**

The `rate_constant` field has **dimensions that depend on reaction order**:
- **Order 0** (spontaneous decay, no order terms): k has units [s⁻¹]
- **Order 1** (2nd-order, one neutral with exponent=1): k has units [m³/s]
- **Order 2** (3rd-order, exponent=2 or two exponent=1 terms): k has units [m⁶/s]

**Always specify k with the correct dimensions for your reaction order!**

---

### Validation Rules (Enforced at Load Time)

ICARION **strictly validates** all order terms:

| **Rule** | **Check** | **Error Example** |
|----------|-----------|-------------------|
| **#1: Exponent range** | `exponent ∈ {0, 1, 2}` | ERROR: `"exponent": 3` → "exponent must be 0, 1, or 2" |
| **#2: Concentration range** | `concentration_m3 ≥ -1.0` | ERROR: `"concentration_m3": -5.0` → "must be ≥ -1.0" |
| **#3: Species exists** | `species ∈ species_db` (if not `"neutral"` or `"M"`) | ERROR: `"species": "XYZ"` → "species 'XYZ' not found" |
| **#4: No duplicate species** | Each `species` appears once | ERROR: Two terms with `"species": "O2"` → "duplicate order term for 'O2'. Use exponent=2 instead." |
| **#5: Max one buffer gas** | At most one term with `concentration_m3 = -1` | ERROR: Two terms with `-1` → "only one term can use buffer gas fallback" |

**Dimensional Consistency Warnings:**

ICARION also warns about likely unit mismatches:

```
⚠  Reaction 'rxn_001': 2nd-order (exponent=1) but k = 1.5e-30 m⁶/s outside typical range [1e-18, 1e-12] m³/s
```

**Typical Rate Constant Ranges:**
- **1st-order (spontaneous):** k ~ 10⁻³ to 10⁶ s⁻¹
- **2nd-order (ion + neutral):** k ~ 10⁻¹⁸ to 10⁻¹² m³/s
- **3rd-order (termolecular):** k ~ 10⁻⁴⁴ to 10⁻³⁸ m⁶/s

**Effective rate calculation:**

The effective rate includes **both temperature and concentration dependencies**:
```
k_eff(T) = k(T) × [S₁]^n₁ × [S₂]^n₂ × ...

where k(T) depends on rate_model:
  Constant:         k(T) = k₀
  Arrhenius:        k(T) = A × exp(-Eₐ/(kB·T))
  ModifiedArrhenius: k(T) = A × (T/T₀)ⁿ × exp(-Eₐ/(kB·T))
```

**Examples:**

- **Arrhenius (rate increases with T):** Proton transfer H₃O⁺ + NH₃ → NH₄⁺ with Eₐ = 0.12 eV
  - At 300 K: k = 1.8×10⁻¹¹ m³/s
  - At 400 K: k = 3.5×10⁻¹¹ m³/s (2× faster!)
  
- **Modified Arrhenius (anti-Arrhenius, rate decreases with T):** Ion-dipole capture H₃O⁺ + H₂O → H₃O⁺·H₂O with n = -0.5
  - At 200 K: k = 2.45×10⁻⁹ m³/s (faster at low T!)
  - At 400 K: k = 1.73×10⁻⁹ m³/s (slower at high T)

### Using Databases in Configs

**Enable species database:**

```json
{
  "species_database": "data/species_database_v1.json",
  "domains": [ ... ]
}
```

**Enable reactions:**

```json
{
  "species_database": "data/species_database_v1.json",
  "reaction_database": "data/reactions_database_v1.json",
  "physics": {
    "enable_reactions": true
  }
}
```

When `enable_reactions` is enabled, ICARION will:

1. **Load the reaction database** from the specified file
2. **Handle reactions stochastically** using Monte Carlo method during simulation
3. **Compute effective rates** at each time step:

   ```text
   k_eff [s⁻¹] = k(T) × ∏ᵢ [Xᵢ]^nᵢ
   ```

   where k(T) is the temperature-dependent rate constant and [Xᵢ] are concentration terms

4. **Calculate reaction probability** for each particle per time step:

   ```text
   P = 1 - exp(-k_eff × dt)
   ```

5. **Select channels** if multiple reactions compete (probability-weighted)
6. **Update species** by replacing reactant particle with product particle

Equilibrium-linked reactions can mark a forward channel with `equilibrium=true` and SI thermochemistry metadata (`delta_r_H_J_mol`, `delta_r_S_J_molK`). The loader then creates a dynamic reverse channel unless one is already explicit; see `schema/reactions.schema.json` for the full schema.

**Equilibrium scope:** v1.1 equilibrium support is a bath kinetics population model for single ion species conversion,
`A+ + nX -> B+ + products`, with optional third-body placeholder `M`. Neutral
products are not explicit state variables. The thermodynamic factor therefore
assumes that concentration order terms represent the relevant stoichiometric
neutral partners. This is appropriate for simple cluster association,
dissociation, and isomerization-style population balances, but not for arbitrary
exchange reactions such as `A+ + X <-> B+ + Y` where the activity or
concentration of `Y` must appear in the reverse direction. Empirical falloff
kinetics or reduced mechanisms where kinetic order and elementary stoichiometry
differ should use explicit forward/reverse rates instead of `equilibrium=true`.
The equilibrium constant is treated as a dimensionless `K_p` with standard
pressure `p0 = 1 bar`; `delta_r_H_J_mol` and `delta_r_S_J_molK` are interpreted
as temperature-independent standard reaction enthalpy and entropy over the
simulated temperature range.

For equilibrium reactions, do not use `species: "neutral"` as a real neutral
partner. The loader rejects that combination because the equilibrium expression
does not know which gas species activity should be included. Use an explicit
neutral species such as `"N2"` or `"H2O"`, or use `"M"` only for a third-body
placeholder that cancels from the population equilibrium.

Reaction events update species, mass, charge, CCS, and mobility, but keep the
ion position and velocity unchanged. Thus `equilibrium=true` enforces chemical
population ratios, not microscopic detailed balance in full phase space; it does
not model reaction enthalpy release, neutral recoil, or product energy
partitioning.

See [Reaction Database Schema](#reaction-database-schema) for details on configuring reactions.

**Path resolution:** Paths can be:

- Absolute: `/absolute/path/to/database.json`
- Relative to config file: `../data/species.json`
- Relative to working directory: `data/species.json`

**Schema validation:** Use provided JSON schemas to validate databases:

```bash
python3 schema/validator.py schema/species.schema.json data/species_database_v1.json
python3 schema/validator.py schema/reactions.schema.json data/reactions_database_v1.json
```

**Example configs:**

- `examples/ims/` - IMS examples
- `examples/reactions/reaction_demo.json` - Ion-molecule reactions

---

## Ion Initialization 

ICARION provides flexible ion initialization with per-species position boundaries and velocity distributions.

### Configuration Syntax

Ion configuration can be specified in two ways:

#### Option 1: Generate from Species List

Each species can have its own spatial region (boundaries) and velocity distribution:

```json
{
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 1000,
        "position": {
          "type": "gaussian",
          "center": [0.0, 0.0, 0.0],
          "std": [0.001, 0.001, 0.001]
        },
        "velocity": {
          "type": "thermal",
          "temperature_K": 300
        }
      }
    ]
  }
}
```

#### Option 2: Load from Ion Cloud File

```json
{
  "ions": {
    "from_file": "ion_clouds/my_cloud.json"
  }
}
```

### Position Distributions

Each species can be created in its own spatial region (boundaries):

#### Point (Single Location)

All ions at same position:

```json
"position": {
  "type": "point",
  "center": [0.0, 0.0, 0.0]  // [m]
}
```

#### Gaussian (Normal Distribution)

3D Gaussian distribution with independent standard deviations per axis:

```json
"position": {
  "type": "gaussian",
  "center": [0.0, 0.0, 0.0],         // Mean position [m]
  "std": [0.001, 0.002, 0.003]       // Std dev per axis [m]
}
```

#### Uniform Sphere

Uniform distribution within a sphere:

```json
"position": {
  "type": "uniform_sphere",
  "center": [0.0, 0.0, 0.0],   // Sphere center [m]
  "radius": 0.005              // Sphere radius [m]
}
```

#### Uniform Cylinder

Uniform distribution within a cylinder (aligned with z-axis):

```json
"position": {
  "type": "uniform_cylinder",
  "center": [0.0, 0.0, 0.0],   // Cylinder center [m]
  "radius": 0.002,             // Radial extent [m]
  "length": 0.01               // Axial extent (z-direction) [m]
}
```

#### Uniform Box

Uniform distribution within a rectangular box:

```json
"position": {
  "type": "uniform_box",
  "min": [-0.001, -0.001, 0.0],    // Min corner [m]
  "max": [0.001, 0.001, 0.01]      // Max corner [m]
}
```

### Velocity Distributions

#### Fixed (All Ions Same Velocity)

```json
"velocity": {
  "type": "fixed",
  "value": [0.0, 0.0, 100.0]  // Velocity [m/s]
}
```

#### Thermal (Maxwell-Boltzmann)

Random directions sampled from Maxwell-Boltzmann distribution. **No directed drift** (always random):

```json
"velocity": {
  "type": "thermal",
  "temperature_K": 300  // Temperature [K]
}
```

**Note:** For thermal distribution, velocity directions are **always random**. Use `kinetic` for directed velocities.

#### Kinetic (Directed with Energy)

Fixed kinetic energy in a specified direction, with optional angular spread:

```json
"velocity": {
  "type": "kinetic",
  "energy_eV": 0.1,              // Kinetic energy [eV]
  "direction": [0, 0, 1],        // Direction vector (normalized automatically)
  "spread_angle_deg": 10         // Optional: cone angle spread [degrees]
}
```

#### Gaussian

Gaussian distribution for each velocity component:

```json
"velocity": {
  "type": "gaussian",
  "mean": [0.0, 0.0, 100.0],     // Mean velocity [m/s]
  "std": [10.0, 10.0, 20.0]      // Std dev per axis [m/s]
}
```

### Multi-Species Example

Different species in different spatial regions with different velocities:

```json
{
  "species_database": "data/species_database_v1.json",
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 500,
        "position": {
          "type": "gaussian",
          "center": [0.0, 0.0, 0.0],
          "std": [0.001, 0.001, 0.001]
        },
        "velocity": {
          "type": "thermal",
          "temperature_K": 300
        }
      },
      {
        "id": "H5O2+",
        "count": 500,
        "position": {
          "type": "uniform_sphere",
          "center": [0.005, 0.0, 0.0],
          "radius": 0.002
        },
        "velocity": {
          "type": "thermal",
          "temperature_K": 350
        }
      },
      {
        "id": "O2-",
        "count": 100,
        "position": {
          "type": "uniform_cylinder",
          "center": [0.0, 0.0, 0.01],
          "radius": 0.001,
          "length": 0.005
        },
        "velocity": {
          "type": "kinetic",
          "energy_eV": 0.1,
          "direction": [0, 0, 1],
          "spread_angle_deg": 10
        }
      }
    ]
  }
}
```

### Ion Cloud File Format

For loading ions from file:

```json
{
  "ions": [
    {
      "species": "H3O+",
      "pos": [0.0, 0.0, 0.0],      // Position [m]
      "vel": [0.0, 0.0, 100.0],    // Velocity [m/s]
      "birth_time": 0.0            // Birth time [s]
    },
    {
      "species": "H5O2+",
      "pos": [0.001, 0.0, 0.0],
      "vel": [10.0, 0.0, 50.0],
      "birth_time": 0.0
    }
  ]
}
```

**Note:** Species must exist in species database. Properties (mass, charge, mobility, CCS) are looked up automatically.

### Validation

Ion configuration is validated during loading:

- Species IDs must exist in species database
- Position/velocity parameters must be physically reasonable
- Warnings for missing species or invalid configurations

Validation errors are collected in `ValidationResult` (not thrown as exceptions).

**Example validation output:**

```
✓ Generated 1100 ions from configuration
  - 500 H3O+ ions (gaussian position, thermal velocity)
  - 500 H5O2+ ions (sphere position, thermal velocity)
  - 100 O2- ions (cylinder position, kinetic velocity)
```

**Example configs:**

- `examples/ion_clouds/default_cloud.json` - Ion cloud file format example

---

## Reaction Database Schema

When `enable_reactions` is set in the `physics` section, you must provide a `reaction_database` JSON file defining all possible reactions.

### Schema Structure

```json
{
  "reactions": [
    {
      "id": "rxn_proton_transfer",
      "reactant": "H3O+",
      "product": "NH4+",
      "rate_model": "Constant",
      "rate_constant": 2.5e-15,
      "order": [
        {
          "species": "NH3",
          "exponent": 1,
          "concentration_m3": -1.0
        }
      ]
    }
  ]
}
```

### Field Descriptions

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | YES | Unique identifier for the reaction |
| `reactant` | string | YES | Species ID of the reactant (must exist in species database) |
| `product` | string | YES | Species ID of the product (must exist in species database) |
| `rate_model` | enum | YES | Temperature model: `"Constant"`, `"Arrhenius"`, `"ModifiedArrhenius"` |
| `rate_constant` | number | YES | Base rate constant (units depend on reaction order, see below) |
| `activation_energy_eV` | number | Conditional | Required if `rate_model` is `"Arrhenius"` or `"ModifiedArrhenius"` |
| `temperature_exponent` | number | Conditional | Required if `rate_model` is `"ModifiedArrhenius"` |
| `reference_temperature_K` | number | Conditional | Required if `rate_model` is `"ModifiedArrhenius"` |
| `order` | array | NO | Concentration terms (if omitted, 1st-order reaction) |

### Order Term Structure

Each entry in the `order` array defines a concentration dependence:

```json
{
  "species": "H2O",
  "exponent": 1,
  "concentration_m3": 2.5e25
}
```

| Field | Type | Values | Description |
|-------|------|--------|-------------|
| `species` | string | Species ID, `"neutral"`, or `"M"` | Which species or placeholder to use for concentration |
| `exponent` | integer | 0, 1, 2 | Power to raise concentration to |
| `concentration_m3` | number | ≥ -1.0 | Explicit concentration [m⁻³] or -1.0 for buffer gas |

**Buffer Gas Fallback:**

- If `concentration_m3 = -1.0` (or omitted), use buffer gas density from `EnvironmentConfig.particle_density_m_3`
- If `concentration_m3 > 0`, use the explicit value
- `species: "neutral"` is allowed for non-equilibrium pseudo-first-order buffer gas kinetics only.
- `species: "M"` is a third-body placeholder. It can be used in equilibrium-linked reactions when the third body cancels from the population equilibrium.

### Rate Constant Units (Order-Dependent)

**Important:** The `rate_constant` field has **order-dependent units**:

| Reaction Order | Rate Constant Unit | Effective Rate Unit | Example |
|----------------|-------------------|---------------------|---------|
| 0th (spontaneous) | [s⁻¹] | [s⁻¹] | Unimolecular decay |
| 2nd (bimolecular) | [m³/s] | [s⁻¹] | A⁺ + X → B⁺ |
| 3rd (termolecular) | [m⁶/s] | [s⁻¹] | A⁺ + X + M → B⁺ |

**Effective Rate Formula:**

```text
k_eff [s⁻¹] = k₀(T) [m³ⁿ⁻³/s] × ∏ᵢ [Xᵢ]^nᵢ
```

Where n = total order = Σ exponents.

### Temperature Models

**1. Constant (no temperature dependence):**

```text
k(T) = k₀
```

**2. Arrhenius:**

```text
k(T) = A × exp(-Eₐ/(kB·T))
```

Where:

- A = `rate_constant` [m³ⁿ⁻³/s]
- Eₐ = `activation_energy_eV` [eV]
- kB = Boltzmann constant

**3. ModifiedArrhenius:**

```text
k(T) = A × (T/T₀)^n × exp(-Eₐ/(kB·T))
```

Where:

- A = `rate_constant` [m³ⁿ⁻³/s]
- T₀ = `reference_temperature_K` [K]
- n = `temperature_exponent` [dimensionless]
- Eₐ = `activation_energy_eV` [eV]

### Examples

**Example 1: 2nd-order reaction with Constant model:**

```json
{
  "id": "rxn_proton_transfer",
  "reactant": "H3O+",
  "product": "NH4+",
  "rate_model": "Constant",
  "rate_constant": 2.5e-15,
  "order": [
    {
      "species": "NH3",
      "exponent": 1,
      "concentration_m3": -1.0
    }
  ]
}
```

**Example 2: 3rd-order reaction with Arrhenius model:**

```json
{
  "id": "rxn_clustering",
  "reactant": "H3O+",
  "product": "H5O2+",
  "rate_model": "Arrhenius",
  "rate_constant": 1.2e-40,
  "activation_energy_eV": 0.05,
  "order": [
    {
      "species": "H2O",
      "exponent": 1,
      "concentration_m3": 2.5e25
    },
    {
      "species": "He",
      "exponent": 1,
      "concentration_m3": -1.0
    }
  ]
}
```

**Example 3: 1st-order unimolecular decay (ModifiedArrhenius):**

```json
{
  "id": "rxn_fragmentation",
  "reactant": "C6H6+",
  "product": "C4H4+",
  "rate_model": "ModifiedArrhenius",
  "rate_constant": 1e5,
  "activation_energy_eV": 1.5,
  "temperature_exponent": 0.5,
  "reference_temperature_K": 300.0
}
```

### Validation

Validate your reaction database with:

```bash
python3 schema/validator.py schema/reactions.schema.json my_reactions.json
```


---

## Configuration Reference

This section provides a quick reference for all configuration sections. For complete field specifications and validation rules, see the JSON schema files in `schema/`.

### Top-Level Structure

```json
{
  "simulation": { /* timing, integrator, GPU settings */ },
  "physics": { /* collision model, reactions, space charge */ },
  "output": { /* HDF5 output settings */ },
  "ions": [ /* initial ion distributions */ ],
  "domains": [ /* simulation regions with fields and boundaries */ ]
}
```

**Schema Files:**
- Complete specification: `schema/icarion-config.schema.json`
- Individual sections: `schema/simulation.schema.json`, `schema/physics.schema.json`, etc.
- Validation: `python3 schema/validator.py schema/icarion-config.schema.json your_config.json`

### Simulation Section

**Required:** `total_time_s`, `dt_s`

```json
"simulation": {
  "total_time_s": 1e-3,
  "dt_s": 1e-9,
  "integrator": "RK4",           // RK4, RK45, Boris
  "write_interval": 100,
  "enable_gpu": false,           // Runtime GPU path falls back to CPU unless CUDA support is available
  "rng_seed": 42
}
```

See `schema/simulation.schema.json` for all options.

**Integrator caveat:** Adaptive RK45 can run with space charge if `ICARION_ADAPTIVE_SC` is not set to `0`. Fields are rebuilt at every RK stage (expensive). Set `ICARION_ADAPTIVE_SC=0` to force the legacy guard (single SC update per macro-step) or use fixed-step RK4 for faster SC runs.

### Physics Section

**Required:** `collision_model`

```json
"physics": {
  "collision_model": "HSS",      // NoCollisions, HSD, HSS, EHSS, Langevin, Friction, InteractionPotentialModel
  "enable_reactions": false,
  "enable_space_charge": false,
  "enable_space_charge_gpu": false, // Prefer GPU P³M (requires CUDA build; falls back automatically)
  "space_charge_model": "auto",     // auto, direct, grid, gpu
  "collision_subcycles_per_step": 1,
  "ipm_orientation_mode": "random"
}
```

**Collision Models:** See [COLLISION_MODELS.md](COLLISION_MODELS.md) for detailed physics and use cases.

**High-pressure stochastic collisions:** `collision_subcycles_per_step > 1` splits each collision application into equal micro-steps and recomputes collision probabilities in each sub-step. `collision_multi_event_mode=true` is a practical approximation for regimes where more than one collision per macro-step is likely; it enforces at least `collision_max_events_per_step` micro-subcycles. Despite the legacy field name, `collision_max_events_per_step` is not a guaranteed upper bound on physical continuous-time collision events. Validate both `lambda * dt_collision` for collision statistics and the global `dt_s` for trajectory accuracy; subcycling does not make a large RK/global step valid if fields, gradients, momentum relaxation, walls, or domain boundaries are under-resolved.

**InteractionPotentialModel controls:** `ipm_orientation_mode` accepts `random` or `fixed`; fixed mode uses `ipm_fixed_orientation_index`. Optional `ipm_vrel_log_prefix` and `ipm_momentum_log_prefix` enable CSV diagnostics. Species should reference offline sample files with the canonical `ipm_samples_file` key.

**Space charge:** Set `enable_space_charge` to true to activate Coulomb coupling. `space_charge_model` selects the backend: `auto` uses Direct for small ion counts and Grid for larger ion counts, `direct` forces the O(N²) direct summation model, `grid` forces the geometry-aware Poisson model, and `gpu` requires an available CUDA-backed model. `enable_space_charge_gpu=true` requests GPU P³M only in `auto` mode; CPU builds or unavailable GPU contexts fall back to Direct/Grid.

See `schema/physics.schema.json` for all options.

### Output Section

```json
"output": {
  "folder": "./results/my_sim",
  "trajectory_file": "output.h5",
  "trajectory_mode": "full",      // full|minimal
  "print_progress": true,
  "buffer_byte_cap": 0          // Optional RAM cap for trajectory buffer (bytes, 0 = unlimited)
}
```

`output.trajectory_mode="minimal"` disables `/trajectory` snapshots and writes compact per-ion final-state data to `/analysis/minimal_transport` instead.

**Output Format:** See [HDF5_OUTPUT_STRUCTURE.md](HDF5_OUTPUT_STRUCTURE.md) for file structure.

See `schema/output.schema.json` for all options.

### Domains Section

**Required:** `name`, `instrument`, `geometry`, `environment`

```json
"domains": [
  {
    "name": "drift_tube",
    "instrument": "IMS",           // IMS, TOF, LQIT, Orbitrap, Quadrupole, FTICR
    "geometry": {
      "radius_m": 0.05,
      "length_m": 0.2,
      "origin_m": [0, 0, 0]
    },
    "environment": {
      "pressure_Pa": 101325,
      "temperature_K": 300,
      "gas_species": "N2"
    },
    "fields": { /* DC, RF, or field arrays */ },
    "boundary": { "type": "Absorption" }
  }
]
```

See `schema/domain.schema.json`, `schema/geometry.schema.json`, `schema/environment.schema.json` for details.

**Performance note:** Domain lookup uses an axial/radial prefilter for cylindrical domains; complex geometries (e.g., Orbitrap) fall back to linear scan. Avoid excessive domain slicing unless necessary.

### Boundary Types

| Type | Description | Use Case |
|------|-------------|----------|
| `Absorption` | Ion deactivated at boundary | Detectors, open boundaries |
| `SpecularReflection` | Mirror-like reflection | Ideal metallic surfaces |
| `DiffuseReflection` | Cosine-distributed reflection | Rough surfaces |
| `ThermalReflection` | Thermalization to wall temperature | Gas-surface collisions |

See `schema/boundary.schema.json` for parameters.

### Field Configuration

**DC Fields:**
```json
"fields": {
  "DC": {
    "axial_V": 100.0,              // Can be constant or waveform
    "EN_Td": 50.0
  }
}
```

**Field Arrays (HDF5):**
```json
"field_array_terms": [
  {
    "file": "field_arrays/my_field.h5",
    "scale_kind": "DC_Axial",     // DC_Axial, RF_Voltage, Constant
    "scale_factor": 1.0
  }
]
```

**Waveforms:**
```json
"axial_V": {
  "type": "sinusoidal",
  "amplitude": 100.0,
  "frequency_Hz": 1e6,
  "offset": 0.0,
  "phase_rad": 0.0
}
```

**Supported waveform types:** `constant`, `linear`, `quadratic`, `exponential`, `sinusoidal`, `pwm`, `pulsed`, `arbitrary`

`environment.pressure_Pa` supports the same value-or-waveform syntax (number, inline waveform object, or `@waveform_ref`).

See `schema/waveform.schema.json` and `schema/fields.schema.json` for complete specifications.

---

## Schema Validation

### JSON Schema Files

All configuration schemas are located in `schema/`:

- **`icarion-config.schema.json`** - Master schema (top-level)
- **`simulation.schema.json`** - Simulation parameters
- **`physics.schema.json`** - Physics options
- **`output.schema.json`** - Output configuration
- **`domain.schema.json`** - Domain definition
- **`geometry.schema.json`** - Geometry specification
- **`environment.schema.json`** - Gas environment
- **`fields.schema.json`** - Electric/magnetic fields
- **`boundary.schema.json`** - Boundary conditions
- **`ions.schema.json`** - Ion initialization
- **`waveform.schema.json`** - Waveform definitions
- **`reactions.schema.json`** - Reaction database
- **`species.schema.json`** - Species database
- **`common-types.schema.json`** - Reusable type definitions

### Validating Your Config

#### Using Python Validator

```bash
# Validate a config file
python3 schema/validator.py schema/icarion-config.schema.json my_config.json

# Validate multiple files
for f in examples/*/*.json; do python3 schema/validator.py schema/icarion-config.schema.json "$f"; done
```

#### Using ICARION CLI

```bash
# Validate config without running simulation
./build/src/icarion_main --validate-config my_config.json

# Dry-run (load config, initialize, but don't simulate)
./build/src/icarion_main --dry-run my_config.json
```

---

## Multi-Gas Configurations

**Status:** Production-ready (v1.1.0)

ICARION supports multi-component gas mixtures for realistic collision and reaction simulations.

### Basic Usage

Define gas mixtures in the domain `environment` section:

```json
"environment": {
  "temperature_K": 300.0,
  "pressure_Pa": 101325.0,
  "gas_mixture": [
    { "species": "N2", "mole_fraction": 0.78 },
    { "species": "O2", "mole_fraction": 0.21 },
    { "species": "Ar", "mole_fraction": 0.01 }
  ]
}
```

**Requirements:**
- Mole fractions must sum to 1.0 (validated at load time)
- All species must exist in species database or have `cross_section_m2` specified

### Collision Cross Sections

**Option 1: Use species database (recommended)**
```json
{ "species": "N2", "mole_fraction": 0.78 }
// CCS loaded from species_database.json
```

**Option 2: Explicit cross section**
```json
{ "species": "N2", "mole_fraction": 0.78, "cross_section_m2": 4.0e-19 }
```

### Backward Compatibility

Single-gas configs still work:
```json
"environment": {
  "gas_species": "He",  // Old format
  "pressure_Pa": 101325.0,
  "temperature_K": 300.0
}
```

See `examples/ims/` for complete multi-gas examples.
## Helper Tools

### Config Creation Script

Located at `scripts/create_config.py`:

```bash
# Interactive mode (guided prompts)
python3 scripts/create_config.py

# Template-based creation
python3 scripts/create_config.py --template ims --output ims_config.json

# Customize template
python3 scripts/create_config.py \
  --template tof \
  --time 5e-4 \
  --timestep 1e-9 \
  --pressure 1e-6 \
  --temperature 300 \
  --output my_tof.json

# List available templates
python3 scripts/create_config.py --list-templates
```

**Available Templates:**

- `ims` - Ion Mobility Spectrometry (drift tube)
- `tof` - Time-of-Flight mass spectrometer
- `lqit` - Linear Quadrupole Ion Trap
- `orbitrap` - Orbitrap mass analyzer
- `minimal` - Minimal valid configuration

---

### Orbitrap Geometry (Hyperlogarithmic Electrodes)

Orbitrap mass analyzers use hyperlogarithmic electrode geometry (Kharchenko et al. 2021, DOI: 10.1007/s13361-011-0325-3). The electrode surfaces follow:

```
z² = 0.5(r² - R²) + R_m² × ln(R/r)
```

Where:
- **z** = axial coordinate (origin at trap center)
- **r** = radial coordinate
- **R** = electrode radius at z=0 (use `radius_in_m` for inner, `radius_out_m` for outer electrode)
- **R_m** = characteristic radius (field shaping parameter, `radius_char_m`)

**Parameter Guidelines:**

1. **Physical Constraint:** `radius_char_m > radius_out_m` (required for stable quadro-logarithmic potential)
2. **Typical Values:**
   - Commercial Orbitrap: R_in ≈ 6mm, R_out ≈ 15mm, R_m ≈ 22mm
3. **Gap Width:** Determined by R_out - R_in (typically 4-10mm)
4. **Field Strength:** Controlled by `k` parameter in DC fields (typically 380-500 V/mm²)

k = 2*V_0/(R_m² × ln(R_out/R_in) - 0.5(R_out² - R_in²)) - see (Kharchenko et al. 2021, DOI: 10.1007/s13361-011-0325-3)

**Boundary Checking:**

ICARION uses bisection to compute r_in(z) and r_out(z) at each ion position, validating:
```
r_in(z) ≤ r ≤ r_out(z)  AND  -length_m/2 ≤ z ≤ length_m/2
```

Ions outside these boundaries are deactivated (electrode collision).

---

## Examples

Complete configuration examples for all supported instruments are available in `examples/`:

- **IMS**: `examples/ims/ims_basic.json` - Ion mobility spectrometry
- **TIMS**: `examples/ims/ims_tims_basic.json` - Trapped IMS axial ramp field with gas flow
- **TOF**: `examples/tof/tof_linear.json` - Time-of-flight mass spectrometry
- **LQIT**: `examples/lqit/lqit_basic.json` - Linear quadrupole ion trap
- **Orbitrap**: `examples/orbitrap/orbitrap_basic.json` - Orbitrap mass analyzer
- **Quadrupole**: `examples/quadrupole/quad_basic.json` - Quadrupole mass filter
- **Field Arrays**: `examples/field_arrays/` - HDF5 field array examples
- **Reactions**: `examples/reactions/` - Chemical reaction simulations

### Quick Example: IMS Drift Tube

```json
{
  "simulation": {
    "total_time_s": 5e-4,
    "dt_s": 1e-8,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "HSS"
  },
  "output": {
    "directory": "results/ims",
    "basename": "ims_output"
  },
  "ions": [
    {
      "species": "H3O+",
      "count": 1000,
      "position": {
        "type": "Gaussian",
        "mean": [0, 0, 0.001],
        "std_dev": [0.001, 0.001, 0.0005]
      }
    }
  ],
  "domains": [
    {
      "name": "drift_region",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0, 0, 0],
        "length_m": 0.05,
        "radius_m": 0.015
      },
      "environment": {
        "pressure_Pa": 200.0,
        "temperature_K": 300.0,
        "gas_species": "He"
      },
      "fields": {
        "dc": {
          "EN_Td": 10.0
        }
      }
    }
  ]
}
```

See `examples/README.md` for detailed documentation of all example configurations.

---

### TIMS Axial Ramp Fields

TIMS can be declared with `instrument: "TIMS"` and a `fields.TIMS` block. The
analytical model evaluates the local axial field as:

```text
E_z(z,t) = (1 - f(t)) * E_initial(z) + f(t) * E_final(z)
```

`axial_field_initial_*` and `axial_field_final_*` can be uniform values or 1D
profiles sampled at `z_positions_m`. The ramp fraction `f(t)` is linear or exponential between
`ramp_start_s` and `ramp_end_s`, or can be supplied directly via
`ramp_fraction`. Standard IMS RF/DC settings still provide radial confinement.

TIMS examples usually also set an axial gas flow:

```json
"env": {
  "temperature_K": 300.0,
  "pressure_Pa": 310.0,
  "gas_species": "N2",
  "flow_model": "axial_uniform",
  "flow_parameters": {
    "axial_flow_velocity_m_s": 134.0
  }
}
```

Axial flow profiles are defined in the local domain coordinate system, matching
the TIMS axial field direction. For rotated domains, ICARION evaluates the
radial profile in local coordinates and rotates the velocity back to global
coordinates before applying damping or stochastic collisions.

See `examples/ims/ims_tims_basic.json` for a complete configuration.

---

## Common Patterns

### Per-Domain Integrator Override

Use the global `simulation.integrator` as default, but override for specific domains:

```json
{
  "simulation": {
    "integrator": "RK4"  // Global default
  },
  "domains": [
    {
      "name": "low_pressure_region",
      "integrator": "RK45",  // Override: use adaptive stepping
      "environment": { "pressure_Pa": 1e-6, ... }
    },
    {
      "name": "high_pressure_region",
      // No integrator specified → uses "RK4" from simulation
      "environment": { "pressure_Pa": 101325, ... }
    }
  ]
}
```

### Collision Model Selection

```json
// No collisions (vacuum, short flight times)
"physics": { "collision_model": "NoCollisions" }

// High pressure (IMS, drift tubes)
"physics": { "collision_model": "Friction" }

// High precision (research simulations)
"physics": { "collision_model": "EHSS" }
```

### Environment cache and boundaries

- Collisions/reactions use the gas environment cached at the start of the macro-timestep (and after explicit domain switches). Mid-step boundary crossings are applied only after integration, so stochastic rates during that step still use the pre-step environment. Keep `timestep` small across sharp pressure/temperature jumps or accept this approximation.

### LQIT Dipolar AC Excitation

Legacy LQIT configs can keep using `fields.AC.voltage_V` and
`fields.AC.frequency_Hz`; that remains a single local x-axis dipolar drive.

For independent radial x/y excitation, use the top-level
`fields.dipolar_excitation` block next to `AC`, `RF`, and `DC`:

```json
"fields": {
  "AC": {
    "voltage_V": 10.0,
    "frequency_Hz": 200000.0,
    "phase_rad": 0.0
  },
  "dipolar_excitation": {
    "x": {
      "enabled": true,
      "amplitude_V": 10.0,
      "frequency_Hz": 200000.0,
      "phase_rad": 0.0
    },
    "y": {
      "enabled": true,
      "amplitude_V": 2.5,
      "frequency_Hz": 200000.0,
      "phase_rad": 1.57079632679,
      "ramp": 1.0
    }
  }
}
```

`amplitude_V`, `frequency_Hz`, and `ramp` use the same value-or-waveform syntax
as other field quantities. `fields.AC.dipolar_excitation` and
`fields.AC.x`/`fields.AC.y` are accepted as compatibility shorthands, but the
top-level `fields.dipolar_excitation` block is the preferred SSOT. The loader
rejects mixed top-level/nested definitions and rejects axes that specify both
`amplitude_V` and its legacy alias `voltage_V`.

---

## Troubleshooting

### Common Errors

**Error: Missing required field 'dt_s' in simulation config**

- Solution: Add `"dt_s": 1e-9` to simulation section

**Error: Missing required field 'total_time_s' in simulation config**

- Solution: Add `"total_time_s": 1e-3` to simulation section

**Error: Domain missing required 'name' field**

- Solution: Add `"name": "my_domain"` to each domain

**Error: Configuration validation failed**

- Solution: Run `--validate-config` to see detailed error messages

**Error: Species 'X+' not found in species database**

- Solution: Add species 'X+' to the specified species database or check for typos.

**Error: Could not open ion cloud JSON file:**

- Solution: Verify the existence of the specified ion cloud and check the file path.

### Best Practices

1. **Always validate configs** before running expensive simulations:

   ```bash
   ./build/src/icarion_main --validate-config my_config.json
   ```

2. **Start with a template** from `examples/` directory

3. **Use appropriate time steps**:

   - IMS (collisional): `dt_s = 1e-10` to `1e-9`
   - TOF (vacuum): `dt_s = 1e-9` to `1e-8`

4. **Check units**: All values in SI units (meters, seconds, Pascal, Kelvin)

5. **Use OpenMP first** for large CPU simulations. Enable GPU only for builds and physics paths that explicitly support it; unsupported GPU paths fall back to CPU.

---

## See Also

- **Schema Documentation**: [`../schema/README.md`](../schema/README.md)
- **CLI Reference**: [`CLI_USAGE.md`](CLI_USAGE.md)
- **Architecture**: [`ARCHITECTURE.md`](ARCHITECTURE.md)
- **Example Configs**: [`../examples/`](../examples/)
- **Python Helper**: `scripts/create_config.py --help`

---

## Version

This guide corresponds to the ICARION v1.1.0 configuration surface.

For schema version history, see the git log for the `schema/` directory.
