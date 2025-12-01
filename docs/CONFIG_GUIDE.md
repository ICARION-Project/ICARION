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
      "env": {
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
- EHSS: uses `CCS_EHSS[gas]` if present, else geometry, else (without geometry) throws.
- Tool: `./ccs_precompute --input species.json --output out.json --species H3O+ --ref-gas He --ref-ccs-A2 110.0 [--model HSS|EHSS] [--override] [--n-orientations 300]`.

---

## Species and Reaction Databases (v1.0)

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
      "rate_constant": 3.5e-9,
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
      "rate_constant": 1.5e-9,
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
      "rate_constant": 2.0e-9,
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
  - `species`: Species ID for concentration dependence (or `"neutral"` for buffer gas)
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
| **#3: Species exists** | `species ∈ species_db` (if not `"neutral"`) | ERROR: `"species": "XYZ"` → "species 'XYZ' not found" |
| **#4: No duplicate species** | Each `species` appears once | ERROR: Two terms with `"species": "O2"` → "duplicate order term for 'O2'. Use exponent=2 instead." |
| **#5: Max one buffer gas** | At most one term with `concentration_m3 = -1` | ERROR: Two terms with `-1` → "only one term can use buffer gas fallback" |

**Dimensional Consistency Warnings:**

ICARION also warns about likely unit mismatches:

```
⚠  Reaction 'rxn_001': 2nd-order (exponent=1) but k = 1.5e-30 m⁶/s outside typical range [1e-12, 1e-6] m³/s
```

**Typical Rate Constant Ranges:**
- **1st-order (spontaneous):** k ~ 10⁻³ to 10⁶ s⁻¹
- **2nd-order (ion + neutral):** k ~ 10⁻¹² to 10⁻⁶ m³/s
- **3rd-order (termolecular):** k ~ 10⁻³⁰ to 10⁻²⁴ m⁶/s

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
      "rate_constant": 2.5e-9,
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
| `species` | string | Species ID or `"neutral"` | Which species to use for concentration |
| `exponent` | integer | 0, 1, 2 | Power to raise concentration to |
| `concentration_m3` | number | ≥ -1.0 | Explicit concentration [m⁻³] or -1.0 for buffer gas |

**Buffer Gas Fallback:**

- If `concentration_m3 = -1.0` (or omitted), use buffer gas density from `EnvironmentConfig.particle_density_m_3`
- If `concentration_m3 > 0`, use the explicit value

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
  "rate_constant": 2.5e-9,
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
  "rate_constant": 1.2e-28,
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

## Required Fields

### Simulation Section

```json
"simulation": {
  "total_time_s": 1e-3,      // REQUIRED: Total simulation time [seconds]
  "dt_s": 1e-9,              // REQUIRED: Time step [seconds]
  "integrator": "RK4",       // Default integrator (can be overridden per-domain)
  "write_interval": 100,     // Steps between trajectory snapshots
  "enable_gpu": false,       // Enable GPU acceleration
  "enable_openmp": false,    // Enable OpenMP threading
  "rng_seed": 42             // Random number generator seed
}
```

**Available Integrators:**

- `"RK4"` - 4th order Runge-Kutta (default, stable)
- `"RK45"` - Runge-Kutta-Fehlberg (adaptive step size)
- `"Boris"` - Boris integrator (for strong magnetic fields)

### Physics Section

```json
"physics": {
  "collision_model": "HSS",            // Collision model (required)
  "enable_reactions": false,           // Enable chemical reactions
  "enable_space_charge": false,        // Enable Coulomb forces between ions
  "enable_ou_thermalization": false,   // Enable Ornstein-Uhlenbeck thermalization
  "force_ou_for_stochastic": false     // Force OU for stochastic models
}
```

**Available Collision Models:**

- `"NoCollisions"` - Free flight (no collisions)

- `"HSD"` - Hard Sphere Deterministic (collision frequency does not depend on velocity)
- `"Friction"` - Frictional-drag model (depending on mobility)
- `"Langevin"` - Langevin dynamics (long-range and velocity-dependent collision frequency)

- `"HSS"` - Hard Sphere Stochastic  
- `"EHSS"` - Exact Hard Sphere Scattering (recommended for low pressures)

**Feature Flags:**

- **`enable_reactions`**: Enable stochastic chemical reactions (requires `reaction_database`)
- **`enable_space_charge`**: Enable Coulomb (space charge) forces between ions
  - **N < 1000**: Direct N-body summation (exact, O(N²))
  - **N ≥ 1000**: Grid-based Poisson solver (fast, O(N log N))
  - Auto-selects method based on ion count
- **`enable_ou_thermalization`**: Apply Ornstein-Uhlenbeck velocity kicks for thermalization
- **`force_ou_for_stochastic`**: Force OU thermalization even for stochastic collision models

### Output Section

```json
"output": {
  "folder": "./results",              // Output directory (relative or absolute)
  "trajectory_file": "trajectories.h5", // HDF5 trajectory file name
  "print_progress": true              // Print progress to console
}
```

### Domains Section

Each domain represents a physical region with specific instrument geometry and conditions:

```json
"domains": [
  {
    "name": "drift_region",           // REQUIRED: Domain identifier
    "instrument": "IMS",              // REQUIRED: Instrument type
    "integrator": "RK4",              // OPTIONAL: Override global integrator
    
    "geometry": {                     // REQUIRED
      "origin_m": [0.0, 0.0, 0.0],    // Starting position [x, y, z] in meters
      "length_m": 0.05,               // Domain length along z-axis
      "radius_m": 0.01,               // Cylindrical radius (for IMS, TOF, LQIT, etc.)
      
      // Orbitrap-specific parameters (hyperlogarithmic electrodes):
      "radius_in_m": 0.006,           // Inner electrode radius at z=0 [m]
      "radius_out_m": 0.015,          // Outer electrode radius at z=0 [m]
      "radius_char_m": 0.022          // Characteristic radius R_m [m], must be > radius_out_m
    },
    
    "env": {                          // REQUIRED
      "pressure_Pa": 101325.0,        // Gas pressure [Pascal]
      "temperature_K": 300.0,         // Gas temperature [Kelvin]
      "gas_species": "He",            // Gas species (He, N2, Ar, etc.) fallback if no mixture
      "gas_mixture": [                // OPTIONAL mixture; overrides gas_species when non-empty
        { "species": "N2", "mole_fraction": 0.78, "cross_section_m2": 4.0e-19 },
        { "species": "O2", "mole_fraction": 0.21 },
        { "species": "H2O", "mole_fraction": 0.01 }
      ],
      "gas_velocity_m_s": [0, 0, 0]   // Gas flow velocity [m/s]
    },
    
    "fields": {                       // OPTIONAL
      "DC": {
        "axial_V": 250.0,             // Axial DC voltage [V]
        "EN_Td": 10.0                 // Reduced field strength [Townsend]
      },
      "RF": {
        "voltage_V": 0.0,             // RF amplitude [V]
        "frequency_Hz": 1e6,          // RF frequency [Hz]
        "phase_rad": 0.0              // RF phase [radians]
      },
      "field_array_terms": [          // OPTIONAL: Load fields from HDF5
        {
          "file": "field_arrays/my_geometry.h5",  // HDF5 file path
          "scale_kind": "DC_Axial",    // Scaling mode (see below)
          "scale_factor": 1.0          // Additional multiplier (default: 1.0)
        }
      ]
    },
    
    "boundary": {                     // OPTIONAL: Ion-boundary interaction
      "type": "Absorption",           // Action type (see below)
      "accommodation_coeff": 1.0,     // Accommodation coefficient [0-1] (for reflections)
      "temperature_K": 300.0          // Wall temperature [K] (for thermal reflections)
    }
  }
]
```

**Available Instrument Types:**

- `"IMS"` - Ion Mobility Spectrometry
- `"TOF"` - Time-of-Flight
- `"LQIT"` - Linear Quadrupole Ion Trap
- `"Orbitrap"` - Orbitrap mass analyzer
- `"Quadrupole"` - Quadrupole mass filter
- `"FTICR"` - Fourier Transform Ion Cyclotron Resonance
- `"NoFixedInstrument"` - Generic region

### Boundary Configuration (Ion-Wall Interactions)

**Status:** Production-ready (v1.0)

ICARION supports multiple boundary action types to model ion-wall interactions realistically.

#### Configuration

```json
"boundary": {
  "type": "Absorption",              // REQUIRED: Boundary action type
  "accommodation_coeff": 1.0,        // OPTIONAL: Accommodation coefficient [0-1]
  "temperature_K": 300.0             // OPTIONAL: Wall temperature [K]
}
```

**Default Behavior:** If `boundary` is omitted, the domain defaults to `Absorption`.

#### Available Boundary Types

| Type | Description | Physics | Required Fields |
|------|-------------|---------|-----------------|
| `"Absorption"` | Ion is deactivated upon boundary collision | Perfectly absorbing wall (no reflection) | `type` only |
| `"SpecularReflection"` | Perfect elastic reflection with angle of incidence = angle of reflection | Mirror-like surface, momentum conserved | `type` only |
| `"DiffuseReflection"` | Random reflection direction (cosine distribution) | Rough surface, no memory of incidence angle | `type`, optionally `accommodation_coeff` |
| `"ThermalReflection"` | Ion thermalized to wall temperature upon reflection | Energy exchange with wall | `type`, `accommodation_coeff`, `temperature_K` |

#### Physics Details

**Absorption:**
- Ion is marked as inactive (lost) when crossing boundary
- No velocity or position update
- Use case: Detectors, absorbing electrodes, open boundaries

**Specular Reflection:**
- Velocity component normal to surface is reversed: v⊥ → -v⊥
- Tangential components preserved: v∥ → v∥
- Total kinetic energy conserved
- Use case: Ideal metallic electrodes, mirror traps

**Diffuse Reflection:**
- Reflection direction sampled from cosine distribution (Lambert's law)
- Accommodation coefficient α controls energy retention:
  - α = 0: Fully elastic (like specular)
  - α = 1: Full thermalization (like thermal)
  - 0 < α < 1: Partial energy exchange
- Use case: Rough surfaces, partial thermalization

**Thermal Reflection:**
- Ion velocity resampled from Maxwell-Boltzmann distribution at wall temperature
- Full thermalization: K.E. = (3/2) k_B T_wall
- Accommodation coefficient α blends thermal and elastic components
- Use case: Gas-surface collisions, thermal accommodation experiments

#### Example Configurations

**1. Absorbing Detector:**
```json
"boundary": {
  "type": "Absorption"
}
```

**2. Ideal Metallic Electrodes:**
```json
"boundary": {
  "type": "SpecularReflection"
}
```

**3. Rough Surface (Partial Thermalization):**
```json
"boundary": {
  "type": "DiffuseReflection",
  "accommodation_coeff": 0.7,     // 70% energy accommodation
  "temperature_K": 300.0
}
```

**4. Gas-Surface Collision (Full Thermalization):**
```json
"boundary": {
  "type": "ThermalReflection",
  "accommodation_coeff": 1.0,     // Full thermalization
  "temperature_K": 300.0          // Wall at room temperature
}
```

#### GPU Acceleration Compatibility (Phase 11)

**GPU Support:**
- `Absorption`: Fully supported on GPU (cylindrical domains only)
- Reflections (`Specular`, `Diffuse`, `Thermal`): CPU fallback required
- Orbitrap geometry: CPU fallback required (hyperlogarithmic surface)

**Automatic Dispatch:**
ICARION automatically selects GPU or CPU boundary checking based on configuration:
- If `boundary.type == "Absorption"` AND `instrument != "Orbitrap"` → GPU used (if enabled)
- Otherwise → CPU fallback (supports all features)

**Performance Note:** Reflection boundary checks have minimal performance impact (<1% of total runtime) compared to integration and collision handling. GPU acceleration is primarily beneficial for integration (Phase 11 complete).

### Field Arrays (HDF5-based Electric Fields)

**Status:** Production-ready

For complex geometries where analytical field solutions are unavailable, ICARION supports loading pre-computed electric field arrays from HDF5 files.

#### Field Array Configuration

```json
"fields": {
  "field_array_terms": [
    {
      "file": "field_arrays/electrode_geometry.h5",
      "scale_kind": "DC_Axial",
      "scale_factor": 1.0
    }
  ]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `file` | string | YES | Path to HDF5 field array file (relative or absolute) |
| `scale_kind` | enum | YES | Voltage scaling mode (see below) |
| `scale_factor` | number | NO | Additional constant multiplier (default: 1.0) |

#### Scaling Modes (`scale_kind`)

Field arrays are stored **normalized to 1 Volt**. At runtime, fields are scaled by applied voltages:

| `scale_kind` | Formula | Use Case | Example |
|--------------|---------|----------|--------|
| `"Constant"` | E = E_array × `scale_factor` | Static scaling | Test fields, uniform fields |
| `"DC_Axial"` | E = E_array × DC_voltage | Drift tubes, TOF | IMS drift field, TOF acceleration |
| `"DC_Quad"` | E = E_array × DC_voltage | Quadrupole fields | Quadrupole mass filter, ion guides |
| `"DC_Radial"` | E = E_array × DC_voltage | Cylindrical electrodes | Cylindrical lenses, apertures |
| `"RF"` | E = E_array × V_rf × cos(2πft+φ) | Time-dependent RF | RF traps, Orbitraps, quadrupoles |

**Example:** For a 50mm drift tube with 100V applied:
- Field array stores: Ez = 20 V/m (= 1V / 0.05m)
- At runtime with `DC_Axial`: E_scaled = 20 V/m × 100 V = 2000 V/m ✓

#### HDF5 File Format

Field arrays must follow this structure:

```
field_array.h5
├── x [1D array]        # x-coordinates [m], size nx
├── y [1D array]        # y-coordinates [m], size ny
├── z [1D array]        # z-coordinates [m], size nz
├── Ex [3D array]       # E-field x-component [V/m], shape (nx, ny, nz)
├── Ey [3D array]       # E-field y-component [V/m], shape (nx, ny, nz)
├── Ez [3D array]       # E-field z-component [V/m], shape (nx, ny, nz)
└── phi [3D array]      # Potential [V] (optional), shape (nx, ny, nz)
```

**Requirements:**
- Coordinate arrays (x, y, z) must be 1D and uniformly spaced
- Field arrays (Ex, Ey, Ez) must be 3D with shape matching (nx, ny, nz)
- Potential (phi) is optional
- All coordinates in meters [m], fields in volts per meter [V/m]

#### Creating Field Arrays with Python

Generate HDF5 field arrays using `h5py`:

```python
import h5py
import numpy as np

# Define 3D grid
x = np.linspace(-5e-3, 5e-3, 10)   # ±5mm, 10 points
y = np.linspace(-5e-3, 5e-3, 10)   # ±5mm, 10 points
z = np.linspace(0, 50e-3, 20)      # 0-50mm, 20 points
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Compute fields (normalized to 1V reference)
# Example: Uniform axial field, 1V over 50mm
Ex = np.zeros_like(X)
Ey = np.zeros_like(Y)
Ez = np.ones_like(Z) * 20.0  # 1V / 0.05m = 20 V/m
phi = -Ez * Z                # Potential φ(z) = -∫E·dz

# Save to HDF5
with h5py.File('my_field.h5', 'w') as f:
    f.create_dataset('x', data=x)
    f.create_dataset('y', data=y)
    f.create_dataset('z', data=z)
    f.create_dataset('Ex', data=Ex)
    f.create_dataset('Ey', data=Ey)
    f.create_dataset('Ez', data=Ez)
    f.create_dataset('phi', data=phi)  # Optional
```

**Example Generator:** `examples/field_arrays/create_example_field_array.py`

#### Field Interpolation

At runtime, ICARION uses **trilinear interpolation** to evaluate fields at arbitrary positions:

- **Method:** 8-corner weighted average
- **Accuracy:** O(h²) convergence with grid spacing h
- **Boundary handling:** Returns zero field outside grid

#### Complete Example

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
      "name": "custom_geometry",
      "instrument": "NoFixedInstrument",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.05,
        "radius_m": 0.01
      },
      "env": {
        "pressure_Pa": 1e-6,
        "temperature_K": 300.0,
        "gas_species": "He"
      },
      "fields": {
        "field_array_terms": [
          {
            "file": "field_arrays/electrode_set_1.h5",
            "scale_kind": "DC_Axial",
            "scale_factor": 1.0
          },
          {
            "file": "field_arrays/electrode_set_2.h5",
            "scale_kind": "Constant",
            "scale_factor": 50.0
          }
        ],
        "DC": {
          "axial_V": 100.0
        }
      }
    }
  ]
}
```

**Note:** Multiple field array terms are summed:
```
E_total = Σᵢ (E_array_i × scale_i) + E_analytical
```

#### Validation

Validate field arrays before running simulations:

```bash
# Test HDF5 file structure
python3 -c "import h5py; f=h5py.File('my_field.h5','r'); print(list(f.keys()))"

# Verify dimensions
python3 -c "
import h5py
with h5py.File('my_field.h5', 'r') as f:
    print(f'x: {f[\"x\"].shape}, Ex: {f[\"Ex\"].shape}')
"

# Run E2E tests
cd build && ctest -R field_array
```

**Test Files:**
- `tests/config/test_field_array_e2e.cpp` - End-to-end validation (source only, not built)
- `tests/config/test_field_array_terms_loader.cpp` - JSON loading tests (executable: `build/tests/config/test_field_array_terms_loader`)

**Recommendation:** Use 1-2mm grid spacing for typical ion optics geometries.

---

### Waveforms (Time-Varying Parameters)

**Status:** Production-ready (v1.0)

ICARION supports time-varying field parameters through a flexible waveform system. Any voltage or frequency parameter can be static OR time-dependent.

#### Quick Examples

**Static Value:**
```json
"fields": {
  "AC": {
    "voltage_V": 100.0,
    "frequency_Hz": 1000000.0
  }
}
```

**Linear Voltage Ramp:**
```json
"fields": {
  "AC": {
    "voltage_V": {
      "type": "linear",
      "start": 0.0,
      "end": 500.0,
      "start_time_s": 0.0,
      "end_time_s": 0.001,
      "clamp": true
    },
    "frequency_Hz": 1000000.0
  }
}
```

**Named Waveforms (Reusable):**
```json
{
  "domains": [
    {
      "name": "domain_1",
      "fields": {
        "waveforms": {
          "ac_voltage_ramp": {
            "type": "linear",
            "start": 0,
            "end": 500,
            "end_time_s": 0.001
          },
          "rf_chirp": {
            "type": "linear",
            "start": 1000000,
            "end": 2000000,
            "end_time_s": 0.01
          }
        },
        "AC": {
          "voltage_V": "@ac_voltage_ramp",
          "frequency_Hz": 1000000
        },
        "RF": {
          "frequency_Hz": "@rf_chirp"
        }
      }
    }
  ]
}
```

#### Waveform Types

ICARION supports 6 waveform types:

| Type | Use Case | Parameters |
|------|----------|------------|
| `constant` | Static value | `value` |
| `linear` | Voltage/frequency sweep | `start`, `end`, `end_time_s`, `start_time_s` (opt), `clamp` (opt) |
| `quadratic` | Acceleration profile | `a`, `b`, `c`, `start_time_s` (opt), `end_time_s` (opt) |
| `sinusoidal` | Amplitude modulation | `offset` (opt), `amplitude`, `frequency_Hz`, `phase_rad` (opt) |
| `pulsed` | Single pulse | `low`, `high`, `pulse_start_s`, `pulse_width_s` |
| `arbitrary` | Custom curve | `times[]`, `values[]`, `interpolation` (opt) |

#### Waveform Type Details

##### 1. Constant Waveform

```json
{
  "type": "constant",
  "value": 123.45
}
```

**Shorthand:** Just use a number: `"voltage_V": 123.45`

##### 2. Linear Waveform

Linearly interpolates from `start` to `end` between `start_time_s` and `end_time_s`:

```json
{
  "type": "linear",
  "start": 0.0,          // Starting value
  "end": 500.0,          // Ending value
  "start_time_s": 0.0,   // Ramp start time [s] (default: 0)
  "end_time_s": 0.001,   // Ramp end time [s] (required)
  "clamp": true          // Hold end value after end_time_s (default: true)
}
```

**Behavior:**
- `t < start_time_s`: Returns `start`
- `start_time_s ≤ t ≤ end_time_s`: Linear interpolation
- `t > end_time_s`:
  - If `clamp=true`: Returns `end` (hold final value)
  - If `clamp=false`: Returns `start` (reset to initial)

**Example Use Cases:**
- IMS voltage ramp: 0V → 500V over 1ms
- TOF acceleration sweep
- Frequency chirp (linear frequency increase)

##### 3. Quadratic Waveform

Quadratic function: **y(t) = a + b×t + c×t²**

```json
{
  "type": "quadratic",
  "a": 10.0,             // Constant term
  "b": 5.0,              // Linear coefficient
  "c": 2.0,              // Quadratic coefficient
  "start_time_s": 0.0,   // Active from this time (default: 0)
  "end_time_s": 1.0      // Active until this time (default: 1e9)
}
```

**Behavior:**
- `t < start_time_s` or `t > end_time_s`: Returns `a` (constant term)
- `start_time_s ≤ t ≤ end_time_s`: Evaluates quadratic

**Example Use Cases:**
- Acceleration profiles (constant acceleration = quadratic position)
- Nonlinear voltage ramps

##### 4. Sinusoidal Waveform

Sinusoidal modulation: **y(t) = offset + amplitude × sin(2πft + φ)**

```json
{
  "type": "sinusoidal",
  "offset": 250.0,       // DC offset (default: 0)
  "amplitude": 50.0,     // Oscillation amplitude
  "frequency_Hz": 100.0, // Modulation frequency [Hz]
  "phase_rad": 0.0       // Phase offset [radians] (default: 0)
}
```

**Example Use Cases:**
- Amplitude modulation of RF voltages
- Oscillating DC bias
- Phase-locked excitation

##### 5. Pulsed Waveform

Single rectangular pulse:

```json
{
  "type": "pulsed",
  "low": 10.0,           // Baseline value
  "high": 100.0,         // Pulse value
  "pulse_start_s": 1.0,  // Pulse start time [s]
  "pulse_width_s": 0.5   // Pulse duration [s]
}
```

**Behavior:**
- `t < pulse_start_s`: Returns `low`
- `pulse_start_s ≤ t < pulse_start_s + pulse_width_s`: Returns `high`
- `t ≥ pulse_start_s + pulse_width_s`: Returns `low`

**Example Use Cases:**
- Pulsed ion injection
- Gating voltages
- Resonant excitation pulses

##### 6. Arbitrary Waveform

Custom waveform defined by time-value pairs with interpolation:

```json
{
  "type": "arbitrary",
  "times": [0.0, 0.001, 0.002, 0.003, 0.004],
  "values": [0, 100, 50, 200, 0],
  "interpolation": "linear"   // "linear" (default), "step", or "cubic"
}
```

**Requirements:**
- `times` and `values` must have same length (≥2 points)
- `times` must be strictly increasing

**Interpolation Modes:**
- `"linear"`: Linear interpolation between points (default)
- `"step"`: Step function (hold previous value)
- `"cubic"`: Cubic Hermite spline (smooth)

**Behavior:**
- `t < times[0]`: Returns `values[0]`
- `times[0] ≤ t ≤ times[n]`: Interpolates
- `t > times[n]`: Returns `values[n]`

**Example Use Cases:**
- Complex voltage profiles from experiments
- Imported waveforms from external tools
- Multi-stage ramps

#### Waveform Library (Named Waveforms)

Define waveforms once and reference them multiple times. Waveforms can be defined at **two levels**:

1. **Global level (top-level `"waveforms"`):** Shared across all domains
2. **Domain level (`"fields.waveforms"`):** Domain-specific overrides

```json
{
  "waveforms": {
    "voltage_ramp": {
      "type": "linear",
      "start": 0,
      "end": 500,
      "start_time_s": 0.0,
      "end_time_s": 0.005,
      "clamp": true
    },
    "freq_chirp": {
      "type": "linear",
      "start": 1000000,
      "end": 2000000,
      "start_time_s": 0.0,
      "end_time_s": 0.01,
      "clamp": true
    }
  },
  "domains": [
    {
      "name": "ims_region_1",
      "fields": {
        "DC": {
          "axial_V": "@voltage_ramp"
        },
        "AC": {
          "voltage_V": 200,
          "frequency_Hz": "@freq_chirp"
        }
      }
    },
    {
      "name": "ims_region_2",
      "fields": {
        "DC": {
          "axial_V": "@voltage_ramp"
        },
        "RF": {
          "voltage_V": 100,
          "frequency_Hz": "@freq_chirp"
        }
      }
    }
  ]
}
```

**Resolution Order:**
1. Check domain-local `fields.waveforms` library first (if exists)
2. Check global top-level `waveforms` library second
3. Error if not found in either

**Benefits:**
- ✅ **SSOT compliance:** Define shared waveforms once at top level
- ✅ **No duplication:** Same waveform used in multiple domains → single definition
- ✅ **Self-documenting:** Named waveforms clarify intent
- ✅ **Override capability:** Domains can override global waveforms locally
- ✅ **Consistency:** Change global waveform → updates all domains that reference it

**Syntax:**
- Waveform reference: `"@waveform_id"`
- Must start with `@` symbol
- ID must exist in global `waveforms` or domain-local `fields.waveforms` library

**Example Use Cases:**
- **Global waveforms:** Voltage ramps, frequency chirps used across multiple domains
- **Domain-local overrides:** Same waveform name but domain-specific parameters
- **Mixed usage:** Some waveforms global (reusable), others local (domain-specific)

#### ValueOrWaveform Pattern

Every field parameter in ICARION can be:

1. **Static value** (number): `"voltage_V": 100.0`
2. **Inline waveform** (object): `"voltage_V": {"type": "linear", ...}`
3. **Waveform reference** (string): `"voltage_V": "@my_waveform"`

**Exactly one** option must be specified.

#### Applicable Fields

Waveforms can be used for:

**DC Fields:**
- `axial_V` - Axial DC voltage
- `quad_V` - Quadrupole DC voltage
- `radial_V` - Radial DC voltage

**RF Fields:**
- `voltage_V` - RF voltage amplitude
- `frequency_Hz` - RF frequency

**AC Fields:**
- `voltage_V` - AC voltage amplitude  
- `frequency_Hz` - AC frequency

#### Static and Dynamic Values

Fields can be specified as either static values or waveforms:

```json
// Static value (simple)
"AC": {
  "voltage_V": 100,
  "frequency_Hz": 1000000
}

// Same as waveform constant:
"AC": {
  "voltage_V": {"constant_value": 100},
  "frequency_Hz": {"constant_value": 1000000}
}
```

**Migration:** Use linear waveforms instead (see Migration Script below).

#### Complete Example

```json
{
  "simulation": {
    "total_time_s": 0.01,
    "dt_s": 1e-9,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "HSS"
  },
  "output": {
    "folder": "./results",
    "trajectory_file": "waveform_test.h5"
  },
  "domains": [
    {
      "name": "ims_region",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0, 0, 0],
        "length_m": 0.05,
        "radius_m": 0.01
      },
      "env": {
        "pressure_Pa": 101325,
        "temperature_K": 300,
        "gas_species": "N2"
      },
      "fields": {
        "waveforms": {
          "voltage_ramp": {
            "type": "linear",
            "start": 0,
            "end": 500,
            "end_time_s": 0.005
          },
          "freq_chirp": {
            "type": "linear",
            "start": 1000000,
            "end": 2000000,
            "end_time_s": 0.01
          },
          "modulation": {
            "type": "sinusoidal",
            "offset": 250,
            "amplitude": 50,
            "frequency_Hz": 100
          }
        },
        "DC": {
          "axial_V": 100.0
        },
        "AC": {
          "voltage_V": "@voltage_ramp",
          "frequency_Hz": "@freq_chirp"
        },
        "RF": {
          "voltage_V": "@modulation",
          "frequency_Hz": 1000000
        }
      }
    }
  ]
}
```

#### Validation

Validate waveform configs:

```bash
# Validate JSON schema
python3 schema/validator.py schema/icarion-config.schema.json my_config.json

# Test waveform evaluation
cd build && ctest -R Waveform
```

**Tests:**
- `test_waveform_types` - 26 tests for waveform evaluation
- `test_waveform_loader` - 22 tests for JSON parsing

#### Performance Notes

- **Waveform evaluation:** Expected to be negligible (simple arithmetic + std::variant dispatch)
- **Memory:** ~128 bytes per ValueOrWaveform field (sizeof estimate)
- **Recommended:** Use waveform library for frequently referenced waveforms
- **Note:** Formal benchmarking pending (see validation/performance/ roadmap)

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

**Status:** Production-ready (v1.1)

ICARION supports multi-component gas mixtures for realistic collision and reaction simulations. This enables modeling of air mixtures, doped gases, and trace contaminants.

### Gas Mixture Basics

Gas mixtures are defined in the domain `env` section using the `gas_mixture` array:

```json
{
  "domains": [{
    "env": {
      "temperature_K": 300.0,
      "pressure_Pa": 101325.0,
      "gas_mixture": [
        {
          "species": "N2",
          "mole_fraction": 0.78,
          "participates_in_collisions": true,
          "participates_in_reactions": true
        },
        {
          "species": "O2",
          "mole_fraction": 0.21,
          "participates_in_collisions": true,
          "participates_in_reactions": true
        },
        {
          "species": "Ar",
          "mole_fraction": 0.01,
          "participates_in_collisions": true,
          "participates_in_reactions": false
        }
      ]
    }
  }]
}
```

**Key Features:**

- **Mole fractions** must sum to 1.0 (validated at load time)
- **Per-gas flags** control whether each gas participates in collisions/reactions
- **Partial densities** computed automatically: `n_i = x_i × n_total`
- **Backward compatible:** Single-gas configs still work (uses `gas_species` field)

### Per-Gas Participation Flags

Use flags to exclude trace species from physics calculations:

```json
{
  "gas_mixture": [
    {
      "species": "N2",
      "mole_fraction": 0.999,
      "participates_in_collisions": true,
      "participates_in_reactions": true
    },
    {
      "species": "VOC_trace",
      "mole_fraction": 0.001,
      "participates_in_collisions": false,  // ✅ Ignore in mobility
      "participates_in_reactions": true     // ✅ But allow reactions
    }
  ]
}
```

**Use cases:**

- **VOC detection:** Trace species (ppm-level) participate in reactions but don't affect ion mobility
- **Dopant gases:** Small amounts of dopant (e.g., 0.1% acetone) for charge transfer without collision effects
- **Inert carriers:** Ar/He carrier gas participates in collisions but not reactions

### Multi-Gas Collision Physics

ICARION computes collision rates per gas component and selects the colliding gas stochastically:

**Algorithm:**

1. Compute collision rate for each gas: `k_i = n_i × σ_i × v_eff`
2. Total collision rate: `k_total = Σ k_i`
3. Collision probability: `P = 1 - exp(-k_total × dt)`
4. Select gas with probability: `P(gas_i) = k_i / k_total`
5. Execute collision with selected gas

**Example:**

For 78% N₂ / 21% O₂ / 1% Ar mixture:

- If σ_N2 = 100 Ų, σ_O2 = 120 Ų, σ_Ar = 80 Ų
- Then P(N₂) ≈ 0.75, P(O₂) ≈ 0.24, P(Ar) ≈ 0.01

**Statistically correct:** Collision frequencies match experimental drift tube measurements.

### Cross-Section Data Requirements

Multi-gas collisions require gas-specific cross-sections (CCS). ICARION provides **automatic fallback** with multiple tiers:

#### Option 1: Precomputed CCS Maps (RECOMMENDED)

Use the `ccs_precompute` tool to generate gas-specific CCS for all common gases:

```bash
./ccs_precompute \
    --input data/species_database_v1.json \
    --output data/species_database_enriched.json \
    --species H3O+ \
    --ref-gas He \
    --ref-ccs-A2 110.0 \
    --model HSS \
    --override
```

**Output:** Adds CCS maps for He, N₂, O₂, Ar, CO₂, Ne, H₂O to the species entry:

```json
{
  "H3O+": {
    "mass_u": 19.0,
    "charge": 1,
    "CCS_m2": 110e-20,
    "ccs_reference_gas": "He",
    "CCS_HSS": {
      "He": 110e-20,
      "N2": 127e-20,
      "O2": 131e-20,
      "Ar": 141e-20,
      "CO2": 155e-20,
      "Ne": 118e-20,
      "H2O": 134e-20
    },
    "CCS_EHSS": {
      "He": 108e-20,
      "N2": 125e-20,
      "O2": 129e-20,
      "Ar": 138e-20,
      "CO2": 152e-20,
      "Ne": 116e-20,
      "H2O": 131e-20
    }
  }
}
```

**Accuracy:** ±5% for spherical ions, ±10% for complex ions (validated against experimental data)

#### Option 2: Automatic CCS Derivation (FALLBACK)

If CCS maps are missing, ICARION automatically derives them using the Hard-Sphere Scattering (HSS) model:

**Formula:**
```
σ_target = π (r_ion + r_target)²
r_ion = sqrt(σ_ref / π) - r_ref
```

**Example:**

If you provide:
- Reference CCS in He: `CCS_m2 = 110e-20`
- Reference gas: `ccs_reference_gas = "He"`

ICARION automatically computes CCS for N₂, O₂, Ar, etc.

**Log output:**
```
[EHSS] Derived CCS for H3O+:N2 = 127.34 Å² (from He reference).
       For better accuracy: ccs_precompute --species H3O+ --ref-gas He --ref-ccs-A2 110.0
```

**Accuracy:** ±10-15% (acceptable for initial simulations, improved with precomputed data)

**Benefits:**
- ✅ No user action required (automatic fallback)
- ✅ Clear warnings guide optimization
- ✅ Maintains SSOT principle (reference CCS is source of truth)

#### Option 3: Reference CCS (NOT RECOMMENDED)

If no reference gas is specified, ICARION falls back to the single reference CCS for all gases:

**Warning:**
```
[EHSS] Using reference CCS (110.0 Å²) for gas O2 - may be inaccurate (17% error)!
       Run: ccs_precompute --species H3O+ --ref-gas He --ref-ccs-A2 110.0
```

**Error:** ~17% error for O₂ collisions (He CCS used incorrectly)

#### CCS Lookup Priority

ICARION uses this 4-tier hierarchy:

1. ✅ **Precomputed CCS map** (`CCS_HSS[gas]` or `CCS_EHSS[gas]`) - BEST
2. ✅ **Runtime geometry** (EHSS only: OAPA projection from atom coordinates) - ACCURATE but slow
3. ✅ **Automatic derivation** (HSS formula from reference CCS) - ACCEPTABLE (~10% error)
4. ⚠️ **Reference CCS fallback** (same CCS for all gases) - INACCURATE (~17% error)
5. ❌ **Error** (no CCS data available)

**Recommendation:** Always run `ccs_precompute` for production simulations.

### Multi-Gas Reaction Rates

Reaction rates are computed per gas component using partial densities:

**Example reaction:**
```json
{
  "reactions": [
    {
      "id": "rxn_with_N2",
      "reactant": "A+",
      "product": "B+",
      "rate_constant": 1e-10,
      "order": [
        {
          "species": "N2",
          "exponent": 1,
          "concentration_m3": -1.0  // ✅ Auto-lookup from gas_mixture
        }
      ]
    },
    {
      "id": "rxn_with_O2",
      "reactant": "A+",
      "product": "C+",
      "rate_constant": 5e-10,
      "order": [
        {
          "species": "O2",
          "exponent": 1,
          "concentration_m3": -1.0  // ✅ Auto-lookup from gas_mixture
        }
      ]
    }
  ]
}
```

**Competing channels:**

For 78% N₂ / 21% O₂ mixture at 300 K, 1 atm:

- `n_N2 = 0.78 × 2.45e25 m⁻³ = 1.91e25 m⁻³`
- `n_O2 = 0.21 × 2.45e25 m⁻³ = 5.15e24 m⁻³`

Reaction rates:
- `k₁ = 1e-10 m³/s × 1.91e25 m⁻³ = 1.91e15 s⁻¹` (N₂ channel)
- `k₂ = 5e-10 m³/s × 5.15e24 m⁻³ = 2.58e15 s⁻¹` (O₂ channel)

Product branching ratio:
- `P(B+) / P(C+) = k₁ / k₂ = 0.74` → **42% B+, 58% C+**

**Gas name normalization:**

ICARION automatically normalizes gas names to prevent mismatches:

- `"N_2"`, `"N2"`, `"n2"` → all normalized to `"n2"` ✅
- `"O_2"`, `"O2"`, `"o2"` → all normalized to `"o2"` ✅

**This prevents reaction lookup failures due to inconsistent naming!**

### Complete Multi-Gas Example

```json
{
  "title": "Air Mixture Collision and Reaction Test",
  
  "simulation": {
    "total_time_s": 0.001,
    "dt_s": 1e-9,
    "integrator": "RK4"
  },
  
  "physics": {
    "collision_model": "HSS",
    "enable_reactions": true
  },
  
  "species_database": "data/species_database_enriched.json",
  "reaction_database": "data/reactions_database_v1.json",
  
  "domains": [{
    "name": "air_region",
    "instrument": "IMS",
    
    "geometry": {
      "origin_m": [0, 0, 0],
      "length_m": 0.05,
      "radius_m": 0.01
    },
    
    "env": {
      "temperature_K": 300.0,
      "pressure_Pa": 101325.0,
      "gas_mixture": [
        {
          "species": "N2",
          "mole_fraction": 0.78,
          "participates_in_collisions": true,
          "participates_in_reactions": true
        },
        {
          "species": "O2",
          "mole_fraction": 0.21,
          "participates_in_collisions": true,
          "participates_in_reactions": true
        },
        {
          "species": "Ar",
          "mole_fraction": 0.01,
          "participates_in_collisions": true,
          "participates_in_reactions": false
        }
      ]
    },
    
    "fields": {
      "DC": {
        "EN_Td": 10.0
      }
    }
  }],
  
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 100,
        "position": {
          "type": "point",
          "value": [0, 0, 0.001]
        }
      }
    ]
  },
  
  "output": {
    "folder": "./results/multi_gas",
    "trajectory_file": "trajectories.h5",
    "sampling": {
      "interval_steps": 100
    }
  }
}
```

### Validation and Testing

**Test multi-gas configurations:**

```bash
# Validate JSON schema
python3 schema/validator.py schema/icarion-config.schema.json examples/ims/ims_multi_gas_air.json

# Run unit tests
cd build && ctest -R MultiGas

# Check collision statistics
./build/src/icarion_main examples/ims/ims_multi_gas_air.json
```

**Expected test results:**
- `MultiGasCollision` - 6 tests for gas-specific CCS lookup and weighted selection
- `MultiGasReaction` - 1 test for competing reaction channels

**Collision statistics output:**
```
=== Collision Statistics ===
Total collisions: 8523
  - N2: 6641 (77.9%)
  - O2: 1795 (21.1%)
  - Ar: 87 (1.0%)
```

**Matches expected mole fractions!**

### Migration from Single-Gas Configs

**Old config (single gas):**
```json
{
  "env": {
    "pressure_Pa": 101325.0,
    "temperature_K": 300.0,
    "gas_species": "N2"
  }
}
```

**New config (equivalent multi-gas):**
```json
{
  "env": {
    "pressure_Pa": 101325.0,
    "temperature_K": 300.0,
    "gas_mixture": [
      {
        "species": "N2",
        "mole_fraction": 1.0
      }
    ]
  }
}
```

**No breaking changes:** Old configs with `gas_species` still work (automatic single-gas fallback).

### Performance Notes

**Multi-gas overhead:**
- CPU: ~5-10% overhead for 3-component mixtures (negligible)
- GPU: Same (gas selection done on GPU)
- Memory: +8 bytes per gas component per timestep (negligible)

**Recommendation:** Use multi-gas freely - performance impact is minimal.

---

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

### Ion Mobility Spectrometry (IMS)

```json
{
  "simulation": {
    "total_time_s": 5e-4,
    "dt_s": 1e-8,
    "integrator": "RK4",
    "enable_gpu": true
  },
  "physics": {
    "collision_model": "HSS"
  },
  "output": {
    "folder": "./results/ims",
    "trajectory_file": "ims_trajectories.h5"
  },
  "domains": [
    {
      "name": "drift_region",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.05,
        "radius_m": 0.015
      },
      "env": {
        "pressure_Pa": 200.0,
        "temperature_K": 300.0,
        "gas_species": "He"
      },
      "fields": {
        "DC": {
          "EN_Td": 10.0
        }
      }
    }
  ]
}
```

### Time-of-Flight (TOF)

```json
{
  "simulation": {
    "total_time_s": 1e-4,
    "dt_s": 1e-10,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "NoCollisions"
  },
  "output": {
    "folder": "./results/tof",
    "trajectory_file": "tof_trajectories.h5"
  },
  "domains": [
    {
      "name": "acceleration_region",
      "instrument": "TOF",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.02,
        "radius_m": 0.01
      },
      "env": {
        "pressure_Pa": 1e-6,
        "temperature_K": 300.0,
        "gas_species": "He"
      },
      "fields": {
        "DC": {
          "axial_V": 5000.0
        }
      }
    },
    {
      "name": "drift_region",
      "instrument": "TOF",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.02],
        "length_m": 0.5,
        "radius_m": 0.01
      },
      "env": {
        "pressure_Pa": 1e-6,
        "temperature_K": 300.0,
        "gas_species": "He"
      }
    }
  ]
}
```

### Multi-Domain Setup

```json
{
  "simulation": {
    "total_time_s": 1e-3,
    "dt_s": 1e-9,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "HSS"
  },
  "output": {
    "folder": "./results/hybrid",
    "trajectory_file": "trajectories.h5"
  },
  "domains": [
    {
      "name": "ims_drift",
      "instrument": "IMS",
      "integrator": "RK4",
      "geometry": { "origin_m": [0, 0, 0], "length_m": 0.05, "radius_m": 0.015 },
      "env": { "pressure_Pa": 200.0, "temperature_K": 300.0, "gas_species": "He" },
      "fields": { "DC": { "EN_Td": 10.0 } }
    },
    {
      "name": "tof_region",
      "instrument": "TOF",
      "integrator": "RK45",
      "geometry": { "origin_m": [0, 0, 0.05], "length_m": 0.5, "radius_m": 0.02 },
      "env": { "pressure_Pa": 1e-6, "temperature_K": 300.0, "gas_species": "He" },
      "fields": { "DC": { "axial_V": 5000.0 } }
    }
  ]
}
```

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
      "env": { "pressure_Pa": 1e-6, ... }
    },
    {
      "name": "high_pressure_region",
      // No integrator specified → uses "RK4" from simulation
      "env": { "pressure_Pa": 101325, ... }
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

5. **Enable GPU** for large simulations (`enable_gpu: true`)

---

## See Also

- **Schema Documentation**: [`../schema/README.md`](../schema/README.md)
- **CLI Reference**: [`CLI_USAGE.md`](CLI_USAGE.md)
- **Architecture**: [`ARCHITECTURE.md`](ARCHITECTURE.md)
- **Example Configs**: [`../examples/`](../examples/)
- **Python Helper**: `scripts/create_config.py --help`

---

## Version

This guide corresponds to ICARION v1.0 configuration schema.

For schema version history, see the git log for the `schema/` directory.
