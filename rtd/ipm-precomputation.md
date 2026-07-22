# Interaction-potential precomputation

## Purpose

`interaction_potential_precompute` integrates classical ion-neutral scattering trajectories offline over molecular orientations, relative velocity bins, impact parameters, and stochastic trajectory samples. At runtime, `InteractionPotentialModel` loads the resulting HDF5 table through the species database field `ipm_samples_file` and samples collision events without reintegrating the molecular interaction potential for every event.

This separates an offline molecular scattering calculation from runtime event-level collision sampling. Neither step by itself is a complete, time-resolved IMS/MS instrument simulation; `icarion` supplies that instrument simulation using the precomputed table.

## Minimal workflow

```bash
cmake --build build --target interaction_potential_precompute

build/src/interaction_potential_precompute \
  --input data/species_database_v1.json \
  --species H3O+ \
  --gas He \
  --output H3O+_He_ipm.h5 \
  --potential lj1264 \
  --seed 12345
```

Reference the result from the resolved ion entry:

```json
{
  "species": {
    "H3O+": {
      "mass_amu": 19.02,
      "charge": 1,
      "ipm_samples_file": "molecules/precomputed_ipm/H3O+_He_ipm.h5"
    }
  }
}
```

Use an explicit seed for publication calculations.

Optional context can be attached without changing numerical identity:

```bash
interaction_potential_precompute \
  --input species.json --species H3O+ --gas He \
  --output H3O+_He_ipm.h5 --note-file ipm_notes.md
```

The exact note bytes, `inline`/`file` source, basename only file provenance, and SHA-256 are stored under `/metadata/annotations`. Notes are limited to 65,536 UTF-8 bytes and may contain identifying text only when supplied by the user. They complement rather than replace resolved inputs, build/RNG provenance, and completion state.

## Numerical format

The runtime contract remains:

```text
format = ipm_offline_samples
version = 1
```

The `/metadata` hierarchy is additive: it neither increments the numerical format version nor changes the existing loader contract.

| Dataset | Purpose |
|---|---|
| `logv_bins` | Logarithmic relative-velocity grid. |
| `orientations_quat` | Molecular orientations as quaternions. |
| `sigma_event_m2` | Event-rate cross section used by the runtime. |
| `sigma_mt_m2` | Diagnostic momentum-transfer cross section. It is not applied again to stored momentum kicks. |
| `b_max_m` | Sampled impact parameter limit. |
| `attempted_trajectories`, `accepted_trajectories`, `rejected_non_asymptotic` | Per-cell trajectory accounting. |
| `max_energy_step_error`, `max_energy_cumulative_error` | Integration diagnostics. |
| `cdf_offsets`, `cdf_counts`, `cdf_values`, `dp_samples` | Full empirical momentum-transfer distributions. |
| `dp_stats` | Lower fidelity compact statistics fallback. |

Full-CDF tables are recommended for scientific production. Compact `dp_stats` files do not preserve the full joint momentum-transfer distribution.

## Reproducibility metadata

```text
/metadata
├── schema
├── software
├── system
├── rng
├── species
├── neutral
├── precompute
├── inputs
│   ├── hashes
│   └── blobs
├── annotations     # optional
└── completion
```

| Group | Contents |
|---|---|
| `/metadata/schema` | Metadata and numerical-format identity. |
| `/metadata/software` | ICARION, source, compiler, build, OpenMP, and CUDA provenance. |
| `/metadata/system` | Non-identifying operating system and hardware context. |
| `/metadata/rng` | Resolved seed and sampling scheme. |
| `/metadata/species` | Resolved ion definition and molecular geometry. |
| `/metadata/neutral` | Final gas model and parameters with source labels. |
| `/metadata/precompute` | Fully resolved numerical and model settings. |
| `/metadata/inputs` | Input filenames, SHA-256 hashes, and immutable embedded snapshots. |
| `/metadata/annotations` | Optional exact-byte user note, source, basename, and SHA-256. |
| `/metadata/completion` | Progress, checkpoint, resume, and timing state. |

### `/metadata/schema`

`metadata_schema_version = 1.0.0`, `data_format = ipm_offline_samples`, and `data_format_version = 1`. The metadata schema may evolve independently of the numerical IPM runtime format.

### `/metadata/software`

This group records the ICARION version, Git hash, Git state (`clean`, `dirty`, or `unknown`), Git state capture mode (`configure_time`), exact CMake build configuration, compiler and version, C++ standard, OpenMP status and effective thread count, CUDA build status, and full or core-only build mode. Git state is captured when CMake configures the build and therefore does not necessarily reflect later changes made without rerunning CMake.

### `/metadata/system`

The operating system, kernel, logical CPU count, CPU model and total memory where available, and metadata timestamp are recorded. IPM publication files intentionally omit username, hostname, and absolute input paths.

### `/metadata/rng`

The file stores the actual `base_seed`, whether it was explicitly supplied, its origin, `std::mt19937_64`, the orientation sampling mode, and the textual `splitmix64` per-cell derivation from base seed, orientation index, and velocity bin index. Automatic seeds use an implementation-dependent `std::hash`; the actual resolved value is nevertheless stored. Use an explicit seed for a strict cross-platform study definition.

### `/metadata/species` and `/metadata/neutral`

These groups describe the resolved calculation, not just source filenames. Species metadata includes ID, charge, mass in u and kg, reference temperature, selected database entry, geometry filename, SHA-256 and embedded content, atom count, and geometry-derived mass. Neutral metadata includes gas ID and model, mass in u and kg, final sigma, epsilon and polarizability, applicable anisotropic N2 polarizabilities, bond and quadrupole parameters. Stored source labels use exactly `built_in`, `parameter_file`, `cli_override`, `derived_isotropic_fallback`, or `not_used`.

### `/metadata/precompute`

This group stores values after defaults and overrides: potential and parameter models, mixing rule, polarization model, orientation grid and count, velocity range and bin count, trials, temperature, impact parameter controls, integration tolerances and maximum steps, scale factors, full-CDF/compact mode, checkpoint interval, requested resume mode, effective OpenMP thread count, a sanitized command line, and `resolved_options_json`. Filenames in the command line are reduced to basenames so absolute paths are not disclosed.

### `/metadata/inputs`

Required text inputs are captured as immutable bytes before numerical work. The same bytes are parsed, hashed with SHA-256, embedded, and reused for every checkpoint and final write. Records can cover the species database, molecular geometry, gas and element parameter files, and a Lebedev grid when used. Only basenames are stored; identical files are embedded once and cross-referenced. A required snapshot failure aborts instead of writing misleading provenance, while unused optional inputs are marked as not used.

### `/metadata/completion`

`success`, `is_checkpoint`, completed and total cells, the original start timestamp, current write/completion timestamp, current process wall time, accumulated wall time across resumed processes, and whether resume was actually used describe file state. A checkpoint is intentionally incomplete and must not be interpreted as a completed production table.

## Checkpoint and resume workflow

Compact output can periodically replace the output with a checkpoint and later resume it:

```bash
build/src/interaction_potential_precompute \
  --input data/species_database_v1.json --species H3O+ --gas He \
  --output H3O+_He_compact.h5 --compact-dp-stats --checkpoint-cells 10 --seed 12345

build/src/interaction_potential_precompute \
  --input data/species_database_v1.json --species H3O+ --gas He \
  --output H3O+_He_compact.h5 --compact-dp-stats --checkpoint-cells 10 --seed 12345 --resume
```

Resume validates numerical settings, model choices, seed, velocity and orientation grids, and every used input hash before accepting completed cells. A mismatch aborts rather than mixing calculations.

Annotations are immutable across checkpoints. Resume without flags preserves an existing annotation. A supplied annotation must match the checkpoint content exactly, while the checkpoint's original source and filename are retained. Adding a note to an unannotated checkpoint, changing or removing one, or encountering malformed/hash-mismatched annotation metadata aborts with exit code 1.

`--stop-after-checkpoint` is only a checkpoint-workflow testing option. It requires compact checkpoint mode, deliberately exits after the first incomplete checkpoint, and returns exit code 2. It is not a successful production completion.

## Inspecting a file

```bash
h5dump -n H3O+_He_ipm.h5
h5dump -d /metadata/rng/base_seed H3O+_He_ipm.h5
h5dump -d /metadata/completion/success H3O+_He_ipm.h5
```

```python
import h5py

with h5py.File("H3O+_He_ipm.h5", "r") as f:
    print(f.attrs["format"], f.attrs["version"])
    print(f["metadata/schema/metadata_schema_version"][()].decode())
    print(f["metadata/rng/base_seed"][()])
    print(f["metadata/inputs/hashes/species_database_sha256"][()].decode())
    print(bool(f["metadata/completion/success"][()]),
          bool(f["metadata/completion/is_checkpoint"][()]))
    print(f["metadata/inputs/blobs/species_database"][()].decode())
```

## Publication checklist

- Use an explicit seed.
- Prefer full-CDF mode unless compact mode is scientifically justified.
- Require acceptable, preferably zero, non-asymptotic rejects.
- Archive the complete `.h5` file unchanged.
- Record the ICARION release or commit.
- Inspect `success` and `is_checkpoint`.
- Retain the generated file instead of rebuilding it from undocumented inputs.
- Report major potential, parameter, gas, sampling, and convergence settings.

See [Collision models](collision-models.md), [Species database](species-database.md), [Output files](output-files.md), [Related and complementary software](related-software.md), and [Validation](validation.md).
