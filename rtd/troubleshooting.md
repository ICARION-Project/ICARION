# Troubleshooting

This page lists common first run problems and the fastest checks.

## Binary not found

If `icarion` is not found, either install a release package or use the build tree
binary:

```bash
./build/src/icarion_main examples/ims/ims_basic.json
```

After installing a Debian package, verify:

```bash
icarion --version
icarion --dump-build-info
```

## Configuration fails to load

Run the checks in increasing depth:

```bash
icarion --validate-schema config.json
icarion --validate-config config.json
icarion --dry-run config.json
```

Typical causes:

- JSON syntax errors, especially trailing commas.
- Misspelled top level keys such as `domains`, `simulation`, or `physics`.
- Lowercase field sections; use `DC`, `RF`, and `AC`.
- Database paths that are correct from your editor but wrong from the current
  working directory.

## Species or reactions are missing

Check the configured database paths:

```json
"species_database": "data/species_database.json",
"reaction_database": "data/reaction_database.json"
```

The species IDs used by `ions[].species_id`, collision data, and reactions must
match the database entries. If reactions are enabled, product species must also
exist in the species database.

For a minimal first test, disable reactions:

```bash
icarion --no-reactions config.json
```

or:

```bash
icarion --set physics.enable_reactions=false config.json
```

## Output file is too large

Trajectory size scales with ion count, number of written steps, and whether full
trajectory datasets are enabled.

Reduce output volume by:

- increasing `simulation.write_interval`
- shortening `simulation.total_time_s`
- reducing the number of ions
- using `output.trajectory_mode: "minimal"` when detailed trajectories are not
  needed
- setting `--buffer-byte-cap` for long runs

See [Output files](output-files.md) before designing large parameter sweeps.

## Expected arrays are missing

If a trajectory dataset is missing, check `output.trajectory_mode`. Minimal mode
intentionally omits some high volume arrays. Also confirm that the simulation
actually reached a write step and did not fail during initialization.

## Ions disappear immediately

Common causes:

- Initial positions are outside the domain.
- The ion cloud is wider than the domain radius or length.
- The timestep is too large for the local field or geometry.
- The boundary is absorptive and the field accelerates ions into a wall.
- The domain origin/length does not match the intended coordinate system.

Debug with fewer ions, a shorter total time, a smaller `simulation.dt_s`, and
verbose logging:

```bash
icarion --verbose --log-file debug.log config.json
```

## Collision results look wrong

Start by confirming the selected global collision model:

```bash
icarion --list-collision-models
icarion --set physics.collision_model=EHSS config.json
```

Important caveats:

- `physics.collision_model` applies globally to all domains.
- HSS/EHSS need reliable species collision parameters or molecular geometries.
- Stochastic collision models require timestep convergence checks.
- Deterministic models and stochastic models answer different modeling
  questions.

See [Collision models](collision-models.md) for model-specific guidance.

## GPU is not used

`simulation.enable_gpu` only requests GPU execution. The binary must be built
with CUDA support, and unsupported paths can fall back to CPU behavior. Check:

```bash
icarion --dump-build-info
```

For v1.0.1 onboarding and validation, prefer the CPU path unless you explicitly
need to evaluate GPU support.
