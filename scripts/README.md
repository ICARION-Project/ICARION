# ICARION Scripts

Utility scripts for configuration generation and analysis.

## `create_config.py`

Generate valid ICARION JSON configuration files from templates.

**Quick start:**
```bash
# List templates
python3 scripts/create_config.py --list-templates

# Generate config
python3 scripts/create_config.py --template ims --output my_config.json

# Run simulation
./build/src/icarion_main my_config.json
```

**Available templates:**
- `minimal` - H₃O⁺, 100 ions, NoCollisions
- `ims` - H₃O⁺, 1000 ions, HSS collisions
- `tof` - ReserpineH⁺, 100 ions, vacuum
- `lqit` - CaffeineH⁺, 500 ions, HSS collisions
- `orbitrap` - ReserpineH⁺, 100 ions, vacuum

All templates are v1.0 schema compliant and tested.

- Simulation time and timestep
- Collision model
- Number of domains
- Field configurations

#### 3. Customize Template Parameters

```bash
python3 scripts/create_config.py \
    --template ims \
    --total-time 1e-3 \
    --timestep 1e-9 \
    --collision-model HSS \
    --output configs/my_config.json
```

### Command-Line Options

```
positional arguments:
## `compute_ccs_maps.py`

Precompute collision cross-section (CCS) maps for temperature/pressure ranges.

**Usage:**
```bash
python3 scripts/compute_ccs_maps.py --species H3O+ --temp-range 100 500 --pressure-range 100 1000
```

## Related Documentation

- `docs/CONFIG_GUIDE.md` - Configuration reference
- `examples/` - Ready-to-run examples
