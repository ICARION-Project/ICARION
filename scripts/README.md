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
- `minimal` - H₃O⁺, 100 ions, NoCollisions, basic testing
- `ims` - H₃O⁺, 1000 ions, HSS collisions, drift tube
- `tof` - ReserpineH⁺, 100 ions, vacuum, time-of-flight
- `lqit` - CaffeineH⁺, 500 ions, HSS collisions, Paul trap
- `orbitrap` - ReserpineH⁺, 100 ions, vacuum, hyperlogarithmic field
- `quadrupole` - CaffeineH⁺, 50 ions, RF/DC fields, mass filter
- `fticr` - ReserpineH⁺, 50 ions, 7T magnetic field, Boris integrator

All templates v1.0.0 compliant and tested.

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

## Related Documentation

- `docs/CONFIG_GUIDE.md` - Configuration reference
- `examples/` - Ready-to-run examples
