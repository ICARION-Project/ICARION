# ICARION Scripts

Utility scripts for ICARION workflow automation.

---

## 📄 `create_config.py`

**Purpose:** Generate valid ICARION configuration files from templates or interactively.

### Quick Start

```bash
# List available templates
python3 scripts/create_config.py --list-templates

# Create config from template
python3 scripts/create_config.py --template ims --output my_config.json

# Validate generated config
python3 src/core/config/schema/validator.py my_config.json
./build/src/icarion_main --validate-config my_config.json
```

### Available Templates

| Template | Description | Use Case |
|----------|-------------|----------|
| `minimal` | Bare minimum configuration | Quick testing, learning |
| `ims` | Ion Mobility Spectrometry | Drift tube IMS experiments |
| `tof` | Time-of-Flight | TOF mass spectrometry |
| `lqit` | Linear Quadrupole Ion Trap | Paul trap simulations |
| `orbitrap` | Orbitrap Mass Analyzer | High-resolution MS |

### Usage Examples

#### 1. Generate IMS Configuration

```bash
python3 scripts/create_config.py --template ims --output configs/my_ims_run.json
```

#### 2. Interactive Mode (Advanced)

```bash
python3 scripts/create_config.py --interactive --output configs/custom.json
```

Follow the prompts to customize:

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
  output                Output JSON file path

optional arguments:
  -h, --help            Show help message
  --template {minimal,ims,tof,lqit,orbitrap}
                        Use predefined template
  --list-templates      List all available templates
  --interactive         Interactive configuration builder
  --validate            Validate generated config automatically
  
Template customization:
  --total-time FLOAT    Simulation time [s]
  --timestep FLOAT      Integration timestep [s]
  --integrator {RK4,RK45,Boris}
                        Integrator choice
  --collision-model {NoCollisions,HSD,HSS,EHSS,SDS}
                        Collision model
  --enable-gpu          Enable GPU acceleration
  --output-folder PATH  Output directory
```

### Generated Config Features

All generated configs include:

- Valid JSON Schema v1.0 format
- Reasonable default parameters
- Inline comments (via `"_comment"` fields)
- Ready for validation and simulation
- Documented field descriptions

### Workflow Integration

**Recommended workflow:**

```bash
# 1. Create config from template
python3 scripts/create_config.py --template ims --output my_sim.json

# 2. Edit manually (optional)
nano my_sim.json

# 3. Validate
python3 src/core/config/schema/validator.py my_sim.json
./build/src/icarion_main --validate-config my_sim.json

# 4. Run simulation
./build/src/icarion_main my_sim.json
```

### Validation Features

The script automatically checks:

- Required fields present
- Valid enum values
- Positive values for physical quantities
- Consistent domain indices
- Field configuration logic

### Error Handling

Common issues and solutions:

| Error | Cause | Solution |
|-------|-------|----------|
| `ValueError: invalid literal` | Wrong value type | Check input types match schema |
| `ValidationError: missing required` | Incomplete config | Use `--validate` to check |
| `FileNotFoundError` | Invalid path | Check file paths are accessible |

### Advanced: Custom Templates

Add your own templates in `create_config.py`:

```python
TEMPLATES = {
    'my_custom': {
        'simulation': {
            'total_time_s': 1e-3,
            'dt_s': 1e-9,
            # ... your parameters
        }
    }
}
```

---

## 🔧 Future Scripts (Planned)

- `convert_legacy_config.py` - Migrate old configs to v1.0
- `compare_configs.py` - Diff two configuration files
- `optimize_params.py` - Parameter sweeps for sensitivity analysis
- `batch_simulate.py` - Run multiple configs in parallel

---

## 📚 Related Documentation

- [CONFIG_GUIDE.md](../docs/CONFIG_GUIDE.md) - Full configuration reference
- [Schema Documentation](../src/core/config/schema/README.md) - JSON Schema details
- [Examples](../examples/) - Production-ready example configs

---

## 🐛 Bug Reports

Found a bug in the script? Report it with:
```bash
python3 scripts/create_config.py --version
python3 --version
uname -a
```

Include your command and error output in the issue.
