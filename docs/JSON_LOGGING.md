# JSON Logging

## Overview

ICARION supports structured logging in JSON format for automated analysis with pandas, performance dashboards, and debugging workflows. Instead of colored console output, you get machine-readable log lines that are easy to filter and analyze.

**Flag:** `--log-format json`  
**Executable:** `icarion_main` (in your build directory, e.g., `./build/src/icarion_main`)

Each log entry is a single-line JSON object with `time`, `level`, `cat` (category), and `msg`:

```json
{"time":"2025-11-21T10:30:52.704","level":"info","cat":"main","msg":"ICARION v1.0.0 starting"}
```

## Generate JSON Logs

```bash
# Console output (JSON format)
./build/src/icarion_main --log-format json examples/ims/ims_basic.json

# Save to file
./build/src/icarion_main --log-format json --log-file simulation.log examples/ims/ims_basic.json

# Traditional text format (default)
./build/src/icarion_main examples/ims/ims_basic.json
```

## Analyze with Python (pandas)

Create your own analysis scripts using pandas:

### Load Logs

```python
import pandas as pd
import json

logs = []
with open('simulation.log') as f:
    for line in f:
        if line.strip().startswith('{'):
            logs.append(json.loads(line))

df = pd.DataFrame(logs)
df['time'] = pd.to_datetime(df['time'])
```

### Filter by Category

```python
# Configuration logs
config_df = df[df['cat'] == 'config']

# Performance logs
perf_df = df[df['cat'] == 'perf']

# All categories
print(df['cat'].unique())
# ['main', 'config', 'integrator', 'hdf5', 'physics', 'perf', 'gpu', 'reactions']
```

### Filter by Level

```python
# Errors only
errors = df[df['level'] == 'error']

# Warnings and errors
issues = df[df['level'].isin(['warn', 'error'])]

# Debug logs
debug = df[df['level'] == 'debug']
```

### Time Analysis & Pattern Matching

```python
# Sort and compute time between logs
df = df.sort_values('time')
df['delta_ms'] = df['time'].diff().dt.total_seconds() * 1000

# Find slowest operations
slowest = df.nlargest(10, 'delta_ms')[['cat', 'msg', 'delta_ms']]

# Pattern search
species_logs = df[df['msg'].str.contains('species', case=False)]
df['particles'] = df['msg'].str.extract(r'(\d+)\s+particles?')[0].astype(float)

# Export
df.to_csv('analysis.csv', index=False)
```

## JSON Format

Each log entry is a single-line JSON object:

```json
{
  "time": "2025-11-21T10:30:52.704",
  "level": "info",
  "cat": "main",
  "msg": "ICARION v1.0.0 starting"
}
```

**Fields:**

- `time`: ISO 8601 timestamp with milliseconds
- `level`: debug | info | warn | error
- `cat`: main | config | integrator | hdf5 | physics | perf | gpu | reactions
- `msg`: Log message (string)

## Common Queries

### Configuration Summary

```python
config = df[df['cat'] == 'config']
print("Species:", config[config['msg'].str.contains('species')]['msg'].values[0])
print("Domains:", config[config['msg'].str.contains('domain')]['msg'].values[0])
```

### Error Analysis & Visualization

```python
# Find errors
errors = df[df['level'] == 'error']
print(f"Found {len(errors)} errors" if len(errors) > 0 else "✓ No errors")

# Plot timeline (optional: requires matplotlib)
df['seconds'] = (df['time'] - df['time'].min()).dt.total_seconds()
for cat in df['cat'].unique():
    cat_df = df[df['cat'] == cat]
    plt.scatter(cat_df['seconds'], [cat]*len(cat_df), alpha=0.6)
plt.savefig('timeline.png')
```

## Requirements

```bash
pip install pandas  # Required
pip install matplotlib  # Optional: for plotting
```

## Tips

- Use `df.head()` and `df.tail()` to preview large logs
- Filter by category/level: `df[df['cat'] == 'perf']`
- Search messages: `df[df['msg'].str.contains('pattern')]`
- Cache parsed data: `df.to_pickle('logs.pkl')`

## See Also

- [`CLI_USAGE.md`](CLI_USAGE.md) - Full command-line reference
- [`../examples/README.md`](../examples/README.md) - Example configurations
- [`CONFIG_GUIDE.md`](CONFIG_GUIDE.md) - Configuration file format
