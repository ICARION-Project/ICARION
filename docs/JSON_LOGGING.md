# JSON Logging

## Overview

ICARION supports structured logging in JSON format for automated analysis with pandas, performance dashboards, and debugging workflows. Instead of colored console output, you get machine-readable log lines that are easy to filter and analyze.

**Flag:** `--log-format json`

Each log entry is a single-line JSON object with `time`, `level`, `cat` (category), and `msg`:

```json
{"time":"2025-11-21T10:30:52.704","level":"info","cat":"main","msg":"ICARION v1.0.0 starting"}
```

## Generate JSON Logs

```bash
# Console output
./icarion_main --log-format json examples/ims_basic.json

# Save to file (workaround for flush issue)
./icarion_main --log-format json examples/ims_basic.json 2>&1 | grep '^{' > simulation.log

# Traditional text format (default)
./icarion_main --log-format text examples/ims_basic.json
```

## Analyze with Python

### Basic Analysis

```bash
python examples/analyze_logs.py simulation.log
```

### Advanced Analysis

```bash
python examples/analyze_logs_advanced.py simulation.log
```

## pandas Quick Snippets

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

### Time Analysis

```python
# Sort by time
df = df.sort_values('time')

# Compute time between logs
df['delta_ms'] = df['time'].diff().dt.total_seconds() * 1000

# Find slowest operations
slowest = df.nlargest(10, 'delta_ms')[['cat', 'msg', 'delta_ms']]
print(slowest)
```

### Pattern Matching

```python
# Find species mentions
species_logs = df[df['msg'].str.contains('species', case=False)]

# Find particle counts
particle_logs = df[df['msg'].str.contains(r'\d+\s+particles?', case=False)]

# Extract numeric values
df['particles'] = df['msg'].str.extract(r'(\d+)\s+particles?')[0].astype(float)
```

### Export Results

```python
# CSV
df.to_csv('analysis.csv', index=False)

# Excel with sheets
with pd.ExcelWriter('analysis.xlsx') as writer:
    df.to_excel(writer, sheet_name='All')
    for cat in df['cat'].unique():
        df[df['cat']==cat].to_excel(writer, sheet_name=cat)

# JSON summary
summary = {
    'total': len(df),
    'categories': df['cat'].value_counts().to_dict(),
    'levels': df['level'].value_counts().to_dict()
}
with open('summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
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

### Error Report

```python
errors = df[df['level'] == 'error']
if len(errors) > 0:
    print(f"Found {len(errors)} errors:")
    for _, row in errors.iterrows():
        print(f"  [{row['time']}] {row['cat']}: {row['msg']}")
else:
    print("✓ No errors")
```

### Category Timeline

```python
import matplotlib.pyplot as plt

df['seconds'] = (df['time'] - df['time'].min()).dt.total_seconds()

for cat in df['cat'].unique():
    cat_df = df[df['cat'] == cat]
    plt.scatter(cat_df['seconds'], [cat]*len(cat_df), label=cat, alpha=0.6)

plt.xlabel('Time (seconds)')
plt.ylabel('Category')
plt.legend()
plt.tight_layout()
plt.savefig('timeline.png')
```

## Requirements

```bash
pip install pandas
pip install openpyxl  # Optional: for Excel export
pip install matplotlib  # Optional: for plotting
```

## Tips

1. **Large logs:** Use `df.head()` and `df.tail()` to preview
2. **Memory:** Process logs in chunks with `pd.read_json(chunksize=1000)`
3. **Search:** Use `df[df['msg'].str.contains('pattern')]` for text search
4. **Debugging:** Filter by category + time range for focused analysis
5. **Performance:** Parse JSON once, cache DataFrame with `df.to_pickle()`

## See Also

- `examples/analyze_logs.py` - Basic analysis script
- `examples/analyze_logs_advanced.py` - Advanced examples
- `docs/JSON_LOGGING_IMPLEMENTATION.md` - Full documentation
- `examples/README.md` - Usage examples
