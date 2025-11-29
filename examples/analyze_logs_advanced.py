"""
ICARION Log Analysis - Advanced Examples
=========================================

Demonstrates advanced log analysis techniques using pandas.

This script shows how to:
- Filter logs by category and level
- Compute time-based statistics
- Extract performance metrics
- Visualize logging patterns

Requirements:
    pip install pandas matplotlib seaborn

Usage:
    python analyze_logs_advanced.py simulation.log
"""

import pandas as pd
import sys
from pathlib import Path
import json

def load_json_logs(log_file: Path) -> pd.DataFrame:
    """Load ICARION JSON logs into DataFrame."""
    logs = []
    with open(log_file, 'r') as f:
        for line in f:
            if line.strip() and line.strip().startswith('{'):
                try:
                    logs.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    
    df = pd.DataFrame(logs)
    df['time'] = pd.to_datetime(df['time'])
    return df


def example_filtering(df: pd.DataFrame):
    """Example 1: Filter logs by category and level."""
    print("\n=== Example 1: Filtering ===")
    
    # Get all configuration logs
    config_logs = df[df['cat'] == 'config']
    print(f"Configuration logs: {len(config_logs)}")
    print(config_logs[['time', 'msg']].to_string(index=False))
    
    # Get all errors and warnings
    issues = df[df['level'].isin(['error', 'warn'])]
    if not issues.empty:
        print(f"\nFound {len(issues)} issues:")
        print(issues[['level', 'cat', 'msg']].to_string(index=False))
    else:
        print("\n✓ No errors or warnings")


def example_time_analysis(df: pd.DataFrame):
    """Example 2: Time-based analysis."""
    print("\n=== Example 2: Time Analysis ===")
    
    # Compute time deltas between log entries
    df = df.sort_values('time').copy()
    df['time_delta'] = df['time'].diff().dt.total_seconds()
    
    # Find slowest operations
    slowest = df.nlargest(5, 'time_delta')[['cat', 'msg', 'time_delta']]
    print("\nSlowest operations:")
    print(slowest.to_string(index=False))
    
    # Group by category and compute stats
    print("\nTime per category:")
    category_times = df.groupby('cat')['time_delta'].agg(['count', 'sum', 'mean'])
    print(category_times)


def example_pattern_matching(df: pd.DataFrame):
    """Example 3: Extract specific patterns from messages."""
    print("\n=== Example 3: Pattern Matching ===")
    
    # Extract species information
    species_logs = df[df['msg'].str.contains('species', case=False, na=False)]
    print(f"\nSpecies-related logs: {len(species_logs)}")
    for msg in species_logs['msg'].head(5):
        print(f"  - {msg}")
    
    # Extract domain information
    domain_logs = df[df['msg'].str.contains('domain', case=False, na=False)]
    print(f"\nDomain-related logs: {len(domain_logs)}")
    for msg in domain_logs['msg'].head(5):
        print(f"  - {msg}")
    
    # Extract numeric values (e.g., particle counts)
    particle_logs = df[df['msg'].str.contains(r'\d+\s+particles?', case=False, na=False)]
    if not particle_logs.empty:
        print(f"\nParticle count logs: {len(particle_logs)}")
        for msg in particle_logs['msg'].head(3):
            print(f"  - {msg}")


def example_category_breakdown(df: pd.DataFrame):
    """Example 4: Detailed category breakdown."""
    print("\n=== Example 4: Category Breakdown ===")
    
    categories = df['cat'].unique()
    for cat in sorted(categories):
        cat_df = df[df['cat'] == cat]
        print(f"\n{cat.upper()} ({len(cat_df)} entries):")
        
        # Show first and last message
        if len(cat_df) > 0:
            print(f"  First: {cat_df.iloc[0]['msg'][:80]}")
            if len(cat_df) > 1:
                print(f"  Last:  {cat_df.iloc[-1]['msg'][:80]}")
        
        # Level distribution
        levels = cat_df['level'].value_counts()
        if len(levels) > 1:
            print(f"  Levels: {dict(levels)}")


def example_export_formats(df: pd.DataFrame, base_path: Path):
    """Example 5: Export to various formats for external tools."""
    print("\n=== Example 5: Export Options ===")
    
    # CSV for spreadsheets
    csv_file = base_path.with_suffix('.csv')
    df.to_csv(csv_file, index=False)
    print(f"✓ CSV: {csv_file}")
    
    # Excel with multiple sheets (if openpyxl available)
    try:
        excel_file = base_path.with_suffix('.xlsx')
        with pd.ExcelWriter(excel_file) as writer:
            df.to_excel(writer, sheet_name='All Logs', index=False)
            for cat in df['cat'].unique():
                cat_df = df[df['cat'] == cat]
                cat_df.to_excel(writer, sheet_name=cat, index=False)
        print(f"✓ Excel: {excel_file}")
    except ImportError:
        print("  (Excel export requires: pip install openpyxl)")
    
    # JSON summary
    summary = {
        'total_entries': len(df),
        'categories': df['cat'].value_counts().to_dict(),
        'levels': df['level'].value_counts().to_dict(),
        'time_range': {
            'start': df['time'].min().isoformat(),
            'end': df['time'].max().isoformat(),
            'duration_seconds': (df['time'].max() - df['time'].min()).total_seconds()
        }
    }
    json_file = base_path.with_name(base_path.stem + '_summary.json')
    with open(json_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"✓ Summary JSON: {json_file}")


def main():
    if len(sys.argv) != 2:
        print("Usage: python analyze_logs_advanced.py <json_log_file>")
        sys.exit(1)
    
    log_file = Path(sys.argv[1])
    if not log_file.exists():
        print(f"Error: File not found: {log_file}")
        sys.exit(1)
    
    print(f"Loading {log_file}...")
    df = load_json_logs(log_file)
    print(f"Loaded {len(df)} log entries\n")
    
    # Run all examples
    example_filtering(df)
    example_time_analysis(df)
    example_pattern_matching(df)
    example_category_breakdown(df)
    example_export_formats(df, log_file)
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)


if __name__ == '__main__':
    main()
