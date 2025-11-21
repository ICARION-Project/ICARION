#!/usr/bin/env python3
"""
ICARION Log Analysis Example
=============================

Parse JSON-formatted ICARION logs and create analysis dashboard.

Usage:
    # Run simulation with JSON logging
    ./icarion_main --log-format json --log-file simulation.log config.json
    
    # Analyze logs
    python analyze_logs.py simulation.log

Requirements:
    pip install pandas matplotlib

Author: ICARION Development Team
Date: 2025-11-21
License: Apache-2.0
"""

import json
import sys
from pathlib import Path
from typing import List, Dict
import pandas as pd

def parse_json_logs(log_file: Path) -> pd.DataFrame:
    """
    Parse ICARION JSON logs into pandas DataFrame.
    
    Args:
        log_file: Path to JSON log file (one JSON object per line)
    
    Returns:
        DataFrame with columns: time, level, category, message
    """
    logs = []
    
    with open(log_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            try:
                log_entry = json.loads(line)
                logs.append({
                    'time': pd.to_datetime(log_entry['time']),
                    'level': log_entry['level'],
                    'category': log_entry['cat'],
                    'message': log_entry['msg']
                })
            except json.JSONDecodeError as e:
                print(f"Warning: Invalid JSON on line {line_num}: {e}", file=sys.stderr)
                continue
            except KeyError as e:
                print(f"Warning: Missing field {e} on line {line_num}", file=sys.stderr)
                continue
    
    if not logs:
        raise ValueError("No valid log entries found")
    
    df = pd.DataFrame(logs)
    df = df.sort_values('time')
    return df


def analyze_logs(df: pd.DataFrame) -> Dict:
    """
    Generate summary statistics from log DataFrame.
    
    Args:
        df: DataFrame from parse_json_logs()
    
    Returns:
        Dictionary with analysis results
    """
    stats = {
        'total_entries': len(df),
        'time_span': (df['time'].max() - df['time'].min()).total_seconds(),
        'level_counts': df['level'].value_counts().to_dict(),
        'category_counts': df['category'].value_counts().to_dict(),
        'errors': df[df['level'] == 'error']['message'].tolist(),
        'warnings': df[df['level'] == 'warn']['message'].tolist(),
    }
    
    # Extract performance metrics (if present)
    perf_logs = df[df['category'] == 'perf']
    if not perf_logs.empty:
        stats['performance_logs'] = len(perf_logs)
    
    return stats


def print_summary(stats: Dict):
    """Print formatted summary of log analysis."""
    print("\n" + "="*60)
    print("ICARION Log Analysis Summary")
    print("="*60)
    
    print(f"\nTotal log entries: {stats['total_entries']}")
    print(f"Time span: {stats['time_span']:.2f} seconds")
    
    print("\n--- Log Levels ---")
    for level, count in stats['level_counts'].items():
        print(f"  {level:8s}: {count:5d} ({100*count/stats['total_entries']:.1f}%)")
    
    print("\n--- Categories ---")
    for cat, count in sorted(stats['category_counts'].items(), key=lambda x: -x[1]):
        print(f"  {cat:12s}: {count:5d}")
    
    if stats['errors']:
        print(f"\n--- Errors ({len(stats['errors'])}) ---")
        for i, err in enumerate(stats['errors'][:5], 1):
            print(f"  {i}. {err[:80]}...")
        if len(stats['errors']) > 5:
            print(f"  ... and {len(stats['errors']) - 5} more")
    
    if stats['warnings']:
        print(f"\n--- Warnings ({len(stats['warnings'])}) ---")
        for i, warn in enumerate(stats['warnings'][:5], 1):
            print(f"  {i}. {warn[:80]}...")
        if len(stats['warnings']) > 5:
            print(f"  ... and {len(stats['warnings']) - 5} more")
    
    if 'performance_logs' in stats:
        print(f"\n--- Performance ---")
        print(f"  Timer measurements: {stats['performance_logs']}")
    
    print("\n" + "="*60)


def main():
    """Main entry point."""
    if len(sys.argv) != 2:
        print("Usage: python analyze_logs.py <json_log_file>", file=sys.stderr)
        print("\nExample:")
        print("  ./icarion_main --log-format json --log-file sim.log config.json")
        print("  python analyze_logs.py sim.log")
        sys.exit(1)
    
    log_file = Path(sys.argv[1])
    
    if not log_file.exists():
        print(f"Error: Log file not found: {log_file}", file=sys.stderr)
        sys.exit(1)
    
    try:
        print(f"Parsing {log_file}...")
        df = parse_json_logs(log_file)
        
        print(f"Loaded {len(df)} log entries")
        print(f"Time range: {df['time'].min()} to {df['time'].max()}")
        
        stats = analyze_logs(df)
        print_summary(stats)
        
        # Export to CSV for further analysis
        csv_file = log_file.with_suffix('.csv')
        df.to_csv(csv_file, index=False)
        print(f"\n✓ Exported to CSV: {csv_file}")
        print("  You can now analyze with: pandas, matplotlib, seaborn, etc.")
        
    except Exception as e:
        print(f"Error analyzing logs: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
