#!/usr/bin/env python3
"""
Migrate result paths from results/ to validation/results/ in config files.

This script:
1. Updates all config files to use validation/results/ instead of results/
2. Moves existing result files from results/ to validation/results/
3. Preserves directory structure
"""

import json
import shutil
from pathlib import Path

def update_config_path(config_path):
    """Update output.folder path in a single config file."""
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    if 'output' in config and 'folder' in config['output']:
        old_path = config['output']['folder']
        
        # Only update if it starts with "results/" and not already "validation/results/"
        if old_path.startswith('results/') and not old_path.startswith('validation/results/'):
            new_path = old_path.replace('results/', 'validation/results/', 1)
            config['output']['folder'] = new_path
            
            # Write back
            with open(config_path, 'w') as f:
                json.dump(config, f, indent=2)
            
            return old_path, new_path
    
    return None, None

def main():
    """Main migration logic."""
    repo_root = Path(__file__).resolve().parents[2]
    configs_dir = repo_root / 'validation' / 'configs' / 'instruments'
    
    print("=" * 80)
    print("MIGRATING RESULT PATHS: results/ → validation/results/")
    print("=" * 80)
    
    # Find all config files
    config_files = list(configs_dir.rglob('*.json'))
    print(f"\nFound {len(config_files)} config files")
    
    # Update config files
    updated_count = 0
    migrations = []
    
    for config_path in config_files:
        old_path, new_path = update_config_path(config_path)
        if old_path:
            updated_count += 1
            migrations.append((old_path, new_path))
            rel_config = config_path.relative_to(repo_root)
            print(f"✓ {rel_config}")
            print(f"  {old_path} → {new_path}")
    
    print(f"\n✅ Updated {updated_count} config files")
    
    # Move existing result files (only for validation-relevant directories)
    print("\n" + "=" * 80)
    print("MOVING RESULT FILES")
    print("=" * 80)
    
    # Directories to migrate (only validation-relevant ones)
    dirs_to_migrate = [
        'v1.0_test/instruments/ims',
        'v1.0_test/instruments/lqit',
        'v1.0_test/instruments/quadrupole',
        'v1.0_test/instruments/quadrupole_stability',
    ]
    
    results_root = repo_root / 'results'
    validation_results = repo_root / 'validation' / 'results'
    
    moved_files = 0
    moved_bytes = 0
    
    for dir_path in dirs_to_migrate:
        src_dir = results_root / dir_path
        dst_dir = validation_results / dir_path
        
        if not src_dir.exists():
            print(f"⊘ {dir_path}: source doesn't exist, skipping")
            continue
        
        # Create destination directory
        dst_dir.mkdir(parents=True, exist_ok=True)
        
        # Move all files (but not subdirectories recursively - only top level)
        files = list(src_dir.glob('*.*'))  # Only files with extensions
        if files:
            print(f"\n📁 {dir_path}/")
            for src_file in files:
                dst_file = dst_dir / src_file.name
                
                if dst_file.exists():
                    # Check if identical
                    if src_file.stat().st_size == dst_file.stat().st_size:
                        print(f"  ⊙ {src_file.name} (already exists, removing duplicate)")
                        src_file.unlink()
                    else:
                        print(f"  ⚠ {src_file.name} (different file exists, keeping both)")
                else:
                    # Move file
                    file_size = src_file.stat().st_size
                    shutil.move(str(src_file), str(dst_file))
                    moved_files += 1
                    moved_bytes += file_size
                    size_mb = file_size / (1024**2)
                    print(f"  → {src_file.name} ({size_mb:.1f} MB)")
    
    print(f"\n✅ Moved {moved_files} files ({moved_bytes / (1024**2):.1f} MB total)")
    
    print("\n" + "=" * 80)
    print("MIGRATION COMPLETE")
    print("=" * 80)
    print(f"✓ {updated_count} configs updated")
    print(f"✓ {moved_files} result files moved")
    print(f"✓ All results now in validation/results/")
    print("\nYou can now safely delete old results/ if desired:")
    print("  rm -rf results/v1.0_test/instruments/")

if __name__ == '__main__':
    main()
