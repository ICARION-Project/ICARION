#!/usr/bin/env python3
"""
Standardize file headers across the ICARION codebase.
Applies consistent MIT License headers to all source files.
"""

import os
import re
from pathlib import Path

# Standard header for all source files
STANDARD_HEADER = """// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors
"""

# Patterns to identify and remove old headers
OLD_HEADER_PATTERNS = [
    r'^// SPDX-License-Identifier:.*?\n',
    r'^// SPDX-FileCopyrightText:.*?\n',
    r'^/\*\*\s*\n\s*\*\s*={5,}\s*\n.*?\*\s*={5,}\s*\n\s*\*/\s*\n',  # Old block style
]

def has_header(content):
    """Check if file already has any ICARION header."""
    first_lines = '\n'.join(content.split('\n')[:5])
    return 'ICARION' in first_lines or 'SPDX' in first_lines or 'MIT License' in first_lines

def remove_old_headers(content):
    """Remove old header styles."""
    # Remove SPDX lines
    content = re.sub(r'^// SPDX-License-Identifier:.*?\n', '', content)
    content = re.sub(r'^// SPDX-FileCopyrightText:.*?\n', '', content)
    
    # Remove old block comment headers (multi-line)
    content = re.sub(
        r'^/\*\*.*?\*\s*={5,}\s*\n\s*\*/\s*\n',
        '',
        content,
        flags=re.DOTALL | re.MULTILINE
    )
    
    return content

def process_file(filepath):
    """Add/update standard header in a source file."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            original = f.read()
        
        # Skip empty files
        if not original.strip():
            return False
        
        # Remove old headers
        content = remove_old_headers(original)
        
        # Remove leading whitespace
        content = content.lstrip()
        
        # Add standard header
        new_content = STANDARD_HEADER + '\n' + content
        
        # Only write if changed
        if new_content != original:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(new_content)
            return True
    
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
    
    return False

def main():
    """Process all source files in the repository."""
    root = Path(__file__).parent.parent
    
    # File extensions to process
    extensions = {'.h', '.cpp', '.cu', '.cuh'}
    
    # Directories to skip
    skip_dirs = {'build', 'cmake', '.git', '__pycache__', 'results', '_deps'}
    
    updated_count = 0
    total_count = 0
    
    # Walk through source directories
    for src_dir in ['src', 'tests', 'tools']:
        src_path = root / src_dir
        if not src_path.exists():
            continue
        
        for filepath in src_path.rglob('*'):
            # Skip directories
            if filepath.is_dir():
                continue
            
            # Skip non-source files
            if filepath.suffix not in extensions:
                continue
            
            # Skip blacklisted directories
            if any(skip in filepath.parts for skip in skip_dirs):
                continue
            
            total_count += 1
            if process_file(filepath):
                updated_count += 1
                print(f"✓ {filepath.relative_to(root)}")
    
    print(f"\nProcessed {total_count} files, updated {updated_count}")

if __name__ == '__main__':
    main()
