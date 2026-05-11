#!/usr/bin/env python3
"""
Bypass parmed angle type 3 incompatibility for gmx_mmpbsa.

This script:
1. Changes angle function type 3 to type 1 in ligand topology,
   adding a comment for tracking. Saves with 'bypass_' prefix.
2. Updates main topology to include the bypass ligand topology,
   saves with 'bypass_' prefix.

Usage:
    python bypass_angle_type3.py <ligand.itp> <main.top>
"""

import argparse
import os
import re
import sys


def modify_ligand_itp(input_file: str, output_file: str) -> int:
    """
    Change angle function type 3 to type 1 in ligand topology.
    
    Returns the number of modifications made.
    """
    modified_count = 0
    in_angles_section = False
    lines_out = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        stripped = line.strip()
        
        if stripped.startswith('[ angles ]'):
            in_angles_section = True
            lines_out.append(line)
            continue
        
        if stripped.startswith('[') and stripped.endswith(']'):
            in_angles_section = False
            lines_out.append(line)
            continue
        
        if in_angles_section and stripped and not stripped.startswith(';'):
            parts = line.split()
            if len(parts) >= 4:
                func_type = parts[3]
                if func_type == '3':
                    parts[3] = '1'
                    
                    original_comment = ''
                    comment_idx = line.find(';')
                    if comment_idx != -1:
                        original_comment = line[comment_idx:].rstrip()
                    
                    indent = line[:len(line) - len(line.lstrip())]
                    data_part = '     '.join(parts[:4]) + '     ' + '     '.join(parts[4:])
                    
                    if 'type 3 to 1' not in original_comment:
                        if original_comment:
                            new_line = indent + data_part + ' type 3 to 1 for parmed bypass\n'
                        else:
                            new_line = indent + data_part + '    ; type 3 to 1 for parmed bypass\n'
                    else:
                        new_line = indent + data_part + '\n'
                    
                    lines_out.append(new_line)
                    modified_count += 1
                    continue
        
        lines_out.append(line)
    
    with open(output_file, 'w') as f:
        f.writelines(lines_out)
    
    return modified_count


def find_ligand_include(lines: list, ligand_filename: str) -> str | None:
    """
    Find ligand include line in main topology matching the ligand filename.
    Returns include filename or None if not found.
    """
    include_pattern = re.compile(r'#include\s+["\']([^"\']+\.itp)["\']')
    ligand_basename = os.path.basename(ligand_filename)
    
    for line in lines:
        match = include_pattern.search(line)
        if match:
            include_file = match.group(1)
            include_basename = os.path.basename(include_file)
            if include_basename == ligand_basename:
                return include_file
    
    return None


def modify_main_top(input_file: str, output_file: str, ligand_filename: str) -> tuple:
    """
    Change topology include of ligand in main topology.
    
    Returns (success, include_filename_used).
    """
    lines_out = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    target_include = find_ligand_include(lines, ligand_filename)
    if target_include is None:
        return False, None
    
    bypass_include = 'bypass_' + os.path.basename(target_include)
    modified = False
    
    for line in lines:
        include_pattern = re.compile(r'(#include\s+["\'])' + re.escape(target_include) + r'(["\'])')
        match = include_pattern.search(line)
        if match:
            new_line = include_pattern.sub(r'\1' + bypass_include + r'\2', line)
            lines_out.append(new_line)
            modified = True
        else:
            lines_out.append(line)
    
    with open(output_file, 'w') as f:
        f.writelines(lines_out)
    
    return modified, target_include


def get_output_path(input_path: str, prefix: str = 'bypass_') -> str:
    """Generate output path with prefix added to filename."""
    dirname = os.path.dirname(input_path)
    basename = os.path.basename(input_path)
    return os.path.join(dirname, prefix + basename)


def main():
    parser = argparse.ArgumentParser(
        description='Bypass parmed angle type 3 incompatibility for gmx_mmpbsa.'
    )
    parser.add_argument('ligand_itp', help='Ligand topology file (.itp)')
    parser.add_argument('main_top', help='Main topology file (.top)')
    parser.add_argument('-o', '--output-dir', help='Output directory (default: same as input)')
    parser.add_argument('--ligand-include', help='Ligand include path in main topology (auto-detected if not specified)')
    
    args = parser.parse_args()
    
    lig_itp = os.path.abspath(args.ligand_itp)
    main_top = os.path.abspath(args.main_top)
    
    if not os.path.exists(lig_itp):
        print(f"Error: Ligand topology file not found: {lig_itp}")
        sys.exit(1)
    if not os.path.exists(main_top):
        print(f"Error: Main topology file not found: {main_top}")
        sys.exit(1)
    
    output_dir = args.output_dir
    if output_dir:
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        bypass_lig_itp = os.path.join(output_dir, 'bypass_' + os.path.basename(lig_itp))
        bypass_main_top = os.path.join(output_dir, 'bypass_' + os.path.basename(main_top))
    else:
        bypass_lig_itp = get_output_path(lig_itp)
        bypass_main_top = get_output_path(main_top)
    
    print(f"Processing ligand topology: {lig_itp}")
    n_modified = modify_ligand_itp(lig_itp, bypass_lig_itp)
    print(f"  Modified {n_modified} angle type 3 -> type 1 entries")
    print(f"  Output: {bypass_lig_itp}")
    
    ligand_include = args.ligand_include if args.ligand_include else os.path.basename(lig_itp)
    print(f"\nProcessing main topology: {main_top}")
    modified, include_used = modify_main_top(main_top, bypass_main_top, ligand_include)
    if modified and include_used:
        print(f"  Updated ligand include: {include_used} -> bypass_{os.path.basename(include_used)}")
        print(f"  Output: {bypass_main_top}")
    else:
        print(f"  Warning: No ligand include found to modify")
    
    print("\nDone.")


if __name__ == '__main__':
    main()
