#!/usr/bin/env python3

import argparse
import os
import re

def find_section(content, section_name):
    """Find the start and end of a section in .itp file"""
    pattern = r'^\s*\[' + re.escape(section_name) + r'\s*\]'
    lines = content.split('\n')
    
    start_idx = None
    end_idx = None
    
    for i, line in enumerate(lines):
        if re.match(pattern, line):
            start_idx = i
            break
    
    if start_idx is not None:
        for i in range(start_idx + 1, len(lines)):
            if lines[i].strip().startswith('['):
                end_idx = i
                break
        if end_idx is None:
            end_idx = len(lines)
    
    return start_idx, end_idx

def extract_receptor_topology(rec_top, rec_itp):
    """Extract receptor topology from receptor topology file to create .itp"""
    with open(rec_top, 'r') as f:
        lines = f.readlines()
    
    start_line = None
    for i, line in enumerate(lines):
        if '[ moleculetype ]' in line:
            start_line = i
            break
    
    if start_line is None:
        raise ValueError(f"Could not find [ moleculetype ] in {rec_top}")
    
    end_idx = None
    for i, line in enumerate(lines):
        if i <= start_line:
            continue
        if line.strip().startswith('#include') and '.ff/' in line:
            end_idx = i
            break
    
    if end_idx is None:
        end_idx = len(lines)
    
    receptor_itp = ''.join(lines[start_line:end_idx])
    
    with open(rec_itp, 'w') as f:
        f.write(receptor_itp)
    
    print(f"Created {rec_itp}")

def extract_water_models(ff_dir):
    """Extract water model names from watermodels.dat in forcefield directory"""
    watermodels_file = os.path.join(ff_dir, 'watermodels.dat')
    water_models = set()
    
    if os.path.exists(watermodels_file):
        with open(watermodels_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and len(line.split()) >= 1:
                    name = line.split()[0]
                    water_models.add(name)
    
    return water_models

def extract_ff_paths_from_top(top_file):
    """Extract forcefield, water model, and ions itp paths from receptor topology"""
    ff_path = None
    water_itp = None
    ions_itp = None
    ff_dir = None
    
    with open(top_file, 'r') as f:
        for line in f:
            if line.strip().startswith('#include') and '.ff/' in line:
                path = line.strip().split('"')[1] if '"' in line else line.strip()
                if 'forcefield.itp' in path:
                    ff_path = path
                    ff_dir = os.path.dirname(path)
                elif 'ions.itp' in path:
                    ions_itp = path
    
    water_models = extract_water_models(ff_dir) if ff_dir else set()
    
    with open(top_file, 'r') as f:
        for line in f:
            if line.strip().startswith('#include') and '.ff/' in line:
                path = line.strip().split('"')[1] if '"' in line else line.strip()
                if path.endswith('.itp') and 'forcefield' not in path and 'ions' not in path:
                    for wm in water_models:
                        if f'{wm}.itp' in path:
                            water_itp = path
                            break
                    if water_itp:
                        break
    
    return ff_path, water_itp, ions_itp

def get_moleculetype_name(top_file):
    """Extract moleculetype name from a topology file (.top or .itp)"""
    with open(top_file, 'r') as f:
        lines = f.readlines()
    
    in_moleculetype = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('[ moleculetype ]') or stripped == '[moleculetype]':
            in_moleculetype = True
            continue
        if in_moleculetype:
            if stripped.startswith('['):
                break
            if stripped and not stripped.startswith(';') and not stripped.startswith('#'):
                return stripped.split()[0]
    
    return None

def get_prefix(filename):
    """Extract prefix from filename (without extension)"""
    basename = os.path.basename(filename)
    if '.' in basename:
        return basename.rsplit('.', 1)[0]
    return basename

def get_atom_count(gro_file):
    """Get atom count from .gro file"""
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    return int(lines[1].strip())

def get_box_line(gro_file):
    """Get box coordinates from .gro file"""
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    return lines[-1].strip()

def combine_coordinates(rec_gro, lig_gro, output_gro):
    """Combine receptor and ligand coordinates into a single .gro file"""
    with open(rec_gro, 'r') as f:
        rec_lines = f.readlines()
    
    with open(lig_gro, 'r') as f:
        lig_lines = f.readlines()
    
    rec_atom_count = int(rec_lines[1].strip())
    lig_atom_count = int(lig_lines[1].strip())
    total_atoms = rec_atom_count + lig_atom_count
    
    box = get_box_line(rec_gro)
    
    combined = []
    combined.append(rec_lines[0].strip() + '\n')
    combined.append(f' {total_atoms}\n')
    
    for line in rec_lines[2:-1]:
        combined.append(line)
    
    for line in lig_lines[2:-1]:
        combined.append(line)
    
    combined.append(box + '\n')
    
    with open(output_gro, 'w') as f:
        f.writelines(combined)
    
    print(f"Created {output_gro} with {total_atoms} atoms")

def clean_itp_for_system(itp_file):
    """Remove [ atomtypes ] section from .itp file if it exists"""
    with open(itp_file, 'r') as f:
        content = f.read()
    
    if '[ atomtypes ]' in content or '[ atomtypes ]\n' in content:
        atomtypes_start = content.find('[ atomtypes ]')
        moleculetype_start = content.find('[ moleculetype ]', atomtypes_start)
        
        if moleculetype_start != -1:
            content = content[:atomtypes_start] + content[moleculetype_start:]
            
            with open(itp_file, 'w') as f:
                f.write(content)
            print(f"Removed [ atomtypes ] section from {itp_file}")
        else:
            print(f"Warning: {itp_file} has [ atomtypes ] but no [ moleculetype ] section")

def create_system_topology(args):
    """Create sys.top combining receptor and ligand topologies"""
    
    for itp in [args.rec_itp, args.lig_itp]:
        if itp and os.path.exists(itp):
            clean_itp_for_system(itp)
    
    sys_top = f"""; Include forcefield parameters
#include "{args.ff_path}"

; Include receptor topology
#include "{args.rec_itp}"

; Include ligand topology
#include "{args.lig_itp}"
"""
    
    if os.path.exists(args.water_itp):
        sys_top += f"""

; Include water topology
#include "{args.water_itp}"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
"""
    
    if os.path.exists(args.ions_itp):
        sys_top += f"""

; Include topology for ions
#include "{args.ions_itp}"
"""
    
    sys_top += f"""

[ system ]
; Name
{args.sys_name}

[ molecules ]
; Compound        #mols
{args.rec_name}             1
{args.lig_name}     1
"""
    
    with open(args.sys_top, 'w') as f:
        f.write(sys_top)
    
    print(f"Created {args.sys_top}")

def main():
    parser = argparse.ArgumentParser(description='Combine receptor and ligand for molecular dynamics')
    
    parser.add_argument('--rec-gro', default='prot.gro', help='Receptor coordinates (.gro)')
    parser.add_argument('--lig-gro', default='lig_redock.gro', help='Ligand coordinates (.gro)')
    parser.add_argument('--rec-top', default='topol.top', help='Receptor topology (.top)')
    parser.add_argument('--rec-itp', default='rec.itp', help='Receptor .itp output')
    parser.add_argument('--lig-itp', default='lig.itp', help='Ligand .itp file')
    
    parser.add_argument('--com-gro', default='com.gro', help='Combined coordinates output')
    parser.add_argument('--sys-top', default='sys.top', help='System topology output')
    
    parser.add_argument('--ff-path', help='Force field path (default: auto-detected from receptor topology)')
    parser.add_argument('--water-itp', help='Water topology (default: auto-detected from receptor topology)')
    parser.add_argument('--ions-itp', help='Ions topology (default: auto-detected from receptor topology)')
    
    parser.add_argument('--sys-name', help='System name (default: auto-generated from receptor and ligand filenames)')
    parser.add_argument('--rec-name', help='Receptor molecule name (default: auto-detected from receptor topology)')
    parser.add_argument('--lig-name', help='Ligand molecule name (default: auto-detected from ligand topology)')
    
    args = parser.parse_args()
    
    ff_path, water_itp, ions_itp = extract_ff_paths_from_top(args.rec_top)
    
    if args.ff_path is None:
        if ff_path:
            args.ff_path = ff_path
            print(f"Auto-detected forcefield: {args.ff_path}")
        else:
            parser.error("Could not auto-detect forcefield path. Please specify --ff-path.")
    
    if args.water_itp is None:
        if water_itp:
            args.water_itp = water_itp
            print(f"Auto-detected water model: {args.water_itp}")
        else:
            args.water_itp = ''
    
    if args.ions_itp is None:
        if ions_itp:
            args.ions_itp = ions_itp
            print(f"Auto-detected ions itp: {args.ions_itp}")
        else:
            args.ions_itp = ''
    
    if args.rec_name is None:
        args.rec_name = get_moleculetype_name(args.rec_top)
        if args.rec_name:
            print(f"Auto-detected receptor name: {args.rec_name}")
        else:
            parser.error("Could not auto-detect receptor molecule name. Please specify --rec-name.")
    
    if args.lig_name is None:
        args.lig_name = get_moleculetype_name(args.lig_itp)
        if args.lig_name:
            print(f"Auto-detected ligand name: {args.lig_name}")
        else:
            parser.error("Could not auto-detect ligand molecule name. Please specify --lig-name.")
    
    if args.sys_name is None:
        rec_prefix = get_prefix(args.rec_gro)
        lig_prefix = get_prefix(args.lig_gro)
        args.sys_name = f"{rec_prefix}_{lig_prefix}"
        print(f"Auto-generated system name: {args.sys_name}")
    
    print("\nExtracting receptor topology...")
    extract_receptor_topology(args.rec_top, args.rec_itp)
    
    print("Combining coordinates...")
    combine_coordinates(args.rec_gro, args.lig_gro, args.com_gro)
    
    print("Creating system topology...")
    create_system_topology(args)
    
    print("\nDock2com completed successfully!")

if __name__ == '__main__':
    main()
