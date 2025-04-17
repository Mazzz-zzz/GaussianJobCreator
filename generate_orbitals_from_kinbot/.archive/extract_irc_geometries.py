#!/usr/bin/env python3
"""
Extract geometries from KinBot IRC log files.

This script:
1. Searches for IRC log files in kinbot-exploration directories
2. Extracts coordinates from IRC calculation steps
3. Creates combined XYZ trajectory files for visualizing the complete reaction path
"""

import os
import sys
import argparse
from pathlib import Path
import re

# Atomic number to element symbol mapping
CODE = {"1" : "H", "2" : "He", "3" : "Li", "4" : "Be", "5" : "B", \
"6"  : "C", "7"  : "N", "8"  : "O", "9" : "F", "10" : "Ne", \
"11" : "Na" , "12" : "Mg" , "13" : "Al" , "14" : "Si" , "15" : "P", \
"16" : "S"  , "17" : "Cl" , "18" : "Ar" , "19" : "K"  , "20" : "Ca", \
"21" : "Sc" , "22" : "Ti" , "23" : "V"  , "24" : "Cr" , "25" : "Mn", \
"26" : "Fe" , "27" : "Co" , "28" : "Ni" , "29" : "Cu" , "30" : "Zn", \
"31" : "Ga" , "32" : "Ge" , "33" : "As" , "34" : "Se" , "35" : "Br", \
"36" : "Kr" , "37" : "Rb" , "38" : "Sr" , "39" : "Y"  , "40" : "Zr", \
"41" : "Nb" , "42" : "Mo" , "43" : "Tc" , "44" : "Ru" , "45" : "Rh", \
"46" : "Pd" , "47" : "Ag" , "48" : "Cd" , "49" : "In" , "50" : "Sn", \
"51" : "Sb" , "52" : "Te" , "53" : "I"  , "54" : "Xe" , "55" : "Cs", \
"56" : "Ba" , "57" : "La" , "58" : "Ce" , "59" : "Pr" , "60" : "Nd", \
"61" : "Pm" , "62" : "Sm" , "63" : "Eu" , "64" : "Gd" , "65" : "Tb", \
"66" : "Dy" , "67" : "Ho" , "68" : "Er" , "69" : "Tm" , "70" : "Yb", \
"71" : "Lu" , "72" : "Hf" , "73" : "Ta" , "74" : "W"  , "75" : "Re", \
"76" : "Os" , "77" : "Ir" , "78" : "Pt" , "79" : "Au" , "80" : "Hg", \
"81" : "Tl" , "82" : "Pb" , "83" : "Bi" , "84" : "Po" , "85" : "At", \
"86" : "Rn" , "87" : "Fr" , "88" : "Ra" , "89" : "Ac" , "90" : "Th", \
"91" : "Pa" , "92" : "U"  , "93" : "Np" , "94" : "Pu" , "95" : "Am", \
"96" : "Cm" , "97" : "Bk" , "98" : "Cf" , "99" : "Es" ,"100" : "Fm", \
"101": "Md" ,"102" : "No" ,"103" : "Lr" ,"104" : "Rf" ,"105" : "Db", \
"106": "Sg" ,"107" : "Bh" ,"108" : "Hs" ,"109" : "Mt" ,"110" : "Ds", \
"111": "Rg" ,"112" : "Uub","113" : "Uut","114" : "Uuq","115" : "Uup", \
"116": "Uuh","117" : "Uus","118" : "Uuo"}


def find_irc_log_files(base_dir="kinbot-exploration"):
    """Find all IRC log files in the kinbot-exploration directory and subdirectories."""
    base_path = Path(base_dir).resolve()
    print(f"Searching for IRC log files in: {base_path}")
    
    irc_log_files = []
    
    # Walk through directories and find IRC log files
    for root, dirs, files in os.walk(base_path):
        for file in files:
            if file.endswith('.log') and ('_IRC_F' in file or '_IRC_R' in file):
                log_path = Path(os.path.join(root, file))
                irc_log_files.append(log_path)
    
    if not irc_log_files:
        print(f"WARNING: No IRC log files found in {base_dir}")
        return []
    
    print(f"Found {len(irc_log_files)} IRC log files")
    
    # Print some of the found files for verification
    if len(irc_log_files) > 0:
        print("\nSample of files found:")
        for file in irc_log_files[:min(5, len(irc_log_files))]:
            print(f"  - {file}")
        if len(irc_log_files) > 5:
            print(f"  ... and {len(irc_log_files) - 5} more")
    
    return irc_log_files


def extract_irc_frames(log_file):
    """Extract molecular geometries from IRC log file."""
    frames = []
    reaction_coords = []
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        # Find all "Standard orientation" or "Input orientation" sections
        orientation_indices = []
        for i, line in enumerate(lines):
            if "Standard orientation:" in line or "Input orientation:" in line:
                orientation_indices.append(i)
        
        # Try to extract IRC step information
        for i, line in enumerate(lines):
            if "IRC-IRC-IRC" in line and "Point Number" in lines[i+1]:
                try:
                    point_line = lines[i+1]
                    point_match = re.search(r'Point Number\s+(\d+)', point_line)
                    if point_match:
                        step_num = int(point_match.group(1))
                        reaction_coords.append(step_num)
                except (IndexError, ValueError):
                    pass
        
        # If we couldn't find IRC points, just number them sequentially
        if not reaction_coords:
            reaction_coords = list(range(len(orientation_indices)))
        
        # Extract coordinates from each orientation section
        for idx, orient_idx in enumerate(orientation_indices):
            try:
                # Get the IRC step number or use sequential numbering
                step_num = reaction_coords[idx] if idx < len(reaction_coords) else idx
                
                # Skip header lines to reach coordinate data
                i = orient_idx + 5
                atoms = []
                
                # Read atom coordinates until we hit the separator
                while i < len(lines) and "-----" not in lines[i]:
                    parts = lines[i].split()
                    if len(parts) >= 6:
                        atom_idx = parts[0]
                        atom_num = parts[1]
                        symbol = CODE.get(atom_num, "X")  # Use X for unknown elements
                        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                        atoms.append((symbol, x, y, z))
                    i += 1
                
                # Only add if we found atoms
                if atoms:
                    # Create timestep label
                    if "_IRC_F" in str(log_file):
                        direction = "Forward"
                    elif "_IRC_R" in str(log_file):
                        direction = "Reverse"
                    else:
                        direction = "Unknown"
                        
                    timestep = f"IRC {direction} Point {step_num}"
                    frames.append((timestep, atoms, step_num, direction))
            except Exception as e:
                print(f"Error extracting frame {idx} from {log_file}: {str(e)}")
        
        return frames
        
    except Exception as e:
        print(f"Error reading {log_file}: {str(e)}")
        return []


def save_xyz_file(frames, output_path):
    """Save a multi-frame XYZ file from extracted geometries."""
    with open(output_path, 'w') as f:
        for timestep, atoms, step_num, direction in frames:
            # Number of atoms
            f.write(f"{len(atoms)}\n")
            # Comment line with step info
            f.write(f"{timestep}\n")
            # Atom coordinates
            for symbol, x, y, z in atoms:
                f.write(f"{symbol} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    
    print(f"  - Saved {len(frames)} frames to {output_path}")
    return output_path


def process_irc_files(irc_files, output_dir="./irc_geometries"):
    """Process IRC log files and extract geometries."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Group IRC files by reaction
    reaction_groups = {}
    for log_file in irc_files:
        filename = log_file.name
        
        # Extract reaction name
        match = re.match(r'(.+?)_(IRC_[FR])', filename)
        if match:
            reaction_name = match.group(1)
            irc_type = match.group(2)
        else:
            # Fallback if pattern is different
            reaction_name = filename.split('_IRC_')[0]
            irc_type = "IRC_F" if "_IRC_F" in filename else "IRC_R"
        
        if reaction_name not in reaction_groups:
            reaction_groups[reaction_name] = {}
        
        reaction_groups[reaction_name][irc_type] = log_file
    
    # Process each reaction
    for reaction_name, files in reaction_groups.items():
        print(f"Processing reaction: {reaction_name}")
        
        # Create directory for this reaction
        reaction_dir = Path(output_dir) / reaction_name
        reaction_dir.mkdir(parents=True, exist_ok=True)
        
        # Collect frames from forward and reverse directions
        forward_frames = []
        reverse_frames = []
        
        # Process forward IRC if available
        if "IRC_F" in files:
            forward_log = files["IRC_F"]
            print(f"  Extracting forward IRC from: {forward_log}")
            forward_frames = extract_irc_frames(forward_log)
            print(f"  - Found {len(forward_frames)} frames in forward IRC")
        
        # Process reverse IRC if available
        if "IRC_R" in files:
            reverse_log = files["IRC_R"]
            print(f"  Extracting reverse IRC from: {reverse_log}")
            reverse_frames = extract_irc_frames(reverse_log)
            print(f"  - Found {len(reverse_frames)} frames in reverse IRC")
        
        # Skip if no frames found
        if not forward_frames and not reverse_frames:
            print(f"  WARNING: No frames extracted for {reaction_name}")
            continue
        
        # Create combined XYZ trajectory
        combined_path = reaction_dir / f"{reaction_name}_combined.xyz"
        with open(combined_path, 'w') as f:
            # Write reverse frames in reverse order if available
            if reverse_frames:
                for i in range(len(reverse_frames)-1, -1, -1):
                    timestep, atoms, step_num, direction = reverse_frames[i]
                    f.write(f"{len(atoms)}\n")
                    f.write(f"{timestep}\n")
                    for symbol, x, y, z in atoms:
                        f.write(f"{symbol} {x:12.6f} {y:12.6f} {z:12.6f}\n")
            
            # Write forward frames if available
            if forward_frames:
                for timestep, atoms, step_num, direction in forward_frames:
                    f.write(f"{len(atoms)}\n")
                    f.write(f"{timestep}\n")
                    for symbol, x, y, z in atoms:
                        f.write(f"{symbol} {x:12.6f} {y:12.6f} {z:12.6f}\n")
        
        total_frames = len(forward_frames) + len(reverse_frames)
        print(f"  - Created combined trajectory with {total_frames} frames: {combined_path}")
        print("----------------------------------------")
    
    return reaction_groups


def main():
    parser = argparse.ArgumentParser(description="Extract geometries from KinBot IRC logs")
    parser.add_argument("--base-dir", default="../kinbot-exploration",
                        help="Base directory containing KinBot results")
    parser.add_argument("--output-dir", default="./irc_geometries",
                        help="Directory to store extracted geometries")
    
    args = parser.parse_args()
    
    print("Searching for IRC log files...")
    irc_files = find_irc_log_files(args.base_dir)
    
    if irc_files:
        reaction_groups = process_irc_files(irc_files, output_dir=args.output_dir)
        
        # Print summary
        print(f"\nSummary:")
        print(f"Processed {len(reaction_groups)} reactions")
        print(f"All combined trajectories have been saved to: {args.output_dir}")
        print(f"\nYou can visualize these XYZ files using tools like VMD, Avogadro, or PyMOL")
    else:
        print("No IRC log files found. No geometries could be extracted.")


if __name__ == "__main__":
    main() 