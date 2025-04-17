#!/usr/bin/env python3
"""
Combined script to extract IRC geometries and generate orbital calculation inputs.

This script:
1. Searches for IRC log files in kinbot-exploration directories
2. Extracts coordinates from IRC calculation steps
3. Creates combined XYZ trajectory files for visualizing reaction paths
4. Generates Gaussian input files for orbital calculations from these XYZ files
"""

import os
import sys
import argparse
import re
import json
from pathlib import Path

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

# ============== IRC GEOMETRY EXTRACTION FUNCTIONS ==============

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


def get_successful_reactions(base_dir="kinbot-exploration"):
    """Find all successful IRC reactions by parsing kinbot_monitor.out files."""
    successful_reactions = {}
    
    # Walk through directories and find monitor files
    for root, dirs, files in os.walk(Path(base_dir)):
        if "kinbot_monitor.out" in files:
            monitor_file = os.path.join(root, "kinbot_monitor.out")
            try:
                with open(monitor_file, 'r') as f:
                    for line in f.readlines():
                        parts = line.strip().split()
                        if len(parts) >= 3 and parts[0] == "-1":  # Successful reaction
                            reaction_name = parts[2]
                            # Check if this is an IRC reaction (either r12_insertion, r13_insertion, etc.)
                            if any(marker in reaction_name for marker in ["_r12_", "_r13_", "_intra_"]):
                                successful_reactions[reaction_name] = True
            except Exception as e:
                print(f"Error reading {monitor_file}: {str(e)}")
    
    print(f"Found {len(successful_reactions)} successful IRC reactions")
    return successful_reactions


def process_irc_files(irc_files, output_dir="./irc_geometries", successful_reactions=None):
    """Process IRC log files and extract geometries."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Store paths to generated XYZ files
    generated_xyz_files = []
    
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
    processed_reactions = 0
    for reaction_name, files in reaction_groups.items():
        # Skip unsuccessful reactions if we have a filter
        if successful_reactions is not None and reaction_name not in successful_reactions:
            print(f"Skipping reaction {reaction_name} as it did not lead to products")
            continue
            
        print(f"Processing reaction: {reaction_name}")
        processed_reactions += 1
        
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
        
        # Add the generated XYZ file to the list
        generated_xyz_files.append(combined_path)
    
    print(f"Successfully processed {processed_reactions} out of {len(reaction_groups)} reactions")
    return generated_xyz_files, reaction_groups

# ============== ORBITAL GENERATION FUNCTIONS ==============

def read_xyz_frames(xyz_file):
    """Extract frames from XYZ trajectory file."""
    frames = []
    try:
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            try:
                n_atoms = int(lines[i].strip())
                timestep = lines[i+1].strip()
                atoms = []
                for j in range(n_atoms):
                    atom_line = lines[i+2+j].strip().split()
                    atoms.append((atom_line[0], float(atom_line[1]), float(atom_line[2]), float(atom_line[3])))
                frames.append((timestep, atoms))
                i += n_atoms + 2
            except (ValueError, IndexError):
                break
    except Exception as e:
        print(f"Error reading {xyz_file}: {str(e)}")
        return []
    
    print(f"  Extracted {len(frames)} frames from {xyz_file}")
    return frames


def create_gaussian_input_file(molecule, reaction_type, timestep, atoms, output_dir, 
                              step_num, method="B3LYP", basis="6-31G(d)"):
    """Create a Gaussian input file for a single frame."""
    # Create base name for files
    base_name = f"{molecule}_{reaction_type}_step{step_num:04d}"
    
    # Ensure the output directory exists
    mol_dir = Path(output_dir) / molecule / reaction_type
    mol_dir.mkdir(parents=True, exist_ok=True)
    
    # Create Gaussian input file
    gjf_file = mol_dir / f"{base_name}.gjf"
    
    # Skip if the file already exists
    if gjf_file.exists():
        print(f"  - Input file already exists: {gjf_file}")
        return str(gjf_file)
    
    with open(gjf_file, 'w') as f:
        f.write(f"%chk={base_name}.chk\n")
        f.write("%mem=8GB\n")
        f.write("%nprocshared=4\n")
        f.write(f"# {method}/{basis} NoSymm pop=full density=current\n\n")
        f.write(f"{molecule} {timestep}\n\n")
        f.write("0 1\n")
        for symbol, x, y, z in atoms:
            f.write(f"{symbol:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
        f.write("\n")
    
    print(f"  - Created input file: {gjf_file}")
    return str(gjf_file)


def process_xyz_files(xyz_files, output_dir="./orbital_inputs", max_frames=0,
                     method="B3LYP", basis="6-31G(d)"):
    """Process XYZ files and create Gaussian input files."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Store information about all generated files
    all_inputs = {}
    
    for i, xyz_file in enumerate(xyz_files):
        xyz_path = Path(xyz_file)
        print(f"Processing file {i+1}/{len(xyz_files)}: {xyz_path}")
        
        # Parse molecule name from file path (for IRC files)
        filename = xyz_path.name
        if "_combined.xyz" in filename:
            molecule = filename.split("_combined.xyz")[0]
            # For IRC files, the "temperature" field is actually the reaction type
            reaction_type = "IRC"
        else:
            # Fallback to simpler naming if pattern doesn't match
            molecule = filename.split('.')[0]
            reaction_type = "unknown"
        
        # Process frames from this XYZ file
        frames = read_xyz_frames(xyz_file)
        
        # Select frames based on max_frames
        if max_frames > 0 and len(frames) > max_frames:
            step = max(1, len(frames) // max_frames)
            selected_frames = frames[::step][:max_frames]
            print(f"  - Processing {len(selected_frames)} selected frames (out of {len(frames)})")
        else:
            selected_frames = frames
            print(f"  - Processing all {len(frames)} frames")
        
        # Store input files for this XYZ file
        molecule_inputs = []
        
        # Create Gaussian input files for each selected frame
        for frame_idx, (timestep, atoms) in enumerate(selected_frames):
            # Use sequential numbering for the combined trajectory
            step_num = frame_idx
            
            # Create input file - use "Combined" as the reaction type instead of Forward/Reverse
            gjf_file = create_gaussian_input_file(
                molecule, "Combined", timestep, atoms, output_dir, 
                step_num, method, basis
            )
            
            # Add to the list of inputs
            molecule_inputs.append({
                "input_file": gjf_file,
                "molecule": molecule,
                "reaction_type": "Combined",
                "step": step_num,
                "original_timestep": timestep  # Preserve the original timestep info
            })
        
        # Add this molecule's inputs to the main dictionary
        if molecule_inputs:
            if molecule not in all_inputs:
                all_inputs[molecule] = {}
            all_inputs[molecule][reaction_type] = molecule_inputs
        
        print(f"  - Created {len(molecule_inputs)} input files")
        print("----------------------------------------")
    
    # Write the summary JSON file
    summary_file = Path(output_dir) / "input_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(all_inputs, f, indent=2)
    
    print(f"\nSummary of generated input files saved to: {summary_file}")
    return all_inputs


def main():
    parser = argparse.ArgumentParser(description="Extract IRC geometries and generate orbital calculation inputs")
    parser.add_argument("--irc-dir", default="../kinbot-exploration",
                        help="Base directory containing KinBot results (default: ../kinbot-exploration)")
    parser.add_argument("--xyz-dir", default="./irc_geometries",
                        help="Directory to store extracted geometries (default: ./irc_geometries)")
    parser.add_argument("--output-dir", default="./orbital_inputs",
                        help="Directory to store orbital input files (default: ./orbital_inputs)")
    parser.add_argument("--max-frames", type=int, default=20,
                        help="Maximum number of frames per trajectory (0 for all frames, default: 20)")
    parser.add_argument("--method", default="B3LYP",
                        help="Computational method to use (default: B3LYP)")
    parser.add_argument("--basis", default="6-31G(d)",
                        help="Basis set to use (default: 6-31G(d))")
    parser.add_argument("--skip-irc", action="store_true",
                        help="Skip IRC extraction (use if XYZ files were already generated)")
    parser.add_argument("--all-reactions", action="store_true",
                        help="Process all reactions, even those without products")
    
    args = parser.parse_args()
    
    # STEP 1: Extract IRC geometries (if not skipped)
    if not args.skip_irc:
        print("=== STEP 1: Extracting IRC geometries ===")
        print("Searching for IRC log files...")
        irc_files = find_irc_log_files(args.irc_dir)
        
        if irc_files:
            # Find successful reactions (only if we're filtering)
            successful_reactions = None
            if not args.all_reactions:
                print("\nIdentifying successful IRC reactions...")
                successful_reactions = get_successful_reactions(args.irc_dir)
                print(f"Found {len(successful_reactions)} successful IRC reactions")
                
            # Process IRC files
            xyz_files, reaction_groups = process_irc_files(
                irc_files, 
                output_dir=args.xyz_dir,
                successful_reactions=successful_reactions
            )
            
            # Print summary
            print(f"\nSummary:")
            print(f"Processed {len(xyz_files)} reactions with successful pathways")
            print(f"All combined trajectories have been saved to: {args.xyz_dir}")
        else:
            print("No IRC log files found. No geometries could be extracted.")
            xyz_files = []
    else:
        print("=== Skipping IRC extraction, using existing XYZ files ===")
        # Find all XYZ files in the XYZ directory
        xyz_path = Path(args.xyz_dir)
        xyz_files = list(xyz_path.glob("**/*.xyz"))
        print(f"Found {len(xyz_files)} existing XYZ files in {args.xyz_dir}")
    
    # STEP 2: Generate orbital input files
    if xyz_files:
        print("\n=== STEP 2: Generating orbital calculation inputs ===")
        all_inputs = process_xyz_files(
            xyz_files,
            output_dir=args.output_dir,
            max_frames=args.max_frames,
            method=args.method,
            basis=args.basis
        )
        
        # Count total inputs
        total_inputs = sum(len(reaction_files) 
                        for molecule in all_inputs.values() 
                        for reaction_files in molecule.values())
        
        print(f"\nNext steps:")
        print(f"1. Generated {total_inputs} Gaussian input files in: {args.output_dir}")
        print(f"2. Use the submit_orbital_calculations.sh script to run the calculations:")
        print(f"   $ sbatch submit_orbital_calculations.sh {args.output_dir}")
        print(f"\nYou can control the number of frames with --max-frames (0 for all frames)")
        print(f"You can also specify different computational methods with --method and --basis")
    else:
        print("No XYZ files available. Cannot generate Gaussian input files.")


if __name__ == "__main__":
    main() 