#!/usr/bin/env python3
"""
Generate orbital cube files from ADMP XYZ trajectory files.

This script:
1. Searches for XYZ files in ADMP results directories
2. Processes frames from XYZ trajectory files
3. Creates input files for Gaussian calculations
4. Organizes the files by molecule and temperature
"""

import os
import argparse
from pathlib import Path
import json

def find_xyz_files(base_dir="../ADMP_decomposition_gaussian/admp_jobs/results"):
    """Find all XYZ trajectory files in the results directory."""
    base_path = Path(base_dir).resolve()
    print(f"Searching for .xyz files in: {base_path}")
    
    # First search method: use glob pattern
    all_xyz_files = list(base_path.glob("**/*.xyz"))
    
    # If no files found with the first method, try alternative search patterns
    if not all_xyz_files:
        print("No .xyz files found, trying broader search...")
        all_xyz_files = list(base_path.glob("**/*ADMP*.xyz"))
    
    # If still no files, try using os.walk which might be more thorough
    if not all_xyz_files:
        print("Still no .xyz files found, trying more thorough search...")
        for root, dirs, files in os.walk(base_path):
            for file in files:
                if file.endswith('.xyz'):
                    all_xyz_files.append(Path(os.path.join(root, file)))
    
    if not all_xyz_files:
        print(f"WARNING: No .xyz files found in {base_dir}")
        print("Please check that:")
        print("  1. The path to your results directory is correct")
        print("  2. ADMP calculations have completed with XYZ trajectory files")
        return []
    
    print(f"Found {len(all_xyz_files)} XYZ trajectory files")
    
    # Print some of the found files for verification
    if len(all_xyz_files) > 0:
        print("\nSample of files found:")
        for file in all_xyz_files[:min(5, len(all_xyz_files))]:
            print(f"  - {file}")
        if len(all_xyz_files) > 5:
            print(f"  ... and {len(all_xyz_files) - 5} more")
    
    return all_xyz_files

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

def create_gaussian_input_file(molecule, temp, timestep, atoms, output_dir, 
                              step_num, method="B3LYP", basis="6-31G(d)"):
    """Create a Gaussian input file for a single frame."""
    # Create base name for files
    base_name = f"{molecule}_{temp}_step{step_num:04d}"
    
    # Ensure the output directory exists
    mol_dir = Path(output_dir) / molecule / temp
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
        f.write(f"# {method}/{basis} pop=full density=current\n\n")
        f.write(f"{molecule} {timestep}\n\n")
        f.write("0 1\n")
        for symbol, x, y, z in atoms:
            f.write(f"{symbol:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
        f.write("\n")
    
    print(f"  - Created input file: {gjf_file}")
    return str(gjf_file)

def process_xyz_files(xyz_files, output_dir="./orbital_inputs", max_frames=10,
                     method="B3LYP", basis="6-31G(d)"):
    """Process XYZ files and create Gaussian input files."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Store information about all generated files
    all_inputs = {}
    
    for i, xyz_file in enumerate(xyz_files):
        xyz_path = Path(xyz_file)
        print(f"Processing file {i+1}/{len(xyz_files)}: {xyz_path}")
        
        # Parse molecule name and temperature from file name
        filename = xyz_path.name
        if '_ADMP_' in filename:
            molecule = filename.split('_ADMP_')[0]
            temp = filename.split('_ADMP_')[1].split('.')[0]
        else:
            # Fallback to directory names if filename pattern is different
            parts = str(xyz_file).split('/')
            molecule_idx = parts.index("results") + 1 if "results" in parts else -3
            temp_idx = molecule_idx + 1 if molecule_idx >= 0 else -2
            
            if molecule_idx >= 0 and temp_idx < len(parts):
                molecule = parts[molecule_idx]
                temp = parts[temp_idx]
            else:
                molecule = filename.split('.')[0]
                temp = "unknown"
        
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
            try:
                step_num = int(timestep.split()[-1])
            except (ValueError, IndexError):
                step_num = frame_idx
            
            # Create input file
            gjf_file = create_gaussian_input_file(
                molecule, temp, timestep, atoms, output_dir, 
                step_num, method, basis
            )
            
            # Add to the list of inputs
            molecule_inputs.append({
                "input_file": gjf_file,
                "molecule": molecule,
                "temperature": temp,
                "step": step_num
            })
        
        # Add this molecule's inputs to the main dictionary
        if molecule_inputs:
            if molecule not in all_inputs:
                all_inputs[molecule] = {}
            all_inputs[molecule][temp] = molecule_inputs
        
        print(f"  - Created {len(molecule_inputs)} input files")
        print("----------------------------------------")
    
    # Write the summary JSON file
    summary_file = Path(output_dir) / "input_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(all_inputs, f, indent=2)
    
    print(f"\nSummary of generated input files saved to: {summary_file}")
    return all_inputs

def main():
    parser = argparse.ArgumentParser(description="Generate Gaussian input files from XYZ trajectories")
    parser.add_argument("--base-dir", default="../ADMP_decomposition_gaussian/admp_jobs/results",
                      help="Base directory containing ADMP results")
    parser.add_argument("--output-dir", default="./orbital_inputs",
                      help="Directory to store generated input files (default: ./orbital_inputs)")
    parser.add_argument("--max-frames", type=int, default=10,
                      help="Maximum number of frames per trajectory (default: 10)")
    parser.add_argument("--method", default="B3LYP",
                      help="Computational method to use (default: B3LYP)")
    parser.add_argument("--basis", default="6-31G(d)",
                      help="Basis set to use (default: 6-31G(d))")
    
    args = parser.parse_args()
    
    print("Searching for XYZ trajectory files...")
    xyz_files = find_xyz_files(args.base_dir)
    
    if xyz_files:
        all_inputs = process_xyz_files(
            xyz_files,
            output_dir=args.output_dir,
            max_frames=args.max_frames,
            method=args.method,
            basis=args.basis
        )
        
        # Count total inputs
        total_inputs = sum(len(temp_files) 
                         for molecule in all_inputs.values() 
                         for temp_files in molecule.values())
        
        print(f"\nNext steps:")
        print(f"1. Generated {total_inputs} Gaussian input files in: {args.output_dir}")
        print(f"2. Run the Gaussian calculations with the submission script:")
        print(f"   $ sbatch generate_orbitals_from_ADMP/submit_orbital_calculations.sh")
        print(f"   (No parameters needed - will use {args.output_dir} by default)")
        print(f"\nYou can control the number of frames with --max-frames")
        print(f"Default is 10 frames per trajectory to keep computation time reasonable.")
        print(f"You can also specify different computational methods with --method and --basis")
    else:
        print("No XYZ files found. Cannot generate Gaussian input files.")

if __name__ == "__main__":
    main() 