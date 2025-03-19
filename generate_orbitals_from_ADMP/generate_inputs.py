#!/usr/bin/env python3
"""
Generate orbital cube files from ADMP formatted checkpoint files (.fchk).

This script:
1. Searches for formatted checkpoint files in ADMP results directories
2. Creates a bash script to generate cube files for HOMO and LUMO orbitals
3. Organizes the cube files by molecule and temperature
"""

import os
import glob
import argparse
import subprocess
from pathlib import Path

def find_fchk_files(base_dir="../ADMP_decomposition_gaussian/admp_jobs/results"):
    """Find all formatted checkpoint files in the results directory."""
    base_path = Path(base_dir).resolve()
    print(f"Searching for .fchk files in: {base_path}")
    
    # First search method: use glob pattern
    all_fchk_files = list(base_path.glob("**/checkpoints/*.fchk"))
    
    # If no files found with the first method, try alternative search patterns
    if not all_fchk_files:
        print("No .fchk files found in checkpoints directories, trying broader search...")
        all_fchk_files = list(base_path.glob("**/*.fchk"))
    
    # If still no files, try using os.walk which might be more thorough
    if not all_fchk_files:
        print("Still no .fchk files found, trying more thorough search...")
        for root, dirs, files in os.walk(base_path):
            for file in files:
                if file.endswith('.fchk'):
                    all_fchk_files.append(Path(os.path.join(root, file)))
    
    if not all_fchk_files:
        print(f"WARNING: No .fchk files found in {base_dir}")
        print("Please check that:")
        print("  1. The path to your results directory is correct")
        print("  2. ADMP calculations have completed with checkpoint files")
        print("  3. Formatted checkpoint files (.fchk) have been created")
        return []
    
    print(f"Found {len(all_fchk_files)} formatted checkpoint files")
    
    # Print some of the found files for verification
    if len(all_fchk_files) > 0:
        print("\nSample of files found:")
        for file in all_fchk_files[:min(5, len(all_fchk_files))]:
            print(f"  - {file}")
        if len(all_fchk_files) > 5:
            print(f"  ... and {len(all_fchk_files) - 5} more")
    
    return all_fchk_files

def create_cubegen_script(fchk_files, output_dir="./cube_files", 
                         grid_size=80, orbitals=None):
    """
    Create a bash script to generate cube files for specified orbitals.
    
    Args:
        fchk_files: List of paths to .fchk files
        output_dir: Directory to store cube files
        grid_size: Resolution of the cube grid (default: 80)
        orbitals: List of orbitals to generate (default: ['HOMO', 'LUMO'])
    """
    if orbitals is None:
        orbitals = ['HOMO', 'LUMO']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create the bash script
    script_path = Path("./generate_cubes.sh")
    
    with open(script_path, 'w') as script:
        script.write("#!/bin/bash\n\n")
        script.write("# This script generates cube files for molecular orbitals from ADMP simulations\n\n")
        script.write("# Load Gaussian module\n")
        script.write("module load Gaussian/g16c01-CUDA-11.7.0\n\n")
        script.write("echo 'Starting cube file generation'\n")
        script.write("echo '----------------------------'\n\n")
        
        for fchk_file in fchk_files:
            # Parse molecule name and temperature from path
            parts = str(fchk_file).split('/')
            molecule_idx = parts.index("results") + 1 if "results" in parts else -3
            temp_idx = molecule_idx + 1 if molecule_idx >= 0 else -2
            
            if molecule_idx >= 0 and temp_idx < len(parts):
                molecule = parts[molecule_idx]
                temperature = parts[temp_idx]
            else:
                # Fallback to filename parsing if directory structure is different
                basename = os.path.basename(fchk_file)
                molecule = basename.split('_ADMP_')[0] if '_ADMP_' in basename else basename.split('.')[0]
                temperature = basename.split('_ADMP_')[1].split('.')[0] if '_ADMP_' in basename else "unknown"
            
            # Create molecule/temperature specific output directory
            mol_dir = f"{output_dir}/{molecule}/{temperature}"
            script.write(f"mkdir -p {mol_dir}\n")
            
            # Get frame number or use filename if no frame info
            base_filename = os.path.basename(fchk_file)
            frame = base_filename.split('.')[0]
            
            # Generate cubegen commands for each orbital
            for orbital in orbitals:
                output_cube = f"{mol_dir}/{orbital.lower()}_{frame}.cube"
                script.write(f"echo 'Generating {orbital} cube file for {molecule} at {temperature}'\n")
                script.write(f"cubegen 0 MO={orbital} {fchk_file} {output_cube} {grid_size}\n")
            
            script.write("\n")
        
        script.write("echo 'All cube files generated successfully'\n")
    
    # Make the script executable
    os.chmod(script_path, 0o755)
    print(f"Created cubegen script at {script_path}")
    return script_path

def main():
    parser = argparse.ArgumentParser(description="Generate cube files from ADMP formatted checkpoint files")
    parser.add_argument("--base-dir", default="../ADMP_decomposition_gaussian/admp_jobs/results", 
                      help="Base directory containing ADMP results")
    parser.add_argument("--output-dir", default="./cube_files",
                      help="Directory to store cube files")
    parser.add_argument("--grid-size", type=int, default=80,
                      help="Resolution of the cube grid (default: 80)")
    parser.add_argument("--orbitals", nargs="+", default=["HOMO", "LUMO"],
                      help="Orbitals to generate (default: HOMO LUMO)")
    
    args = parser.parse_args()
    
    print("Searching for formatted checkpoint files...")
    fchk_files = find_fchk_files(args.base_dir)
    
    if fchk_files:
        script_path = create_cubegen_script(
            fchk_files, 
            output_dir=args.output_dir,
            grid_size=args.grid_size,
            orbitals=args.orbitals
        )
        
        print(f"\nNext steps:")
        print(f"1. Review the generated script: {script_path}")
        print(f"2. Run the script to generate cube files: bash {script_path}")
        print(f"3. Cube files will be created in: {args.output_dir}")
        print("\nYou can also specify additional orbitals with the --orbitals option.")
        print("Example: python generate_inputs.py --orbitals HOMO LUMO HOMO-1 LUMO+1")

if __name__ == "__main__":
    main() 