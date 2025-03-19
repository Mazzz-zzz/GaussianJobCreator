#!/usr/bin/env python3
"""
Generate orbital cube files from ADMP formatted checkpoint files (.fchk).

This script:
1. Searches for formatted checkpoint files in ADMP results directories
2. Creates a SLURM submission script to run cubegen for HOMO and LUMO orbitals
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

def create_slurm_script(fchk_files, output_dir="./cube_files", 
                        grid_size=80, orbitals=None, max_time="12:00:00"):
    """
    Create a SLURM submission script to process all checkpoint files with cubegen.
    
    Args:
        fchk_files: List of paths to .fchk files
        output_dir: Directory to store cube files
        grid_size: Resolution of the cube grid (default: 80)
        orbitals: List of orbitals to generate (default: ['HOMO', 'LUMO'])
        max_time: Maximum time allowed for the SLURM job (default: "12:00:00")
    """
    if orbitals is None:
        orbitals = ['HOMO', 'LUMO']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create the SLURM submission script
    script_path = Path("./submit_orbital_generation.s")
    
    with open(script_path, 'w') as script:
        # Write SLURM header
        script.write("""#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time={0}
#SBATCH --mem=5G
#SBATCH --partition=sapphire
#SBATCH --job-name=ADMP_orbitals
#SBATCH --output=ADMP_orbitals_%j.log

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0

echo "Starting orbital visualization generation from ADMP checkpoint files"
echo "-------------------------------------------------------------------"

""".format(max_time))

        # Function to validate and prepare output directories
        script.write("""# Function to check for existing cube files
check_existing_cube() {
    local cube_file="$1"
    local fchk_file="$2"
    local orbital="$3"
    
    if [ -f "$cube_file" ]; then
        echo "  - Cube file for $orbital already exists: $cube_file"
        return 0
    else
        echo "  - Generating $orbital cube file from: $fchk_file"
        return 1
    fi
}

""")

        # Process all checkpoint files
        script.write("echo \"Processing " + str(len(fchk_files)) + " checkpoint files\"\n")
        script.write("echo \"----------------------------------------\"\n\n")
        
        for i, fchk_file in enumerate(fchk_files):
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
            script.write(f"echo \"Processing file {i+1}/{len(fchk_files)}: {fchk_file}\"\n")
            script.write(f"mkdir -p \"{mol_dir}\"\n")
            
            # Get frame number or use filename if no frame info
            base_filename = os.path.basename(fchk_file)
            frame = base_filename.split('.')[0]
            
            # Generate cubegen commands for each orbital
            for orbital in orbitals:
                output_cube = f"{mol_dir}/{orbital.lower()}_{frame}.cube"
                
                # Add check to skip existing files
                script.write(f"\nif ! check_existing_cube \"{output_cube}\" \"{fchk_file}\" \"{orbital}\"; then\n")
                script.write(f"    echo \"  Running: cubegen 0 MO={orbital} {fchk_file} {output_cube} {grid_size} h\"\n")
                script.write(f"    cubegen 0 MO={orbital} {fchk_file} {output_cube} {grid_size} h\n")
                script.write(f"    \n")
                script.write(f"    if [ -f \"{output_cube}\" ]; then\n")
                script.write(f"        echo \"  ✓ Successfully created {orbital} cube file\"\n")
                script.write(f"    else\n")
                script.write(f"        echo \"  ✗ ERROR: Failed to create {orbital} cube file\"\n")
                script.write(f"    fi\n")
                script.write(f"fi\n")
            
            script.write("\necho \"----------------------------------------\"\n")
        
        # Add summary section
        script.write("""
echo "Orbital visualization generation completed"
echo "Checking results..."

# Count generated cube files
CUBE_COUNT=$(find "${0}" -name "*.cube" | wc -l)
echo "Number of cube files generated: $CUBE_COUNT"

if [ "$CUBE_COUNT" -gt 0 ]; then
    echo "Cube files were generated in the following locations:"
    find "${0}" -type d -not -path "${0}" | sort | head -n 10
    DIR_COUNT=$(find "${0}" -type d -not -path "${0}" | wc -l)
    if [ "$DIR_COUNT" -gt 10 ]; then
        echo "... and $(($DIR_COUNT - 10)) more directories"
    fi
else
    echo "WARNING: No cube files were generated. Check the log for errors."
fi

echo "Job completed"

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s
""".format(output_dir))
    
    # Make the script executable
    os.chmod(script_path, 0o755)
    print(f"Created SLURM submission script at {script_path}")
    
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
    parser.add_argument("--max-time", default="12:00:00",
                      help="Maximum time for SLURM job (default: 12:00:00)")
    
    args = parser.parse_args()
    
    print("Searching for formatted checkpoint files...")
    fchk_files = find_fchk_files(args.base_dir)
    
    if fchk_files:
        script_path = create_slurm_script(
            fchk_files, 
            output_dir=args.output_dir,
            grid_size=args.grid_size,
            orbitals=args.orbitals,
            max_time=args.max_time
        )
        
        print(f"\nNext steps:")
        print(f"1. Submit the job to the cluster with: sbatch {script_path}")
        print(f"2. Cube files will be created in: {args.output_dir}")
        print("\nYou can also specify additional orbitals with the --orbitals option.")
        print("Example: python generate_inputs.py --orbitals HOMO LUMO HOMO-1 LUMO+1")
    else:
        print("No formatted checkpoint files found. Cannot generate orbital visualization.")

if __name__ == "__main__":
    main() 