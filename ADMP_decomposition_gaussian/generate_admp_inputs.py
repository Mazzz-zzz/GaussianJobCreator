#!/usr/bin/env python

"""
Script to generate Gaussian input files for ADMP (Atom-centered Density Matrix Propagation)
calculations at various temperatures using optimized geometries from the gaussian_projects folder.

ADMP is an ab initio molecular dynamics method that can be used to study thermal decomposition
pathways of molecules by simulating their behavior at elevated temperatures.
"""

import os
import glob
import re
from pathlib import Path

# Configuration
DEFAULT_CONFIG = {
    'input_dir': "../geom_optimise_guassian/gaussian_projects",
    'output_dir': "admp_jobs",
    'temperatures': [800],  # in Kelvin
    'max_points': 2000,
    'delta_t': 0.5,  # femtoseconds
    'method': "B3LYP",
    'basis': "6-31G(d)",
    'mem': "8GB",
    'nproc': 8,
    'rstf': 10  # Save checkpoint every n steps (ADMP RSTF parameter)
}

def extract_geometry(file_path):
    """Extract geometry, charge, and multiplicity from a Gaussian input file."""
    geometry_lines = []
    charge = 0
    multiplicity = 1
    atom_coordinate_pattern = re.compile(r'^([A-Za-z]+|\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
            # Find the charge and multiplicity line
            matches = re.findall(r'\n\s*(\d+)\s+(\d+)\s*\n', content)
            if matches:
                charge = int(matches[0][0])
                multiplicity = int(matches[0][1])
                print(f"Found charge {charge} and multiplicity {multiplicity} from {os.path.basename(file_path)}")
            else:
                print(f"Warning: Could not identify charge/multiplicity in {file_path}")
                return [], charge, multiplicity
                
            # The geometry section follows the charge/multiplicity line
            sections = content.split(f"{charge} {multiplicity}")
            if len(sections) < 2:
                print(f"Warning: Could not split content properly for {file_path}")
                return [], charge, multiplicity
                
            # Geometry is in the section after the charge/multiplicity
            geometry_section = sections[1].strip()
            
            # Process line by line to get valid atom coordinates
            for line in geometry_section.split('\n'):
                line = line.strip()
                if not line:
                    continue  # Skip empty lines
                    
                # Check if line matches atom coordinate pattern
                match = atom_coordinate_pattern.match(line)
                if match:
                    geometry_lines.append(line)
                elif len(geometry_lines) > 0:
                    # If we hit a line that's not an atom coordinate after 
                    # we've started collecting geometries, we're done
                    break
                    
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return [], charge, multiplicity
        
    return geometry_lines, charge, multiplicity


def create_admp_input(molecule_path, temp, output_dir, 
                     max_points=2000, delta_t=0.5, full_step=100,
                     method='B3LYP', basis='6-31G(d)', mem='8GB', nproc=8, rstf=10):
    """Create Gaussian input file for ADMP calculation at specified temperature.
    Note: Temperature control in ADMP is achieved through initial velocities,
    not through direct parameters to the ADMP keyword."""
    
    molecule_name = os.path.basename(molecule_path).replace('.gjf', '')
    output_path = Path(output_dir) / f"{molecule_name}_ADMP_{temp}K.gjf"
    
    # Extract geometry, charge and multiplicity
    geometry, charge, multiplicity = extract_geometry(molecule_path)
    
    if not geometry:
        print(f"WARNING: No geometry found in {molecule_path}, skipping.")
        return False
        
    with open(output_path, 'w') as f:
        # Header with proper ADMP specifications based on Gaussian documentation
        f.write(f"%mem={mem}\n")
        f.write(f"%nprocshared={nproc}\n")
        # Write route section with ADMP keyword - note that temperature is controlled via initial velocities
        #f.write(f"# {method}/{basis} ADMP int=ultrafine\n\n")
        f.write(f"# {method}/{basis} ADMP int=ultrafine Temperature={temp}\n\n")
        
        # Title
        f.write(f"{molecule_name} ADMP thermal decomposition simulation targeting {temp}K\n\n")
        
        # Charge and multiplicity
        f.write(f"{charge} {multiplicity}\n")
        
        # Write geometry
        for line in geometry:
            f.write(f"{line}\n")
            
        # Note: Initial velocities would need to be added here for temperature control
        # This requires additional implementation to generate appropriate velocities
        # based on the desired temperature
        
        # End file with newline
        f.write("\n")
        
    # Verify the file content to ensure no syntax errors
    try:
        with open(output_path, 'r') as f:
            content = f.read()
            if "'" in content or '"' in content:
                print(f"WARNING: Unexpected quote characters found in {output_path}")
                # Remove any quote characters
                content = content.replace("'", "").replace('"', "")
                with open(output_path, 'w') as f_clean:
                    f_clean.write(content)
                print(f"Cleaned quote characters from {output_path}")
    except Exception as e:
        print(f"Error verifying file {output_path}: {e}")
        
    print(f"Created {output_path}")
    return True


def main():
    config = DEFAULT_CONFIG
    
    # Create output directory
    output_dir = Path(config['output_dir'])
    output_dir.mkdir(exist_ok=True, parents=True)
    print(f"Output will be saved to {output_dir}")
    
    # Find input files
    potential_paths = [
        Path(config['input_dir']),
        Path("../geom_optimise_guassian/gaussian_projects"),
        Path("geom_optimise_guassian/gaussian_projects"),
        Path("../geom_optimise_gaussian/gaussian_projects"),  # Note the correct spelling
        Path("geom_optimise_gaussian/gaussian_projects")      # Note the correct spelling
    ]
    
    input_dir = None
    for path in potential_paths:
        if path.exists():
            input_dir = path
            print(f"Found input directory at {input_dir}")
            break
    
    if not input_dir:
        print("ERROR: Cannot find input directory for optimized geometries!")
        print("Searched the following paths:")
        for path in potential_paths:
            print(f"  - {path}")
        return
    
    # Find all .gjf and .com files
    input_files = list(input_dir.glob("*.gjf")) + list(input_dir.glob("*.com"))
    
    if not input_files:
        print(f"No input files found in {input_dir}")
        return
        
    print(f"Found {len(input_files)} input files.")
    
    # Create directory structure for temperatures
    for temp in config['temperatures']:
        temp_dir = output_dir / f"{temp}K"
        temp_dir.mkdir(exist_ok=True)
        
        print(f"\nGenerating ADMP inputs for {temp}K:")
        
        # Generate ADMP inputs for each molecule at this temperature
        for molecule_path in input_files:
            create_admp_input(
                molecule_path=molecule_path,
                temp=temp,
                output_dir=temp_dir,
                max_points=config['max_points'],
                delta_t=config['delta_t'],
                method=config['method'],
                basis=config['basis'],
                mem=config['mem'],
                nproc=config['nproc'],
                rstf=config['rstf']
            )
            
    print(f"\nGenerated ADMP input files for {len(input_files)} molecules at {len(config['temperatures'])} temperatures.")
    print("Run these Gaussian calculations to simulate thermal decomposition processes.")
    print("\nTo analyze results:")
    print("1. Extract snapshots from ADMP trajectories at points where bonds break")
    print("2. Use these geometries as starting points for transition state searches")
    print("3. Perform IRC calculations to confirm decomposition pathways")


if __name__ == "__main__":
    main() 