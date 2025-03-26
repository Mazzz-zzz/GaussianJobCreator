#!/usr/bin/env python

"""
Script to generate Gaussian input files for barrier energy calculations using
optimized geometries from geom_optimise_gaussian folder.
"""

import os
from pathlib import Path
import shutil
import re

def create_ts_input(reactant_path, product_path, ts_name, output_dir):
    """Create Gaussian input file for transition state search using QST3"""
    output_path = Path(output_dir) / f"{ts_name}.gjf"
    
    # Extract charge and multiplicity from reactant file
    reactant_charge = 0
    reactant_multiplicity = 1
    try:
        with open(reactant_path, 'r') as r:
            lines = r.readlines()
            for i, line in enumerate(lines):
                if line.strip() and not line.startswith("#") and not line.startswith("%"):
                    # Skip title line
                    if i > 0 and lines[i-1].strip() and not lines[i-1].startswith("#") and not lines[i-1].startswith("%"):
                        # This should be the charge and multiplicity line
                        parts = line.strip().split()
                        if len(parts) >= 2 and all(p.isdigit() for p in parts[:2]):
                            reactant_charge = int(parts[0])
                            reactant_multiplicity = int(parts[1])
                        break
    except Exception as e:
        print(f"Warning: Could not extract charge and multiplicity from reactant file: {e}")
    
    # Extract geometry data
    def extract_geometry(file_path):
        """Extract only valid atomic coordinates from a Gaussian input file."""
        geometry_lines = []
        atom_coordinate_pattern = re.compile(r'^([A-Za-z]+|\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
        
        with open(file_path, 'r') as f:
            content = f.read()
            
            # Find the charge and multiplicity line
            matches = re.findall(r'\n\s*(\d+\s+\d+)\s*\n', content)
            if not matches:
                print(f"Warning: Could not identify charge/multiplicity section in {file_path}")
                return []
                
            # The geometry section follows the charge/multiplicity line
            sections = content.split(matches[0])
            if len(sections) < 2:
                print(f"Warning: Could not split content properly for {file_path}")
                return []
                
            # Geometry is in the section after the charge/multiplicity
            geometry_section = sections[1].strip()
            
            # Process line by line to get valid atom coordinates
            for line in geometry_section.split('\n'):
                line = line.strip()
                if not line:
                    break  # Empty line marks end of coordinates
                    
                # Check if line matches atom coordinate pattern
                match = atom_coordinate_pattern.match(line)
                if match:
                    geometry_lines.append(line + '\n')
                else:
                    # If we hit a line that's not an atom coordinate, we're done
                    # with the geometry section
                    break
                    
        return geometry_lines
    
    # Get geometries
    reactant_geom = extract_geometry(reactant_path)
    product_geom = extract_geometry(product_path)
    
    # Debug output - print what we've extracted
    print(f"Extracted {len(reactant_geom)} atom lines from reactant")
    print(f"Extracted {len(product_geom)} atom lines from product")
    
    if not reactant_geom or not product_geom:
        print("WARNING: Failed to extract geometry data!")
        print(f"Reactant file: {reactant_path}")
        print(f"Product file: {product_path}")
    
    with open(output_path, 'w') as f:
        # Header with TS search specifications
        f.write("%mem=3GB\n")
        f.write("%nprocshared=4\n")
        f.write("# opt=(ts,calcfc,noeigen,qst3) freq m062x/def2tzvp int=ultrafine scf=(tight,xqc)\n\n")
        
        # Title
        f.write(f"{ts_name} transition state search\n\n")
        
        # Charge and multiplicity for reactant
        f.write(f"{reactant_charge} {reactant_multiplicity}\n")
        
        # Write reactant geometry
        for line in reactant_geom:
            f.write(line)
        f.write("\n")
        
        # Separator for product section
        f.write("--Link1--\n\n")
        
        # Add header again for product section
        f.write("%mem=3GB\n")
        f.write("%nprocshared=4\n")
        f.write("# opt=(ts,calcfc,noeigen,qst3) freq m062x/def2tzvp int=ultrafine scf=(tight,xqc)\n\n")
        f.write(f"{ts_name} product geometry\n\n")
        
        # Charge and multiplicity for product
        f.write(f"{reactant_charge} {reactant_multiplicity}\n")
        
        # Write product geometry
        for line in product_geom:
            f.write(line)
        f.write("\n")
        
        # Separator for TS guess section
        f.write("--Link1--\n\n")
        
        # Add header again for TS guess section
        f.write("%mem=3GB\n")
        f.write("%nprocshared=4\n")
        f.write("# opt=(ts,calcfc,noeigen,qst3) freq m062x/def2tzvp int=ultrafine scf=(tight,xqc)\n\n")
        f.write(f"{ts_name} TS guess geometry\n\n")
        
        # Charge and multiplicity for TS guess
        f.write(f"{reactant_charge} {reactant_multiplicity}\n")
        
        # For TS guess, we'll use the reactant geometry
        for line in reactant_geom:
            f.write(line)
        f.write("\n")
        
    # Debug - print the total size of the file
    file_size = os.path.getsize(output_path)
    print(f"Created {output_path}, size: {file_size} bytes")
    
    # If file is suspiciously small, alert
    if file_size < 500:
        print("WARNING: Output file is very small, something might be wrong!")

def setup_reaction_paths():
    """Define reaction pathways and their components"""
    return {
        "PFMS_TS1": {
            "reactant": "PFMS",
            "ts_name": "TS1M",
            "products": ["TS1M_Product1", "HF"]
        },
        "PFMS_TS2": {
            "reactant": "PFMS",
            "ts_name": "TS2M",
            "products": ["HCF3", "SO3"]
        },
        #"PFMS_TS3": {
        #    "reactant": "PFMS",
        #    "ts_name": "TS3M",
        #    "products": ["TS3M_Product1", "CF3_Radical"]
        #}
        #"PFMS_TS4": {
        #    "reactant": "PFMS",
        #    "ts_name": "TS4M",
        #    "products": ["TS4M_Product1", "F_TS4M"]
        #}
        #"PFMS_TS5": {
        #    "reactant": "PFMS",
        #    "ts_name": "TS5M",
        #    "products": ["TS5M_Product1", "F_TS5M"]
        #}
        #"PFMS_TS6": {
        #    "reactant": "PFMS",
        #    "ts_name": "TS6M",
        #    "products": ["F_TS6M", "CF2O", "SO2"]
        #}
    }

def main():
    # Create directory structure
    transition_state_dir = Path("transition_state_job")
    transition_state_dir.mkdir(exist_ok=True)
    
    base_dir = transition_state_dir 
    
    # Get optimized geometries directory
    geom_opt_dir = Path("../geom_optimise_guassian/gaussian_projects")
    
    # Check if input directory exists
    if not geom_opt_dir.exists():
        print(f"WARNING: Input directory {geom_opt_dir} does not exist!")
        print("Trying alternate path...")
        geom_opt_dir = Path("geom_optimise_guassian/gaussian_projects")
        if not geom_opt_dir.exists():
            print(f"ERROR: Cannot find input directory at {geom_opt_dir} either!")
            print("Please check path to optimized geometry files.")
            return
        else:
            print(f"Found input directory at {geom_opt_dir}")
    
    # Setup reaction paths
    reaction_paths = setup_reaction_paths()
    
    # Generate TS input files for each reaction path
    for rxn_name, components in reaction_paths.items():
        reactant_file = geom_opt_dir / f"{components['reactant']}.gjf"
        product_files = [geom_opt_dir / f"{p}.gjf" for p in components['products']]
        
        # Check if input files exist
        if not reactant_file.exists():
            print(f"ERROR: Reactant file {reactant_file} does not exist!")
            continue
            
        if not product_files[0].exists():
            print(f"ERROR: Product file {product_files[0]} does not exist!")
            continue
        
        print(f"Creating TS input for {rxn_name}...")
        print(f"  Reactant: {reactant_file}")
        print(f"  Product: {product_files[0]}")
        
        # Create TS input
        create_ts_input(
            reactant_file,
            product_files[0],  # Using first product for initial TS guess
            components['ts_name'],
            base_dir
        )
    
    print(f"Generated transition state search files in {base_dir}/")
    print("Run Gaussian calculations for all .gjf files in the transition_state_job directory")

if __name__ == "__main__":
    main()