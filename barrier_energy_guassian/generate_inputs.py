#!/usr/bin/env python

"""
Script to generate Gaussian input files for barrier energy calculations using
optimized geometries from geom_optimise_gaussian folder.
"""

import os
from pathlib import Path
import shutil

def create_ts_input(reactant_path, product_path, ts_name, output_dir):
    """Create Gaussian input file for transition state search using QST3"""
    output_path = Path(output_dir) / f"{ts_name}.gjf"
    
    with open(output_path, 'w') as f:
        # Header with TS search specifications
        f.write("%mem=3GB\n")
        f.write("%nprocshared=4\n")
        f.write("# opt=(ts,calcfc,noeigen,qst3) freq m062x/def2tzvp geom=connectivity int=ultrafine scf=(tight,xqc)\n\n")
        
        # Title
        f.write(f"{ts_name} transition state search\n\n")
        
        # Copy reactant geometry
        f.write("--Link1--\n")
        with open(reactant_path, 'r') as r:
            for line in r:
                if "opt=" not in line and "#" not in line:
                    f.write(line)
        
        # Copy product geometry
        f.write("--Link1--\n")
        with open(product_path, 'r') as p:
            for line in p:
                if "opt=" not in line and "#" not in line:
                    f.write(line)
        
        # Initial TS guess (average of reactant and product)
        f.write("--Link1--\n")
        # ... (will be generated from reactant/product interpolation)

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

def create_barrier_calculation_script(output_dir):
    """Create a Python script to calculate barrier energies from Gaussian outputs"""
    script_path = Path(output_dir) / "calculate_barriers.py"
    
    with open(script_path, 'w') as f:
        f.write("""#!/usr/bin/env python

import os
from pathlib import Path

def extract_energy(log_file):
    \"\"\"Extract final energy from Gaussian output file\"\"\"
    energy = None
    with open(log_file, 'r') as f:
        for line in f:
            if "SCF Done:" in line:
                energy = float(line.split()[4])
    return energy

def calculate_barrier(reactant_log, ts_log, product_logs):
    \"\"\"Calculate barrier energy and reaction energy\"\"\"
    reactant_e = extract_energy(reactant_log)
    ts_e = extract_energy(ts_log)
    product_e = sum(extract_energy(p) for p in product_logs)
    
    barrier = (ts_e - reactant_e) * 627.509  # Convert Hartrees to kcal/mol
    rxn_energy = (product_e - reactant_e) * 627.509
    
    return barrier, rxn_energy

def main():
    reaction_paths = {
        # Add your reaction paths here
    }
    
    results = {}
    for rxn, paths in reaction_paths.items():
        barrier, rxn_energy = calculate_barrier(
            f"{paths['reactant']}.log",
            f"{paths['ts_name']}.log",
            [f"{p}.log" for p in paths['products']]
        )
        results[rxn] = {
            'barrier': barrier,
            'reaction_energy': rxn_energy
        }
    
    # Print results
    print("\\nBarrier Energy Results (kcal/mol):")
    print("-" * 50)
    for rxn, data in results.items():
        print(f"{rxn}:")
        print(f"  Barrier Energy: {data['barrier']:.2f}")
        print(f"  Reaction Energy: {data['reaction_energy']:.2f}")

if __name__ == "__main__":
    main()
""")

def main():
    # Create directory structure
    base_dir = Path("barrier_energy_gaussian")
    base_dir.mkdir(exist_ok=True)
    
    # Get optimized geometries directory
    geom_opt_dir = Path("../geom_optimise_guassian/gaussian_projects")
    
    # Setup reaction paths
    reaction_paths = setup_reaction_paths()
    
    # Generate TS input files for each reaction path
    for rxn_name, components in reaction_paths.items():
        reactant_file = geom_opt_dir / f"{components['reactant']}.gjf"
        product_files = [geom_opt_dir / f"{p}.gjf" for p in components['products']]
        
        # Create TS input
        create_ts_input(
            reactant_file,
            product_files[0],  # Using first product for initial TS guess
            components['ts_name'],
            base_dir
        )
        
        # Copy reactant and product input files
        for file in [reactant_file] + product_files:
            shutil.copy2(file, base_dir)
    
    # Create barrier calculation script
    create_barrier_calculation_script(base_dir)
    
    print("Generated barrier energy calculation files in barrier_energy_gaussian/")
    print("1. Run Gaussian calculations for all .gjf files")
    print("2. Run calculate_barriers.py to get barrier energies")

if __name__ == "__main__":
    main()