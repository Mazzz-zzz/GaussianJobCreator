#!/usr/bin/env python

"""
Script to generate Gaussian QST3 transition state search input files
using existing optimized geometries from .gjf files and generating an initial 
transition state guess using linear synchronous transit (LST) method.
"""

import os
from pathlib import Path
import re

# Import from our custom modules
from align_molecules import align_coordinates, ensure_matching_atoms
from make_ts_guess import generate_transition_state_guess

# Configuration
DEFAULT_CONFIG = {
    # Directory paths
    'input_dir': "../geom_optimise_guassian/gaussian_projects", 
    'output_dir': "qst3_jobs",
    
    # Gaussian calculation parameters
    'level_of_theory': "b3lyp/6-31G(d)",
    'mem': "8GB",
    'nproc': 8,
    'add_freq': True,    # Add frequency calculation
    'geom_connect': False, # Use connectivity information
    'redundant': False,   # Use redundant internal coordinates
    
    # Calculation options
    'temperature': None,  # Optional temperature for calculations
    'solvent': None,     # Optional solvent model
    
    # Execution mode
    'mode': 'all',       # 'all', 'reaction', or 'list'
    'reaction': None,    # Specific reaction to process (if mode='reaction')
    
    # Output control
    'verbose': True      # Print detailed progress information
}

def ensure_directory(directory):
    """Create directory if it doesn't exist"""
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created directory: {directory}")

def create_qst3_input(reactant_file, product_file, output_dir, config=None):
    """
    Create a Gaussian QST3 input file using geometries from reactant and product files
    and generating a transition state guess.
    
    Parameters:
    -----------
    reactant_file : str
        Path to the reactant .gjf file
    product_file : str
        Path to the product .gjf file
    output_dir : str
        Directory to save the QST3 input file
    config : dict
        Configuration dictionary with calculation parameters
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    # Ensure reactant and product have matching atoms
    reactant_data, product_data = ensure_matching_atoms(reactant_file, product_file)
    
    # Get base names for title
    reactant_name = os.path.basename(reactant_file).replace('.gjf', '')
    product_name = os.path.basename(product_file).replace('.gjf', '')
    
    # Create output filename
    output_file = os.path.join(output_dir, f"TS_{reactant_name}_to_{product_name}_QST3.gjf")
    
    # Check if charge and multiplicity match
    if reactant_data['charge'] != product_data['charge'] or reactant_data['multiplicity'] != product_data['multiplicity']:
        print(f"WARNING: Charge/multiplicity mismatch between reactant and product:")
        print(f"  Reactant: {reactant_data['charge']} {reactant_data['multiplicity']}")
        print(f"  Product:  {product_data['charge']} {product_data['multiplicity']}")
        print(f"QST3 requires same charge and multiplicity for both structures.")
    
    # Generate transition state guess using the specified method
    ts_guess_coordinates = generate_transition_state_guess(
        reactant_data, 
        product_data, 
        method="lst"
    )
    
    with open(output_file, 'w') as f:
        # Write header and route
        f.write(f"%mem={config['mem']}\n")
        f.write(f"%nprocshared={config['nproc']}\n")
        f.write(f"%chk=TS_{reactant_name}_to_{product_name}_QST3.chk\n")
        
        # Build route line
        if config['redundant']:
            route = f"#p {config['level_of_theory']} opt=(QST3,redundant)"
        else:
            route = f"#p {config['level_of_theory']} opt=QST3"
            
        if config['add_freq']:
            route += " Freq=NoRaman"
        if config['geom_connect']:
            route += " Geom=Connect"
        if config['temperature']:
            route += f" Temperature={config['temperature']}"
        if config['solvent']:
            route += f" SCRF=(Solvent={config['solvent']})"
        
        f.write(f"{route}\n\n")
        
        # Title card for reagents (reactants)
        f.write(f"Reactants: {reactant_name}\n\n")
        
        # Write charge/multiplicity and coordinates for reactant
        f.write(f"{int(reactant_data['charge'])} {int(reactant_data['multiplicity'])}\n")
        reactant_aligned = align_coordinates(reactant_data['coordinates'])
        for line in reactant_aligned:
            f.write(line)
        
        # Add blank line after reactant coordinates
        f.write("\n")
        
        # Title card for products
        f.write(f"Products: {product_name}\n\n")
        
        # Align product coordinates to reactant and write charge/multiplicity and coordinates
        from align_molecules import align_product_to_reactant
        product_aligned, _, _, _ = align_product_to_reactant(
            reactant_data['coordinates'], 
            product_data['coordinates']
        )
        f.write(f"{int(product_data['charge'])} {int(product_data['multiplicity'])}\n")
        for line in product_aligned:
            f.write(line)
        
        # Add blank line after product coordinates
        f.write("\n")
        
        # Title card for transition state guess
        f.write(f"Transition State Guess: {reactant_name} to {product_name}\n\n")
        
        # Write TS guess structure with same charge/multiplicity
        f.write(f"{int(reactant_data['charge'])} {int(reactant_data['multiplicity'])}\n")
        for line in ts_guess_coordinates:
            f.write(line)
        
        # Ensure a blank line at the end of the file
        f.write("\n")
    
    print(f"Created QST3 input file: {output_file}")
    return output_file

def list_available_molecules(input_dir):
    """List all available molecules from .gjf files in the input directory"""
    gjf_files = list(Path(input_dir).glob("*.gjf"))
    print("Available molecules:")
    for i, file in enumerate(gjf_files):
        print(f"{i+1}. {file.stem}")
    return gjf_files

def setup_reaction_pairs():
    """Dictionary of common reaction pairs based on molecules in generate_inputs.py"""
    return {
        "SN2_halogen_exchange_F_Cl": {
            "reactants": ["Na", "CH3F"],
            "products": ["CH3Cl", "Na"]
        },
        "SN2_halogen_exchange_F_Br": {
            "reactants": ["Na", "CH3F"],
            "products": ["CH3Br", "Na"]
        },
        "TS1M_PFMS_Pathway": {
            "reactants": ["PFMS"],
            "products": ["TS1M_Product1", "HF"]
        },
        "TS2M_PFMS_Pathway": {
            "reactants": ["PFMS"],
            "products": ["HCF3", "SO3"]
        },
        "TS3M_PFMS_Pathway": {
            "reactants": ["ISOPFMS"],
            "products": ["TS3M_Product1", "CF3_Radical"]
        },
        "TS4M_PFMS_Pathway": {
            "reactants": ["PFMS"],
            "products": ["TS4M_Product1", "F_TS4M"]
        },
        "TS5M_PFMS_Pathway": {
            "reactants": ["PFMS"],
            "products": ["TS5M_Product1", "F_TS5M"]
        },
        "TS6M_PFMS_Pathway": {
            "reactants": ["PFMS"],
            "products": ["F_TS6M", "CF2O", "SO2"]
        }
    }

def main(config=None):
    """
    Main function to generate QST3 input files based on configuration
    
    Parameters:
    -----------
    config : dict, optional
        Configuration dictionary. If None, uses DEFAULT_CONFIG
    """
    # Use default configuration if none provided
    if config is None:
        config = DEFAULT_CONFIG.copy()
    
    # Create output directory
    ensure_directory(config['output_dir'])
    
    # Get reaction pairs
    reaction_pairs = setup_reaction_pairs()
    
    # Handle different execution modes
    if config['mode'] == 'list':
        print("Available predefined reactions:")
        for i, name in enumerate(reaction_pairs.keys()):
            reactants = ", ".join(reaction_pairs[name]["reactants"])
            products = ", ".join(reaction_pairs[name]["products"])
            print(f"{i+1}. {name}: {reactants} -> {products}")
        return
    
    # Find input directory
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
            if config['verbose']:
                print(f"Found input directory at {input_dir}")
            break
    
    if not input_dir:
        print("ERROR: Cannot find input directory for optimized geometries!")
        print("Searched the following paths:")
        for path in potential_paths:
            print(f"  - {path}")
        return
        
    config['input_dir'] = str(input_dir)
    
    # Process specific reaction or all predefined reactions
    reactions_to_process = []
    if config['mode'] == 'reaction':
        if config['reaction'] not in reaction_pairs:
            print(f"Error: Reaction '{config['reaction']}' not found in predefined reactions")
            print("Use mode='list' to see available reactions")
            return
        reactions_to_process = [config['reaction']]
    else:
        # Process all predefined reactions
        reactions_to_process = list(reaction_pairs.keys())
    
    for reaction_name in reactions_to_process:
        if config['verbose']:
            print(f"\nProcessing reaction: {reaction_name}")
        reaction = reaction_pairs[reaction_name]
        
        # For reactions with multiple reactants/products, we need to create a combined structure
        # This is just a placeholder - in practice, you would need to align these molecules properly
        # For now, we'll just use the first reactant and first product
        reactant_file = os.path.join(config['input_dir'], f"{reaction['reactants'][0]}.gjf")
        product_file = os.path.join(config['input_dir'], f"{reaction['products'][0]}.gjf")
        
        if not os.path.exists(reactant_file):
            print(f"Warning: Reactant file {reactant_file} not found")
            continue
            
        if not os.path.exists(product_file):
            print(f"Warning: Product file {product_file} not found")
            continue
        
        # Always use LST method
        if config['verbose']:
            print("\nGenerating transition state guess using Linear Synchronous Transit (LST)...")
        create_qst3_input(
            reactant_file=reactant_file,
            product_file=product_file,
            output_dir=config['output_dir'],
            config=config
        )

    if config['verbose']:
        print(f"\nQST3 input generation complete. Files saved to {config['output_dir']}")
        print("Run these Gaussian calculations to locate transition states.")
        print("\nTips for verifying transition states:")
        print("1. Check that the optimization converged to a stationary point")
        print("2. Verify that the structure has exactly ONE imaginary frequency")
        print("3. Examine the imaginary mode to confirm it corresponds to the reaction coordinate")
        print("4. For confirmation, perform IRC calculations from the transition state")

if __name__ == "__main__":
    main() 