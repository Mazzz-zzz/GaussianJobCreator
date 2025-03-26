#!/usr/bin/env python

"""
Script to generate Gaussian QST2 transition state search input files
using existing optimized geometries from .gjf files
"""

import os
from pathlib import Path
import argparse
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
import numpy as np

# Configuration
DEFAULT_CONFIG = {
    # Directory paths
    'input_dir': "../geom_optimise_guassian/gaussian_projects", 
    'output_dir': "qst2_jobs",
    
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

def parse_gaussian_input(gjf_file):
    """Parse a Gaussian .gjf input file to extract geometry and other parameters"""
    with open(gjf_file, 'r') as f:
        lines = f.readlines()
    
    # Find charge and multiplicity line
    charge_mult_line = None
    for i, line in enumerate(lines):
        if re.match(r'^[+-]?\d+\s+\d+\s*$', line.strip()):
            charge_mult_line = i
            break
    
    if charge_mult_line is None:
        raise ValueError(f"Could not find charge/multiplicity line in {gjf_file}")
    
    charge, multiplicity = map(int, lines[charge_mult_line].split())
    
    # Extract coordinates - start from charge/multiplicity line and continue until blank line
    coordinates = []
    i = charge_mult_line + 1
    while i < len(lines) and lines[i].strip():
        coordinates.append(lines[i])
        i += 1
    
    # Look for connectivity information (if present)
    connectivity = []
    while i < len(lines):
        if lines[i].strip() and lines[i].strip()[0].isdigit():
            connectivity.append(lines[i])
        i += 1
        if i < len(lines) and not lines[i].strip():
            break
    
    return {
        'charge': charge,
        'multiplicity': multiplicity,
        'coordinates': coordinates,
        'connectivity': connectivity
    }

def align_coordinates(coordinate_lines):
    symbols = []
    positions = []
    for line in coordinate_lines:
        tokens = line.split()
        if len(tokens) < 4:
            continue
        symbols.append(tokens[0])
        positions.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
    atoms = Atoms(symbols=symbols, positions=positions)
    com = atoms.get_center_of_mass()
    atoms.translate(-com)
    aligned_lines = []
    for s, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        aligned_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    return aligned_lines

def kabsch_align(ref_coords, coords):
    # Ensure the arrays have compatible shapes
    if ref_coords.shape != coords.shape:
        raise ValueError(f"Incompatible shapes: reactant {ref_coords.shape} vs product {coords.shape}")
    
    if len(ref_coords) < 3:
        print("WARNING: Kabsch alignment needs at least 3 atoms for proper 3D alignment")
        print("Using identity rotation matrix instead")
        return np.eye(3)
    
    # Compute covariance matrix
    H = np.dot(coords.T, ref_coords)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = np.dot(Vt.T, U.T)
    return R

def align_product_to_reactant(reactant_lines, product_lines):
    # Parse reactant coordinates
    ref_symbols = []
    ref_positions = []
    for line in reactant_lines:
        tokens = line.split()
        if len(tokens) < 4:
            continue
        ref_symbols.append(tokens[0])
        ref_positions.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
    ref_positions = np.array(ref_positions)
    ref_com = np.mean(ref_positions, axis=0)
    ref_centered = ref_positions - ref_com

    # Parse product coordinates
    prod_symbols = []
    prod_positions = []
    for line in product_lines:
        tokens = line.split()
        if len(tokens) < 4:
            continue
        prod_symbols.append(tokens[0])
        prod_positions.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
    prod_positions = np.array(prod_positions)
    
    # Check if reactant and product have the same number of atoms
    if len(ref_symbols) != len(prod_symbols):
        print(f"WARNING: Reactant has {len(ref_symbols)} atoms but product has {len(prod_symbols)} atoms")
        print("QST2 requires the same atoms in the same order for both reactant and product")
        print("Attempting to use just the first min(len(reactant), len(product)) atoms...")
        
        # Use only the minimum number of atoms from both structures
        min_atoms = min(len(ref_symbols), len(prod_symbols))
        ref_symbols = ref_symbols[:min_atoms]
        ref_positions = ref_positions[:min_atoms]
        ref_centered = ref_centered[:min_atoms]
        prod_symbols = prod_symbols[:min_atoms]
        prod_positions = prod_positions[:min_atoms]
    
    # Check if the atom types match in order
    for i, (ref, prod) in enumerate(zip(ref_symbols, prod_symbols)):
        if ref != prod:
            print(f"WARNING: Atom mismatch at position {i+1}: reactant has {ref} but product has {prod}")
            print("QST2 requires identical atom ordering between reactant and product")
            print("This may affect the quality of the transition state search")
    
    # Center the product coordinates
    prod_com = np.mean(prod_positions, axis=0)
    prod_centered = prod_positions - prod_com
    
    try:
        # Compute rotation matrix to align product to reactant
        R = kabsch_align(ref_centered, prod_centered)
        prod_aligned = np.dot(prod_centered, R)
        
        # Format aligned product coordinates
        aligned_lines = []
        for s, pos in zip(prod_symbols, prod_aligned):
            aligned_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    except Exception as e:
        print(f"WARNING: Could not align product to reactant: {e}")
        print("Using centered but unaligned product coordinates instead")
        
        # Use centered but unaligned coordinates if alignment fails
        aligned_lines = []
        for s, pos in zip(prod_symbols, prod_centered):
            aligned_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    
    return aligned_lines

def ensure_matching_atoms(reactant_file, product_file, temp_dir=None):
    """
    Ensure that reactant and product have the same atoms in the same order.
    Instead of creating temp files, returns matched data structures.
    
    Parameters:
    -----------
    reactant_file : str
        Path to the reactant .gjf file
    product_file : str
        Path to the product .gjf file
    temp_dir : str
        Unused parameter, kept for backward compatibility
    
    Returns:
    --------
    tuple:
        (reactant_data, product_data) with matched atoms
    """
    # Parse input files
    reactant_data = parse_gaussian_input(reactant_file)
    product_data = parse_gaussian_input(product_file)
    
    # Extract atom symbols from coordinates
    reactant_symbols = []
    for line in reactant_data['coordinates']:
        tokens = line.split()
        if len(tokens) >= 1:
            reactant_symbols.append(tokens[0])
    
    product_symbols = []
    for line in product_data['coordinates']:
        tokens = line.split()
        if len(tokens) >= 1:
            product_symbols.append(tokens[0])
    
    # Check if the atom lists match
    if len(reactant_symbols) == len(product_symbols):
        # Check if atom types match in the same order
        matching = True
        for i, (r, p) in enumerate(zip(reactant_symbols, product_symbols)):
            if r != p:
                print(f"WARNING: Atom mismatch at position {i+1}: reactant has {r} but product has {p}")
                matching = False
        
        if matching:
            print("Reactant and product have matching atoms in the same order - good for QST2")
            return (reactant_data, product_data)  # Return original data if matching
    else:
        print(f"WARNING: Reactant has {len(reactant_symbols)} atoms but product has {len(product_symbols)} atoms")
        print("QST2 requires the same atoms in the same order for both reactant and product")
    
    # If we're here, the atoms don't match - we need to match them in memory
    print("Matching atoms for QST2 compatibility...")
    
    # For simplicity, we'll use the atoms from the reactant for both files
    # Keep the reactant data as is
    matched_product_data = product_data.copy()
    
    # Create matched product coordinates
    matched_product_coords = []
    min_atoms = min(len(reactant_symbols), len(product_symbols))
    
    for i in range(min_atoms):
        product_line = product_data['coordinates'][i]
        tokens = product_line.split()
        if len(tokens) >= 4:
            coord_part = " ".join(tokens[1:])
            matched_product_coords.append(f"{reactant_symbols[i]} {coord_part}\n")
    
    # If reactant has more atoms, add them with product's last atom coordinates
    if len(reactant_symbols) > min_atoms:
        last_product_coords = product_data['coordinates'][min_atoms-1].split()[1:]
        for i in range(min_atoms, len(reactant_symbols)):
            matched_product_coords.append(f"{reactant_symbols[i]} {' '.join(last_product_coords)}\n")
    
    # Update the product data with matched coordinates
    matched_product_data['coordinates'] = matched_product_coords
    
    print("Atoms matched for QST2 compatibility")
    
    return (reactant_data, matched_product_data)

def create_qst2_input(reactant_file, product_file, output_dir, config=None):
    """
    Create a Gaussian QST2 input file using geometries from reactant and product files
    
    Parameters:
    -----------
    reactant_file : str
        Path to the reactant .gjf file
    product_file : str
        Path to the product .gjf file
    output_dir : str
        Directory to save the QST2 input file
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
    output_file = os.path.join(output_dir, f"TS_{reactant_name}_to_{product_name}.gjf")
    
    # Check if charge and multiplicity match
    if reactant_data['charge'] != product_data['charge'] or reactant_data['multiplicity'] != product_data['multiplicity']:
        print(f"WARNING: Charge/multiplicity mismatch between reactant and product:")
        print(f"  Reactant: {reactant_data['charge']} {reactant_data['multiplicity']}")
        print(f"  Product:  {product_data['charge']} {product_data['multiplicity']}")
        print(f"QST2 requires same charge and multiplicity for both structures.")
    
    with open(output_file, 'w') as f:
        # Write header and route
        f.write(f"%mem={config['mem']}\n")
        f.write(f"%nprocshared={config['nproc']}\n")
        f.write(f"%chk=TS_{reactant_name}_to_{product_name}.chk\n")
        
        # Build route line
        if config['redundant']:
            route = f"#p {config['level_of_theory']} opt=(QST2,redundant)"
        else:
            route = f"#p {config['level_of_theory']} opt=QST2"
            
        if config['add_freq']:
            route += " Freq=NoRaman"
        if config['geom_connect']:
            route += " Geom=Connect"
        if config['temperature']:
            route += f" Temperature={config['temperature']}"
        if config['solvent']:
            route += f" SCRF=(Solvent={config['solvent']})"
        
        f.write(f"{route}\n\n")
        
        # First title section for reactants
        f.write(f"Reactants: {reactant_name}\n\n")
        
        # Write charge/multiplicity and coordinates for reactant
        f.write(f"{int(reactant_data['charge'])} {int(reactant_data['multiplicity'])}\n")
        reactant_aligned = align_coordinates(reactant_data['coordinates'])
        for line in reactant_aligned:
            f.write(line)
        
        # Add blank line after reactant coordinates
        f.write("\n")
        
        # Second title section for products
        f.write(f"Products: {product_name}\n\n")
        
        # Align product coordinates to reactant and write charge/multiplicity and coordinates
        product_aligned = align_product_to_reactant(reactant_data['coordinates'], product_data['coordinates'])
        f.write(f"{int(product_data['charge'])} {int(product_data['multiplicity'])}\n")
        for line in product_aligned:
            f.write(line)
        
        # Ensure a blank line at the end of the file
        f.write("\n")
    
    print(f"Created QST2 input file: {output_file}")
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
    Main function to generate QST2 input files based on configuration
    
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
        
        create_qst2_input(
            reactant_file=reactant_file,
            product_file=product_file,
            output_dir=config['output_dir'],
            config=config
        )

    if config['verbose']:
        print(f"\nQST2 input generation complete. Files saved to {config['output_dir']}")
        print("Run these Gaussian calculations to locate transition states.")
        print("\nTips for verifying transition states:")
        print("1. Check that the optimization converged to a stationary point")
        print("2. Verify that the structure has exactly ONE imaginary frequency")
        print("3. Examine the imaginary mode to confirm it corresponds to the reaction coordinate")
        print("4. For confirmation, perform IRC calculations from the transition state")

if __name__ == "__main__":
    main()
