#!/usr/bin/env python

"""
Script to generate Gaussian QST3 transition state search input files
using existing optimized geometries from .gjf files and generating an initial 
transition state guess using various methods.
"""

import os
from pathlib import Path
import argparse
import re
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.geometry import get_distances
from scipy.interpolate import CubicSpline
import itertools

# Configuration
DEFAULT_CONFIG = {
    'input_dir': "../geom_optimise_guassian/gaussian_projects", 
    'output_dir': "qst3_jobs",
    'level_of_theory': "b3lyp/6-31G(d)",
    'mem': "8GB",
    'nproc': 8,
    'add_freq': True,    # Add frequency calculation
    'geom_connect': True, # Use connectivity information
    'ts_guess_method': 'lst', # Default transition state guess method
    'temperature': None,  # Optional temperature for calculations
    'solvent': None,      # Optional solvent model
    'redundant': False    # Use redundant internal coordinates
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

def get_atoms_from_coordinates(coordinate_lines):
    """Extract atoms and coordinates from coordinate lines"""
    symbols = []
    positions = []
    for line in coordinate_lines:
        tokens = line.split()
        if len(tokens) < 4:
            continue
        symbols.append(tokens[0])
        positions.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
    atoms = Atoms(symbols=symbols, positions=positions)
    return atoms

def align_coordinates(coordinate_lines):
    """Align molecule coordinates to center of mass and return formatted lines"""
    atoms = get_atoms_from_coordinates(coordinate_lines)
    com = atoms.get_center_of_mass()
    atoms.translate(-com)
    aligned_lines = []
    for s, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        aligned_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    return aligned_lines

def kabsch_align(ref_coords, coords):
    """Align coordinates using the Kabsch algorithm"""
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
    """Align product coordinates to reactant for proper comparison"""
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
        print("QST3 requires the same atoms in the same order for both reactant and product")
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
            print("QST3 requires identical atom ordering between reactant and product")
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
    
    return aligned_lines, ref_symbols, ref_positions, prod_aligned

def identify_key_bond_changes(reactant_data, product_data, threshold=0.4):
    """
    Identify key bonds that are changing during the reaction
    Returns a list of atom pairs (indices) that show significant distance changes
    """
    # Get atom coordinates from both structures
    reactant_atoms = get_atoms_from_coordinates(reactant_data['coordinates'])
    product_atoms = get_atoms_from_coordinates(product_data['coordinates'])
    
    # Ensure they have the same number of atoms
    if len(reactant_atoms) != len(product_atoms):
        min_atoms = min(len(reactant_atoms), len(product_atoms))
        reactant_atoms = Atoms(
            symbols=reactant_atoms.get_chemical_symbols()[:min_atoms],
            positions=reactant_atoms.get_positions()[:min_atoms]
        )
        product_atoms = Atoms(
            symbols=product_atoms.get_chemical_symbols()[:min_atoms],
            positions=product_atoms.get_positions()[:min_atoms]
        )
    
    # Get all pairwise distances in both structures
    r_distances = reactant_atoms.get_all_distances()
    p_distances = product_atoms.get_all_distances()
    
    # Find pairs with significant distance changes
    key_pairs = []
    n_atoms = len(reactant_atoms)
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            dist_diff = abs(r_distances[i, j] - p_distances[i, j])
            if dist_diff > threshold:
                key_pairs.append((i, j))
    
    return key_pairs, reactant_atoms, product_atoms

def generate_ts_guess_lst(reactant_data, product_data):
    """
    Generate a transition state guess using Linear Synchronous Transit (LST)
    This is a simple interpolation halfway between reactant and product
    """
    # First align the product to the reactant
    aligned_results = align_product_to_reactant(
        reactant_data['coordinates'], 
        product_data['coordinates']
    )
    
    aligned_prod_lines, ref_symbols, ref_positions, prod_aligned_positions = aligned_results
    
    # Simple linear interpolation at halfway point (t=0.5)
    t = 0.5  # halfway point
    ts_positions = ref_positions + t * (prod_aligned_positions - ref_positions)
    
    # Format the resulting TS guess coordinates
    ts_lines = []
    for s, pos in zip(ref_symbols, ts_positions):
        ts_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    
    return ts_lines

def generate_ts_guess_qst(reactant_data, product_data, num_points=20):
    """
    Generate a transition state guess using Quadratic Synchronous Transit (QST)
    This performs a more sophisticated interpolation to estimate the highest energy point
    
    In a real implementation, we would compute the energy at each point using quantum chemistry.
    For this example, we'll use a heuristic to approximate where the maximum might be.
    """
    # First align the product to the reactant
    aligned_results = align_product_to_reactant(
        reactant_data['coordinates'], 
        product_data['coordinates']
    )
    
    aligned_prod_lines, ref_symbols, ref_positions, prod_aligned_positions = aligned_results
    
    # Identify key bond changes between reactant and product
    key_pairs, r_atoms, p_atoms = identify_key_bond_changes(reactant_data, product_data)
    
    if not key_pairs:
        print("WARNING: No significant bond changes detected, using simple LST interpolation")
        return generate_ts_guess_lst(reactant_data, product_data)
    
    # Choose the pair with the largest distance change as our reaction coordinate
    max_change = 0
    key_pair = key_pairs[0]
    
    for i, j in key_pairs:
        r_dist = np.linalg.norm(ref_positions[i] - ref_positions[j])
        p_dist = np.linalg.norm(prod_aligned_positions[i] - prod_aligned_positions[j])
        change = abs(r_dist - p_dist)
        
        if change > max_change:
            max_change = change
            key_pair = (i, j)
    
    print(f"Identified key bond change between atoms {key_pair[0]+1} and {key_pair[1]+1}")
    print(f"Distance changes from {np.linalg.norm(ref_positions[key_pair[0]] - ref_positions[key_pair[1]]):.3f} Å to {np.linalg.norm(prod_aligned_positions[key_pair[0]] - prod_aligned_positions[key_pair[1]]):.3f} Å")
    
    # Create a simple quadratic reaction path model - in practice would need energy calculations
    # This is a heuristic approximation
    
    # For bond breaking/forming, the TS is often closer to the reactants than products
    # for endothermic reactions, and closer to products for exothermic reactions
    # Here we'll simply place it at 60% of the path assuming most reactions are endothermic
    t_value = 0.6
    
    # Generate the TS guess at this point along the path
    ts_positions = ref_positions + t_value * (prod_aligned_positions - ref_positions)
    
    # For the key bond(s), make some additional adjustments to better approximate TS
    for i, j in key_pairs[:3]:  # Consider up to 3 most changed bonds
        r_dist = np.linalg.norm(ref_positions[i] - ref_positions[j])
        p_dist = np.linalg.norm(prod_aligned_positions[i] - prod_aligned_positions[j])
        
        # If bond is breaking (getting longer)
        if p_dist > r_dist:
            # Transition states often have longer bonds than simple interpolation would predict
            adjustment = (p_dist - r_dist) * 0.2  # 20% extra stretching
            direction = (ts_positions[j] - ts_positions[i]) / np.linalg.norm(ts_positions[j] - ts_positions[i])
            ts_positions[i] -= direction * adjustment / 2
            ts_positions[j] += direction * adjustment / 2
        # If bond is forming (getting shorter)
        else:
            # Transition states often have longer bonds than simple interpolation would predict
            # So we don't compress as much as the linear model
            adjustment = (r_dist - p_dist) * 0.2  # 20% less compression
            direction = (ts_positions[j] - ts_positions[i]) / np.linalg.norm(ts_positions[j] - ts_positions[i])
            ts_positions[i] += direction * adjustment / 2
            ts_positions[j] -= direction * adjustment / 2
    
    # Format the resulting TS guess coordinates
    ts_lines = []
    for s, pos in zip(ref_symbols, ts_positions):
        ts_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    
    return ts_lines

def generate_ts_guess_chemical_intuition(reactant_data, product_data):
    """
    Generate a transition state guess using chemical intuition
    This approach identifies key bond changes and adjusts them based on known
    properties of transition states
    """
    # First align the product to the reactant
    aligned_results = align_product_to_reactant(
        reactant_data['coordinates'], 
        product_data['coordinates']
    )
    
    aligned_prod_lines, ref_symbols, ref_positions, prod_aligned_positions = aligned_results
    
    # Identify key bond changes between reactant and product
    key_pairs, r_atoms, p_atoms = identify_key_bond_changes(reactant_data, product_data)
    
    if not key_pairs:
        print("WARNING: No significant bond changes detected, using simple LST interpolation")
        return generate_ts_guess_lst(reactant_data, product_data)
    
    # Start with reactant structure
    ts_positions = ref_positions.copy()
    
    print("Applying chemical intuition to adjust key bond distances:")
    
    # Look at each key bond change
    for i, j in key_pairs:
        r_dist = np.linalg.norm(ref_positions[i] - ref_positions[j])
        p_dist = np.linalg.norm(prod_aligned_positions[i] - prod_aligned_positions[j])
        print(f"  Atoms {i+1}-{j+1}: {ref_symbols[i]}-{ref_symbols[j]} distance changes from {r_dist:.3f} Å to {p_dist:.3f} Å")
        
        # Determine if bond is breaking or forming
        if p_dist > r_dist:
            # Bond is breaking - extend it by ~50% of the total change
            adjustment = (p_dist - r_dist) * 0.5
            new_dist = r_dist + adjustment
            print(f"    Bond breaking: adjusting to {new_dist:.3f} Å")
        else:
            # Bond is forming - move atoms ~40% closer
            adjustment = (r_dist - p_dist) * 0.4
            new_dist = r_dist - adjustment
            print(f"    Bond forming: adjusting to {new_dist:.3f} Å")
        
        # Adjust the positions
        direction = ts_positions[j] - ts_positions[i]
        current_dist = np.linalg.norm(direction)
        direction = direction / current_dist  # normalize
        
        # Scale to the new distance
        ts_positions[i] -= direction * (current_dist - new_dist) / 2
        ts_positions[j] += direction * (current_dist - new_dist) / 2
    
    # Format the resulting TS guess coordinates
    ts_lines = []
    for s, pos in zip(ref_symbols, ts_positions):
        ts_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    
    return ts_lines

def ensure_matching_atoms(reactant_file, product_file):
    """
    Ensure that reactant and product have the same atoms in the same order.
    Returns matched data structures.
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
            print("Reactant and product have matching atoms in the same order - good for QST3")
            return (reactant_data, product_data)  # Return original data if matching
    else:
        print(f"WARNING: Reactant has {len(reactant_symbols)} atoms but product has {len(product_symbols)} atoms")
        print("QST3 requires the same atoms in the same order for both reactant and product")
    
    # If we're here, the atoms don't match - we need to match them in memory
    print("Matching atoms for QST3 compatibility...")
    
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
    
    print("Atoms matched for QST3 compatibility")
    
    return (reactant_data, matched_product_data)

def generate_transition_state_guess(reactant_data, product_data, method='lst'):
    """Generate transition state guess using specified method"""
    print(f"Generating transition state guess using {method} method...")
    
    if method == 'lst':
        return generate_ts_guess_lst(reactant_data, product_data)
    elif method == 'qst':
        return generate_ts_guess_qst(reactant_data, product_data)
    elif method == 'intuition':
        return generate_ts_guess_chemical_intuition(reactant_data, product_data)
    else:
        print(f"Unknown method '{method}', falling back to LST")
        return generate_ts_guess_lst(reactant_data, product_data)

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
        method=config['ts_guess_method']
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

def main():
    # Use default configuration
    config = DEFAULT_CONFIG.copy()
    
    # Setup command-line arguments - these will override defaults if specified
    parser = argparse.ArgumentParser(description="Generate Gaussian QST3 transition state search input files with TS guessing")
    parser.add_argument("--input_dir", help=f"Directory containing optimized .gjf files (default: {config['input_dir']})")
    parser.add_argument("--output_dir", help=f"Directory to save QST3 input files (default: {config['output_dir']})")
    parser.add_argument("--level", help=f"Level of theory for QST3 calculation (default: {config['level_of_theory']})")
    parser.add_argument("--mem", help=f"Memory allocation for Gaussian (default: {config['mem']})")
    parser.add_argument("--nproc", type=int, help=f"Number of processors to use (default: {config['nproc']})")
    parser.add_argument("--no_freq", action="store_true", help="Disable frequency calculation after optimization")
    parser.add_argument("--no_connect", action="store_true", help="Disable use of connectivity specification")
    parser.add_argument("--ts_method", choices=["lst", "qst", "intuition"], 
                        help=f"Method for generating transition state guess (default: {config['ts_guess_method']})")
    parser.add_argument("--reaction", help="Specific reaction to set up (see available reactions with --list_reactions)")
    parser.add_argument("--list_reactions", action="store_true", help="List available predefined reactions")
    parser.add_argument("--manual", action="store_true", help="Manually select reactant and product files")
    parser.add_argument("--temperature", help="Temperature for calculation (e.g., 298.15)")
    parser.add_argument("--solvent", help="Solvent for SCRF calculation (e.g., Water)")
    parser.add_argument("--redundant", action="store_true", help="Use redundant internal coordinates in optimization")
    args = parser.parse_args()
    
    # Update config with any command-line options
    if args.input_dir:
        config['input_dir'] = args.input_dir
    if args.output_dir:
        config['output_dir'] = args.output_dir
    if args.level:
        config['level_of_theory'] = args.level
    if args.mem:
        config['mem'] = args.mem
    if args.nproc:
        config['nproc'] = args.nproc
    if args.no_freq:
        config['add_freq'] = False
    if args.no_connect:
        config['geom_connect'] = False
    if args.ts_method:
        config['ts_guess_method'] = args.ts_method
    if args.temperature:
        config['temperature'] = args.temperature
    if args.solvent:
        config['solvent'] = args.solvent
    if args.redundant:
        config['redundant'] = True
    
    # Create output directory
    ensure_directory(config['output_dir'])
    
    # Handle listing reactions
    reaction_pairs = setup_reaction_pairs()
    if args.list_reactions:
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
            print(f"Found input directory at {input_dir}")
            break
    
    if not input_dir:
        print("ERROR: Cannot find input directory for optimized geometries!")
        print("Searched the following paths:")
        for path in potential_paths:
            print(f"  - {path}")
        return
        
    config['input_dir'] = str(input_dir)
    
    # Manual selection mode
    if args.manual:
        gjf_files = list_available_molecules(config['input_dir'])
        
        if not gjf_files:
            print(f"No .gjf files found in {config['input_dir']}")
            return
            
        print("\nSelect reactant file number:")
        reactant_idx = int(input("> ")) - 1
        reactant_file = gjf_files[reactant_idx]
        
        print("\nSelect product file number:")
        product_idx = int(input("> ")) - 1
        product_file = gjf_files[product_idx]
        
        print("\nSelect TS guess method:")
        print("1. Linear Synchronous Transit (LST)")
        print("2. Quadratic Synchronous Transit (QST)")
        print("3. Chemical Intuition")
        ts_method_choice = int(input("> "))
        
        ts_methods = {1: "lst", 2: "qst", 3: "intuition"}
        config['ts_guess_method'] = ts_methods.get(ts_method_choice, "lst")
        
        create_qst3_input(
            reactant_file=str(reactant_file),
            product_file=str(product_file),
            output_dir=config['output_dir'],
            config=config
        )
        return
    
    # Process specific reaction or all predefined reactions
    reactions_to_process = []
    if args.reaction:
        if args.reaction not in reaction_pairs:
            print(f"Error: Reaction '{args.reaction}' not found in predefined reactions")
            print("Use --list_reactions to see available reactions")
            return
        reactions_to_process = [args.reaction]
    else:
        # Process all predefined reactions
        reactions_to_process = list(reaction_pairs.keys())
    
    for reaction_name in reactions_to_process:
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
        
        # Try each TS guessing method
        for method in ["lst", "qst", "intuition"]:
            print(f"\nTrying {method} method for TS guessing...")
            config['ts_guess_method'] = method
            create_qst3_input(
                reactant_file=reactant_file,
                product_file=product_file,
                output_dir=config['output_dir'],
                config=config
            )

    print(f"\nQST3 input generation complete. Files saved to {config['output_dir']}")
    print("Run these Gaussian calculations to locate transition states.")
    print("\nTips for verifying transition states:")
    print("1. Check that the optimization converged to a stationary point")
    print("2. Verify that the structure has exactly ONE imaginary frequency")
    print("3. Examine the imaginary mode to confirm it corresponds to the reaction coordinate")
    print("4. For confirmation, perform IRC calculations from the transition state")

if __name__ == "__main__":
    main() 