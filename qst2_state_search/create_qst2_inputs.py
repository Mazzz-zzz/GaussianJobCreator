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
    'input_dir': "../geom_optimise_gaussian/gaussian_projects", 
    'output_dir': "qst2_jobs",
    
    # Gaussian calculation parameters
    'level_of_theory': "b3lyp/6-31G(d)",
    'mem': "8GB",
    'nproc': 8,
    'add_freq': True,    # Add frequency calculation
    'geom_connect': False, # Use connectivity information
    'use_zmatrix': True,   # Use Z-matrix format (True) or Cartesian (False)
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

def parse_gaussian_input(com_file):
    """
    Parse a Gaussian .com input file to extract geometry and other parameters.
    Now supports both Cartesian coordinates and Z-matrix format.
    """
    with open(com_file, 'r') as f:
        lines = f.readlines()
    
    # Find charge and multiplicity line
    charge_mult_line = None
    coord_start = None
    
    # Skip header lines (memory, processors, etc.)
    header_end = 0
    for i, line in enumerate(lines):
        if line.strip().startswith('#'):
            # Skip until we find a blank line after the route section
            j = i + 1
            while j < len(lines) and lines[j].strip():
                j += 1
            header_end = j + 1  # Skip the blank line too
            break
    
    # Skip title and blank line
    title_line = header_end
    if title_line < len(lines) and lines[title_line].strip():
        title_line += 1  # Skip title
    if title_line < len(lines) and not lines[title_line].strip():
        title_line += 1  # Skip blank line after title
    
    # Find charge and multiplicity
    for i in range(title_line, len(lines)):
        line = lines[i].strip()
        # Skip empty lines
        if not line:
            continue
        
        # Charge and multiplicity line (either comma or space separated)
        if re.match(r'^\s*[+-]?\d+\s*[,\s]\s*\d+\s*$', line):
            charge_mult_line = i
            coord_start = i + 1
            break
    
    if charge_mult_line is None:
        raise ValueError(f"Could not find charge/multiplicity line in {com_file}")
    
    # Parse charge and multiplicity (handle both space and comma separated)
    charge_mult = lines[charge_mult_line].strip().replace(',', ' ').split()
    charge, multiplicity = map(int, charge_mult)
    
    # Extract coordinates - start from after charge/multiplicity line
    coordinates = []
    i = coord_start
    is_zmatrix = False
    
    # Try to detect if it's a Z-matrix
    # In a Z-matrix, after the first atom, lines typically have:
    # atom1
    # atom2 bond_to_atom index
    # atom3 bond_to_atom index angle_with_atom index
    # atom4+ bond_to_atom index angle_with_atom index dihedral_with_atom index
    
    # Let's check for Z-matrix pattern in the first few lines
    zmatrix_pattern = True
    first_line_tokens = None
    
    if i < len(lines) and lines[i].strip():
        first_line_tokens = lines[i].strip().split()
        # First line in Z-matrix should normally just have an atom symbol
        if len(first_line_tokens) > 1 and len(first_line_tokens) != 4:  # Not Cartesian (4 tokens: atom x y z)
            zmatrix_pattern = True
        elif len(first_line_tokens) == 4:
            # Could be Cartesian if the tokens after the first look like numbers
            try:
                float(first_line_tokens[1])
                float(first_line_tokens[2])
                float(first_line_tokens[3])
                zmatrix_pattern = False  # It's Cartesian if they're all numbers
            except ValueError:
                zmatrix_pattern = True   # Not Cartesian if not numbers
    
    # For filenames with _zmat, assume Z-matrix format regardless of content pattern
    if "_zmat" in com_file.lower():
        is_zmatrix = True
    else:
        is_zmatrix = zmatrix_pattern
    
    # Read coordinates until blank line
    while i < len(lines) and lines[i].strip():
        coordinates.append(lines[i])
        i += 1
    
    # Skip any blank lines
    while i < len(lines) and not lines[i].strip():
        i += 1
    
    # Look for connectivity information
    connectivity = []
    while i < len(lines):
        line = lines[i].strip()
        if not line:  # Stop at blank line
            break
        if line[0].isdigit():  # Only include lines that start with a number
            connectivity.append(line)
        i += 1
    
    if not coordinates:
        raise ValueError(f"No valid coordinates found in {com_file}")
    
    return {
        'charge': charge,
        'multiplicity': multiplicity,
        'coordinates': coordinates,
        'connectivity': connectivity,
        'is_zmatrix': is_zmatrix
    }

def align_coordinates(coordinate_lines):
    """
    Align Cartesian coordinates around their center of mass.
    Only used for Cartesian coordinates, not Z-matrix.
    """
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
    """
    Align product Cartesian coordinates to match reactant orientation.
    Only used for Cartesian coordinates, not Z-matrix.
    """
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
        Path to the reactant .com file
    product_file : str
        Path to the product .com file
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
    
    # Check if both files use the same format (Z-matrix or Cartesian)
    if reactant_data['is_zmatrix'] != product_data['is_zmatrix']:
        print(f"WARNING: Format mismatch between reactant and product:")
        print(f"  Reactant: {'Z-matrix' if reactant_data['is_zmatrix'] else 'Cartesian'}")
        print(f"  Product:  {'Z-matrix' if product_data['is_zmatrix'] else 'Cartesian'}")
        print(f"QST2 works best when both structures use the same format.")
    
    # If using Z-matrix, we don't need to match atoms - just return the data as is
    if reactant_data['is_zmatrix'] and product_data['is_zmatrix']:
        print("Both reactant and product use Z-matrix format - keeping original structures")
        return (reactant_data, product_data)
    
    # For Cartesian coordinates, we need to ensure atom matching
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
    
    # Update connectivity if present
    if matched_product_data['connectivity']:
        # Keep only connectivity for atoms that exist in both structures
        new_connectivity = []
        for line in matched_product_data['connectivity']:
            parts = line.strip().split()
            if not parts or not parts[0].isdigit():
                continue
            atom_idx = int(parts[0])
            if atom_idx <= min_atoms:
                new_connectivity.append(line)
        matched_product_data['connectivity'] = new_connectivity
    
    print("Atoms matched for QST2 compatibility")
    
    return (reactant_data, matched_product_data)

def combine_molecules(molecule_files, is_zmatrix=False, spacing=4.0):
    """
    Combine multiple molecule files into a single structure with proper spacing.
    
    Parameters:
    -----------
    molecule_files : list
        List of paths to molecule files
    is_zmatrix : bool
        Whether the files are in Z-matrix format
    spacing : float
        Spacing between molecules in Angstroms (only used for Cartesian coordinates)
        
    Returns:
    --------
    dict:
        Combined molecule data with charge, multiplicity, coordinates, connectivity
    """
    if not molecule_files:
        raise ValueError("No molecule files provided to combine")
    
    # If only one molecule, just parse and return it
    if len(molecule_files) == 1:
        return parse_gaussian_input(molecule_files[0])
    
    # For Z-matrix files, we need special handling
    if is_zmatrix:
        print("WARNING: Combining multiple Z-matrix files requires special handling")
        print("Creating a combined structure by concatenating Z-matrices with variable renaming")
        
        combined_data = {
            'charge': 0,
            'multiplicity': 1,
            'coordinates': [],
            'connectivity': [],
            'is_zmatrix': True
        }
        
        var_counter = 1
        atom_counter = 0
        variables = {}  # To store all variables and their values
        
        # First parse all files to get total charge and max multiplicity
        for file_path in molecule_files:
            mol_data = parse_gaussian_input(file_path)
            combined_data['charge'] += mol_data['charge']
            combined_data['multiplicity'] = max(combined_data['multiplicity'], mol_data['multiplicity'])
        
        # Now process each molecule
        for file_idx, file_path in enumerate(molecule_files):
            mol_data = parse_gaussian_input(file_path)
            
            # For first molecule, keep the coordinates as is
            if file_idx == 0:
                for line in mol_data['coordinates']:
                    combined_data['coordinates'].append(line)
                atom_counter += len(mol_data['coordinates'])
                continue
            
            # For subsequent molecules, we need to carefully add them
            # with references to atoms in the first molecule
            first_atoms_added = False
            
            for i, line in enumerate(mol_data['coordinates']):
                line = line.strip()
                if not line:
                    continue
                
                tokens = line.split()
                
                # First atom of each additional molecule - reference to first atom of first molecule
                if i == 0:
                    # Format: atomSymbol atom1 distance
                    combined_data['coordinates'].append(f"{tokens[0]} 1 R{var_counter}\n")
                    variables[f"R{var_counter}"] = spacing * (file_idx + 1)  # Arbitrary distance
                    var_counter += 1
                    first_atoms_added = True
                # Second atom - reference first atom of this molecule and first atom of first molecule
                elif i == 1 and first_atoms_added:
                    # Format: atomSymbol atom1 distance atom2 angle
                    combined_data['coordinates'].append(f"{tokens[0]} {atom_counter} R{var_counter} 1 A{var_counter}\n")
                    variables[f"R{var_counter}"] = 1.5  # Arbitrary reasonable bond length
                    variables[f"A{var_counter}"] = 90.0  # Arbitrary angle in degrees
                    var_counter += 1
                # Third atom - reference first two atoms of this molecule
                elif i == 2 and first_atoms_added:
                    # Format: atomSymbol atom1 distance atom2 angle atom3 dihedral
                    combined_data['coordinates'].append(
                        f"{tokens[0]} {atom_counter} R{var_counter} {atom_counter+1} A{var_counter} 1 D{var_counter}\n"
                    )
                    variables[f"R{var_counter}"] = 1.5  # Arbitrary bond length
                    variables[f"A{var_counter}"] = 90.0  # Arbitrary angle
                    variables[f"D{var_counter}"] = 180.0  # Arbitrary dihedral
                    var_counter += 1
                # Remaining atoms - reference to the first three atoms of this molecule
                else:
                    # Format: atomSymbol atom1 distance atom2 angle atom3 dihedral
                    combined_data['coordinates'].append(
                        f"{tokens[0]} {atom_counter} R{var_counter} {atom_counter+1} A{var_counter} {atom_counter+2} D{var_counter}\n"
                    )
                    variables[f"R{var_counter}"] = 1.5  # Arbitrary bond length
                    variables[f"A{var_counter}"] = 90.0  # Arbitrary angle
                    variables[f"D{var_counter}"] = 180.0  # Arbitrary dihedral
                    var_counter += 1
                
                atom_counter += 1
        
        # Add variables section at the end
        if variables:
            combined_data['coordinates'].append("\n")  # Blank line before variables
            for var, value in variables.items():
                combined_data['coordinates'].append(f"{var}={value:.6f}\n")
        
        # Don't need connectivity for Z-matrix
        return combined_data
    
    # For Cartesian coordinates, we'll combine them with spatial separation
    combined_data = {
        'charge': 0,
        'multiplicity': 1,
        'coordinates': [],
        'connectivity': [],
        'is_zmatrix': False
    }
    
    # Track current position for molecule placement
    current_position = np.array([0.0, 0.0, 0.0])
    atom_counter = 0
    
    for file_path in molecule_files:
        mol_data = parse_gaussian_input(file_path)
        
        # Add charges
        combined_data['charge'] += mol_data['charge']
        
        # For multiplicity, use the maximum value
        combined_data['multiplicity'] = max(combined_data['multiplicity'], mol_data['multiplicity'])
        
        # Extract Cartesian coordinates
        atoms = []
        for line in mol_data['coordinates']:
            tokens = line.split()
            if len(tokens) >= 4:
                try:
                    symbol = tokens[0]
                    position = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])
                    atoms.append((symbol, position))
                except ValueError:
                    print(f"Warning: Skipping non-Cartesian line: {line}")
        
        if not atoms:
            print(f"Warning: No Cartesian coordinates found in {file_path}")
            continue
        
        # Calculate center of mass
        positions = np.array([pos for _, pos in atoms])
        center_of_mass = np.mean(positions, axis=0)
        
        # Translate molecule to current position
        for symbol, position in atoms:
            # Adjust position: center the molecule and then offset by current_position
            new_position = position - center_of_mass + current_position
            combined_data['coordinates'].append(
                f"{symbol:2s}    {new_position[0]: .4f}   {new_position[1]: .4f}   {new_position[2]: .4f}\n"
            )
        
        # Update connectivity if present
        if mol_data['connectivity']:
            for line in mol_data['connectivity']:
                parts = line.strip().split()
                if parts and parts[0].isdigit():
                    # Adjust atom indices
                    new_line = f"{int(parts[0]) + atom_counter}"
                    for i in range(1, len(parts)):
                        if parts[i].isdigit():
                            new_line += f" {int(parts[i]) + atom_counter}"
                        else:
                            new_line += f" {parts[i]}"
                    combined_data['connectivity'].append(new_line)
        
        # Update atom counter for connectivity adjustment
        atom_counter += len(atoms)
        
        # Move current position for next molecule
        # We'll place molecules along the x-axis with specified spacing
        current_position[0] += spacing + np.max(positions[:, 0]) - np.min(positions[:, 0])
    
    return combined_data

def create_qst2_input(reactant_file, product_files, output_dir, config=None):
    """
    Create a Gaussian QST2 input file using geometries from reactant and product files
    
    Parameters:
    -----------
    reactant_file : str
        Path to the reactant .com file
    product_files : list
        List of paths to product .com files
    output_dir : str
        Directory to save the QST2 input file
    config : dict
        Configuration dictionary with calculation parameters
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    # Parse reactant data
    reactant_data = parse_gaussian_input(reactant_file)
    
    # For multiple product files, combine them with spatial separation
    if isinstance(product_files, list) and len(product_files) > 1:
        print(f"Combining {len(product_files)} product files: {[os.path.basename(f) for f in product_files]}")
        product_data = combine_molecules(product_files, is_zmatrix=reactant_data['is_zmatrix'])
    else:
        # For a single product file, parse it normally
        product_file = product_files[0] if isinstance(product_files, list) else product_files
        product_data = parse_gaussian_input(product_file)
    
    # Get base names for title
    reactant_name = os.path.basename(reactant_file).replace('_zmat.com', '')
    if isinstance(product_files, list):
        product_names = [os.path.basename(f).replace('_zmat.com', '') for f in product_files]
        product_name = "_and_".join(product_names)
    else:
        product_name = os.path.basename(product_files).replace('_zmat.com', '')
    
    # Create output filename (using .gjf extension)
    output_file = os.path.join(output_dir, f"TS_{reactant_name}_to_{product_name}.gjf")
    
    # Check if charge and multiplicity match
    if reactant_data['charge'] != product_data['charge']:
        print(f"WARNING: Charge mismatch between reactant and product:")
        print(f"  Reactant: {reactant_data['charge']}")
        print(f"  Product:  {product_data['charge']}")
        print(f"QST2 requires same charge for both structures.")
    
    # Check if products have the same number of atoms as reactants
    reactant_atoms = len(reactant_data['coordinates'])
    product_atoms = len(product_data['coordinates'])
    
    if reactant_atoms != product_atoms:
        print(f"WARNING: Atom count mismatch between reactant ({reactant_atoms}) and product ({product_atoms})")
        print("QST2 requires same number of atoms in both structures.")
        
        if reactant_atoms > product_atoms:
            print(f"Need to add {reactant_atoms - product_atoms} atoms to product structure.")
            # If using Z-matrix, we can't easily add atoms - warn and proceed with caution
            if product_data['is_zmatrix']:
                print("WARNING: Cannot automatically add atoms to Z-matrix format")
                print("The QST2 calculation may fail with unequal atom counts")
            else:
                # For Cartesian, duplicate the last atom at slightly different positions
                last_atom_line = product_data['coordinates'][-1]
                tokens = last_atom_line.split()
                if len(tokens) >= 4:
                    symbol = tokens[0]
                    coords = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])
                    
                    for i in range(reactant_atoms - product_atoms):
                        # Add a small offset to avoid exact overlaps
                        offset = np.array([0.1 * (i + 1), 0.1 * (i + 1), 0.1 * (i + 1)])
                        new_coords = coords + offset
                        new_line = f"{symbol:2s}    {new_coords[0]: .4f}   {new_coords[1]: .4f}   {new_coords[2]: .4f}\n"
                        product_data['coordinates'].append(new_line)
                    
                    print(f"Added {reactant_atoms - product_atoms} {symbol} atoms to match reactant count")
        else:
            print(f"Need to remove {product_atoms - reactant_atoms} atoms from product structure.")
            # Truncate product data to match reactant count
            product_data['coordinates'] = product_data['coordinates'][:reactant_atoms]
            if product_data['connectivity']:
                # Keep only connectivity for atoms that exist
                new_connectivity = []
                for line in product_data['connectivity']:
                    parts = line.strip().split()
                    if not parts or not parts[0].isdigit():
                        continue
                    atom_idx = int(parts[0])
                    if atom_idx <= reactant_atoms:
                        new_connectivity.append(line)
                product_data['connectivity'] = new_connectivity
            
            print(f"Removed {product_atoms - reactant_atoms} atoms from product to match reactant count")
    
    with open(output_file, 'w') as f:
        # Write header and route
        f.write(f"%mem={config['mem']}\n")
        f.write(f"%nprocshared={config['nproc']}\n")
        f.write(f"%chk=TS_{reactant_name}_to_{product_name}.chk\n")
        
        # Build route line
        route = f"# opt=(QST2,maxcyc=999,noeigen) freq m062x/def2tzvp"
        
        # Add geometry specification based on config and file format
        if config['geom_connect'] and not reactant_data['is_zmatrix']:
            route += " geom=connectivity"
        elif reactant_data['is_zmatrix'] or config['use_zmatrix']:
            route += " geom=z-matrix"
            
        route += " int=ultrafine scf=(tight,xqc)"
        
        f.write(f"{route}\n\n")
        
        # First title section for reactants
        f.write(f"Reactants: {reactant_name}\n\n")
        
        # Write charge/multiplicity and coordinates for reactant
        f.write(f"{int(reactant_data['charge'])} {int(reactant_data['multiplicity'])}\n")  # Space-separated for .gjf
        
        # For Z-matrix, use original coordinates; for Cartesian, align around center of mass
        if reactant_data['is_zmatrix']:
            for line in reactant_data['coordinates']:
                f.write(line if line.endswith('\n') else line + '\n')
        else:
            reactant_aligned = align_coordinates(reactant_data['coordinates'])
            for line in reactant_aligned:
                f.write(line)
        
        # Add connectivity for reactant only if using Cartesian with connectivity
        if config['geom_connect'] and not reactant_data['is_zmatrix'] and reactant_data['connectivity']:
            f.write("\n")
            # Write connectivity with a blank line at the end
            for line in reactant_data['connectivity']:
                f.write(line + "\n")
            f.write("\n")
        else:
            if reactant_data['is_zmatrix']:
                print(f"Using Z-matrix format for {reactant_file}")
            elif config['geom_connect'] and not reactant_data['connectivity']:
                print(f"Warning: No connectivity information found in {reactant_file}")
            f.write("\n")
        
        # Second title section for products
        f.write(f"Products: {product_name}\n\n")
        
        # Write charge/multiplicity for product
        f.write(f"{int(product_data['charge'])} {int(product_data['multiplicity'])}\n")  # Space-separated for .gjf
        
        # For Z-matrix, use original coordinates; for Cartesian, align to reactant
        if product_data['is_zmatrix']:
            for line in product_data['coordinates']:
                f.write(line if line.endswith('\n') else line + '\n')
        else:
            # For combined products, we've already built proper coordinates
            if isinstance(product_files, list) and len(product_files) > 1:
                for line in product_data['coordinates']:
                    f.write(line)
            else:
                # For single product, align to reactant
                product_aligned = align_product_to_reactant(reactant_data['coordinates'], product_data['coordinates'])
                for line in product_aligned:
                    f.write(line)
        
        # Add connectivity for product only if using Cartesian with connectivity
        if config['geom_connect'] and not product_data['is_zmatrix'] and product_data['connectivity']:
            f.write("\n")
            # Write connectivity with a blank line at the end
            for line in product_data['connectivity']:
                f.write(line + "\n")
            f.write("\n")
        else:
            if product_data['is_zmatrix']:
                print(f"Using Z-matrix format for product files")
            elif config['geom_connect'] and not product_data['connectivity']:
                print(f"Warning: No connectivity information found in product files")
            f.write("\n")
        
        # Ensure a blank line at the end of the file
        f.write("\n")
    
    print(f"Created QST2 input file: {output_file}")
    print(f"  - Reactant: {reactant_name} ({reactant_atoms} atoms)")
    
    if isinstance(product_files, list):
        print(f"  - Products: {len(product_files)} molecules ({product_atoms} total atoms)")
        for prod in product_files:
            print(f"    * {os.path.basename(prod)}")
    else:
        print(f"  - Product: {product_name} ({product_atoms} atoms)")
    
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
            "reactants": ["Na_zmat", "CH3F_zmat"],
            "products": ["CH3Cl_zmat", "Na_zmat"]
        },
        "SN2_halogen_exchange_F_Br": {
            "reactants": ["Na_zmat", "CH3F_zmat"],
            "products": ["CH3Br_zmat", "Na_zmat"]
        },
        "TS1M_PFMS_Pathway": {
            "reactants": ["PFMS_zmat"],
            "products": ["TS1M_Product1_zmat", "HF_zmat"]
        },
        "TS2M_PFMS_Pathway": {
            "reactants": ["PFMS_zmat"],
            "products": ["HCF3_zmat", "SO3_zmat"]
        },
        "TS3M_PFMS_Pathway": {
            "reactants": ["ISOPFMS_zmat"],
            "products": ["TS3M_Product1_zmat", "CF3_Radical_zmat"]
        },
        "TS4M_PFMS_Pathway": {
            "reactants": ["PFMS_zmat"],
            "products": ["TS4M_Product1_zmat", "F_TS4M_zmat"]
        },
        "TS5M_PFMS_Pathway": {
            "reactants": ["PFMS_zmat"],
            "products": ["TS5M_Product1_zmat", "F_TS5M_zmat"]
        },
        "TS6M_PFMS_Pathway": {
            "reactants": ["PFMS_zmat"],
            "products": ["F_TS6M_zmat", "CF2O_zmat", "SO2_zmat"]
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
        
        # Get the first reactant file
        if len(reaction['reactants']) > 0:
            reactant_file = os.path.join(config['input_dir'], f"{reaction['reactants'][0]}.com")
            if not os.path.exists(reactant_file):
                print(f"Warning: Reactant file {reactant_file} not found")
                continue
        else:
            print(f"Warning: No reactants defined for {reaction_name}")
            continue
        
        # Get all product files
        product_files = []
        for product in reaction['products']:
            product_file = os.path.join(config['input_dir'], f"{product}.com")
            if os.path.exists(product_file):
                product_files.append(product_file)
            else:
                print(f"Warning: Product file {product_file} not found")
        
        if not product_files:
            print(f"Warning: No valid product files found for {reaction_name}")
            continue
        
        # Create QST2 input file with all products
        create_qst2_input(
            reactant_file=reactant_file,
            product_files=product_files,
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
