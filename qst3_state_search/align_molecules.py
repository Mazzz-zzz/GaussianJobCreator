#!/usr/bin/env python

"""
Functions for aligning molecules and manipulating coordinates for QST3 calculations.
"""

import numpy as np
from ase import Atoms

def parse_gaussian_input(gjf_file):
    """Parse a Gaussian .gjf input file to extract geometry and other parameters"""
    import re
    
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
    """
    Align product coordinates to reactant using a more robust approach.
    Handles cases with different numbers of atoms and atom type mismatches.
    
    Args:
        reactant_lines: List of coordinate lines from reactant structure
        product_lines: List of coordinate lines from product structure
        
    Returns:
        tuple: (aligned_product_lines, ref_symbols, ref_centered, prod_aligned)
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
    
    # Handle different numbers of atoms
    if len(ref_symbols) != len(prod_symbols):
        print(f"WARNING: Reactant has {len(ref_symbols)} atoms but product has {len(prod_symbols)} atoms")
        print("QST3 requires the same atoms in the same order for both reactant and product")
        print("Creating a matched product structure...")
        
        # Create a new product structure with the same atoms as the reactant
        new_prod_symbols = ref_symbols.copy()
        new_prod_positions = np.zeros_like(ref_positions)
        
        # Copy over the positions we have from the product
        min_atoms = min(len(ref_symbols), len(prod_symbols))
        new_prod_positions[:min_atoms] = prod_positions[:min_atoms]
        
        # For any remaining atoms, use positions from reactant but displaced
        if min_atoms < len(ref_symbols):
            new_prod_positions[min_atoms:] = ref_positions[min_atoms:] + np.array([3.0, 0.0, 0.0])
        
        prod_symbols = new_prod_symbols
        prod_positions = new_prod_positions
    
    # Check and fix atom type mismatches
    for i, (ref, prod) in enumerate(zip(ref_symbols, prod_symbols)):
        if ref != prod:
            print(f"WARNING: Atom mismatch at position {i+1}: reactant has {ref} but product has {prod}")
            print(f"Setting product atom type to match reactant ({ref})")
            prod_symbols[i] = ref
    
    # Center both structures at their geometric centers
    ref_com = np.mean(ref_positions, axis=0)
    prod_com = np.mean(prod_positions, axis=0)
    ref_centered = ref_positions - ref_com
    prod_centered = prod_positions - prod_com
    
    # Calculate the maximum extent of the reactant molecule
    max_extent = np.max(np.abs(ref_centered))
    
    # Place product at a distance proportional to the molecule size, but at least 3Å
    SEPARATION = max(3.0, max_extent * 1.5)  # Angstroms
    
    # Align product along principal axes
    prod_aligned = prod_centered.copy()
    
    # Move product in x direction by separation distance
    prod_aligned[:, 0] += SEPARATION
    
    # Check for invalid structures (atoms too close together)
    MIN_DISTANCE = 0.7  # Minimum allowed distance between atoms in Angstroms
    for i in range(len(prod_aligned)):
        for j in range(i + 1, len(prod_aligned)):
            dist = np.linalg.norm(prod_aligned[i] - prod_aligned[j])
            if dist < MIN_DISTANCE:
                print(f"WARNING: Very close atoms detected in product ({dist:.2f}Å), adjusting positions...")
                # Move atoms apart slightly in random direction
                direction = np.random.rand(3) - 0.5
                direction = direction / np.linalg.norm(direction)
                prod_aligned[j] = prod_aligned[i] + direction * MIN_DISTANCE * 1.2
    
    # Format aligned product coordinates
    aligned_lines = []
    for s, pos in zip(prod_symbols, prod_aligned):
        aligned_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    
    return aligned_lines, ref_symbols, ref_centered, prod_aligned

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