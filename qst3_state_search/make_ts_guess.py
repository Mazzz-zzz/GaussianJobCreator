#!/usr/bin/env python

"""
Functions for generating transition state guesses for QST3 calculations
using Linear Synchronous Transit (LST) method.
"""

import numpy as np
from align_molecules import align_product_to_reactant

def check_valid_structure(positions, symbols):
    """
    Check if a molecular structure has valid coordinates
    Returns True if the structure seems valid, False otherwise
    """
    if len(positions) < 2:
        return False
    
    # Check for duplicate coordinates
    unique_positions = set()
    for pos in positions:
        pos_tuple = tuple(np.round(pos, 3))  # Round to 3 decimal places
        unique_positions.add(pos_tuple)
    
    # If more than half the coordinates are duplicates, structure is suspicious
    if len(unique_positions) < len(positions) / 2:
        print("WARNING: Structure has many duplicate atom positions - may be invalid")
        return False
    
    # Check for unreasonable bond lengths (too short < 0.7Å or too long > 5Å)
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist < 0.7 and symbols[i] != 'H' and symbols[j] != 'H':
                print(f"WARNING: Very short bond ({dist:.2f} Å) between {symbols[i]}{i+1}-{symbols[j]}{j+1}")
                return False
    
    return True

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
    
    # Check if product structure is valid
    product_valid = check_valid_structure(prod_aligned_positions, ref_symbols)
    if not product_valid:
        print("WARNING: Product structure appears invalid. Using modified reactant structure for TS guess.")
        # Instead of interpolation, just slightly modify the reactant structure
        ts_positions = ref_positions.copy()
        # Apply a small random perturbation to break symmetry
        np.random.seed(42)  # For reproducibility
        perturbation = np.random.uniform(-0.2, 0.2, ts_positions.shape)
        ts_positions += perturbation
    else:
        # Simple linear interpolation at halfway point (t=0.5)
        t = 0.5  # halfway point
        ts_positions = ref_positions + t * (prod_aligned_positions - ref_positions)
    
    # Format the resulting TS guess coordinates
    ts_lines = []
    for s, pos in zip(ref_symbols, ts_positions):
        ts_lines.append(f"{s:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
    
    return ts_lines

def generate_transition_state_guess(reactant_data, product_data, method='lst'):
    """Generate transition state guess using LST method"""
    print("Generating transition state guess using Linear Synchronous Transit (LST) method...")
    
    # Check if reactant structure is valid
    reactant_atoms = []
    reactant_positions = []
    for line in reactant_data['coordinates']:
        tokens = line.split()
        if len(tokens) >= 4:
            reactant_atoms.append(tokens[0])
            reactant_positions.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
    
    reactant_valid = check_valid_structure(np.array(reactant_positions), reactant_atoms)
    if not reactant_valid:
        print("ERROR: Reactant structure appears invalid. Cannot generate valid TS guess.")
        # Return a minimal valid structure to avoid crashes
        ts_lines = []
        for atom, pos in zip(reactant_atoms, reactant_positions):
            ts_lines.append(f"{atom:2s}    {pos[0]: .4f}   {pos[1]: .4f}   {pos[2]: .4f}\n")
        return ts_lines
    
    return generate_ts_guess_lst(reactant_data, product_data) 