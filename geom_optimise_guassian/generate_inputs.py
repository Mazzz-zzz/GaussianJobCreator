#!/usr/bin/env python

"""
Script to generate Gaussian input files for reaction optimization study using SMILES
"""

import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

def get_connectivity_matrix(mol):
    """Generate connectivity matrix from RDKit molecule"""
    connectivity = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx() + 1  # 1-based indexing for Gaussian
        bonds = []
        for bond in atom.GetBonds():
            other_idx = bond.GetOtherAtomIdx(atom.GetIdx()) + 1
            if other_idx > idx:  # Only include forward connections
                bonds.append(f"{other_idx} 1.0")
        if bonds:
            connectivity.append(f"{idx} " + " ".join(bonds))
        else:
            connectivity.append(f"{idx}")
    return connectivity

def create_gaussian_input(mol, name, charge, multiplicity, output_dir="gaussian_projects"):
    """Create a Gaussian input file with parameters"""
    output_path = Path(output_dir) / f"{name}.gjf"
    
    with open(output_path, 'w') as f:
        # Header
        f.write("%mem=3GB\n")
        f.write("%nprocshared=4\n")
        f.write("# opt=(maxcyc=999,noeigen) freq m062x/def2tzvp geom=connectivity int=ultrafine scf=(tight,xqc)\n\n")
        
        # Title
        f.write(f"{name} optimization\n\n")
        
        # Charge and multiplicity
        f.write(f"{charge} {multiplicity}\n")
        
        # Coordinates
        conf = mol.GetConformer()
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            f.write(f"{atom.GetSymbol():2s}    {pos.x:10.6f}    {pos.y:10.6f}    {pos.z:10.6f}\n")
        
        # Connectivity matrix
        f.write("\n")
        connectivity = get_connectivity_matrix(mol)
        for line in connectivity:
            f.write(line + "\n")
        
        # End with blank line
        f.write("\n")

# Dictionary of molecules with their specifications
molecules = {
    # Reactants
    'Na': {
        'smiles': '[Na]',
        'charge': 0,
        'multiplicity': 2
    },
    'CH3F': {
        'smiles': 'CF',
        'charge': 0,
        'multiplicity': 1
    },
    'CH3Cl': {
        'smiles': 'CCl',
        'charge': 0,
        'multiplicity': 1
    },
    'CH3Br': {
        'smiles': 'CBr',
        'charge': 0,
        'multiplicity': 1
    },
    'HCl': {
        'smiles': 'Cl',  # Will be converted to HCl
        'charge': 0,
        'multiplicity': 1
    },
    'HBr': {
        'smiles': 'Br',  # Will be converted to HBr
        'charge': 0,
        'multiplicity': 1
    },
    # Products
    'CH3': {
        'smiles': '[CH3]',
        'charge': 0,
        'multiplicity': 2
    },
    'H': {
        'smiles': '[H]',
        'charge': 0,
        'multiplicity': 2
    },

    # PFMS and all its reaction pathways
    'PFMS': {
        'smiles': 'C(F)(F)(F)S(=O)(=O)O',  # Reactant
        'charge': 0,
        'multiplicity': 1
    },
    
    # TS1M Products
    'TS1M_Product1': {
        'smiles': 'FC1(F)OS(=O)(=O)1',
        'charge': 0,
        'multiplicity': 1
    },
    'HF': {
        'smiles': 'F',  # Will be converted to HF
        'charge': 0,
        'multiplicity': 1
    },
    
    # TS2M Products
    'HCF3': {
        'smiles': 'C(F)(F)F',
        'charge': 0,
        'multiplicity': 1
    },
    'SO3': {
        'smiles': 'S(=O)(=O)=O',
        'charge': 0,
        'multiplicity': 1
    },
    
    # ISOPFMS Product
    'ISOPFMS': {
        'smiles': 'C(F)(F)(F)OS(O)(=O)O',
        'charge': 0,
        'multiplicity': 1
    },
    
    # TS3M Products
    'TS3M_Product1': {
        'smiles': '[O]S(O[H])=O',  # Radical
        'charge': 0,
        'multiplicity': 2
    },
    'CF3_Radical': {
        'smiles': '[*](F)(F)F',  # Radical
        'charge': 0,
        'multiplicity': 2
    },
    
    # TS4M Products
    'TS4M_Product1': {
        'smiles': 'FC(F)(F)OS(=O)(=O)1',
        'charge': 0,
        'multiplicity': 1
    },
    'F_TS4M': {  # Adding suffix to distinguish from other F products
        'smiles': '[F]',
        'charge': 0,
        'multiplicity': 2
    },
    
    # TS5M Products (same products as TS4M but different pathway)
    'TS5M_Product1': {
        'smiles': 'FC(F)(F)OS(=O)(=O)1',
        'charge': 0,
        'multiplicity': 1
    },
    'F_TS5M': {
        'smiles': '[F]',
        'charge': 0,
        'multiplicity': 2
    },
    
    # TS6M Products
    'F_TS6M': {
        'smiles': '[F]',
        'charge': 0,
        'multiplicity': 2
    },
    'CF2O': {
        'smiles': 'FC(=O)F',
        'charge': 0,
        'multiplicity': 1
    },
    'SO2': {
        'smiles': 'O=S=O',
        'charge': 0,
        'multiplicity': 1
    }
}

def main():
    # Create output directory if it doesn't exist
    os.makedirs("gaussian_projects", exist_ok=True)
    
    # Generate input files for each molecule
    for name, specs in molecules.items():
        try:
            # Create RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(specs['smiles'])
            if mol is None:
                raise ValueError(f"Failed to create molecule from SMILES: {specs['smiles']}")
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # Optimize with MMFF94
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Create Gaussian input file
            create_gaussian_input(
                mol=mol,
                name=name,
                charge=specs['charge'],
                multiplicity=specs['multiplicity']
            )
            print(f"Created input file for {name}")
            
        except Exception as e:
            print(f"Error processing {name}: {str(e)}")

if __name__ == "__main__":
    main() 