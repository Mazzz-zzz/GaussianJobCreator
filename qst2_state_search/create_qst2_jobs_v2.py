#!/usr/bin/env python3
"""
Generate QST2 Gaussian input files from Z-matrix formatted .com files.

This script expects the input folder to contain reactant and product files
specified in the REACTION_PAIRS dictionary in the CONFIG SECTION.
For each reaction pair, the script will:
  - Read the reactant and product .com files.
  - Extract the charge/multiplicity, geometry (Z-matrix block) and connectivity.
  - Translate the product geometry by a fixed amount along the x-axis.
  - Create a new Gaussian input file for a QST2 job that combines both sections.

All configuration settings are self-contained in the CONFIG SECTION below.
To modify any settings, just edit the values in that section.
  
Usage:
    python create_qst2_jobs_v2.py
"""

import os
from pathlib import Path

# ============================================================================
# CONFIG SECTION - All configurable parameters for QST2 job generation
# ============================================================================

# Default configuration
DEFAULT_CONFIG = {
    # Input/output directories
    'input_dir': "../geom_optimise_guassian/gaussian_projects",
    'output_dir': "qst2_jobs",
    
    # Geometry settings
    'separation': 3.0,  # Separation distance (in Å) for product geometry
    
    # Gaussian calculation settings
    'level_of_theory': "m062x/def2tzvp",
    'mem': "8GB",
    'nproc': 8,
    'add_freq': True,    # Add frequency calculation
    'geom_connect': True, # Use connectivity information
    'int_ultrafine': True, # Use ultrafine integration grid
    'scf_tight': True,   # Use tight SCF convergence
    'scf_xqc': True,     # Use XQC for SCF convergence
    
    # Execution mode
    'mode': 'all',       # 'all', 'reaction', or 'list'
    'reaction': None,    # Specific reaction to process (if mode='reaction')
    
    # Output control
    'verbose': True      # Print detailed progress information
}

# File naming patterns (for file extensions)
FILE_EXTENSION = "_zmat.com"
OUTPUT_EXTENSION = ".gjf"

# ============================================================================
# END CONFIG SECTION
# ============================================================================

def setup_reaction_pairs():
    """Dictionary of common reaction pairs based on molecules in generate_inputs.py"""
    return {
        "SN2_halogen_exchange_F_Cl": {
            "type": "decomposition",
            "reactants": ["Na", "CH3F"],
            "products": ["CH3Cl", "Na"]
        },
        "SN2_halogen_exchange_F_Br": {
            "type": "decomposition",
            "reactants": ["Na", "CH3F"],
            "products": ["CH3Br", "Na"]
        },
        "TS1M_PFMS_Pathway": {
            "type": "decomposition",
            "reactants": ["PFMS"],
            "products": ["TS1M_Product1", "HF"]
        },
        "TS2M_PFMS_Pathway": {
            "type": "decomposition",
            "reactants": ["PFMS"],
            "products": ["HCF3", "SO3"]
        },
        "TS3M_PFMS_Pathway": {
            "type": "decomposition",
            "reactants": ["ISOPFMS"],
            "products": ["TS1M_Product1", "CF3_Radical"]
        },
        "TS4M_PFMS_Pathway": {
            "type": "decomposition",
            "reactants": ["PFMS"],
            "products": ["TS4M_Product1", "F_TS4M"]
        },
        "TS5M_PFMS_Pathway": {
            "type": "decomposition",
            "reactants": ["PFMS"],
            "products": ["TS5M_Product1", "F_TS5M"]
        },
        "TS6M_PFMS_Pathway": {
            "type": "decomposition",
            "reactants": ["PFMS"],
            "products": ["F_TS6M", "CF2O", "SO2"]
        }
    }

def ensure_directory(directory):
    """Create directory if it doesn't exist"""
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created directory: {directory}")

def parse_zmatrix_file(file_path):
    """
    Parse a Gaussian .com file in Z-matrix format.
    
    Assumes the file structure:
      - Optional route/header lines (starting with '#' or '%')
      - Blank line(s)
      - Title line
      - Blank line(s)
      - Charge/multiplicity line (e.g., " 0,1" or "0 1")
      - Z-matrix geometry lines (atom symbol + numbers)
      - Blank line(s)
      - Connectivity lines (if any)
    
    Returns a dictionary with keys:
      'title': Title line,
      'charge_line': Charge/multiplicity line (as a string),
      'geom_lines': List of geometry lines,
      'connectivity_lines': List of connectivity lines.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    i = 0
    # Skip any header/route lines starting with '#' or '%'
    while i < len(lines) and lines[i].strip().startswith(("#", "%")):
        i += 1
    # Skip blank lines
    while i < len(lines) and not lines[i].strip():
        i += 1
    # Get title line
    title = lines[i].strip() if i < len(lines) else ""
    i += 1
    # Skip blank lines
    while i < len(lines) and not lines[i].strip():
        i += 1
    # Charge/multiplicity line
    charge_line = lines[i].strip()
    i += 1
    # Read geometry lines until a blank line is encountered
    geom_lines = []
    while i < len(lines) and lines[i].strip():
        geom_lines.append(lines[i])
        i += 1
    # Skip blank lines
    while i < len(lines) and not lines[i].strip():
        i += 1
    # The rest are connectivity lines (if any)
    connectivity_lines = []
    while i < len(lines) and lines[i].strip():
        connectivity_lines.append(lines[i])
        i += 1

    return {
        'title': title,
        'charge_line': charge_line,
        'geom_lines': geom_lines,
        'connectivity_lines': connectivity_lines
    }

def translate_geometry(geom_lines, dx=0.0, dy=0.0, dz=0.0):
    """
    Translate the geometry lines by adding dx, dy, dz to each coordinate.
    Assumes each geometry line is in the format:
        AtomSymbol   x   y   z
    Returns a list of modified geometry lines (with similar formatting).
    """
    new_lines = []
    for line in geom_lines:
        tokens = line.split()
        # Expect at least 4 tokens: symbol, x, y, z.
        if len(tokens) < 4:
            new_lines.append(line)
            continue
        symbol = tokens[0]
        try:
            x = float(tokens[1])
            y = float(tokens[2])
            z = float(tokens[3])
        except ValueError:
            new_lines.append(line)
            continue
        new_x = x + dx
        new_y = y + dy
        new_z = z + dz
        # Format similar to original with fixed-width formatting.
        new_line = f" {symbol:<2s}    {new_x: .12f}     {new_y: .12f}     {new_z: .12f}\n"
        new_lines.append(new_line)
    return new_lines

def build_route_line(config):
    """Build the Gaussian route line based on configuration settings"""
    route = f"# opt=(maxcyc=999,noeigen) {config['level_of_theory']}"
    
    if config['add_freq']:
        route += " freq"
    if config['geom_connect']:
        route += " geom=connectivity"
    if config['int_ultrafine']:
        route += " int=ultrafine"
    
    scf_options = []
    if config['scf_tight']:
        scf_options.append("tight")
    if config['scf_xqc']:
        scf_options.append("xqc")
    
    if scf_options:
        route += f" scf=({','.join(scf_options)})"
    
    return route

def write_qst2_job(reactant_data, product_data, output_file, config, reaction_name):
    """
    Write a Gaussian QST2 input file that combines the reactant and product geometries.
    
    The product geometry is translated by 'separation' along the x-axis.
    """
    # Build route line from config
    route_line = build_route_line(config)
    
    # Use route line from config
    header = f"%mem={config['mem']}\n"
    header += f"%nprocshared={config['nproc']}\n"
    header += f"%chk=QST2_{reaction_name}.chk\n"
    header += f"{route_line}\n\n"
    
    job_title = f"QST2 job for {reaction_name}: {reactant_data['title']} to {product_data['title']}\n\n"
    
    # Get charge/multiplicity lines
    reactant_charge = reactant_data['charge_line'] + "\n"
    product_charge = product_data['charge_line'] + "\n"
    
    # Translate the product geometry
    translated_product_geom = translate_geometry(product_data['geom_lines'], dx=config['separation'])
    
    with open(output_file, 'w') as f:
        f.write(header)
        f.write(job_title)
        # Write reactant section
        f.write("Reactant geometry:\n\n")
        f.write(reactant_charge)
        for line in reactant_data['geom_lines']:
            f.write(line)
        f.write("\n")
        if reactant_data['connectivity_lines']:
            for line in reactant_data['connectivity_lines']:
                f.write(line)
            f.write("\n")
        # Blank line between sections
        f.write("\n")
        # Write product section
        f.write("Product geometry:\n\n")
        f.write(product_charge)
        for line in translated_product_geom:
            f.write(line)
        f.write("\n")
        if product_data['connectivity_lines']:
            for line in product_data['connectivity_lines']:
                f.write(line)
            f.write("\n")
        f.write("\n")
    print(f"Wrote QST2 input to {output_file}")

def list_available_molecules(input_dir):
    """List all available molecules from .com files in the input directory"""
    com_files = list(Path(input_dir).glob(f"*{FILE_EXTENSION}"))
    print("Available molecules:")
    for i, file in enumerate(com_files):
        print(f"{i+1}. {file.stem.replace(FILE_EXTENSION, '')}")
    return com_files

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
            reaction_type = reaction_pairs[name]["type"]
            print(f"{i+1}. {name} ({reaction_type}): {reactants} -> {products}")
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
        reaction_type = reaction["type"]
        
        if config['verbose']:
            print(f"Reaction type: {reaction_type}")
        
        # For reactions with multiple reactants/products, we need to create a combined structure
        # This is just a placeholder - in practice, you would need to align these molecules properly
        # For now, we'll just use the first reactant and first product
        reactant_file = os.path.join(config['input_dir'], f"{reaction['reactants'][0]}{FILE_EXTENSION}")
        product_file = os.path.join(config['input_dir'], f"{reaction['products'][0]}{FILE_EXTENSION}")
        
        if not os.path.exists(reactant_file):
            print(f"Warning: Reactant file {reactant_file} not found")
            continue
            
        if not os.path.exists(product_file):
            print(f"Warning: Product file {product_file} not found")
            continue
        
        # Parse files and create QST2 job
        reactant_data = parse_zmatrix_file(reactant_file)
        product_data = parse_zmatrix_file(product_file)
        
        # For decomposition reactions, we need to ensure the product atoms are positioned
        # relative to the reactant atoms. This is handled by the translate_geometry function
        # which moves the product geometry by the specified separation.
        if reaction_type == "decomposition":
            if config['verbose']:
                print(f"Handling decomposition reaction: Moving product atoms relative to reactant atoms")
        
        # Output QST2 file
        output_file = os.path.join(config['output_dir'], f"QST2_{reaction_name}{OUTPUT_EXTENSION}")
        write_qst2_job(reactant_data, product_data, output_file, config, reaction_name)

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
