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

# Input/output directories
INPUT_FOLDER = "../geom_optimise_gaussian/gaussian_projects"
OUTPUT_FOLDER = "qst2_jobs"

# Geometry settings
SEPARATION = 3.0  # Separation distance (in Å) for product geometry

# Gaussian calculation settings
GAUSSIAN_ROUTE_LINE = "# opt=(maxcyc=999,noeigen) freq m062x/def2tzvp geom=connectivity int=ultrafine scf=(tight,xqc)"

# File naming patterns (for file extensions)
FILE_EXTENSION = ".com"
OUTPUT_EXTENSION = ".gjf"

# Reaction pairs to process
# Format: {
#   "reaction_name": {
#     "reactant": "path/to/reactant_file",  # Relative to INPUT_FOLDER
#     "product": "path/to/product_file"     # Relative to INPUT_FOLDER
#   }
# }
REACTION_PAIRS = {
    "SN2_halogen_exchange_F_Cl": {
        "reactant": "Na_CH3F_reactant",
        "product": "Na_CH3Cl_product"
    },
    "SN2_halogen_exchange_F_Br": {
        "reactant": "Na_CH3F_reactant",
        "product": "Na_CH3Br_product"
    },
    "TS1M_PFMS_Pathway": {
        "reactant": "PFMS_reactant",
        "product": "TS1M_Product1_HF_product"
    },
    "TS2M_PFMS_Pathway": {
        "reactant": "PFMS_reactant",
        "product": "HCF3_SO3_product"
    }
}

# ============================================================================
# END CONFIG SECTION
# ============================================================================

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

def write_qst2_job(reactant_data, product_data, output_file, separation, reaction_name):
    """
    Write a Gaussian QST2 input file that combines the reactant and product geometries.
    
    The product geometry is translated by 'separation' along the x-axis.
    """
    # Use route line from config
    header = f"{GAUSSIAN_ROUTE_LINE}\n\n"
    job_title = f"QST2 job for {reaction_name}: {reactant_data['title']} to {product_data['title']}\n\n"
    
    # Get charge/multiplicity lines
    reactant_charge = reactant_data['charge_line'] + "\n"
    product_charge = product_data['charge_line'] + "\n"
    
    # Translate the product geometry
    translated_product_geom = translate_geometry(product_data['geom_lines'], dx=separation)
    
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

def main():
    # Create output directory if it does not exist
    output_dir = Path(OUTPUT_FOLDER)
    output_dir.mkdir(exist_ok=True)
    
    # Base path for input files
    input_dir = Path(INPUT_FOLDER)
    
    # Process each reaction pair from config
    for reaction_name, files in REACTION_PAIRS.items():
        print(f"Processing reaction: {reaction_name}")
        
        # Get full paths for reactant and product files
        reactant_file = input_dir / f"{files['reactant']}{FILE_EXTENSION}"
        product_file = input_dir / f"{files['product']}{FILE_EXTENSION}"
        
        # Check if both files exist
        if not reactant_file.exists():
            print(f"Warning: Reactant file not found: {reactant_file}")
            continue
        if not product_file.exists():
            print(f"Warning: Product file not found: {product_file}")
            continue
            
        # Parse files and create QST2 job
        reactant_data = parse_zmatrix_file(reactant_file)
        product_data = parse_zmatrix_file(product_file)
        
        # Output QST2 file
        output_file = output_dir / f"QST2_{reaction_name}{OUTPUT_EXTENSION}"
        write_qst2_job(reactant_data, product_data, output_file, SEPARATION, reaction_name)
        
if __name__ == "__main__":
    main()
