#!/usr/bin/env python
# python Extract_Optimized_From_Gaussian.py filename

from __future__ import print_function
import sys, os

def extract_all(text, target):
    linenums = []
    # Start count at 1 because files start at line 1 not 0
    count = 1
    for line in text:
        if (line.find(target)) > -1:
            linenums.append(count)
        count += 1
    return linenums

code = {"1" : "H", "2" : "He", "3" : "Li", "4" : "Be", "5" : "B", \
"6"  : "C", "7"  : "N", "8"  : "O", "9" : "F", "10" : "Ne", \
"11" : "Na" , "12" : "Mg" , "13" : "Al" , "14" : "Si" , "15" : "P", \
"16" : "S"  , "17" : "Cl" , "18" : "Ar" , "19" : "K"  , "20" : "Ca", \
"21" : "Sc" , "22" : "Ti" , "23" : "V"  , "24" : "Cr" , "25" : "Mn", \
"26" : "Fe" , "27" : "Co" , "28" : "Ni" , "29" : "Cu" , "30" : "Zn", \
"31" : "Ga" , "32" : "Ge" , "33" : "As" , "34" : "Se" , "35" : "Br", \
"36" : "Kr" , "37" : "Rb" , "38" : "Sr" , "39" : "Y"  , "40" : "Zr", \
"41" : "Nb" , "42" : "Mo" , "43" : "Tc" , "44" : "Ru" , "45" : "Rh", \
"46" : "Pd" , "47" : "Ag" , "48" : "Cd" , "49" : "In" , "50" : "Sn", \
"51" : "Sb" , "52" : "Te" , "53" : "I"  , "54" : "Xe" , "55" : "Cs", \
"56" : "Ba" , "57" : "La" , "58" : "Ce" , "59" : "Pr" , "60" : "Nd", \
"61" : "Pm" , "62" : "Sm" , "63" : "Eu" , "64" : "Gd" , "65" : "Tb", \
"66" : "Dy" , "67" : "Ho" , "68" : "Er" , "69" : "Tm" , "70" : "Yb", \
"71" : "Lu" , "72" : "Hf" , "73" : "Ta" , "74" : "W"  , "75" : "Re", \
"76" : "Os" , "77" : "Ir" , "78" : "Pt" , "79" : "Au" , "80" : "Hg", \
"81" : "Tl" , "82" : "Pb" , "83" : "Bi" , "84" : "Po" , "85" : "At", \
"86" : "Rn" , "87" : "Fr" , "88" : "Ra" , "89" : "Ac" , "90" : "Th", \
"91" : "Pa" , "92" : "U"  , "93" : "Np" , "94" : "Pu" , "95" : "Am", \
"96" : "Cm" , "97" : "Bk" , "98" : "Cf" , "99" : "Es" ,"100" : "Fm", \
"101": "Md" ,"102" : "No" ,"103" : "Lr" ,"104" : "Rf" ,"105" : "Db", \
"106": "Sg" ,"107" : "Bh" ,"108" : "Hs" ,"109" : "Mt" ,"110" : "Ds", \
"111": "Rg" ,"112" : "Uub","113" : "Uut","114" : "Uuq","115" : "Uup", \
"116": "Uuh","117" : "Uus","118" : "Uuo"}

def process_log_file(logfile_fn):
    """Process a single ADMP log file and create corresponding xyz file"""
    logfile_bn = os.path.splitext(os.path.basename(logfile_fn))[0]
    
    try:
        with open(logfile_fn, 'r') as logfile_fh:
            text = logfile_fh.readlines()
        
        # Create output file in same directory as input
        output_dir = os.path.dirname(logfile_fn)
        outfile_path = os.path.join(output_dir, logfile_bn + '.xyz')
        
        # Find all "Input orientation:" headers
        input_orient_indices = []
        for i, line in enumerate(text):
            if "Input orientation:" in line:
                input_orient_indices.append(i)
        
        with open(outfile_path, 'w') as outfile:
            for step, idx in enumerate(input_orient_indices):
                # Skip header lines to reach coordinates
                coord_idx = idx + 5
                
                # Collect atoms until we hit the separator
                atoms = []
                while coord_idx < len(text) and "---" not in text[coord_idx]:
                    line = text[coord_idx].strip()
                    if line:  # Skip empty lines
                        parts = line.split()
                        if len(parts) >= 6 and parts[0].isdigit() and parts[1].isdigit():
                            atoms.append((parts[1], parts[3], parts[4], parts[5]))
                    coord_idx += 1
                
                # Write the atoms to file
                outfile.write(f"{len(atoms)}\n")
                outfile.write(f"Time step {step}\n")
                
                for atom in atoms:
                    atom_num, x, y, z = atom
                    atom_symbol = code[atom_num]
                    outfile.write(f"{atom_symbol} {x} {y} {z}\n")
        
        print(f"Successfully processed: {logfile_fn}")
        return True
    
    except Exception as e:
        print(f"Error processing {logfile_fn}: {str(e)}")
        return False

def main():
    # Walk through the admp_jobs/results directory
    results_dir = "./admp_jobs/results"
    if not os.path.exists(results_dir):
        print(f"Error: Results directory '{results_dir}' not found!")
        sys.exit(1)
    
    processed_count = 0
    error_count = 0
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk(results_dir):
        for file in files:
            if file.endswith('.log'):
                log_path = os.path.join(root, file)
                if process_log_file(log_path):
                    processed_count += 1
                else:
                    error_count += 1
    
    print(f"\nProcessing complete!")
    print(f"Successfully processed: {processed_count} files")
    print(f"Errors encountered: {error_count} files")

if __name__ == "__main__":
    main()