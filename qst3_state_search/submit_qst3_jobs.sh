#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --mem=8G
#SBATCH --partition=sapphire
#SBATCH --job-name=QST3_TS
#SBATCH --output=QST3_TS_%j.log

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0
module load Python/3.10.4

export GAUSS_PDEF=${SLURM_CPUS_PER_TASK}

echo "Starting QST3 Transition State calculations"
echo "----------------------------------------"

# Function to validate and fix Gaussian input file
validate_input_file() {
    local input_file="$1"
    echo "Validating input file: $input_file"
    
    # Add %chk directive if not present
    if ! grep -q "^%chk=" "$input_file"; then
        local basename_no_ext=$(basename "$input_file" .gjf)
        echo "%chk=${basename_no_ext}.chk" | cat - "$input_file" > temp && mv temp "$input_file"
        echo "Added %chk directive to input file."
    fi
    
    # Verify three molecular specifications exist (required for QST3)
    local charge_mult_count=$(grep -c -E '^\s*-?[0-9]+\s+[0-9]+\s*$' "$input_file")
    if [ "$charge_mult_count" -lt 3 ]; then
        echo "ERROR: QST3 requires three molecular specifications (reactant, TS guess, and product). Only found $charge_mult_count."
        return 1
    fi
    
    # Check for title cards for each structure
    if ! grep -q "Reactants:" "$input_file" || ! grep -q "Products:" "$input_file" || ! grep -q "Transition State" "$input_file"; then
        echo "WARNING: Proper title cards for reactants, products, and transition state may be missing."
        echo "This might affect how Gaussian interprets the input."
    fi
    
    echo "Input file validation complete."
    return 0
}

# Process all QST3 input files in qst3_jobs directory
find ./qst3_jobs -name "TS_*_QST3.gjf" | sort | while read -r file; do
    echo "Processing: $file"
    echo "----------------------------------------"
    
    # Create results directory based on file path
    filename=$(basename "$file")
    ts_name=${filename%.gjf}
    
    results_dir="./qst3_jobs/results/$ts_name"
    mkdir -p "$results_dir"
    
    # Copy input file to results directory
    cp "$file" "$results_dir/"
    
    # Validate the input file
    if ! validate_input_file "$results_dir/$(basename "$file")"; then
        echo "SKIPPING due to validation errors: $file"
        echo "----------------------------------------"
        continue
    fi
    
    # Run Gaussian on the file
    cd "$results_dir"
    echo "Running Gaussian for $ts_name"
    g16 "$(basename "$file")"
    
    # Check completion status
    if grep -q "Normal termination" "$(basename "${file%.gjf}.log")"; then
        echo "QST3 calculation for $ts_name completed successfully."
        
        # Process checkpoint file
        base_name=$(basename "$file" .gjf)
        if [ -f "${base_name}.chk" ]; then
            echo "Converting checkpoint file to formatted checkpoint..."
            formchk "${base_name}.chk" "${base_name}.fchk"
            
            # Check if a true transition state was found (one imaginary frequency)
            if grep -q "1 imaginary frequencies" "$(basename "${file%.gjf}.log")"; then
                echo "Verified as valid transition state with one imaginary frequency!"
                
                # Extract the imaginary frequency value
                imag_freq=$(grep -A3 "Frequencies" "$(basename "${file%.gjf}.log")" | head -1 | awk '{print $3}')
                echo "Imaginary frequency: $imag_freq cm⁻¹"
                
                # Recommend IRC calculation based on the TS
                echo "Consider running IRC calculation to confirm reaction path:"
                echo "  Create an input file with:"
                echo "    # ${LEVEL_THEORY} IRC Geom=CheckPoint"
                echo "  And use the checkpoint file from this calculation."
            else
                # Extract frequency information to check manually
                echo "Frequency information:"
                grep -A3 "Frequencies" "$(basename "${file%.gjf}.log")" | head -4
                echo "WARNING: Verification of transition state status needs attention."
                echo "Check the log file and verify the structure has exactly ONE imaginary frequency."
                
                # Additional diagnostics
                if grep -q "0 imaginary frequencies" "$(basename "${file%.gjf}.log")"; then
                    echo "NOTE: This structure has 0 imaginary frequencies (minimum, not a TS)."
                    echo "Consider using a different initial TS guess or method."
                elif grep -q "[2-9] imaginary frequencies" "$(basename "${file%.gjf}.log")"; then
                    echo "NOTE: This structure has multiple imaginary frequencies (higher-order saddle point)."
                    echo "The optimization has not found a proper transition state."
                fi
            fi
        else
            echo "WARNING: Checkpoint file not found"
        fi
    else
        echo "WARNING: QST3 calculation for $ts_name may not have completed successfully."
        # Examine error message
        echo "Error details:"
        grep -A5 "Error termination" "$(basename "${file%.gjf}.log")" || echo "No specific error message found."
        
        # Check for common QST3-specific errors
        if grep -q "Coordinate system is invalid" "$(basename "${file%.gjf}.log")"; then
            echo "ERROR: Invalid coordinate system. Check molecular specifications."
        elif grep -q "GradGradGrad failed" "$(basename "${file%.gjf}.log")"; then
            echo "ERROR: Gradient calculation failed. Try a different TS guess or level of theory."
        elif grep -q "Convergence failure" "$(basename "${file%.gjf}.log")"; then
            echo "ERROR: Convergence failure. Try using a better initial guess or adjustment steps."
            echo "Suggestion: Try QST2 first, then use the resulting structure as input for a standard TS optimization."
        fi
    fi
    
    echo "----------------------------------------"
    cd - > /dev/null
done

echo "All QST3 transition state jobs completed."
echo "----------------------------------------"
echo "Next steps:"
echo "1. Verify each transition state has exactly ONE imaginary frequency"
echo "2. Examine the imaginary mode to confirm it follows the reaction coordinate"
echo "3. Run IRC calculations to confirm the pathway connects reactants and products"
echo "4. For structures with multiple imaginary frequencies, try refining the TS guess"
echo "5. For failed convergence, try different levels of theory or optimization settings"
echo "----------------------------------------"

# Compare QST3 results with any available QST2 results
if [ -d "../qst2_state_search/qst2_jobs/results" ]; then
    echo "QST2 results directory found. You may want to compare results from both methods."
    echo "QST3 often provides better convergence when a good initial guess is available."
fi

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s 