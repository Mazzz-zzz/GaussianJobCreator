#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --mem=8G
#SBATCH --partition=sapphire
#SBATCH --job-name=QST2_TS
#SBATCH --output=QST2_TS_%j.log

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0
module load Python/3.10.4

export GAUSS_PDEF=${SLURM_CPUS_PER_TASK}

echo "Starting QST2 Transition State calculations"
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

    
    # Verify two molecular specifications exist (required for QST2)
    local charge_mult_count=$(grep -c -E '^\s*-?[0-9]+\s+[0-9]+\s*$' "$input_file")
    if [ "$charge_mult_count" -lt 2 ]; then
        echo "ERROR: QST2 requires two molecular specifications (reactant and product). Only found $charge_mult_count."
        return 1
    fi
    
    echo "Input file validation complete."
    return 0
}

# Process all QST2 input files in qst2_jobs directory
find ./qst2_jobs -name "TS_*.gjf" | sort | while read -r file; do
    echo "Processing: $file"
    echo "----------------------------------------"
    
    # Create results directory based on file path
    filename=$(basename "$file")
    ts_name=${filename%.gjf}
    
    results_dir="./qst2_jobs/results/$ts_name"
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
        echo "QST2 calculation for $ts_name completed successfully."
        
        # Process checkpoint file
        base_name=$(basename "$file" .gjf)
        if [ -f "${base_name}.chk" ]; then
            echo "Converting checkpoint file to formatted checkpoint..."
            formchk "${base_name}.chk" "${base_name}.fchk"
            
            # Check if a true transition state was found (one imaginary frequency)
            if grep -q "OptimizedTSHasOneNegativeFreq" "$(basename "${file%.gjf}.log")"; then
                echo "Verified as valid transition state with one imaginary frequency!"
            else
                # Extract frequency information to check manually
                grep -A3 "Frequencies" "$(basename "${file%.gjf}.log")" | head -4
                echo "WARNING: Verification of transition state status could not be automated."
                echo "Please check the log file and verify the structure has exactly ONE imaginary frequency."
            fi
        else
            echo "WARNING: Checkpoint file not found"
        fi
    else
        echo "WARNING: QST2 calculation for $ts_name may not have completed successfully."
        # Examine error message
        echo "Error details:"
        grep -A5 "Error termination" "$(basename "${file%.gjf}.log")" || echo "No specific error message found."
    fi
    
    echo "----------------------------------------"
    cd - > /dev/null
done

echo "All QST2 transition state jobs completed."
echo "----------------------------------------"
echo "Next steps:"
echo "1. Verify each transition state has exactly ONE imaginary frequency"
echo "2. Examine the imaginary mode to confirm it follows the reaction coordinate"
echo "3. Consider running IRC calculations to confirm the pathway connects reactants and products"
echo "4. For challenging convergence cases, try QST3 or relaxed scan approaches"
echo "----------------------------------------"

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s 