#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=5G
#SBATCH --partition=sapphire
#SBATCH --job-name=ADMP_decomp
#SBATCH --output=ADMP_decomp_%j.log

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0

export GAUSS_PDEF=${SLURM_CPUS_PER_TASK}

echo "Starting ADMP decomposition calculations"
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
    
    
    echo "Input file validation complete."
    return 0
}

# Process all ADMP input files in admp_jobs directory and its subdirectories
find ./admp_jobs -name "*_ADMP_*.gjf" | sort | while read -r file; do
    echo "Processing: $file"
    echo "----------------------------------------"
    
    # Create results directory based on file path
    filename=$(basename "$file")
    molecule=${filename%_ADMP_*}
    temp=${filename#*_ADMP_}
    temp=${temp%.gjf}
    
    results_dir="./admp_jobs/results/$molecule/$temp"
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
    echo "Running Gaussian for $molecule at $temp"
    g16 "$(basename "$file")"
    
    # Check completion status
    if grep -q "Normal termination" "$(basename "${file%.gjf}.log")"; then
        echo "ADMP calculation for $molecule at $temp completed successfully."
        
        # Process checkpoint file
        base_name=$(basename "$file" .gjf)
        if [ -f "${base_name}.chk" ]; then
            echo "Converting checkpoint file to formatted checkpoint..."
            formchk "${base_name}.chk" "${base_name}.fchk"

        else
            echo "WARNING: Checkpoint file not found"
        fi
    else
        echo "WARNING: ADMP calculation for $molecule at $temp may not have completed successfully."
        # Examine error message
        echo "Error details:"
        grep -A5 "Error termination" "$(basename "${file%.gjf}.log")" || echo "No specific error message found."
    fi
    
    echo "----------------------------------------"
    cd - > /dev/null
done

echo "All ADMP jobs completed."

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s