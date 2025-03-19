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
    
    # Check for syntax issues in ADMP parameters
    if grep -q "'" "$input_file" || grep -q '"' "$input_file"; then
        echo "WARNING: Found quote characters in input file, removing them..."
        sed -i "s/'//g; s/\"//g" "$input_file"
        echo "Cleaned input file."
    fi
    
    # Check for incomplete ADMP parameters
    if grep -q "Rtemp=\s*$" "$input_file"; then
        echo "ERROR: Incomplete ADMP parameter (Rtemp=) detected."
        local temp=$(echo "$input_file" | grep -o "[0-9]\+K" | grep -o "[0-9]\+")
        if [ -n "$temp" ]; then
            echo "Fixing Rtemp parameter with temperature $temp"
            sed -i "s/Rtemp=\s*$/Rtemp=$temp/" "$input_file"
            echo "Fixed Rtemp parameter."
        else
            echo "Could not determine temperature from filename, skipping this file."
            return 1
        fi
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
    chk_dir="$results_dir/checkpoints"
    mkdir -p "$results_dir" "$chk_dir"
    
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
        
        # Process checkpoint files (main and trajectory)
        base_name=$(basename "$file" .gjf)
        
        # Process main checkpoint file if it exists
        if [ -f "${base_name}.chk" ]; then
            echo "Converting main checkpoint file to formatted checkpoint..."
            formchk "${base_name}.chk" "${base_name}.fchk"
            mv "${base_name}.chk" "${base_name}.fchk" "$chk_dir/"
            echo "Main checkpoint files moved to: $chk_dir"
        else
            echo "WARNING: Main checkpoint file not found"
        fi
        
        # Process trajectory checkpoint files for each frame
        echo "Checking for trajectory checkpoint files..."
        chk_count=0
        for traj_chk in ${base_name}.chk.*; do
            if [ -f "$traj_chk" ] && [ "$traj_chk" != "${base_name}.chk.*" ]; then
                chk_count=$((chk_count + 1))
                frame_num="${traj_chk##*.}"
                frame_name="${base_name}_frame${frame_num}"
                echo "Converting trajectory checkpoint file $traj_chk (frame $frame_num)..."
                formchk "$traj_chk" "${frame_name}.fchk"
                cp "$traj_chk" "$chk_dir/${frame_name}.chk"
                mv "${frame_name}.fchk" "$chk_dir/"
            fi
        done
        if [ $chk_count -gt 0 ]; then
            echo "Processed $chk_count trajectory checkpoint files"
            rm -f ${base_name}.chk.*
            rm -f ${base_name}.chk
        else
            echo "No trajectory checkpoint files found"
        fi
    else
        echo "WARNING: ADMP calculation for $molecule at $temp may not have completed successfully."
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