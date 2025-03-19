#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --partition=sapphire
#SBATCH --job-name=XYZ_orbitals
#SBATCH --output=XYZ_orbitals_%j.log

# Help message
usage() {
    echo "Usage: $0 [OPTIONS] [INPUT_DIR]"
    echo "Run Gaussian calculations on input files"
    echo
    echo "Options:"
    echo "  -t, --time HOURS      Set time limit in hours (default: 24)"
    echo "  -m, --memory GB       Set memory limit in GB (default: 8)"
    echo "  -o, --output DIR      Set output directory (default: ./orbital_results)"
    echo "  -i, --input DIR       Set input directory (default: ./orbital_inputs)"
    echo "  -h, --help            Show this help message"
    echo
    echo "If no INPUT_DIR is specified, the default ./orbital_inputs will be used."
    echo
    exit 1
}

# Default directories
INPUT_DIR="./orbital_inputs"
OUTPUT_DIR="./orbital_results"
TIME_LIMIT="24:00:00"
MEMORY="8G"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--time)
            TIME_LIMIT="$2:00:00"
            shift 2
            ;;
        -m|--memory)
            MEMORY="${2}G"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            # If an argument is provided without a flag, assume it's the input directory
            INPUT_DIR="$1"
            shift
            ;;
    esac
done

# Check that input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "WARNING: Input directory does not exist: $INPUT_DIR"
    echo "Creating it now. Make sure to place your .gjf files there before running calculations."
    mkdir -p "$INPUT_DIR"
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0

export GAUSS_PDEF=4

echo "Starting orbital calculations"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "------------------------------------------------"

# Function to check for existing calculations
check_existing_calc() {
    local log_file="$1"
    local base_name="$2"
    
    if [ -f "$log_file" ]; then
        echo "  - Calculation already exists: $log_file"
        return 0
    else
        echo "  - Setting up calculation for: $base_name"
        return 1
    fi
}

# Find all Gaussian input files
GJF_FILES=$(find "$INPUT_DIR" -name "*.gjf")
GJF_COUNT=$(echo "$GJF_FILES" | wc -l)

if [ "$GJF_COUNT" -eq 0 ]; then
    echo "ERROR: No .gjf files found in $INPUT_DIR"
    echo "Please run the generate_orbitals_from_xyz.py script first to create input files,"
    echo "or specify a different input directory with the --input flag."
    exit 1
fi

echo "Found $GJF_COUNT Gaussian input files to process"
echo "------------------------------------------------"

# Process each input file
COUNTER=0
SUCCESS=0
FAILED=0
SKIPPED=0

# Get a sorted list of all input files
mapfile -t ALL_GJF_FILES < <(find "$INPUT_DIR" -name "*.gjf" | sort)

for gjf_file in "${ALL_GJF_FILES[@]}"; do
    COUNTER=$((COUNTER + 1))
    
    # Get base name and directory
    base_name=$(basename "$gjf_file" .gjf)
    input_dir=$(dirname "$gjf_file")
    
    # Extract molecule and temperature from directory structure
    rel_path=${input_dir#$INPUT_DIR/}
    molecule=$(echo "$rel_path" | cut -d '/' -f1)
    temp=$(echo "$rel_path" | cut -d '/' -f2)
    
    # Create corresponding output directory - ensure it exists
    output_subdir="$OUTPUT_DIR/$molecule/$temp"
    mkdir -p "$output_subdir"
    
    # Output files
    log_file="$output_subdir/${base_name}.log"
    
    echo "[$COUNTER/$GJF_COUNT] Processing: $base_name"
    echo "  - Input file: $gjf_file"
    echo "  - Output directory: $output_subdir"
    
    # Check if calculation already exists
    if ! check_existing_calc "$log_file" "$base_name"; then
        # Verify input file exists
        if [ ! -f "$gjf_file" ]; then
            echo "  ✗ ERROR: Input file does not exist: $gjf_file"
            FAILED=$((FAILED + 1))
            echo "------------------------------------------------"
            continue
        fi
        
        # Copy input file to output directory
        cp "$gjf_file" "$output_subdir/"
        
        # Run Gaussian calculation
        echo "  - Running Gaussian calculation..."
        
        # Ensure we're in the output directory before running Gaussian
        if ! cd "$output_subdir"; then
            echo "  ✗ ERROR: Cannot change to output directory: $output_subdir"
            FAILED=$((FAILED + 1))
            echo "------------------------------------------------"
            continue
        fi
        
        # Run Gaussian and capture output
        g16 "${base_name}.gjf" > "${base_name}_g16.out" 2>&1
        G16_STATUS=$?
        
        # Check if calculation succeeded
        if [ $G16_STATUS -eq 0 ]; then
            echo "  - Gaussian calculation completed successfully"
            
            # Generate cube files
            if [ -f "${base_name}.chk" ]; then
                echo "  - Generating formatted checkpoint and cube files"
                formchk "${base_name}.chk"
                cubegen 0 MO=HOMO "${base_name}.fchk" "${base_name}_homo.cube" 80 h
                cubegen 0 MO=LUMO "${base_name}.fchk" "${base_name}_lumo.cube" 80 h
                cubegen 0 density "${base_name}.fchk" "${base_name}_density.cube" 80 h
                ## Adding potential too!
                cubegen 0 Potential=scf "${base_name}.fchk" "${base_name}_pot.cube" 80 h
                echo "  ✓ Successfully created cube files"
                SUCCESS=$((SUCCESS + 1))
            else
                echo "  ✗ ERROR: Checkpoint file not found"
                FAILED=$((FAILED + 1))
            fi
        else
            echo "  ✗ ERROR: Gaussian calculation failed with status $G16_STATUS"
            echo "  - See ${output_subdir}/${base_name}_g16.out for details"
            FAILED=$((FAILED + 1))
        fi
        
        # Return to original directory
        cd "$SLURM_SUBMIT_DIR" || cd /tmp
    else
        SKIPPED=$((SKIPPED + 1))
    fi
    
    echo "------------------------------------------------"
done

echo "Orbital calculations completed"
echo "Summary:"
echo "  - Total input files: $GJF_COUNT"
echo "  - Successfully processed: $SUCCESS"
echo "  - Failed: $FAILED"
echo "  - Skipped (already existed): $SKIPPED"
echo

# Count generated files
LOG_COUNT=$(find "$OUTPUT_DIR" -name "*.log" | wc -l)
CUBE_COUNT=$(find "$OUTPUT_DIR" -name "*.cube" | wc -l)
echo "Number of completed calculations: $LOG_COUNT"
echo "Number of cube files generated: $CUBE_COUNT"

if [ "$CUBE_COUNT" -gt 0 ]; then
    echo "Results were generated in the following locations:"
    find "$OUTPUT_DIR" -type d -not -path "$OUTPUT_DIR" | sort | head -n 10
    DIR_COUNT=$(find "$OUTPUT_DIR" -type d -not -path "$OUTPUT_DIR" | wc -l)
    if [ "$DIR_COUNT" -gt 10 ]; then
        echo "... and $(($DIR_COUNT - 10)) more directories"
    fi
else
    echo "WARNING: No cube files were generated. Check the log for errors."
fi

echo "Job completed"

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s 