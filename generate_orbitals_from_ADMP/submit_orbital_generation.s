#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --mem=5G
#SBATCH --partition=sapphire
#SBATCH --job-name=orbital_gen
#SBATCH --output=orbital_gen_%j.log

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0
module load Python/3.10.4

echo "Starting orbital visualization generation"
echo "----------------------------------------"

# Generate the cubegen script
cd $SLURM_SUBMIT_DIR
python generate_inputs.py

# Check if script generation was successful
if [ -f "generate_cubes.sh" ]; then
    echo "Executing cubegen commands to generate orbital visualization files"
    bash generate_cubes.sh
    
    # Check for success
    if [ -d "./cube_files" ]; then
        echo "Orbital visualization completed successfully"
        echo "Number of cube files generated:"
        find ./cube_files -name "*.cube" | wc -l
    else
        echo "ERROR: Cube file generation may have failed"
    fi
else
    echo "ERROR: Failed to generate cubegen script"
    exit 1
fi

echo "----------------------------------------"
echo "Job completed"

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s 