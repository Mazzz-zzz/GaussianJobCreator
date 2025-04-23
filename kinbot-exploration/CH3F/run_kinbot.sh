#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --partition=sapphire
#SBATCH --job-name=kinbot_PFMS
#SBATCH --output=kinbot_%j.log
#SBATCH --error=kinbot_%j.err

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0
module load Python/3.10.4
module load Anaconda3/2024.02-1

# Initialize conda
eval "$(conda shell.bash hook)"

# Activate KinBot environment
conda activate kinbot-dev

# Run KinBot
kinbot PFOS-V1.json

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s 