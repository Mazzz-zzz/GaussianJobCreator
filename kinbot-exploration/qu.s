#! /bin/bash -f

#SBATCH -N 1
#SBATCH -c {ppn}
#SBATCH -q {queue_name}
#SBATCH -o {errdir}/{name}.stdout
#SBATCH -e {errdir}/{name}.err
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
{slurm_feature}

# Load required modules
module purge
module load Gaussian/g16c01-CUDA-11.7.0
module load Python/3.10.4
module load Anaconda3/2024.02-1

# Initialize conda
eval "$(conda shell.bash hook)"

# Activate KinBot environment
conda activate kinbot-env

# Set up environment variables
export OMP_NUM_THREADS={ppn}
export GAUSS_SCRDIR=/scratch/akhalilov

# Make sure we're in the right directory
cd $SLURM_SUBMIT_DIR

# Run the calculation
{command}