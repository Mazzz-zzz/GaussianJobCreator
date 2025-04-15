#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --partition=sapphire
#SBATCH --job-name=convert_chk
#SBATCH --output=convert_chk_%j.log
#SBATCH --error=convert_chk_%j.err

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0

# Convert the checkpoint file
formchk 1502984803620600000001/1502984803620600000001_r12_insertion_R_6_5_1_IRC_F.chk 1502984803620600000001/1502984803620600000001_r12_insertion_R_6_5_1_IRC_F.fchk

echo "Conversion complete!" 