#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=5G
#SBATCH --partition=sapphire
#SBATCH --job-name=marshal_paper_opt
#SBATCH --output=marshal_paper_opt_%j.log

module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0

export GAUSS_PDEF=${SLURM_CPUS_PER_TASK}


# Process all .gjf files in the directory
for file in *.gjf; do
    echo "Processing $file..."
    g16 "$file"
    echo "Completed $file"
    echo "----------------------------------------"
done


##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s