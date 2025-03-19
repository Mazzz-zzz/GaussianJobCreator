#!/bin/bash
#SBATCH --account="punim0131"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=5G
#SBATCH --partition=sapphire
#SBATCH --job-name=ADMP_orbitals
#SBATCH --output=ADMP_orbitals_%j.log

# Load required modules
module purge
module load NVHPC/22.11-CUDA-11.7.0
module load Gaussian/g16c01-CUDA-11.7.0

echo "Starting orbital visualization generation from ADMP checkpoint files"
echo "-------------------------------------------------------------------"

# Function to check for existing cube files
check_existing_cube() {
    local cube_file="$1"
    local fchk_file="$2"
    local orbital="$3"
    
    if [ -f "$cube_file" ]; then
        echo "  - Cube file for $orbital already exists: $cube_file"
        return 0
    else
        echo "  - Generating $orbital cube file from: $fchk_file"
        return 1
    fi
}

echo "Processing 14 checkpoint files"
echo "----------------------------------------"

echo "Processing CF2O at 600K - 1 frames"
mkdir -p "./cube_files/CF2O/600K"
echo "\nProcessing file 1/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF2O/600K/CF2O_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/CF2O/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF2O/600K/CF2O_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF2O/600K/CF2O_ADMP_600K.fchk ./cube_files/CF2O/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF2O/600K/CF2O_ADMP_600K.fchk ./cube_files/CF2O/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/CF2O/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/CF2O/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF2O/600K/CF2O_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF2O/600K/CF2O_ADMP_600K.fchk ./cube_files/CF2O/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF2O/600K/CF2O_ADMP_600K.fchk ./cube_files/CF2O/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/CF2O/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing CF3 at Radical_600K - 1 frames"
mkdir -p "./cube_files/CF3/Radical_600K"
echo "\nProcessing file 2/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF3_Radical/600K/CF3_Radical_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/CF3/Radical_600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF3_Radical/600K/CF3_Radical_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF3_Radical/600K/CF3_Radical_ADMP_600K.fchk ./cube_files/CF3/Radical_600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF3_Radical/600K/CF3_Radical_ADMP_600K.fchk ./cube_files/CF3/Radical_600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/CF3/Radical_600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/CF3/Radical_600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF3_Radical/600K/CF3_Radical_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF3_Radical/600K/CF3_Radical_ADMP_600K.fchk ./cube_files/CF3/Radical_600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CF3_Radical/600K/CF3_Radical_ADMP_600K.fchk ./cube_files/CF3/Radical_600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/CF3/Radical_600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing CH3 at 600K - 1 frames"
mkdir -p "./cube_files/CH3/600K"
echo "\nProcessing file 3/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3/600K/CH3_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/CH3/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3/600K/CH3_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3/600K/CH3_ADMP_600K.fchk ./cube_files/CH3/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3/600K/CH3_ADMP_600K.fchk ./cube_files/CH3/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/CH3/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/CH3/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3/600K/CH3_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3/600K/CH3_ADMP_600K.fchk ./cube_files/CH3/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3/600K/CH3_ADMP_600K.fchk ./cube_files/CH3/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/CH3/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing CH3Cl at 600K - 1 frames"
mkdir -p "./cube_files/CH3Cl/600K"
echo "\nProcessing file 4/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3Cl/600K/CH3Cl_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/CH3Cl/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3Cl/600K/CH3Cl_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3Cl/600K/CH3Cl_ADMP_600K.fchk ./cube_files/CH3Cl/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3Cl/600K/CH3Cl_ADMP_600K.fchk ./cube_files/CH3Cl/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/CH3Cl/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/CH3Cl/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3Cl/600K/CH3Cl_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3Cl/600K/CH3Cl_ADMP_600K.fchk ./cube_files/CH3Cl/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3Cl/600K/CH3Cl_ADMP_600K.fchk ./cube_files/CH3Cl/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/CH3Cl/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing CH3F at 600K - 1 frames"
mkdir -p "./cube_files/CH3F/600K"
echo "\nProcessing file 5/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3F/600K/CH3F_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/CH3F/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3F/600K/CH3F_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3F/600K/CH3F_ADMP_600K.fchk ./cube_files/CH3F/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3F/600K/CH3F_ADMP_600K.fchk ./cube_files/CH3F/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/CH3F/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/CH3F/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3F/600K/CH3F_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3F/600K/CH3F_ADMP_600K.fchk ./cube_files/CH3F/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/CH3F/600K/CH3F_ADMP_600K.fchk ./cube_files/CH3F/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/CH3F/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing HCF3 at 600K - 1 frames"
mkdir -p "./cube_files/HCF3/600K"
echo "\nProcessing file 6/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCF3/600K/HCF3_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/HCF3/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCF3/600K/HCF3_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCF3/600K/HCF3_ADMP_600K.fchk ./cube_files/HCF3/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCF3/600K/HCF3_ADMP_600K.fchk ./cube_files/HCF3/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/HCF3/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/HCF3/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCF3/600K/HCF3_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCF3/600K/HCF3_ADMP_600K.fchk ./cube_files/HCF3/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCF3/600K/HCF3_ADMP_600K.fchk ./cube_files/HCF3/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/HCF3/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing HCl at 600K - 1 frames"
mkdir -p "./cube_files/HCl/600K"
echo "\nProcessing file 7/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCl/600K/HCl_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/HCl/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCl/600K/HCl_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCl/600K/HCl_ADMP_600K.fchk ./cube_files/HCl/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCl/600K/HCl_ADMP_600K.fchk ./cube_files/HCl/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/HCl/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/HCl/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCl/600K/HCl_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCl/600K/HCl_ADMP_600K.fchk ./cube_files/HCl/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HCl/600K/HCl_ADMP_600K.fchk ./cube_files/HCl/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/HCl/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing HF at 600K - 1 frames"
mkdir -p "./cube_files/HF/600K"
echo "\nProcessing file 8/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HF/600K/HF_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/HF/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HF/600K/HF_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HF/600K/HF_ADMP_600K.fchk ./cube_files/HF/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HF/600K/HF_ADMP_600K.fchk ./cube_files/HF/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/HF/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/HF/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HF/600K/HF_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HF/600K/HF_ADMP_600K.fchk ./cube_files/HF/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/HF/600K/HF_ADMP_600K.fchk ./cube_files/HF/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/HF/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing ISOPFMS at 600K - 1 frames"
mkdir -p "./cube_files/ISOPFMS/600K"
echo "\nProcessing file 9/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/ISOPFMS/600K/ISOPFMS_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/ISOPFMS/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/ISOPFMS/600K/ISOPFMS_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/ISOPFMS/600K/ISOPFMS_ADMP_600K.fchk ./cube_files/ISOPFMS/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/ISOPFMS/600K/ISOPFMS_ADMP_600K.fchk ./cube_files/ISOPFMS/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/ISOPFMS/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/ISOPFMS/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/ISOPFMS/600K/ISOPFMS_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/ISOPFMS/600K/ISOPFMS_ADMP_600K.fchk ./cube_files/ISOPFMS/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/ISOPFMS/600K/ISOPFMS_ADMP_600K.fchk ./cube_files/ISOPFMS/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/ISOPFMS/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing PFMS at 600K - 1 frames"
mkdir -p "./cube_files/PFMS/600K"
echo "\nProcessing file 10/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/PFMS/600K/PFMS_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/PFMS/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/PFMS/600K/PFMS_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/PFMS/600K/PFMS_ADMP_600K.fchk ./cube_files/PFMS/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/PFMS/600K/PFMS_ADMP_600K.fchk ./cube_files/PFMS/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/PFMS/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/PFMS/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/PFMS/600K/PFMS_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/PFMS/600K/PFMS_ADMP_600K.fchk ./cube_files/PFMS/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/PFMS/600K/PFMS_ADMP_600K.fchk ./cube_files/PFMS/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/PFMS/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing SO2 at 600K - 1 frames"
mkdir -p "./cube_files/SO2/600K"
echo "\nProcessing file 11/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO2/600K/SO2_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/SO2/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO2/600K/SO2_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO2/600K/SO2_ADMP_600K.fchk ./cube_files/SO2/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO2/600K/SO2_ADMP_600K.fchk ./cube_files/SO2/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/SO2/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/SO2/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO2/600K/SO2_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO2/600K/SO2_ADMP_600K.fchk ./cube_files/SO2/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO2/600K/SO2_ADMP_600K.fchk ./cube_files/SO2/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/SO2/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing SO3 at 600K - 1 frames"
mkdir -p "./cube_files/SO3/600K"
echo "\nProcessing file 12/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO3/600K/SO3_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/SO3/600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO3/600K/SO3_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO3/600K/SO3_ADMP_600K.fchk ./cube_files/SO3/600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO3/600K/SO3_ADMP_600K.fchk ./cube_files/SO3/600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/SO3/600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/SO3/600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO3/600K/SO3_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO3/600K/SO3_ADMP_600K.fchk ./cube_files/SO3/600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/SO3/600K/SO3_ADMP_600K.fchk ./cube_files/SO3/600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/SO3/600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing TS1M at Product1_600K - 1 frames"
mkdir -p "./cube_files/TS1M/Product1_600K"
echo "\nProcessing file 13/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS1M_Product1/600K/TS1M_Product1_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/TS1M/Product1_600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS1M_Product1/600K/TS1M_Product1_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS1M_Product1/600K/TS1M_Product1_ADMP_600K.fchk ./cube_files/TS1M/Product1_600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS1M_Product1/600K/TS1M_Product1_ADMP_600K.fchk ./cube_files/TS1M/Product1_600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/TS1M/Product1_600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/TS1M/Product1_600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS1M_Product1/600K/TS1M_Product1_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS1M_Product1/600K/TS1M_Product1_ADMP_600K.fchk ./cube_files/TS1M/Product1_600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS1M_Product1/600K/TS1M_Product1_ADMP_600K.fchk ./cube_files/TS1M/Product1_600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/TS1M/Product1_600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"
echo "Processing TS3M at Product1_600K - 1 frames"
mkdir -p "./cube_files/TS3M/Product1_600K"
echo "\nProcessing file 14/14: /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS3M_Product1/600K/TS3M_Product1_ADMP_600K.fchk (Frame 600)"

if ! check_existing_cube "./cube_files/TS3M/Product1_600K/homo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS3M_Product1/600K/TS3M_Product1_ADMP_600K.fchk" "HOMO"; then
    echo "  Running: cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS3M_Product1/600K/TS3M_Product1_ADMP_600K.fchk ./cube_files/TS3M/Product1_600K/homo_frame600.cube 80 h"
    cubegen 0 MO=HOMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS3M_Product1/600K/TS3M_Product1_ADMP_600K.fchk ./cube_files/TS3M/Product1_600K/homo_frame600.cube 80 h
    
    if [ -f "./cube_files/TS3M/Product1_600K/homo_frame600.cube" ]; then
        echo "  ✓ Successfully created HOMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create HOMO cube file for frame 600"
    fi
fi

if ! check_existing_cube "./cube_files/TS3M/Product1_600K/lumo_frame600.cube" "/home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS3M_Product1/600K/TS3M_Product1_ADMP_600K.fchk" "LUMO"; then
    echo "  Running: cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS3M_Product1/600K/TS3M_Product1_ADMP_600K.fchk ./cube_files/TS3M/Product1_600K/lumo_frame600.cube 80 h"
    cubegen 0 MO=LUMO /home/akhalilov/GaussianJobCreator/ADMP_decomposition_gaussian/admp_jobs/results/TS3M_Product1/600K/TS3M_Product1_ADMP_600K.fchk ./cube_files/TS3M/Product1_600K/lumo_frame600.cube 80 h
    
    if [ -f "./cube_files/TS3M/Product1_600K/lumo_frame600.cube" ]; then
        echo "  ✓ Successfully created LUMO cube file for frame 600"
    else
        echo "  ✗ ERROR: Failed to create LUMO cube file for frame 600"
    fi
fi

echo "----------------------------------------"

echo "Orbital visualization generation completed"
echo "Checking results..."

# Count generated cube files
CUBE_COUNT=$(find "$./cube_files" -name "*.cube" | wc -l)
echo "Number of cube files generated: $CUBE_COUNT"

if [ "$CUBE_COUNT" -gt 0 ]; then
    echo "Cube files were generated in the following locations:"
    find "$./cube_files" -type d -not -path "$./cube_files" | sort | head -n 10
    DIR_COUNT=$(find "$./cube_files" -type d -not -path "$./cube_files" | wc -l)
    if [ "$DIR_COUNT" -gt 10 ]; then
        echo "... and $(($DIR_COUNT - 10)) more directories"
    fi
    
    # Count frames processed
    echo "Frame summary by molecule/temperature:"
    for dir in $(find "$./cube_files" -mindepth 2 -type d | sort); do
        frame_count=$(find "$dir" -name "*.cube" | sort | wc -l)
        if [ "$frame_count" -gt 0 ]; then
            orbital_count=$(find "$dir" -name "*.cube" | sort | grep -o "_[^_]*\.cube" | sort | uniq | wc -l)
            echo "  - $dir: $(($frame_count / $orbital_count)) frames, $orbital_count orbitals per frame"
        fi
    done
else
    echo "WARNING: No cube files were generated. Check the log for errors."
fi

echo "Job completed"

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s
