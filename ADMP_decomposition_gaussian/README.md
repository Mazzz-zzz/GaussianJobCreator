# ADMP Decomposition Calculator for Gaussian

This directory contains scripts to set up and run ADMP (Atom-centered Density Matrix Propagation) calculations in Gaussian for studying thermal decomposition pathways of molecules at various temperatures.

Reference: https://wanglab.hosted.uark.edu/g03guide/G03Guide/www.gaussian.com/g_ur/k_admp.htm

## What is ADMP?

ADMP is an ab initio molecular dynamics method that allows you to simulate the time evolution of molecular systems with quantum mechanical accuracy. It's particularly useful for:

- Studying thermal decomposition pathways
- Exploring reaction mechanisms without predefined reaction coordinates
- Finding unanticipated chemical transformations
- Observing structural changes at elevated temperatures

### ADMP Syntax in Gaussian

According to Gaussian documentation, the correct format for ADMP calculations is:

```
# Method/Basis ADMP [options]
```

Where options include:
- `MaxStep=n`: Maximum number of steps (default: 1000)
- `DeltaT=x`: Time step in femtoseconds (default: 0.1)
- `Temp=n`: Initial temperature in Kelvin
- `FullStep=n`: Frequency of full output (default: 1)

Example:
```
# B3LYP/6-31G(d) ADMP MaxStep=2000 DeltaT=0.5 Temp=1000 FullStep=100 
```

**Note**: Do not use the format `ADMP=(MaxPoints=n,DeltaT=x,...)` as this causes syntax errors in Gaussian.

## Files in this Directory

- `generate_admp_inputs.py`: Main script to generate Gaussian input files for ADMP calculations
- `submit_admp_jobs.sh`: SLURM submission script for running ADMP calculations on a cluster
- `admp_jobs/`: Directory containing the generated input files, organized by temperature

## How to Use

### 1. Generate ADMP Input Files

```bash
# Basic usage with default settings
python generate_admp_inputs.py

# Specify custom temperatures
python generate_admp_inputs.py --temperatures="600,800,1000,1200,1500"

# Customize other ADMP parameters
python generate_admp_inputs.py --max_points=1000 --delta_t=0.5 --method="M062X" --basis="def2svp"
```

Available options:
- `--input_dir`: Directory containing optimized molecule files
- `--output_dir`: Directory to store ADMP input files (default: "admp_jobs")
- `--temperatures`: Comma-separated list of temperatures in Kelvin (default: "600,800,1000,1200,1500")
- `--max_points`: Maximum number of steps for ADMP simulation (default: 2000)
- `--delta_t`: Time step for ADMP simulation in femtoseconds (default: 0.5)
- `--method`: Computational method for ADMP (default: "B3LYP")
- `--basis`: Basis set for ADMP (default: "6-31G(d)")
- `--mem`: Memory allocation for Gaussian (default: "8GB")
- `--nproc`: Number of processors for calculation (default: 8)

### 2. Submit ADMP Jobs to SLURM

```bash
# Generate the input files first
python generate_admp_inputs.py

# Submit job to process all ADMP input files
sbatch submit_admp_jobs.s
```

The script will:
- Find all ADMP input files (ending with *_ADMP_*.gjf) in the admp_jobs directory
- Process them one by one with Gaussian
- Create an organized results directory structure
- Display progress and results directly in the SLURM output file

This simplified approach uses a single command to process all files without requiring manual specification of which molecules or temperatures to run.

### 3. Analyzing Results

After ADMP calculations complete, analyze the trajectories:

1. Look for frames where bond breaking occurs in the `.log` files
2. Extract those geometries as potential transition states
3. Perform geometry optimizations on the resulting fragments
4. Run IRC calculations from transition states to confirm decomposition pathways

## Customizing the SLURM Script

Edit `submit_admp_jobs.sh` to:
- Change the list of molecules to study
- Modify time limits or resource requests
- Add different temperatures
- Adjust Gaussian module loading for your HPC system

## Best Practices

- Start with shorter simulations (lower `max_points`) to test
- Use higher temperatures (1000K+) to observe decomposition in a reasonable timeframe
- Watch for energy conservation issues in the output
- For production runs, use at least 2000 steps with a 0.5 fs time step
- Higher levels of theory will be more accurate but significantly more expensive 