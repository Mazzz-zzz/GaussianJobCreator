# KinBot Setup

This directory contains the setup for [KinBot](https://github.com/zadorlab/KinBot), a tool for automated reaction mechanism generation.

## Prerequisites

KinBot requires the following software:

- Python >=3.8
- A PBS or SLURM workload manager
- Gaussian (tested with G09 revision D01)
  - Note: The `formchk` utility needs to be available on the head node where KinBot is running

## Dependencies

KinBot automatically installs the following Python libraries:
- numpy>=1.17.0
- ase>=3.19
- networkx
- rmsd

## Optional Dependencies

For additional functionality, you may want to install:
- MESS or MESMER (master equation solvers for rate coefficient calculations)
- PESViewer (for PES visualization)
- Pybel (Python wrapper of OpenBabel)
- RDKit (including Python bindings)

## Installation

1. Load Anaconda module:
   ```bash
   module load Anaconda3/2024.02-1
   ```

2. Initialize conda for your shell:
   ```bash
   eval "$(conda shell.bash hook)"
   ```

3. Activate the KinBot environment:
   ```bash
   conda activate kinbot-env
   ```

4. Install required packages:
   ```bash
   conda install -c conda-forge openbabel
   pip install kinbot
   ```

## Configuration

After installation, you'll need to:
1. Set up your computational environment
2. Configure your batch system (SLURM/PBS)
3. Set up your quantum chemistry software (Gaussian)

## Usage

Basic usage:
```bash
# First activate the environment
module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate kinbot-env

# Then run KinBot
kinbot input.py
```

For more detailed instructions, refer to the [KinBot documentation](https://github.com/zadorlab/KinBot/wiki). 