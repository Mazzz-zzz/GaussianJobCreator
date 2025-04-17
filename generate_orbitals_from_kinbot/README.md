# IRC to Orbitals

This workflow extracts molecular geometries from IRC (Intrinsic Reaction Coordinate) calculations performed with KinBot and generates orbital visualization files using Gaussian.

## Overview

The workflow consists of two main steps:
1. Extract geometries from IRC log files and create XYZ trajectory files
2. Generate Gaussian input files for orbital calculations from these trajectories

## Usage

### Combined Workflow

The easiest way to run the complete workflow is using the combined script `irc_to_orbitals.py`:

```bash
python irc_to_orbitals.py --irc-dir ../kinbot-exploration --output-dir ./orbital_inputs
```

This will:
1. Find all IRC log files in the KinBot exploration directory
2. Extract geometries and create XYZ trajectory files
3. Generate Gaussian input files for orbital calculations

### Running Calculations

Once you have generated the input files, you can run the calculations using the provided submission script:

```bash
sbatch submit_orbital_calculations.sh
```

You can also use the submission script to run the entire workflow in one step:

```bash
sbatch submit_orbital_calculations.sh --run-prep --kinbot ../kinbot-exploration
```

### Command Line Options

#### irc_to_orbitals.py

```
usage: irc_to_orbitals.py [-h] [--irc-dir IRC_DIR] [--xyz-dir XYZ_DIR]
                          [--output-dir OUTPUT_DIR] [--max-frames MAX_FRAMES]
                          [--method METHOD] [--basis BASIS] [--skip-irc]

Extract IRC geometries and generate orbital calculation inputs

optional arguments:
  -h, --help            show this help message and exit
  --irc-dir IRC_DIR     Base directory containing KinBot results
                        (default: ../kinbot-exploration)
  --xyz-dir XYZ_DIR     Directory to store extracted geometries
                        (default: ./irc_geometries)
  --output-dir OUTPUT_DIR
                        Directory to store orbital input files
                        (default: ./orbital_inputs)
  --max-frames MAX_FRAMES
                        Maximum number of frames per trajectory
                        (0 for all frames, default: 20)
  --method METHOD       Computational method to use (default: B3LYP)
  --basis BASIS         Basis set to use (default: 6-31G(d))
  --skip-irc            Skip IRC extraction (use if XYZ files
                        were already generated)
```

#### submit_orbital_calculations.sh

```
Usage: submit_orbital_calculations.sh [OPTIONS] [INPUT_DIR]
Run orbital calculations on geometries extracted from IRC paths

Options:
  -t, --time HOURS      Set time limit in hours (default: 24)
  -m, --memory GB       Set memory limit in GB (default: 8)
  -o, --output DIR      Set output directory for calculations (default: ./orbital_results)
  -i, --input DIR       Set input directory with Gaussian input files (default: ./orbital_inputs)
  -r, --run-prep        Run IRC to orbitals preparation before calculations
  -k, --kinbot DIR      Directory containing KinBot IRC log files (with -r flag)
  -x, --xyz DIR         Directory to store/read IRC geometries (with -r flag)
  -f, --frames N        Max frames to use per trajectory (with -r flag, default: 20)
  -q, --qm-method STR   Quantum method to use (with -r flag, default: B3LYP)
  -b, --basis STR       Basis set to use (with -r flag, default: 6-31G(d))
  -s, --skip-irc        Skip IRC extraction when preparing inputs (with -r flag)
  -h, --help            Show this help message
```

## Output

The workflow produces:

1. **XYZ trajectory files** in the `irc_geometries` directory, organized by reaction
2. **Gaussian input files** in the `orbital_inputs` directory
3. **Calculation results and cube files** in the `orbital_results` directory

The cube files can be visualized using programs like VMD, Avogadro, or PyMOL.

## Example

```bash
# Run the complete workflow with default settings
python irc_to_orbitals.py

# Run calculations with SLURM
sbatch submit_orbital_calculations.sh

# Run everything in one step with custom settings
sbatch submit_orbital_calculations.sh --run-prep --kinbot /path/to/kinbot --frames 30 --qm-method M062X --basis "6-311+G(d,p)"
``` 