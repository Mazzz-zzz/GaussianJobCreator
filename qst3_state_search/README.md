# QST3 Transition State Search

This tool generates Gaussian QST3 input files for transition state searches. Unlike QST2, QST3 requires an initial guess for the transition state structure, which this script can generate using several methods.

## Features

- Creates Gaussian QST3 input files from reactant and product structures
- Automatically generates transition state guesses using various methods:
  1. **Linear Synchronous Transit (LST)**: Simple linear interpolation between reactant and product
  2. **Quadratic Synchronous Transit (QST)**: Enhanced interpolation with bond-specific adjustments
  3. **Chemical Intuition**: Intelligently modifies bonds based on chemical knowledge of transition states
- Handles atom matching and alignment between reactant and product
- Identifies key bonds involved in the reaction
- Supports predefined reactions or custom molecule pairs
- Optional solvent effects (SCRF) and temperature specification
- Proper title cards for each structure as required by Gaussian

## Dependencies

- Python 3.6+
- NumPy
- SciPy
- RDKit
- ASE (Atomic Simulation Environment)

## Usage

### Basic Usage

```bash
python create_qst3_inputs.py --manual
```

This will:
1. Display available molecules from the input directory
2. Let you select reactant and product files
3. Let you select a transition state guessing method
4. Generate a QST3 input file

### Command Line Options

```bash
python create_qst3_inputs.py [options]
```

Options:
- `--input_dir DIR`: Directory containing optimized .gjf files
- `--output_dir DIR`: Directory to save QST3 input files
- `--level THEORY`: Level of theory for calculation (e.g., "b3lyp/6-31G(d)")
- `--mem MEMORY`: Memory allocation (e.g., "8GB")
- `--nproc N`: Number of processors to use
- `--no_freq`: Disable frequency calculation after optimization
- `--no_connect`: Disable use of connectivity specification
- `--ts_method {lst,qst,intuition}`: Method for generating transition state guess
- `--reaction NAME`: Specific reaction to set up (use --list_reactions to see options)
- `--list_reactions`: List available predefined reactions
- `--manual`: Manually select reactant and product files
- `--temperature VALUE`: Specify temperature for calculation (e.g., "298.15")
- `--solvent NAME`: Specify solvent for SCRF model (e.g., "Water")
- `--redundant`: Use redundant internal coordinates for optimization

### Examples

Generate QST3 inputs for all predefined reactions:
```bash
python create_qst3_inputs.py
```

Generate QST3 input for a specific reaction:
```bash
python create_qst3_inputs.py --reaction SN2_halogen_exchange_F_Cl
```

Generate QST3 input with specific settings:
```bash
python create_qst3_inputs.py --ts_method intuition --level "m062x/6-31++G(d,p)" --temperature 373.15 --solvent Water --redundant
```

## Transition State Guessing Methods

### 1. Linear Synchronous Transit (LST)
- The simplest method
- Creates a straight-line interpolation between reactant and product
- Uses a halfway point (t=0.5) as the TS guess

### 2. Quadratic Synchronous Transit (QST)
- More sophisticated than LST
- Identifies key bond changes in the reaction
- Places the TS guess at t=0.6 along the reaction path (biased towards products)
- Makes additional adjustments to key bonds based on bond forming/breaking rules

### 3. Chemical Intuition
- Most chemically informed method
- Starts with the reactant structure
- For breaking bonds: extends them by 50% of the total change
- For forming bonds: moves atoms 40% closer
- Particularly useful for reactions with clear bond changes

## QST3 Input File Format

The script generates QST3 input files with the following format:

```
%mem=8GB
%nprocshared=8
%chk=TS_reactant_to_product_QST3.chk

#p b3lyp/6-31G(d) Opt=QST3 Freq=NoRaman

Reactants: reactant_name

0 1
C    0.0000   0.0000   0.0000
...

Products: product_name

0 1
C    1.0000   0.0000   0.0000
...

Transition State Guess: reactant_name to product_name

0 1
C    0.5000   0.0000   0.0000
...
```

This format includes proper title cards for reactants, products, and the TS guess, following Gaussian's requirements for QST3 calculations.

## Running QST3 Jobs on HPC

This package includes a SLURM batch script for running QST3 jobs on HPC systems. The script handles validation, job submission, and basic error checking.

### Batch Processing with SLURM

1. Generate QST3 input files using the methods described above
2. Submit the jobs using the provided batch script:

```bash
sbatch submit_qst3_jobs.sh
```

The batch script will:
- Process all QST3 input files in the qst3_jobs directory
- Validate input files (checking for proper QST3 format)
- Run Gaussian calculations
- Convert checkpoint files to formatted checkpoints
- Verify transition state properties (imaginary frequencies)
- Provide diagnostics and error handling

### Customizing the Batch Script

You may need to modify the SLURM directives at the top of the script to match your HPC environment:

```bash
#SBATCH --account="your_account"    # Your account/project ID
#SBATCH --cpus-per-task=8           # Number of CPUs per task
#SBATCH --time=3-00:00:00           # Requested runtime (3 days)
#SBATCH --mem=8G                    # Memory request
#SBATCH --partition=your_partition  # Partition/queue name
```

### Comparing QST2 vs QST3 Results

The batch script includes functionality to compare QST3 results with QST2 results if available. QST3 often provides better convergence when a good initial guess is available, but QST2 can be simpler to set up and may work well for straightforward reactions.

## Tips for Verifying Transition States

After running the Gaussian calculations:
1. Check that the optimization converged to a stationary point
2. Verify that the structure has exactly ONE imaginary frequency
3. Examine the imaginary mode to confirm it corresponds to the reaction coordinate
4. For confirmation, perform IRC calculations from the transition state 