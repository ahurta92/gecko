#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
input_file="optimize.dal"
mol_file="H2O_d-aug-cc-pVTZ.mol"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=dalton_optimize_H2O_d-aug-cc-pVTZ
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=96
#SBATCH --time=48:00:00
#SBATCH --output=optimize_H2O_d-aug-cc-pVTZ.out
#SBATCH --error=optimize_H2O_d-aug-cc-pVTZ.err
#SBATCH --partition=hbm-long-96core

module purge
source ~/load_xeonmax.sh

echo "Running Dalton job: optimize_H2O_d-aug-cc-pVTZ"
echo "Input file: ${input_file}"
echo "Molecule file: ${mol_file}"
echo "--------------------------------------"

dalton -N \$SLURM_NTASKS -gb 10 -dal \${input_file} -mol \${mol_file}

echo "Job completed successfully."
EOF
