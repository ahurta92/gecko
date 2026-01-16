#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
input_file="mad.raman.in"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=mad.raman
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=8
#SBATCH --time=08:00:00
#SBATCH --output=mad.raman.out
#SBATCH --error=mad.raman.err
#SBATCH --partition=hbm-large-96core

module purge
source ~/load_xeonmax.sh

export MAD_NUM_THREADS=20

echo "Running MADNESS job: ${input_file}"
echo "--------------------------------------"

mpirun --map-by numa madqc --wf=response --input=${input_file}

echo "Job completed successfully."
EOF
