#!/bin/bash
# Usage: ./run_input.sh <input_file>
# Example: ./run_input.sh h2o_raman.in

# Exit on error and print commands
set -euo pipefail

# --- check argument ---
if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file=$1
base_name=$(basename "$input_file" .in)

# --- sbatch section ---
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${base_name}
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=2
#SBATCH --time=08:00:00
#SBATCH --output=${base_name}.out
#SBATCH --error=${base_name}.err
#SBATCH --partition=hbm-large-96core

module purge
source ~/load_xeonmax.sh

export MAD_NUM_THREADS=20

mpirun --map-by numa madqc --wf=response --input=${input_file}
EOF

