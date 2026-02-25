#!/bin/bash
# Usage: ./run_madness.sh <input_file.in> [madqc_executable]
# Example: ./run_madness.sh H2O_raman.in /path/to/madqc

set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
  echo "Usage: $0 <input_file.in> [madqc_executable]"
  exit 1
fi

input_file=$1
madqc_bin=${2:-${MADQC_BIN:-madqc}}
base_name=$(basename "$input_file" .in)

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${base_name}
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --time=08:00:00
#SBATCH --output=${base_name}.out
#SBATCH --error=${base_name}.err
#SBATCH --partition=hbm-long-96core

module purge
source ~/load_xeonmax.sh

export MAD_NUM_THREADS=10

mpirun --map-by numa "${madqc_bin}" --wf=response --input=${input_file}
EOF
