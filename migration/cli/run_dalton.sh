#!/bin/bash
# Usage: ./run_dalton.sh <input_file> <mol_file>
# Example: ./run_dalton.sh raman.dal H2O-d-aug-cc-pVQZ.mol

set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <input_file> <mol_file>"
  exit 1
fi

input_file=$1
mol_file=$2
base_name=$(basename "$input_file" .dal)

# --- submit the job dynamically ---
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=dalton_${base_name}
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=96
#SBATCH --time=48:00:00
#SBATCH --output=${base_name}.out
#SBATCH --error=${base_name}.err
#SBATCH --partition=hbm-long-96core

module purge
source ~/load_xeonmax.sh

echo "Running Dalton job: ${base_name}"
echo "Input file: ${input_file}"
echo "Molecule file: ${mol_file}"
echo "--------------------------------------"

dalton -N \$SLURM_NTASKS -gb 10 -dal ${input_file} -mol ${mol_file}

echo "Job completed successfully."
EOF

