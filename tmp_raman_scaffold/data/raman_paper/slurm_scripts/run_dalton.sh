#!/bin/bash
# Usage: ./run_dalton.sh <input_file.dal> <molecule_file.mol>
# Examples:
#   ./run_dalton.sh optimize.dal H2O_aug-cc-pVDZ.mol
#   ./run_dalton.sh raman.dal opt_H2O_aug-cc-pVDZ.mol

set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <input_file.dal> <molecule_file.mol>"
  exit 1
fi

input_file=$1
mol_file=$2
base_name=$(basename "$input_file" .dal)

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=dalton_${base_name}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --time=48:00:00
#SBATCH --output=${base_name}.out
#SBATCH --error=${base_name}.err
#SBATCH --partition=hbm-1tb-long-96core

module purge
source ~/load_xeonmax.sh

echo "Running Dalton job: ${base_name}"
echo "Input file: ${input_file}"
echo "Molecule file: ${mol_file}"
echo "--------------------------------------"

dalton -N \$SLURM_NTASKS -gb 10 -dal ${input_file} -mol ${mol_file}

echo "Job completed successfully."
EOF
