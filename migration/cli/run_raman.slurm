#!/bin/bash
#SBATCH --job-name=dalton_nlo      # Job name
#SBATCH --nodes=6                 # Number of nodes
#SBATCH --ntasks-per-node=96      # Number of tasks (cores) per node
#SBATCH --time=048:00:00            # Maximum runtime (dd:hh:mm:ss)
#SBATCH --output=dalton_nlo_%j.out # Output file name (includes job ID)
#SBATCH --error=dalton_nlo_%j.err  # Error file name (includes job ID)
#SBATCH --partition=hbm-long-96core



input=raman.dal
mol_file=H2O-d-aug-cc-pVQZ.mol

# Load required modules (adjust paths as needed)
module purge
source ~/load_xeonmax.sh

dalton -N $SLURM_NTASKS -gb 1 -dal ${input} -mol ${mol_file}

echo "Job completed successfully."
