#!/bin/bash
# Usage: ./run_madness_hbm96_medium.sh <input_file.in> [madqc_executable]

set -euo pipefail

die() {
  echo "ERROR: $*" >&2
  exit 1
}

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
  echo "Usage: $0 <input_file.in> [madqc_executable]"
  exit 1
fi

input_file=$1
madqc_bin=${2:-${MADQC_BIN:-madqc}}
input_abs=$(readlink -f "$input_file")
input_basename=$(basename "$input_file")
input_lower=$(echo "$input_basename" | tr '[:upper:]' '[:lower:]')
base_name=$(basename "$input_lower" .in)

partition=${SEAWULF_PARTITION:-hbm-medium-96core}
nodes=${SEAWULF_NODES:-8}
ntasks_per_node=${SEAWULF_TASKS_PER_NODE:-8}
mad_threads=${SEAWULF_MAD_THREADS:-10}
env_script=${SEAWULF_ENV_SCRIPT:-$HOME/load_xeonmax.sh}
launcher=${SEAWULF_LAUNCHER:-"mpirun --map-by numa numactl --preferred-many=8-15"}
log_root_default=/gpfs/projects/rjh/adrian/development/gecko/analysis/slurm_logs
log_root=${SLURM_LOG_ROOT:-$log_root_default}
log_dir=${log_root}/${partition}

bin_family="unknown"
case "$madqc_bin" in
  */40core/*) bin_family="core40" ;;
  */96core/*) bin_family="core96" ;;
esac

mpi_family="unknown"
if command -v ldd >/dev/null 2>&1 && [[ -x "$madqc_bin" ]]; then
  if ldd "$madqc_bin" 2>/dev/null | grep -qi "mvapich"; then
    mpi_family="mvapich"
  elif ldd "$madqc_bin" 2>/dev/null | grep -qi "openmpi"; then
    mpi_family="openmpi"
  fi
fi

if [[ "$bin_family" != "unknown" && "$bin_family" != "core96" ]]; then
  die "Binary appears to be for ${bin_family}; expected core96 build: ${madqc_bin}"
fi
if [[ "$mpi_family" == "mvapich" ]]; then
  die "Detected MVAPICH-linked binary on 96-core OpenMPI run path: ${madqc_bin}"
fi

mkdir -p "$log_dir"

log_out=${log_dir}/${base_name}-%j.out
log_err=${log_dir}/${base_name}-%j.err

copy_input=1
if [[ "$input_abs" == "$(readlink -f "$input_lower" 2>/dev/null || true)" ]]; then
  copy_input=0
fi

sbatch <<EOF2
#!/bin/bash
#SBATCH --job-name=${base_name}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${ntasks_per_node}
#SBATCH --time=08:00:00
#SBATCH --output=${log_out}
#SBATCH --error=${log_err}
#SBATCH --partition=${partition}

set -euo pipefail
module purge

# oneAPI setvars can trip under nounset; relax only while sourcing environment.
set +u
source "${env_script}"
set -u

export MAD_NUM_THREADS=${mad_threads}
export __mkl_tmp_TARGET_ARCH_DIR="${__mkl_tmp_TARGET_ARCH_DIR:-}"

if [[ "${copy_input}" == "1" ]]; then
  cp "${input_abs}" "${input_lower}"
fi

${launcher} "${madqc_bin}" --wf=response --input="${input_lower}"
EOF2

echo "Submitted ${base_name}"
echo "Slurm logs: ${log_dir}"
