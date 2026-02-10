#!/usr/bin/env bash
# Source this script to configure Gecko env vars and activate Gecko's venv.
#
# Usage:
#   source slurm_scripts/use_gecko_env.sh [/path/to/gecko]
#
# Resolution order for Gecko root:
#   1) first positional arg
#   2) existing GECKO_ROOT env var
#   3) GECKO_ROOT_DEFAULT env var
#   4) built-in default below

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "This script must be sourced, not executed." >&2
  echo "Example: source ${BASH_SOURCE[0]} /gpfs/projects/rjh/adrian/development/gecko" >&2
  exit 1
fi

DEFAULT_GECKO_ROOT="/gpfs/projects/rjh/adrian/development/gecko"

if [[ $# -gt 1 ]]; then
  echo "Usage: source ${BASH_SOURCE[0]} [/path/to/gecko]" >&2
  return 1
fi

_gecko_root="${1:-${GECKO_ROOT:-${GECKO_ROOT_DEFAULT:-${DEFAULT_GECKO_ROOT}}}}"
_gecko_root="${_gecko_root%/}"

if [[ ! -d "${_gecko_root}" ]]; then
  echo "[gecko-env] Gecko root not found: ${_gecko_root}" >&2
  return 1
fi

export GECKO_ROOT="${_gecko_root}"
export GECKO_SRC="${GECKO_ROOT}/src"

if [[ -d "${GECKO_SRC}" ]]; then
  case ":${PYTHONPATH:-}:" in
    *":${GECKO_SRC}:"*) ;;
    *) export PYTHONPATH="${GECKO_SRC}${PYTHONPATH:+:${PYTHONPATH}}" ;;
  esac
else
  echo "[gecko-env] warning: src directory missing at ${GECKO_SRC}" >&2
fi

if [[ -f "${GECKO_ROOT}/.venv/bin/activate" ]]; then
  # shellcheck disable=SC1090
  source "${GECKO_ROOT}/.venv/bin/activate"
  export GECKO_VENV="${GECKO_ROOT}/.venv"
else
  echo "[gecko-env] warning: venv not found at ${GECKO_ROOT}/.venv" >&2
  echo "[gecko-env] continuing without venv activation" >&2
fi

echo "[gecko-env] GECKO_ROOT=${GECKO_ROOT}"
echo "[gecko-env] GECKO_SRC=${GECKO_SRC}"

if command -v python >/dev/null 2>&1; then
  python - <<'PY'
import importlib.util
spec = importlib.util.find_spec("gecko")
if spec is None:
    print("[gecko-env] gecko import check: NOT FOUND")
else:
    print("[gecko-env] gecko import check: OK")
    print(f"[gecko-env] gecko module path: {spec.origin}")
PY
fi
