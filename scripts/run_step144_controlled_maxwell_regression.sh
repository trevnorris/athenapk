#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
INPUT="${REPO_ROOT}/inputs/em4d_pulse.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step144_controlled_maxwell_regression_local"

NX1="64"
NX2="16"
NX3="16"
MB1="32"
MB2="16"
MB3="16"
TLIM="0.35"
NLIM="20000"
HISTORY_DT="0.002"

SMOKE_TLIM="0.05"
SMOKE_NLIM="200"

usage() {
  cat <<'USAGE'
Run Step-144 controlled-limit Maxwell regression suite.

Goal:
  Execute the GPT-14/14.5 recommended controlled-limit Maxwell checks before
  further Harris sweeps:
    1) controlled parity smoke on canonical hydro benchmarks
    2) analytic KK Coulomb+Yukawa consistency check
    3) strict controlled Nw=1 vacuum pulse in multiple quadrature settings
       to verify boundedness + quadrature invariance of key diagnostics

Usage:
  ./scripts/run_step144_controlled_maxwell_regression.sh [options]

Options:
  --binary PATH        Path to athenaPK binary
  --workdir DIR        Working directory
  --input PATH         Input deck (default: inputs/em4d_pulse.in)
  --output-root DIR    Output directory root
  --nx1 INT            Mesh nx1 (default: 64)
  --nx2 INT            Mesh nx2 (default: 16)
  --nx3 INT            Mesh nx3 (default: 16)
  --mb1 INT            Meshblock nx1 (default: 32)
  --mb2 INT            Meshblock nx2 (default: 16)
  --mb3 INT            Meshblock nx3 (default: 16)
  --tlim FLOAT         Runtime tlim for vacuum cases (default: 0.35)
  --nlim INT           Runtime nlim for vacuum cases (default: 20000)
  --history-dt FLOAT   HST dt for vacuum cases (default: 0.002)
  --smoke-tlim FLOAT   tlim for parity smoke (default: 0.05)
  --smoke-nlim INT     nlim for parity smoke (default: 200)
  -h, --help           Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --input) INPUT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --nx1) NX1="$2"; shift 2 ;;
    --nx2) NX2="$2"; shift 2 ;;
    --nx3) NX3="$2"; shift 2 ;;
    --mb1) MB1="$2"; shift 2 ;;
    --mb2) MB2="$2"; shift 2 ;;
    --mb3) MB3="$2"; shift 2 ;;
    --tlim) TLIM="$2"; shift 2 ;;
    --nlim) NLIM="$2"; shift 2 ;;
    --history-dt) HISTORY_DT="$2"; shift 2 ;;
    --smoke-tlim) SMOKE_TLIM="$2"; shift 2 ;;
    --smoke-nlim) SMOKE_NLIM="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if command -v realpath >/dev/null 2>&1; then
  BINARY="$(realpath "${BINARY}")"
  WORKDIR="$(realpath "${WORKDIR}")"
  INPUT="$(realpath "${INPUT}")"
  OUTPUT_ROOT="$(realpath -m "${OUTPUT_ROOT}")"
fi

if [[ ! -x "${BINARY}" ]]; then
  echo "Binary not found or not executable: ${BINARY}" >&2
  exit 1
fi
if [[ ! -f "${INPUT}" ]]; then
  echo "Input deck not found: ${INPUT}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_ROOT}"

CHECKS_CSV="${OUTPUT_ROOT}/checks.csv"
echo "check,status,log" > "${CHECKS_CSV}"

run_check() {
  local name="$1"
  shift
  local log_path="${OUTPUT_ROOT}/${name}.log"
  local status="PASS"
  if ! "$@" >"${log_path}" 2>&1; then
    status="FAIL"
  fi
  echo "${name},${status},${log_path}" >> "${CHECKS_CSV}"
  echo "${name},${status},${log_path}"
}

COMMON_ARGS=(
  "parthenon/mesh/nx1=${NX1}"
  "parthenon/mesh/nx2=${NX2}"
  "parthenon/mesh/nx3=${NX3}"
  "parthenon/meshblock/nx1=${MB1}"
  "parthenon/meshblock/nx2=${MB2}"
  "parthenon/meshblock/nx3=${MB3}"
  "parthenon/time/tlim=${TLIM}"
  "parthenon/time/nlim=${NLIM}"
  "parthenon/output0/dt=${HISTORY_DT}"
  "parthenon/output1/dt=-1"
  "modes4d/n_modes=1"
  "modes4d/response_w0_enable=false"
  "modes4d/response_lambda_enable=false"
  "modes4d/response_w0_dynamic_mixing_enable=false"
  "modes4d/response_w0_local_enable=false"
  "modes4d/response_w0_local_gradient_mixing_enable=false"
  "modes4d/response_w0_geometry_shift_enable=false"
  "modes4d/em_source_current_gain=0.0"
  "modes4d/em_source_timelike_gain=0.0"
  "modes4d/plasma_force_source_gain=0.0"
  "modes4d/plasma_momw_source_gain=0.0"
  "modes4d/plasma_momw_pressure_source_gain=0.0"
  "modes4d/plasma_energy_source_gain=0.0"
  "modes4d/plasma_w_flux_source_gain=0.0"
  "modes4d/plasma_pressure_transport_gain=0.0"
  "modes4d/plasma_pressure_rusanov_gain=0.0"
  "modes4d/plasma_momw_damping=1.0"
  "modes4d/hard_controlled_limit_enable=true"
  "modes4d/hard_controlled_limit_identity_mode0=true"
  "modes4d/hard_controlled_limit_strict_solver_path=true"
)

run_case() {
  local case_name="$1"
  local nquad="$2"
  shift 2

  local case_dir="${OUTPUT_ROOT}/${case_name}"
  local log_path="${case_dir}/run.log"
  mkdir -p "${case_dir}"

  rm -f "${case_dir}/parthenon.out0.hst" "${case_dir}/parthenon.out1.hst"

  local cmd=("${BINARY}" "-i" "${INPUT}")
  for arg in "${COMMON_ARGS[@]}"; do
    cmd+=("${arg}")
  done
  cmd+=("modes4d/n_quadrature=${nquad}")
  for arg in "$@"; do
    cmd+=("${arg}")
  done

  echo "=== ${case_name} ===" | tee "${log_path}"
  echo "command,${cmd[*]}" | tee -a "${log_path}"
  set +e
  (cd "${case_dir}" && "${cmd[@]}") 2>&1 | tee -a "${log_path}"
  local rc=${PIPESTATUS[0]}
  set -e
  echo "${rc}" > "${case_dir}/exit_code.txt"
  echo "case_status,${case_name},rc=${rc}" | tee -a "${log_path}"
}

run_check "controlled_parity_smoke" \
  python3 "${REPO_ROOT}/scripts/modes4d_controlled_parity_regression.py" \
    --binary "${BINARY}" \
    --workdir "${WORKDIR}" \
    --output-dir "${OUTPUT_ROOT}/controlled_parity_smoke" \
    --smoke-tlim "${SMOKE_TLIM}" \
    --smoke-nlim "${SMOKE_NLIM}"

run_check "kk_coulomb_yukawa_analytic" \
  python3 "${REPO_ROOT}/scripts/kk_coulomb_yukawa_regression.py" \
    --lambdas "0.5,1.0,2.0" \
    --n-terms 64

run_case "maxwell_strict_q1" "1"
run_case "maxwell_strict_q4" "4"
run_case "maxwell_strict_q8" "8"

python3 "${REPO_ROOT}/scripts/summarize_step144_controlled_maxwell_regression.py" \
  --glob "${OUTPUT_ROOT}/maxwell_strict_q*" \
  --baseline-case "maxwell_strict_q1" \
  --checks-csv "${CHECKS_CSV}" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo ""
echo "Step-144 controlled Maxwell regression suite complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
