#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
INPUT="${REPO_ROOT}/inputs/em4d_pulse.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step136_laplacian_isolation_vacuum_local"

NX1="64"
NX2="16"
NX3="16"
MB1="32"
MB2="16"
MB3="16"
TLIM="1.00"
NLIM="120000"
HISTORY_DT="0.002"

usage() {
  cat <<'USAGE'
Run Step-136 Laplacian-path isolation in the fast vacuum harness.

Goal:
  Isolate whether EM Laplacian-path scaling/sign dominates controlled/vacuum
  growth, then cross-check interaction with damping.

Usage:
  ./scripts/run_step136_laplacian_isolation_vacuum.sh [options]

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
  --tlim FLOAT         Runtime tlim (default: 1.00)
  --nlim INT           Runtime nlim (default: 120000)
  --history-dt FLOAT   HST dt (default: 0.002)
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
  "modes4d/n_quadrature=1"
  "modes4d/response_w0_enable=false"
  "modes4d/response_lambda_enable=false"
  "modes4d/response_w0_dynamic_mixing_enable=false"
  "modes4d/response_w0_local_enable=false"
  "modes4d/response_w0_local_gradient_mixing_enable=false"
  "modes4d/response_w0_geometry_shift_enable=false"
  "modes4d/hard_controlled_limit_enable=false"
  "modes4d/hard_controlled_limit_identity_mode0=true"
  "modes4d/em_source_current_gain=0.0"
  "modes4d/em_source_timelike_gain=0.0"
  "modes4d/em_source_mass_gain=1.0"
  "modes4d/em_source_gauge_gain=0.0"
  "modes4d/plasma_force_source_gain=0.0"
  "modes4d/plasma_momw_source_gain=0.0"
  "modes4d/plasma_momw_pressure_source_gain=0.0"
  "modes4d/plasma_energy_source_gain=0.0"
  "modes4d/plasma_w_flux_source_gain=0.0"
  "modes4d/plasma_pressure_transport_gain=0.0"
  "modes4d/plasma_pressure_rusanov_gain=0.0"
  "modes4d/plasma_momw_damping=1.0"
)

run_case() {
  local case_name="$1"
  shift
  local case_dir="${OUTPUT_ROOT}/${case_name}"
  local log_path="${case_dir}/run.log"
  mkdir -p "${case_dir}"

  rm -f "${case_dir}/parthenon.out0.hst" "${case_dir}/parthenon.out1.hst"

  local cmd=("${BINARY}" "-i" "${INPUT}")
  for arg in "${COMMON_ARGS[@]}"; do
    cmd+=("${arg}")
  done
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

run_case "lap_baseline" \
  "modes4d/em_damping=0.0" \
  "modes4d/em_source_damping_gain=1.0" \
  "modes4d/em_source_laplacian_gain=1.0"

run_case "lap_off" \
  "modes4d/em_damping=0.0" \
  "modes4d/em_source_damping_gain=1.0" \
  "modes4d/em_source_laplacian_gain=0.0"

run_case "lap_half" \
  "modes4d/em_damping=0.0" \
  "modes4d/em_source_damping_gain=1.0" \
  "modes4d/em_source_laplacian_gain=0.5"

run_case "lap_flip" \
  "modes4d/em_damping=0.0" \
  "modes4d/em_source_damping_gain=1.0" \
  "modes4d/em_source_laplacian_gain=-1.0"

run_case "lap_baseline_damping_on" \
  "modes4d/em_damping=0.10" \
  "modes4d/em_source_damping_gain=1.0" \
  "modes4d/em_source_laplacian_gain=1.0"

run_case "lap_off_damping_on" \
  "modes4d/em_damping=0.10" \
  "modes4d/em_source_damping_gain=1.0" \
  "modes4d/em_source_laplacian_gain=0.0"

python3 "${REPO_ROOT}/scripts/summarize_step135_em_term_isolation.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --baseline-case "lap_baseline" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo ""
echo "Step-136 Laplacian isolation (vacuum) complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
