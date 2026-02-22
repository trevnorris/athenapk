#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step91_lowgain_guarded_$(date +%Y%m%d_%H%M%S)"

NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
TLIM="0.35"
NLIM="10000"
HISTORY_DT="0.002"

usage() {
  cat <<'USAGE'
Run guarded low-gain bridge-response A/B cases on the Harris setup.

Usage:
  ./scripts/run_step91_lowgain_guarded.sh [options]

Options:
  --binary PATH            Path to athenaPK binary
  --workdir DIR            Working directory for runs
  --controlled-input PATH  Controlled input deck
  --full-input PATH        Full input deck
  --output-root DIR        Output root for all cases
  --nx1 INT                Mesh nx1 (default: 128)
  --nx2 INT                Mesh nx2 (default: 64)
  --nx3 INT                Mesh nx3 (default: 64)
  --mb1 INT                Meshblock nx1 (default: 32)
  --mb2 INT                Meshblock nx2 (default: 16)
  --mb3 INT                Meshblock nx3 (default: 16)
  --tlim FLOAT             Time limit (default: 0.35)
  --nlim INT               Cycle limit (default: 10000)
  --history-dt FLOAT       History dt (default: 0.002)
  -h, --help               Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --controlled-input) CONTROLLED_INPUT="$2"; shift 2 ;;
    --full-input) FULL_INPUT="$2"; shift 2 ;;
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

if [[ ! -x "${BINARY}" ]]; then
  echo "Binary not found or not executable: ${BINARY}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_ROOT}"

COMMON_ATHENA_ARGS=(
  "parthenon/mesh/nx1=${NX1}"
  "parthenon/mesh/nx2=${NX2}"
  "parthenon/mesh/nx3=${NX3}"
  "parthenon/meshblock/nx1=${MB1}"
  "parthenon/meshblock/nx2=${MB2}"
  "parthenon/meshblock/nx3=${MB3}"
  "parthenon/time/tlim=${TLIM}"
  "parthenon/time/nlim=${NLIM}"
  "parthenon/output0/dt=-1"
  "parthenon/output1/dt=${HISTORY_DT}"
  "diffusion/resistivity=ohmic"
  "diffusion/resistivity_coeff=fixed"
  "diffusion/integrator=rkl2"
  "diffusion/rkl2_max_dt_ratio=3.0"
  "diffusion/ohm_diff_coeff_code=1.788854e-02"
)

COMMON_FULL_ARGS=(
  "modes4d/em_conservative_transport=true"
  "problem/harris_4d/drift_current_scale=8.0e-2"
  "problem/harris_4d/aw_mode2_amp=2.0e-2"
  "problem/harris_4d/a0_mode2_amp=1.0e-2"
  "problem/harris_4d/piw_mode2_amp=1.0e-2"
)

run_case() {
  local case_name="$1"
  shift
  local case_dir="${OUTPUT_ROOT}/${case_name}"
  local log_file="${case_dir}/run.log"
  mkdir -p "${case_dir}"

  local cmd=(
    python3 "${REPO_ROOT}/scripts/harris_scan_matrix.py"
    --binary "${BINARY}"
    --workdir "${WORKDIR}"
    --controlled-input "${CONTROLLED_INPUT}"
    --full-input "${FULL_INPUT}"
    --output-dir "${case_dir}"
    --fail-on-check
    --fail-on-nonfinite-hst
    --check-history-nonfinite
    --check-transport-closure
    --check-em-bulk-ledger
  )

  for arg in "${COMMON_ATHENA_ARGS[@]}"; do
    cmd+=(--athena-arg "${arg}")
  done
  for arg in "${COMMON_FULL_ARGS[@]}"; do
    cmd+=(--full-arg "${arg}")
  done
  for arg in "$@"; do
    cmd+=(--full-arg "${arg}")
  done

  echo "=== ${case_name} ===" | tee "${log_file}"
  echo "command,${cmd[*]}" >> "${log_file}"

  set +e
  "${cmd[@]}" | tee -a "${log_file}"
  local rc=${PIPESTATUS[0]}
  set -e

  echo "${rc}" > "${case_dir}/exit_code.txt"
  if [[ ${rc} -eq 0 ]]; then
    echo "case_status,${case_name},PASS" | tee -a "${log_file}"
  else
    echo "case_status,${case_name},FAIL,rc=${rc}" | tee -a "${log_file}"
  fi
}

# Baseline (no response/bridge drive)
run_case "baseline_bridge0"

# Low-gain bridge-response routing checks
run_case "bridge_primary_g0p30" \
  "modes4d/em_source_mode0_even_bridge_gain=0.30" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_from_bridge_gain=0.30" \
  "modes4d/response_w0_drive_from_bridge_channel=primary" \
  "modes4d/response_w0_drive_from_bridge_abs=false"

run_case "bridge_primary_g1p00" \
  "modes4d/em_source_mode0_even_bridge_gain=1.00" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_from_bridge_gain=1.00" \
  "modes4d/response_w0_drive_from_bridge_channel=primary" \
  "modes4d/response_w0_drive_from_bridge_abs=false"

run_case "bridge_primary_abs_g1p00" \
  "modes4d/em_source_mode0_even_bridge_gain=1.00" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_from_bridge_gain=1.00" \
  "modes4d/response_w0_drive_from_bridge_channel=primary" \
  "modes4d/response_w0_drive_from_bridge_abs=true"

echo ""
echo "Step-91 guarded low-gain campaign complete."
echo "output_root,${OUTPUT_ROOT}"
