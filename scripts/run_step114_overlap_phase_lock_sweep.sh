#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step114_overlap_phase_lock_$(date +%Y%m%d_%H%M%S)"

NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
HISTORY_DT="0.002"
TLIM="0.60"
NLIM="60000"

usage() {
  cat <<'USAGE'
Run Step-114 overlap/phase-lock sweep.

Goal:
  Following Step-113 (channel components active but JwEw dormant), sweep phase
  alignment and drift sign/magnitude to maximize overlap in the JwEw work path.

Usage:
  ./scripts/run_step114_overlap_phase_lock_sweep.sh [options]

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
  --history-dt FLOAT       History dt (default: 0.002)
  --tlim FLOAT             Run tlim (default: 0.60)
  --nlim INT               Run nlim (default: 60000)
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
    --history-dt) HISTORY_DT="$2"; shift 2 ;;
    --tlim) TLIM="$2"; shift 2 ;;
    --nlim) NLIM="$2"; shift 2 ;;
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

COMMON_MESH_ARGS=(
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
  "modes4d/lambda=1.0"
  "modes4d/em_conservative_transport=true"
  "problem/harris_4d/drift_current_scale=8.0e-2"
  # Response/mixing off: pure overlap tuning.
  "modes4d/response_w0_enable=false"
  "modes4d/response_lambda_enable=false"
  "modes4d/response_w0_dynamic_mixing_enable=false"
  "modes4d/response_lambda_dynamic_mixing_enable=false"
  "modes4d/response_w0_local_enable=false"
  "modes4d/response_w0_projection_enable=false"
  "modes4d/response_w0_geometry_shift_enable=false"
  "modes4d/response_w0_lambda_shift_enable=false"
  # Default combined seed amplitudes.
  "problem/harris_4d/aw_mode1_amp=4.0e-1"
  "problem/harris_4d/piw_mode1_amp=2.0e-1"
  "problem/harris_4d/ay_mode1_amp=2.0e-1"
  "problem/harris_4d/piy_mode1_amp=1.0e-1"
  "problem/harris_4d/drift_momw_scale=5.0e-2"
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
    --fail-on-nonfinite-hst
    --check-history-nonfinite
    --check-transport-closure
    --check-em-bulk-ledger
    --check-topology
    --check-mechanism
  )

  for arg in "${COMMON_MESH_ARGS[@]}"; do
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

# Reference from Step-113 combined-phase settings.
run_case "ref_combined_phase" \
  "problem/harris_4d/aw_mode1_phase_x=3.141593e-1" \
  "problem/harris_4d/aw_mode1_phase_z=1.099558e+0" \
  "problem/harris_4d/piw_mode1_phase_x=2.513274e-1" \
  "problem/harris_4d/piw_mode1_phase_z=8.482300e-1" \
  "problem/harris_4d/ay_mode1_phase_x=6.283185e-1" \
  "problem/harris_4d/ay_mode1_phase_z=4.398230e-1" \
  "problem/harris_4d/piy_mode1_phase_x=5.654867e-1" \
  "problem/harris_4d/piy_mode1_phase_z=3.141593e-1"

# Phase-lock families.
run_case "lock_phi0_drift_p005" \
  "problem/harris_4d/drift_momw_scale=5.0e-2" \
  "problem/harris_4d/aw_mode1_phase_x=0.0" \
  "problem/harris_4d/aw_mode1_phase_z=0.0" \
  "problem/harris_4d/piw_mode1_phase_x=0.0" \
  "problem/harris_4d/piw_mode1_phase_z=0.0" \
  "problem/harris_4d/ay_mode1_phase_x=0.0" \
  "problem/harris_4d/ay_mode1_phase_z=0.0" \
  "problem/harris_4d/piy_mode1_phase_x=0.0" \
  "problem/harris_4d/piy_mode1_phase_z=0.0"

run_case "lock_phi0_drift_m005" \
  "problem/harris_4d/drift_momw_scale=-5.0e-2" \
  "problem/harris_4d/aw_mode1_phase_x=0.0" \
  "problem/harris_4d/aw_mode1_phase_z=0.0" \
  "problem/harris_4d/piw_mode1_phase_x=0.0" \
  "problem/harris_4d/piw_mode1_phase_z=0.0" \
  "problem/harris_4d/ay_mode1_phase_x=0.0" \
  "problem/harris_4d/ay_mode1_phase_z=0.0" \
  "problem/harris_4d/piy_mode1_phase_x=0.0" \
  "problem/harris_4d/piy_mode1_phase_z=0.0"

run_case "lock_phi_halfpi_drift_p005" \
  "problem/harris_4d/drift_momw_scale=5.0e-2" \
  "problem/harris_4d/aw_mode1_phase_x=1.570796e+0" \
  "problem/harris_4d/aw_mode1_phase_z=1.570796e+0" \
  "problem/harris_4d/piw_mode1_phase_x=1.570796e+0" \
  "problem/harris_4d/piw_mode1_phase_z=1.570796e+0" \
  "problem/harris_4d/ay_mode1_phase_x=1.570796e+0" \
  "problem/harris_4d/ay_mode1_phase_z=1.570796e+0" \
  "problem/harris_4d/piy_mode1_phase_x=1.570796e+0" \
  "problem/harris_4d/piy_mode1_phase_z=1.570796e+0"

run_case "lock_phi_pi_drift_p005" \
  "problem/harris_4d/drift_momw_scale=5.0e-2" \
  "problem/harris_4d/aw_mode1_phase_x=3.141593e+0" \
  "problem/harris_4d/aw_mode1_phase_z=3.141593e+0" \
  "problem/harris_4d/piw_mode1_phase_x=3.141593e+0" \
  "problem/harris_4d/piw_mode1_phase_z=3.141593e+0" \
  "problem/harris_4d/ay_mode1_phase_x=3.141593e+0" \
  "problem/harris_4d/ay_mode1_phase_z=3.141593e+0" \
  "problem/harris_4d/piy_mode1_phase_x=3.141593e+0" \
  "problem/harris_4d/piy_mode1_phase_z=3.141593e+0"

# Magnitude variants for overlap stress.
run_case "lock_phi0_drift_p010" \
  "problem/harris_4d/drift_momw_scale=1.0e-1" \
  "problem/harris_4d/aw_mode1_phase_x=0.0" \
  "problem/harris_4d/aw_mode1_phase_z=0.0" \
  "problem/harris_4d/piw_mode1_phase_x=0.0" \
  "problem/harris_4d/piw_mode1_phase_z=0.0" \
  "problem/harris_4d/ay_mode1_phase_x=0.0" \
  "problem/harris_4d/ay_mode1_phase_z=0.0" \
  "problem/harris_4d/piy_mode1_phase_x=0.0" \
  "problem/harris_4d/piy_mode1_phase_z=0.0"

run_case "lock_phi0_drift_m010" \
  "problem/harris_4d/drift_momw_scale=-1.0e-1" \
  "problem/harris_4d/aw_mode1_phase_x=0.0" \
  "problem/harris_4d/aw_mode1_phase_z=0.0" \
  "problem/harris_4d/piw_mode1_phase_x=0.0" \
  "problem/harris_4d/piw_mode1_phase_z=0.0" \
  "problem/harris_4d/ay_mode1_phase_x=0.0" \
  "problem/harris_4d/ay_mode1_phase_z=0.0" \
  "problem/harris_4d/piy_mode1_phase_x=0.0" \
  "problem/harris_4d/piy_mode1_phase_z=0.0"

python3 "${REPO_ROOT}/scripts/summarize_transverse_channel_activation.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

python3 "${REPO_ROOT}/scripts/summarize_projected_closure.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_projected_closure.csv" \
  --out-md "${OUTPUT_ROOT}_projected_closure.md"

python3 "${REPO_ROOT}/scripts/summarize_step114_overlap_phase_lock.py" \
  --summary-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-csv "${OUTPUT_ROOT}_overlap_ranked.csv" \
  --out-md "${OUTPUT_ROOT}_overlap_ranked.md"

echo ""
echo "Step-114 overlap/phase-lock sweep complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "projected_closure_md,${OUTPUT_ROOT}_projected_closure.md"
echo "overlap_ranked_md,${OUTPUT_ROOT}_overlap_ranked.md"
