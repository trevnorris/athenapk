#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step119_transverse_pinning_audit_$(date +%Y%m%d_%H%M%S)"

NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
HISTORY_DT="0.002"
TLIM="0.80"
NLIM="80000"

usage() {
  cat <<'USAGE'
Run Step-119 transverse-channel pinning audit.

Goal:
  Determine whether runs are still pinned near controlled-limit channels
  (Ew/C/Jw inactive), and identify cases that activate physical transfer
  pathways (JwEw / leakage / topology EMF).

Usage:
  ./scripts/run_step119_transverse_pinning_audit.sh [options]

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
  --tlim FLOAT             Run tlim (default: 0.80)
  --nlim INT               Run nlim (default: 80000)
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

# Step-118-aligned baseline + local warp gradient channel.
COMMON_FULL_ARGS=(
  "modes4d/lambda=1.0"
  "modes4d/em_conservative_transport=true"
  "problem/harris_4d/drift_current_scale=8.0e-2"
  "problem/harris_4d/aw_mode2_amp=2.0e-2"
  "problem/harris_4d/a0_mode2_amp=1.0e-2"
  "problem/harris_4d/piw_mode2_amp=1.0e-2"
  "problem/harris_4d/aw_mode1_amp=4.0e-1"
  "problem/harris_4d/piw_mode1_amp=2.0e-1"
  "problem/harris_4d/ay_mode1_amp=2.0e-1"
  "problem/harris_4d/piy_mode1_amp=1.0e-1"
  "problem/harris_4d/drift_momw_scale=5.0e-2"
  "problem/harris_4d/aw_mode1_phase_x=3.141593e-1"
  "problem/harris_4d/aw_mode1_phase_z=1.099558e+0"
  "problem/harris_4d/piw_mode1_phase_x=2.513274e-1"
  "problem/harris_4d/piw_mode1_phase_z=8.482300e-1"
  "problem/harris_4d/ay_mode1_phase_x=6.283185e-1"
  "problem/harris_4d/ay_mode1_phase_z=4.398230e-1"
  "problem/harris_4d/piy_mode1_phase_x=5.654867e-1"
  "problem/harris_4d/piy_mode1_phase_z=3.141593e-1"
  "modes4d/em_source_current_gain=4.0"
  "modes4d/em_source_timelike_gain=4.0"
  "modes4d/plasma_momw_source_gain=4.0"
  "modes4d/plasma_momw_pressure_source_gain=4.0"
  "modes4d/em_source_damping_gain=0.25"
  "modes4d/response_lambda_enable=false"
  "modes4d/response_w0_enable=true"
  "modes4d/response_w0_dynamic_mixing_enable=true"
  "modes4d/response_w0_dynamic_mixing_gain=2.0"
  "modes4d/response_w0_local_enable=true"
  "modes4d/response_w0_local_mode1_amp=5.0e-2"
  "modes4d/response_w0_local_mode1_kx=3.141592653589793"
  "modes4d/response_w0_local_mode1_kz=6.283185307179586"
  "modes4d/response_w0_local_mode1_omega=0.0"
  "modes4d/response_w0_local_mode2_amp=3.0e-2"
  "modes4d/response_w0_local_mode2_kx=6.283185307179586"
  "modes4d/response_w0_local_mode2_kz=12.566370614359172"
  "modes4d/response_w0_local_mode2_phase=1.0471975511965976"
  "modes4d/response_w0_local_mode2_omega=0.0"
  "modes4d/response_w0_local_gradient_mixing_enable=true"
  "modes4d/response_w0_local_gradient_mixing_gain=4.0"
  "modes4d/response_w0_geometry_shift_enable=true"
  "modes4d/response_w0_projection_enable=false"
  "modes4d/response_w0_projection_gain=1.0"
  "modes4d/response_w0_projection_max_abs=0.0"
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
  "${cmd[@]}" 2>&1 | tee -a "${log_file}"
  local rc=${PIPESTATUS[0]}
  set -e

  echo "${rc}" > "${case_dir}/exit_code.txt"
  if [[ ${rc} -eq 0 ]]; then
    echo "case_status,${case_name},PASS" | tee -a "${log_file}"
  else
    echo "case_status,${case_name},FAIL,rc=${rc}" | tee -a "${log_file}"
  fi
}

# 0) Step-118-like baseline.
run_case "baseline_step118_fixed"

# 1) Stronger direct J^w kick.
run_case "jw_kick_hi" \
  "problem/harris_4d/drift_momw_scale=1.5e-1"

# 2) Stronger mixed-sector seeds (Aw/Piw + Ay/Piy + mode2 support).
run_case "mixed_seed_hi" \
  "problem/harris_4d/aw_mode1_amp=6.0e-1" \
  "problem/harris_4d/piw_mode1_amp=3.0e-1" \
  "problem/harris_4d/ay_mode1_amp=3.0e-1" \
  "problem/harris_4d/piy_mode1_amp=1.5e-1" \
  "problem/harris_4d/aw_mode2_amp=4.0e-2" \
  "problem/harris_4d/a0_mode2_amp=2.0e-2" \
  "problem/harris_4d/piw_mode2_amp=2.0e-2"

# 3) Combined forced activation (seeds + source gains + stronger J^w kick).
run_case "combined_forced_hi" \
  "problem/harris_4d/drift_momw_scale=1.5e-1" \
  "problem/harris_4d/aw_mode1_amp=6.0e-1" \
  "problem/harris_4d/piw_mode1_amp=3.0e-1" \
  "problem/harris_4d/ay_mode1_amp=3.0e-1" \
  "problem/harris_4d/piy_mode1_amp=1.5e-1" \
  "problem/harris_4d/aw_mode2_amp=4.0e-2" \
  "problem/harris_4d/a0_mode2_amp=2.0e-2" \
  "problem/harris_4d/piw_mode2_amp=2.0e-2" \
  "modes4d/em_source_current_gain=8.0" \
  "modes4d/em_source_timelike_gain=8.0" \
  "modes4d/plasma_momw_source_gain=8.0" \
  "modes4d/plasma_momw_pressure_source_gain=8.0" \
  "modes4d/em_source_damping_gain=0.20"

# 4) Sign-flipped J^w kick symmetry check.
run_case "combined_forced_hi_momw_neg" \
  "problem/harris_4d/drift_momw_scale=-1.5e-1" \
  "problem/harris_4d/aw_mode1_amp=6.0e-1" \
  "problem/harris_4d/piw_mode1_amp=3.0e-1" \
  "problem/harris_4d/ay_mode1_amp=3.0e-1" \
  "problem/harris_4d/piy_mode1_amp=1.5e-1" \
  "problem/harris_4d/aw_mode2_amp=4.0e-2" \
  "problem/harris_4d/a0_mode2_amp=2.0e-2" \
  "problem/harris_4d/piw_mode2_amp=2.0e-2" \
  "modes4d/em_source_current_gain=8.0" \
  "modes4d/em_source_timelike_gain=8.0" \
  "modes4d/plasma_momw_source_gain=8.0" \
  "modes4d/plasma_momw_pressure_source_gain=8.0" \
  "modes4d/em_source_damping_gain=0.20"

python3 "${REPO_ROOT}/scripts/summarize_transverse_channel_activation.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

python3 "${REPO_ROOT}/scripts/summarize_projected_closure.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_projected_closure.csv" \
  --out-md "${OUTPUT_ROOT}_projected_closure.md"

python3 "${REPO_ROOT}/scripts/summarize_step119_transverse_pinning.py" \
  --summary-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-csv "${OUTPUT_ROOT}_pinning.csv" \
  --out-md "${OUTPUT_ROOT}_pinning.md"

echo ""
echo "Step-119 transverse-channel pinning audit complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "projected_closure_md,${OUTPUT_ROOT}_projected_closure.md"
echo "pinning_md,${OUTPUT_ROOT}_pinning.md"
