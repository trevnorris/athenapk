#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step116_local_warp_gradient_probe_$(date +%Y%m%d_%H%M%S)"

NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
HISTORY_DT="0.002"
TLIM="0.35"
NLIM="30000"

usage() {
  cat <<'USAGE'
Run Step-116 local warp gradient-mixing probe.

Goal:
  Add local-warp spatial-gradient couplings (∇w0) on top of the existing
  w0-dot conservative mixing and measure whether physically anchored channels
  move off prior floors.

Usage:
  ./scripts/run_step116_local_warp_gradient_probe.sh [options]

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
  --tlim FLOAT             Run tlim (default: 0.35)
  --nlim INT               Run nlim (default: 30000)
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
  "modes4d/response_w0_init=0.0"
  "modes4d/response_w0_dot_init=0.0"
  "modes4d/response_w0_drive_gain=0.0"
  "modes4d/response_w0_drive_from_bridge_gain=0.0"
  "modes4d/response_w0_drive_from_work_gain=0.0"
  "modes4d/response_w0_drive_from_leak_gain=0.0"
  "modes4d/response_w0_stiffness=0.0"
  "modes4d/response_w0_damping=0.0"
  "modes4d/response_w0_max_abs=1.0"
  "modes4d/response_w0_projection_enable=false"
  "modes4d/response_w0_geometry_shift_enable=false"
  "modes4d/response_w0_lambda_shift_enable=false"
  "modes4d/response_w0_dynamic_mixing_enable=true"
  "modes4d/response_w0_dynamic_mixing_gain=2.0"
  "modes4d/response_w0_local_enable=true"
  "modes4d/response_w0_local_mode1_amp=5.0e-2"
  "modes4d/response_w0_local_mode1_kx=3.141592653589793"
  "modes4d/response_w0_local_mode1_kz=6.283185307179586"
  "modes4d/response_w0_local_mode1_phase=0.0"
  "modes4d/response_w0_local_mode1_omega=4.0"
  "modes4d/response_w0_local_mode2_amp=3.0e-2"
  "modes4d/response_w0_local_mode2_kx=6.283185307179586"
  "modes4d/response_w0_local_mode2_kz=12.566370614359172"
  "modes4d/response_w0_local_mode2_phase=1.0471975511965976"
  "modes4d/response_w0_local_mode2_omega=2.0"
  "modes4d/response_w0_local_gradient_mixing_enable=false"
  "modes4d/response_w0_local_gradient_mixing_gain=0.0"
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

# Dot-only local warp (reference)
run_case "local_dual_dotmix_only" \
  "modes4d/response_w0_local_gradient_mixing_enable=false" \
  "modes4d/response_w0_local_gradient_mixing_gain=0.0"

# Dot + gradient local warp couplings
run_case "local_dual_dotmix_grad_g1" \
  "modes4d/response_w0_local_gradient_mixing_enable=true" \
  "modes4d/response_w0_local_gradient_mixing_gain=1.0"

run_case "local_dual_dotmix_grad_g4" \
  "modes4d/response_w0_local_gradient_mixing_enable=true" \
  "modes4d/response_w0_local_gradient_mixing_gain=4.0"

# Gradient-only local warp coupling (static local warp, no local dot drive)
run_case "local_dual_static_grad_g4" \
  "modes4d/response_w0_local_mode1_omega=0.0" \
  "modes4d/response_w0_local_mode2_omega=0.0" \
  "modes4d/response_w0_local_gradient_mixing_enable=true" \
  "modes4d/response_w0_local_gradient_mixing_gain=4.0"

python3 "${REPO_ROOT}/scripts/summarize_transverse_channel_activation.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

python3 "${REPO_ROOT}/scripts/summarize_projected_closure.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_projected_closure.csv" \
  --out-md "${OUTPUT_ROOT}_projected_closure.md"

python3 "${REPO_ROOT}/scripts/summarize_local_warp_gradient_probe.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_local_warp_gradient.csv" \
  --out-md "${OUTPUT_ROOT}_local_warp_gradient.md"

echo ""
echo "Step-116 local warp gradient-mixing probe complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "projected_closure_md,${OUTPUT_ROOT}_projected_closure.md"
echo "local_warp_gradient_md,${OUTPUT_ROOT}_local_warp_gradient.md"
