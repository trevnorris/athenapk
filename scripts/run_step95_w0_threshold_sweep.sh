#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step95_w0_threshold_sweep_$(date +%Y%m%d_%H%M%S)"

NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
HISTORY_DT="0.002"

SHORT_TLIM="0.60"
SHORT_NLIM="60000"
LONG_TLIM="1.20"
LONG_NLIM="120000"

usage() {
  cat <<'USAGE'
Run Step-95 dynamic w0 threshold sweep (conservative, direct bridge forcing off).

Cases:
  - baseline_fixed_center_t0p6
  - dyn_g1_geom0p3_t0p6
  - dyn_g3_geom1p0_t0p6
  - dyn_g6_geom2p0_t0p6
  - dyn_g3_geom1p0_t1p2
  - dyn_g6_geom2p0_t1p2

Usage:
  ./scripts/run_step95_w0_threshold_sweep.sh [options]

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
  --short-tlim FLOAT       Short-run tlim (default: 0.60)
  --short-nlim INT         Short-run nlim (default: 60000)
  --long-tlim FLOAT        Long-run tlim (default: 1.20)
  --long-nlim INT          Long-run nlim (default: 120000)
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
    --short-tlim) SHORT_TLIM="$2"; shift 2 ;;
    --short-nlim) SHORT_NLIM="$2"; shift 2 ;;
    --long-tlim) LONG_TLIM="$2"; shift 2 ;;
    --long-nlim) LONG_NLIM="$2"; shift 2 ;;
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
  "parthenon/output0/dt=-1"
  "parthenon/output1/dt=${HISTORY_DT}"
  "diffusion/resistivity=ohmic"
  "diffusion/resistivity_coeff=fixed"
  "diffusion/integrator=rkl2"
  "diffusion/rkl2_max_dt_ratio=3.0"
  "diffusion/ohm_diff_coeff_code=1.788854e-02"
)

# Conservative + reservoir seed defaults.
COMMON_FULL_ARGS=(
  "modes4d/em_conservative_transport=true"
  "problem/harris_4d/drift_current_scale=8.0e-2"
  "problem/harris_4d/aw_mode2_amp=2.0e-2"
  "problem/harris_4d/a0_mode2_amp=1.0e-2"
  "problem/harris_4d/piw_mode2_amp=1.0e-2"
  # Keep direct bridge forcing off to isolate dynamic geometry response.
  "modes4d/em_source_mode0_even_bridge_gain=0.0"
  "modes4d/response_w0_drive_from_bridge_gain=0.0"
)

run_case() {
  local case_name="$1"
  shift
  local tlim="$1"
  shift
  local nlim="$1"
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
    --athena-arg "parthenon/time/tlim=${tlim}"
    --athena-arg "parthenon/time/nlim=${nlim}"
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

  echo "=== ${case_name} (tlim=${tlim}, nlim=${nlim}) ===" | tee "${log_file}"
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

# Baseline
run_case "baseline_fixed_center_t0p6" "${SHORT_TLIM}" "${SHORT_NLIM}"

# Short-horizon threshold sweep
run_case "dyn_g1_geom0p3_t0p6" "${SHORT_TLIM}" "${SHORT_NLIM}" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=1.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=0.30" \
  "modes4d/response_w0_geometry_shift_include_constant=true"

run_case "dyn_g3_geom1p0_t0p6" "${SHORT_TLIM}" "${SHORT_NLIM}" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=3.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=1.00" \
  "modes4d/response_w0_geometry_shift_include_constant=true"

run_case "dyn_g6_geom2p0_t0p6" "${SHORT_TLIM}" "${SHORT_NLIM}" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=6.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=2.00" \
  "modes4d/response_w0_geometry_shift_include_constant=true"

# Long-horizon escalation
run_case "dyn_g3_geom1p0_t1p2" "${LONG_TLIM}" "${LONG_NLIM}" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=3.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=1.00" \
  "modes4d/response_w0_geometry_shift_include_constant=true"

run_case "dyn_g6_geom2p0_t1p2" "${LONG_TLIM}" "${LONG_NLIM}" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=6.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=2.00" \
  "modes4d/response_w0_geometry_shift_include_constant=true"

# Standard post-processing
python3 "${REPO_ROOT}/scripts/summarize_lowgain_outputs.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

python3 "${REPO_ROOT}/scripts/analyze_windowed_residuals.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_windowed.csv" \
  --out-md "${OUTPUT_ROOT}_windowed.md"

python3 "${REPO_ROOT}/scripts/analyze_ledger_term_decomposition.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_decomp.csv" \
  --out-md "${OUTPUT_ROOT}_decomp.md"

echo ""
echo "Step-95 dynamic w0 threshold sweep complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "windowed_md,${OUTPUT_ROOT}_windowed.md"
echo "decomp_md,${OUTPUT_ROOT}_decomp.md"
