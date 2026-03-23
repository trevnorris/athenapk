#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step218_postpatch_proxy_calibration_local"

NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
HISTORY_DT="0.001"
TLIM="1.20"
NLIM="120000"
RESUME="true"
ZBLOCK_COEFF="8.304215e-01"
MODE1_OMEGA="6.283185307179586"
MODE2_OMEGA="12.566370614359172"
GATE_T_ON="3.0e-1"
GATE_T_OFF="9.0e-1"
GATE_PERIOD="1.5e-1"
GATE_DUTY="4.0e-1"
GATE_RAMP="1.5e-2"
MODE1_AMP_FORCED="1.25e-4"
MODE2_AMP_FORCED="7.5e-5"

usage() {
  cat <<'USAGE'
Run Step-218 post-patch proxy calibration A/B.
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
    --resume) RESUME="$2"; shift 2 ;;
    --zblock-coeff) ZBLOCK_COEFF="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
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

CONTROLLED_STRICT_BASE_ARGS=(
  "modes4d/em_conservative_transport=true"
  "modes4d/hard_controlled_limit_enable=true"
  "modes4d/hard_controlled_limit_identity_mode0=true"
  "modes4d/hard_controlled_limit_strict_solver_path=true"
  "modes4d/n_quadrature=1"
  "modes4d/response_w0_enable=false"
  "modes4d/response_lambda_enable=false"
  "modes4d/response_w0_dynamic_mixing_enable=false"
  "modes4d/response_w0_local_enable=false"
  "modes4d/response_w0_local_gradient_mixing_enable=false"
  "modes4d/response_w0_geometry_shift_enable=false"
  "modes4d/continuity_consistent_current_projection_enable=false"
)

FULL_BASE_ARGS=(
  "modes4d/lambda=1.0"
  "modes4d/em_conservative_transport=true"
  "modes4d/continuity_consistent_current_projection_enable=true"
  "modes4d/bulk_energy_zblock_correction_coeff=${ZBLOCK_COEFF}"
  "problem/harris_4d/drift_current_scale=8.0e-2"
  "problem/harris_4d/drift_momw_scale=1.5e-1"
  "problem/harris_4d/aw_mode1_amp=6.0e-1"
  "problem/harris_4d/piw_mode1_amp=3.0e-1"
  "problem/harris_4d/ay_mode1_amp=3.0e-1"
  "problem/harris_4d/piy_mode1_amp=1.5e-1"
  "problem/harris_4d/aw_mode2_amp=4.0e-2"
  "problem/harris_4d/a0_mode2_amp=2.0e-2"
  "problem/harris_4d/piw_mode2_amp=2.0e-2"
  "problem/harris_4d/aw_mode1_phase_x=3.141593e-1"
  "problem/harris_4d/aw_mode1_phase_z=1.099558e+0"
  "problem/harris_4d/piw_mode1_phase_x=2.513274e-1"
  "problem/harris_4d/piw_mode1_phase_z=8.482300e-1"
  "problem/harris_4d/ay_mode1_phase_x=6.283185e-1"
  "problem/harris_4d/ay_mode1_phase_z=4.398230e-1"
  "problem/harris_4d/piy_mode1_phase_x=5.654867e-1"
  "problem/harris_4d/piy_mode1_phase_z=3.141593e-1"
  "modes4d/em_source_mass_gain=1.0"
  "modes4d/em_source_laplacian_gain=1.0"
  "modes4d/em_source_current_gain=14.0"
  "modes4d/em_source_timelike_gain=14.0"
  "modes4d/plasma_force_source_gain=1.0"
  "modes4d/plasma_momw_source_gain=14.0"
  "modes4d/plasma_momw_pressure_source_gain=14.0"
  "modes4d/plasma_energy_source_gain=1.0"
  "modes4d/plasma_w_flux_source_gain=0.0"
  "modes4d/plasma_rho_divj_gain=1.0"
  "modes4d/plasma_pressure_transport_gain=1.0"
  "modes4d/plasma_pressure_rusanov_gain=1.0"
  "modes4d/em_source_damping_gain=2.0e-1"
  "modes4d/em_source_gauge_gain=0.0"
  "modes4d/response_lambda_enable=false"
  "modes4d/response_w0_enable=true"
  "modes4d/response_w0_dynamic_mixing_enable=true"
  "modes4d/response_w0_dynamic_mixing_gain=3.5"
  "modes4d/response_w0_local_enable=true"
  "modes4d/response_w0_local_mode1_kx=3.141592653589793"
  "modes4d/response_w0_local_mode1_kz=6.283185307179586"
  "modes4d/response_w0_local_mode1_omega=${MODE1_OMEGA}"
  "modes4d/response_w0_local_mode2_kx=6.283185307179586"
  "modes4d/response_w0_local_mode2_kz=12.566370614359172"
  "modes4d/response_w0_local_mode2_phase=1.0471975511965976"
  "modes4d/response_w0_local_mode2_omega=${MODE2_OMEGA}"
  "modes4d/response_w0_local_gradient_mixing_enable=true"
  "modes4d/response_w0_local_gradient_mixing_gain=7.0"
  "modes4d/response_w0_local_gradient_norm_form=quadrature"
  "modes4d/response_w0_geometry_shift_enable=true"
  "modes4d/response_w0_projection_enable=false"
  "modes4d/response_w0_projection_gain=1.0"
  "modes4d/response_w0_projection_max_abs=0.0"
  "modes4d/response_w0_local_gate_t_on=${GATE_T_ON}"
  "modes4d/response_w0_local_gate_t_off=${GATE_T_OFF}"
)

run_case() {
  local case_name="$1"
  local mode1_amp="$2"
  local mode2_amp="$3"
  local gate_period="$4"
  local gate_duty="$5"
  local gate_ramp="$6"
  local case_dir="${OUTPUT_ROOT}/${case_name}"
  local log_file="${case_dir}/run.log"
  local exit_file="${case_dir}/exit_code.txt"
  mkdir -p "${case_dir}"

  if [[ "${RESUME}" == "true" ]] && [[ -f "${exit_file}" ]]; then
    local prior_rc
    prior_rc="$(<"${exit_file}")"
    if [[ "${prior_rc}" == "0" ]]; then
      echo "=== ${case_name} ===" | tee -a "${log_file}"
      echo "case_status,${case_name},SKIP,resume_hit" | tee -a "${log_file}"
      return 0
    else
      rm -f "${exit_file}"
    fi
  fi

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
  for arg in "${COMMON_MESH_ARGS[@]}"; do cmd+=(--athena-arg "${arg}"); done
  for arg in "${CONTROLLED_STRICT_BASE_ARGS[@]}"; do cmd+=(--controlled-arg "${arg}"); done
  for arg in "${FULL_BASE_ARGS[@]}"; do cmd+=(--full-arg "${arg}"); done
  cmd+=(--full-arg "modes4d/response_w0_local_mode1_amp=${mode1_amp}")
  cmd+=(--full-arg "modes4d/response_w0_local_mode2_amp=${mode2_amp}")
  cmd+=(--full-arg "modes4d/response_w0_local_gate_period=${gate_period}")
  cmd+=(--full-arg "modes4d/response_w0_local_gate_duty=${gate_duty}")
  cmd+=(--full-arg "modes4d/response_w0_local_gate_ramp=${gate_ramp}")

  echo "=== ${case_name} ===" | tee "${log_file}"
  echo "command,${cmd[*]}" >> "${log_file}"
  set +e
  "${cmd[@]}" 2>&1 | tee -a "${log_file}"
  local rc=${PIPESTATUS[0]}
  set -e
  echo "${rc}" > "${exit_file}"
  if [[ ${rc} -eq 0 ]]; then
    echo "case_status,${case_name},PASS" | tee -a "${log_file}"
  else
    echo "case_status,${case_name},FAIL,rc=${rc}" | tee -a "${log_file}"
    return ${rc}
  fi
}

run_case "patched_env_ref" "0.0" "0.0" "0.0" "1.0" "0.0"
run_case "patched_env_train_quad_amp0010" "${MODE1_AMP_FORCED}" "${MODE2_AMP_FORCED}" "${GATE_PERIOD}" "${GATE_DUTY}" "${GATE_RAMP}"

python3 "${REPO_ROOT}/scripts/summarize_transverse_channel_activation.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"
python3 "${REPO_ROOT}/scripts/summarize_projected_closure.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_projected_closure.csv" \
  --out-md "${OUTPUT_ROOT}_projected_closure.md"
python3 "${REPO_ROOT}/scripts/summarize_step121_phase2_predictive.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_phase2.csv" \
  --out-md "${OUTPUT_ROOT}_phase2.md"
python3 "${REPO_ROOT}/scripts/summarize_step123_controlled_bounds.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_controlled_bounds.csv" \
  --out-md "${OUTPUT_ROOT}_controlled_bounds.md"
python3 "${REPO_ROOT}/scripts/summarize_step124_tradeoff.py" \
  --phase2-csv "${OUTPUT_ROOT}_phase2.csv" \
  --controlled-bounds-csv "${OUTPUT_ROOT}_controlled_bounds.csv" \
  --summary-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-csv "${OUTPUT_ROOT}_tradeoff.csv" \
  --out-md "${OUTPUT_ROOT}_tradeoff.md"
python3 "${REPO_ROOT}/scripts/summarize_step132_solver_equivalence.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_solver_equivalence.csv" \
  --out-md "${OUTPUT_ROOT}_solver_equivalence.md"
python3 "${REPO_ROOT}/scripts/summarize_step140_mode0_term_isolation.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_mode0_terms.csv" \
  --out-md "${OUTPUT_ROOT}_mode0_terms.md"
python3 "${REPO_ROOT}/scripts/summarize_step149_transport_components.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --summary-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-csv "${OUTPUT_ROOT}_transport_ledger.csv" \
  --out-md "${OUTPUT_ROOT}_transport_ledger.md"
python3 "${REPO_ROOT}/scripts/summarize_step218_postpatch_proxy_calibration.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --summary-csv "${OUTPUT_ROOT}_summary.csv" \
  --phase2-csv "${OUTPUT_ROOT}_phase2.csv" \
  --gate-on "${GATE_T_ON}" \
  --gate-off "${GATE_T_OFF}" \
  --out-csv "${OUTPUT_ROOT}_proxy_calibration.csv" \
  --out-md "${OUTPUT_ROOT}_proxy_calibration.md"
