#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/weekend_geomshift_$(date +%Y%m%d_%H%M%S)"

NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
HISTORY_DT="0.002"

usage() {
  cat <<'EOF'
Run an intensive weekend geometry-response campaign for 4D spillback tests.

Usage:
  ./scripts/run_weekend_geomshift_campaign.sh [options]

Options:
  --binary PATH            Path to athenaPK binary
                           (default: ./build-gpu-cuda124-ada89/bin/athenaPK)
  --workdir DIR            Working directory for runs (default: repo root)
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
  -h, --help               Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary)
      BINARY="$2"
      shift 2
      ;;
    --workdir)
      WORKDIR="$2"
      shift 2
      ;;
    --controlled-input)
      CONTROLLED_INPUT="$2"
      shift 2
      ;;
    --full-input)
      FULL_INPUT="$2"
      shift 2
      ;;
    --output-root)
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    --nx1)
      NX1="$2"
      shift 2
      ;;
    --nx2)
      NX2="$2"
      shift 2
      ;;
    --nx3)
      NX3="$2"
      shift 2
      ;;
    --mb1)
      MB1="$2"
      shift 2
      ;;
    --mb2)
      MB2="$2"
      shift 2
      ;;
    --mb3)
      MB3="$2"
      shift 2
      ;;
    --history-dt)
      HISTORY_DT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n\n' "$1" >&2
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
    --fail-on-nonfinite-hst
    --athena-arg "parthenon/time/tlim=${tlim}"
    --athena-arg "parthenon/time/nlim=${nlim}"
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

# Phase A: screening runs (moderate horizon)
run_case "screen_baseline_t0p8" "0.80" "80000"
run_case "screen_geom_g1_t0p8" "0.80" "80000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.05" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.20" \
  "modes4d/response_w0_drive_gain=0.5" \
  "modes4d/response_w0_drive_from_bridge_gain=50.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=1.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true"
run_case "screen_geom_g3_t0p8" "0.80" "80000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.05" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.20" \
  "modes4d/response_w0_drive_gain=1.5" \
  "modes4d/response_w0_drive_from_bridge_gain=150.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=3.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true"
run_case "screen_geom_g6_t0p8" "0.80" "80000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=3.0" \
  "modes4d/response_w0_drive_from_bridge_gain=300.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=6.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true"
run_case "screen_geom_g6_noconst_t0p8" "0.80" "80000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=3.0" \
  "modes4d/response_w0_drive_from_bridge_gain=300.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=6.0" \
  "modes4d/response_w0_geometry_shift_include_constant=false"
run_case "screen_bridge2000_geom8_t0p8" "0.80" "80000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.08" \
  "modes4d/response_w0_stiffness=0.03" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.30" \
  "modes4d/response_w0_drive_gain=0.0" \
  "modes4d/response_w0_drive_from_bridge_gain=2000.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=8.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true"
run_case "screen_bridge2000_geom8_proj_t0p8" "0.80" "80000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.08" \
  "modes4d/response_w0_stiffness=0.03" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.30" \
  "modes4d/response_w0_drive_gain=0.0" \
  "modes4d/response_w0_drive_from_bridge_gain=2000.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=8.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true" \
  "modes4d/response_w0_projection_enable=true" \
  "modes4d/response_w0_projection_gain=1.0" \
  "modes4d/response_w0_projection_max_abs=0.30"

# Phase B: long-horizon follow-up on strongest settings
run_case "long_baseline_t2p0" "2.00" "220000"
run_case "long_geom_g6_t2p0" "2.00" "220000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.10" \
  "modes4d/response_w0_stiffness=0.04" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.25" \
  "modes4d/response_w0_drive_gain=3.0" \
  "modes4d/response_w0_drive_from_bridge_gain=300.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=6.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true"
run_case "long_bridge2000_geom8_t2p0" "2.00" "220000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.08" \
  "modes4d/response_w0_stiffness=0.03" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.30" \
  "modes4d/response_w0_drive_gain=0.0" \
  "modes4d/response_w0_drive_from_bridge_gain=2000.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=8.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true"
run_case "long_bridge2000_geom8_proj_t2p0" "2.00" "220000" \
  "modes4d/response_w0_enable=true" \
  "modes4d/response_w0_mass=0.08" \
  "modes4d/response_w0_stiffness=0.03" \
  "modes4d/response_w0_damping=0.02" \
  "modes4d/response_w0_max_abs=0.30" \
  "modes4d/response_w0_drive_gain=0.0" \
  "modes4d/response_w0_drive_from_bridge_gain=2000.0" \
  "modes4d/response_w0_drive_from_parity_gain=0.0" \
  "modes4d/response_w0_geometry_shift_enable=true" \
  "modes4d/response_w0_geometry_shift_gain=8.0" \
  "modes4d/response_w0_geometry_shift_include_constant=true" \
  "modes4d/response_w0_projection_enable=true" \
  "modes4d/response_w0_projection_gain=1.0" \
  "modes4d/response_w0_projection_max_abs=0.30"

SUMMARY_CSV="${OUTPUT_ROOT}/campaign_summary.csv"
SUMMARY_MD="${OUTPUT_ROOT}/campaign_summary.md"

python3 - "${OUTPUT_ROOT}" "${REPO_ROOT}" "${SUMMARY_CSV}" "${SUMMARY_MD}" <<'PY'
import csv
import math
import sys
from pathlib import Path

out_root = Path(sys.argv[1])
repo_root = Path(sys.argv[2])
summary_csv = Path(sys.argv[3])
summary_md = Path(sys.argv[4])

sys.path.insert(0, str(repo_root / "scripts"))
import harris_scan_matrix as hsm  # type: ignore

analysis_args = (5.0e-2, 1.0e-8, 1.0e-8, 1.0, 1.0, 1.0e-6, 1.0e-12, 1.0e-3, 1.0e-12)

metrics = [
    "final_psi0_span",
    "final_psi_proj_span",
    "final_psi_w0_span",
    "final_psi_w0_minus_psi_proj_span",
    "final_response_w0",
    "final_response_w0_geometry_center",
    "final_response_w0_drive",
    "final_int_src_em_geometry_shift_abs",
    "final_int_geometry_shift_power_mode0",
    "final_int_response_w0_power_net",
    "final_int_response_bridge_power_sum",
    "final_int_src_mode0_total_abs",
    "final_emf_vw_c_abs",
]

status_metrics = [
    "closure_status",
    "closure_local_mode0_status",
    "transport_closure_status",
    "em_bulk_ledger_status",
]

def ratio(full_v, ctrl_v):
    if (
        isinstance(full_v, (int, float))
        and isinstance(ctrl_v, (int, float))
        and math.isfinite(float(full_v))
        and math.isfinite(float(ctrl_v))
        and abs(float(ctrl_v)) > 0.0
    ):
        return float(full_v) / float(ctrl_v)
    return math.nan

rows = []
for case_dir in sorted(p for p in out_root.iterdir() if p.is_dir()):
    row = {"case": case_dir.name}
    rc_path = case_dir / "exit_code.txt"
    row["exit_code"] = int(rc_path.read_text().strip()) if rc_path.exists() else 999
    c_hst = case_dir / "harris_controlled.out1.hst"
    f_hst = case_dir / "harris_full.out1.hst"
    if c_hst.exists() and f_hst.exists():
        try:
            c_cols = hsm.parse_hst(c_hst)
            f_cols = hsm.parse_hst(f_hst)
            c = hsm.analyze_case("controlled", c_cols, *analysis_args)
            f = hsm.analyze_case("full", f_cols, *analysis_args)
            c_nonfinite = hsm.detect_nonfinite_columns(c_cols)
            f_nonfinite = hsm.detect_nonfinite_columns(f_cols)
            row["history_rows_controlled"] = len(c_cols.get("time", []))
            row["history_rows_full"] = len(f_cols.get("time", []))
            row["final_time_controlled"] = c_cols.get("time", [math.nan])[-1]
            row["final_time_full"] = f_cols.get("time", [math.nan])[-1]
            row["controlled_history_nonfinite_status"] = (
                "FAIL" if c_nonfinite["has_nonfinite"] else "PASS"
            )
            row["full_history_nonfinite_status"] = (
                "FAIL" if f_nonfinite["has_nonfinite"] else "PASS"
            )
            row["controlled_history_nonfinite_count"] = c_nonfinite["total_count"]
            row["full_history_nonfinite_count"] = f_nonfinite["total_count"]
            row["controlled_history_nonfinite_columns"] = c_nonfinite["columns_csv"]
            row["full_history_nonfinite_columns"] = f_nonfinite["columns_csv"]
            for m in metrics:
                cv = c.get(m, math.nan)
                fv = f.get(m, math.nan)
                row[f"controlled_{m}"] = cv
                row[f"full_{m}"] = fv
                row[f"ratio_{m}"] = ratio(fv, cv)
            for m in status_metrics:
                row[f"controlled_{m}"] = c.get(m, "N/A")
                row[f"full_{m}"] = f.get(m, "N/A")
        except Exception as exc:
            row["parse_error"] = str(exc)
    rows.append(row)

fieldnames = sorted({k for r in rows for k in r.keys()})
with summary_csv.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for r in rows:
        writer.writerow(r)

def ffmt(v):
    if isinstance(v, str):
        return v
    if not isinstance(v, (int, float)) or not math.isfinite(float(v)):
        return "nan"
    return f"{float(v):.6e}"

lines = []
lines.append("# Weekend Geometry-Shift Campaign Summary")
lines.append("")
lines.append(f"- output root: `{out_root}`")
lines.append(f"- summary csv: `{summary_csv}`")
lines.append("")
lines.append("| case | rc | ctrl nonfinite | full nonfinite | psi0 ratio | psi_proj ratio | psi_w0 ratio | w0-psiproj(full) | response_w0(full) | geom_shift_src_abs(full) | response_power_net(full) | closure_local(full) | ledger(full) |")
lines.append("|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|")
for r in rows:
    lines.append(
        "| {case} | {rc} | {cnf} | {fnf} | {psi0} | {psip} | {psiw} | {dw} | {rw0} | {gsrc} | {rpow} | {cloc} | {led} |".format(
            case=r.get("case", ""),
            rc=r.get("exit_code", 999),
            cnf=r.get("controlled_history_nonfinite_status", "N/A"),
            fnf=r.get("full_history_nonfinite_status", "N/A"),
            psi0=ffmt(r.get("ratio_final_psi0_span", math.nan)),
            psip=ffmt(r.get("ratio_final_psi_proj_span", math.nan)),
            psiw=ffmt(r.get("ratio_final_psi_w0_span", math.nan)),
            dw=ffmt(r.get("full_final_psi_w0_minus_psi_proj_span", math.nan)),
            rw0=ffmt(r.get("full_final_response_w0", math.nan)),
            gsrc=ffmt(r.get("full_final_int_src_em_geometry_shift_abs", math.nan)),
            rpow=ffmt(r.get("full_final_int_response_w0_power_net", math.nan)),
            cloc=r.get("full_closure_local_mode0_status", "N/A"),
            led=r.get("full_em_bulk_ledger_status", "N/A"),
        )
    )

summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"campaign_summary_csv,{summary_csv}")
print(f"campaign_summary_md,{summary_md}")
PY

echo ""
echo "Weekend campaign complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_csv,${SUMMARY_CSV}"
echo "summary_md,${SUMMARY_MD}"
