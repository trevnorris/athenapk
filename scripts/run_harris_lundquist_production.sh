#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-baseline/bin/athenaPK"
WORKDIR="${REPO_ROOT}/build-baseline/modes4d_harris_lundquist_production"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full.in"
OUTPUT_DIR="${WORKDIR}/outputs"
S_VALUES="250,500,1000,2000"
TLIM="0.20"
NLIM="2000"
OUTPUT_DT="0.02"
CLOSURE_LOCAL_MODE0_ABS_RATE_TOL="3.0e-8"
FULL_MIN_JW_EW_ABS="1.0e-14"
FULL_MIN_S_LEAK_ABS="1.0e-13"
FULL_MIN_MIXED_EW2="1.0e-2"
FULL_MIN_JW_MODE_ACTIVITY_PROXY="1.0e-14"
FULL_MIN_ABS_EM_LEAK_W="1.0e-4"
FULL_MIN_ABS_HELICITY_SUB="1.0e-11"
FULL_MIN_ABS_EDOTB_SUB="1.0e-10"
PSI0_ENHANCEMENT_MIN="1.0"
MAX_ABS_DPSI0_ENHANCEMENT_MIN="1.0"
CONTROLLED_MAX_JW_EW_ABS="1.0e-12"
CONTROLLED_MAX_S_LEAK_ABS="1.0e-12"
CONTROLLED_MAX_JW_MODE_ACTIVITY_PROXY="1.0e-12"
SCAN_MIN_CORR_SPAN_PSI0="1.0e-5"
SCAN_MIN_CORR_SPAN_JW_EW_ABS="1.0e-15"
SCAN_MIN_CORR_SPAN_S_LEAK_ABS="5.0e-15"
SCAN_MIN_CORR_SPAN_JW_MODE_ACTIVITY_PROXY="1.0e-15"
PLOT_OUTPUT_DIR=""
SKIP_PLOTS=0
REPORT_PATH=""
SKIP_REPORT=0

usage() {
  cat <<'EOF'
Run a longer-horizon controlled/full Harris Lundquist scan with strict gates.

Usage:
  ./scripts/run_harris_lundquist_production.sh [options]

Options:
  --binary PATH            athenaPK binary (default: ./build-baseline/bin/athenaPK)
  --workdir DIR            Working directory (default: ./build-baseline/modes4d_harris_lundquist_production)
  --controlled-input PATH  Controlled deck (default: ./inputs/harris_4d_controlled.in)
  --full-input PATH        Full deck (default: ./inputs/harris_4d_full.in)
  --output-dir DIR         Output directory (default: <workdir>/outputs)
  --s-values CSV           Lundquist values (default: 250,500,1000,2000)
  --tlim FLOAT             Runtime limit per run (default: 0.20)
  --nlim INT               Cycle cap per run (default: 2000)
  --output-dt FLOAT        History-stream dt override (default: 0.02)
  --closure-local-mode0-abs-rate-tol FLOAT
                          Local mode-0 closure abs-rate tolerance forwarded to
                          harris_scan_matrix.py (default: 3.0e-8)
  --full-min-jw-ew-abs FLOAT
                          Full-case minimum |J^wE_w| activity floor (default: 1.0e-14)
  --full-min-s-leak-abs FLOAT
                          Full-case minimum S_leak_abs activity floor (default: 1.0e-13)
  --full-min-mixed-ew2 FLOAT
                          Full-case minimum mixed_ew2 activity floor (default: 1.0e-2)
  --full-min-jw-mode-activity-proxy FLOAT
                          Full-case minimum max(jw_mode_l2_1, |J^wE_w|) floor (default: 1.0e-14)
  --full-min-abs-em-leak-w FLOAT
                          Full-case minimum |EM leak_w| activity floor (default: 1.0e-4)
  --full-min-abs-helicity-sub FLOAT
                          Full-case minimum |helicity_sub| activity floor (default: 1.0e-11)
  --full-min-abs-edotb-sub FLOAT
                          Full-case minimum |edotb_sub| activity floor (default: 1.0e-10)
  --psi0-enhancement-min FLOAT
                          Minimum min(full/control psi0 span ratio) report gate (default: 1.0)
  --max-abs-dpsi0-enhancement-min FLOAT
                          Minimum min(full/control max|dpsi0/dt| ratio) report gate (default: 1.0)
  --controlled-max-jw-ew-abs FLOAT
                          Maximum max controlled |J^wE_w| report gate (default: 1.0e-12)
  --controlled-max-s-leak-abs FLOAT
                          Maximum max controlled S_leak_abs report gate (default: 1.0e-12)
  --controlled-max-jw-mode-activity-proxy FLOAT
                          Maximum max controlled jw_mode_activity_proxy report gate (default: 1.0e-12)
  --scan-min-corr-span-psi0 FLOAT
                          Minimum span for corr_logS_full_final_psi0_span signal (default: 1.0e-5)
  --scan-min-corr-span-jw-ew-abs FLOAT
                          Minimum span for corr_logS_full_final_jw_ew_abs signal (default: 1.0e-15)
  --scan-min-corr-span-s-leak-abs FLOAT
                          Minimum span for corr_logS_full_final_s_leak_abs signal (default: 5.0e-15)
  --scan-min-corr-span-jw-mode-activity-proxy FLOAT
                          Minimum span for corr_logS_full_final_jw_mode_activity_proxy (default: 1.0e-15)
  --plot-output-dir DIR   Plot output directory (default: <output-dir>/plots)
  --skip-plots            Skip post-scan plotting
  --report-path PATH      Markdown report path (default: <output-dir>/production_report.md)
  --skip-report           Skip markdown report generation
  -h, --help               Show this help
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
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --s-values)
      S_VALUES="$2"
      shift 2
      ;;
    --tlim)
      TLIM="$2"
      shift 2
      ;;
    --nlim)
      NLIM="$2"
      shift 2
      ;;
    --output-dt)
      OUTPUT_DT="$2"
      shift 2
      ;;
    --closure-local-mode0-abs-rate-tol)
      CLOSURE_LOCAL_MODE0_ABS_RATE_TOL="$2"
      shift 2
      ;;
    --full-min-jw-ew-abs)
      FULL_MIN_JW_EW_ABS="$2"
      shift 2
      ;;
    --full-min-s-leak-abs)
      FULL_MIN_S_LEAK_ABS="$2"
      shift 2
      ;;
    --full-min-mixed-ew2)
      FULL_MIN_MIXED_EW2="$2"
      shift 2
      ;;
    --full-min-jw-mode-activity-proxy)
      FULL_MIN_JW_MODE_ACTIVITY_PROXY="$2"
      shift 2
      ;;
    --full-min-abs-em-leak-w)
      FULL_MIN_ABS_EM_LEAK_W="$2"
      shift 2
      ;;
    --full-min-abs-helicity-sub)
      FULL_MIN_ABS_HELICITY_SUB="$2"
      shift 2
      ;;
    --full-min-abs-edotb-sub)
      FULL_MIN_ABS_EDOTB_SUB="$2"
      shift 2
      ;;
    --psi0-enhancement-min)
      PSI0_ENHANCEMENT_MIN="$2"
      shift 2
      ;;
    --max-abs-dpsi0-enhancement-min)
      MAX_ABS_DPSI0_ENHANCEMENT_MIN="$2"
      shift 2
      ;;
    --controlled-max-jw-ew-abs)
      CONTROLLED_MAX_JW_EW_ABS="$2"
      shift 2
      ;;
    --controlled-max-s-leak-abs)
      CONTROLLED_MAX_S_LEAK_ABS="$2"
      shift 2
      ;;
    --controlled-max-jw-mode-activity-proxy)
      CONTROLLED_MAX_JW_MODE_ACTIVITY_PROXY="$2"
      shift 2
      ;;
    --scan-min-corr-span-psi0)
      SCAN_MIN_CORR_SPAN_PSI0="$2"
      shift 2
      ;;
    --scan-min-corr-span-jw-ew-abs)
      SCAN_MIN_CORR_SPAN_JW_EW_ABS="$2"
      shift 2
      ;;
    --scan-min-corr-span-s-leak-abs)
      SCAN_MIN_CORR_SPAN_S_LEAK_ABS="$2"
      shift 2
      ;;
    --scan-min-corr-span-jw-mode-activity-proxy)
      SCAN_MIN_CORR_SPAN_JW_MODE_ACTIVITY_PROXY="$2"
      shift 2
      ;;
    --plot-output-dir)
      PLOT_OUTPUT_DIR="$2"
      shift 2
      ;;
    --skip-plots)
      SKIP_PLOTS=1
      shift
      ;;
    --report-path)
      REPORT_PATH="$2"
      shift 2
      ;;
    --skip-report)
      SKIP_REPORT=1
      shift
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

mkdir -p "${WORKDIR}" "${OUTPUT_DIR}"

python3 "${REPO_ROOT}/scripts/harris_lundquist_scan.py" \
  --binary "${BINARY}" \
  --workdir "${WORKDIR}" \
  --controlled-input "${CONTROLLED_INPUT}" \
  --full-input "${FULL_INPUT}" \
  --output-dir "${OUTPUT_DIR}" \
  --s-values "${S_VALUES}" \
  --fail-on-check \
  --check-transport-closure \
  --fail-on-scan-check \
  --scan-arg=--athena-arg \
  --scan-arg="parthenon/time/tlim=${TLIM}" \
  --scan-arg=--athena-arg \
  --scan-arg="parthenon/time/nlim=${NLIM}" \
  --scan-arg=--athena-arg \
  --scan-arg="parthenon/output0/dt=-1" \
  --scan-arg=--athena-arg \
  --scan-arg="parthenon/output1/dt=${OUTPUT_DT}" \
  --scan-arg=--closure-local-mode0-abs-rate-tol \
  --scan-arg="${CLOSURE_LOCAL_MODE0_ABS_RATE_TOL}" \
  --scan-arg=--full-min-jw-ew-abs \
  --scan-arg="${FULL_MIN_JW_EW_ABS}" \
  --scan-arg=--full-min-s-leak-abs \
  --scan-arg="${FULL_MIN_S_LEAK_ABS}" \
  --scan-arg=--full-min-mixed-ew2 \
  --scan-arg="${FULL_MIN_MIXED_EW2}" \
  --scan-arg=--full-min-jw-mode-activity-proxy \
  --scan-arg="${FULL_MIN_JW_MODE_ACTIVITY_PROXY}" \
  --scan-arg=--full-min-abs-em-leak-w \
  --scan-arg="${FULL_MIN_ABS_EM_LEAK_W}" \
  --scan-arg=--full-min-abs-helicity-sub \
  --scan-arg="${FULL_MIN_ABS_HELICITY_SUB}" \
  --scan-arg=--full-min-abs-edotb-sub \
  --scan-arg="${FULL_MIN_ABS_EDOTB_SUB}" \
  --scan-min-abs-corr corr_logS_full_final_psi0_span=7.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_s_leak_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_jw_ew_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_mixed_ew2=7.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_em_leak_w_abs=2.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_helicity_sub_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_edotb_sub_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_jw_mode_activity_proxy=5.0e-1 \
  --scan-min-metric min_full_final_jw_ew_abs="${FULL_MIN_JW_EW_ABS}" \
  --scan-min-metric min_full_final_s_leak_abs="${FULL_MIN_S_LEAK_ABS}" \
  --scan-min-metric min_full_final_mixed_ew2="${FULL_MIN_MIXED_EW2}" \
  --scan-min-metric min_full_final_jw_mode_activity_proxy="${FULL_MIN_JW_MODE_ACTIVITY_PROXY}" \
  --scan-min-metric min_full_final_em_leak_w_abs="${FULL_MIN_ABS_EM_LEAK_W}" \
  --scan-min-metric min_full_final_helicity_sub_abs="${FULL_MIN_ABS_HELICITY_SUB}" \
  --scan-min-metric min_full_final_edotb_sub_abs="${FULL_MIN_ABS_EDOTB_SUB}" \
  --scan-min-corr-span corr_logS_full_final_psi0_span="${SCAN_MIN_CORR_SPAN_PSI0}" \
  --scan-min-corr-span corr_logS_full_final_jw_ew_abs="${SCAN_MIN_CORR_SPAN_JW_EW_ABS}" \
  --scan-min-corr-span corr_logS_full_final_s_leak_abs="${SCAN_MIN_CORR_SPAN_S_LEAK_ABS}" \
  --scan-min-corr-span corr_logS_full_final_jw_mode_activity_proxy="${SCAN_MIN_CORR_SPAN_JW_MODE_ACTIVITY_PROXY}" \
  --scan-max-corr corr_logS_full_final_helicity_sub_abs=-2.0e-1 \
  --scan-max-corr corr_logS_full_final_edotb_sub_abs=-2.0e-1

SUMMARY_CSV="${OUTPUT_DIR}/lundquist_scan_summary.csv"
if [[ ! -f "${SUMMARY_CSV}" ]]; then
  printf 'Expected summary CSV missing: %s\n' "${SUMMARY_CSV}" >&2
  exit 1
fi

python3 - "${SUMMARY_CSV}" <<'PY'
import csv
import pathlib
import sys

summary_csv = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader(summary_csv.open("r", encoding="utf-8")))

print("production_summary_csv," + str(summary_csv))
print(
    "S,psi0_ctrl,psi0_full,psi0_ratio,max_abs_dpsi0_dt_ctrl,max_abs_dpsi0_dt_full,max_abs_dpsi0_dt_ratio,"
    "jw_ew_ctrl,jw_ew_full,s_leak_ctrl,s_leak_full,mixed_ew2_ctrl,mixed_ew2_full,mixed_c2_ctrl,mixed_c2_full,"
    "jw_mode_l2_1,jw_mode_activity_proxy,"
    "closure,transport,correlation,activity"
)
for row in rows:
    print(
        ",".join(
            [
                row.get("S", ""),
                row.get("controlled_final_psi0_span", ""),
                row.get("full_final_psi0_span", ""),
                row.get("ratio_full_over_controlled_psi0_span", ""),
                row.get("controlled_max_abs_dpsi0_dt", ""),
                row.get("full_max_abs_dpsi0_dt", ""),
                row.get("ratio_full_over_controlled_max_abs_dpsi0_dt", ""),
                row.get("controlled_final_jw_ew", ""),
                row.get("full_final_jw_ew", ""),
                row.get("controlled_final_s_leak_abs", ""),
                row.get("full_final_s_leak_abs", ""),
                row.get("controlled_final_mixed_ew2", ""),
                row.get("full_final_mixed_ew2", ""),
                row.get("controlled_final_mixed_c2", ""),
                row.get("full_final_mixed_c2", ""),
                row.get("full_final_jw_mode_l2_1", ""),
                row.get("full_final_jw_mode_activity_proxy", ""),
                row.get("full_closure_status", ""),
                row.get("full_transport_closure_status", ""),
                row.get("full_correlation_status", ""),
                row.get("full_activity_status", ""),
            ]
        )
    )
PY

if [[ "${SKIP_PLOTS}" -eq 0 ]]; then
  if [[ -z "${PLOT_OUTPUT_DIR}" ]]; then
    PLOT_OUTPUT_DIR="${OUTPUT_DIR}/plots"
  fi
  python3 "${REPO_ROOT}/scripts/plot_harris_lundquist_summary.py" \
    --summary-csv "${SUMMARY_CSV}" \
    --output-dir "${PLOT_OUTPUT_DIR}"
fi

if [[ "${SKIP_REPORT}" -eq 0 ]]; then
  if [[ -z "${PLOT_OUTPUT_DIR}" ]]; then
    PLOT_OUTPUT_DIR="${OUTPUT_DIR}/plots"
  fi
  if [[ -z "${REPORT_PATH}" ]]; then
    REPORT_PATH="${OUTPUT_DIR}/production_report.md"
  fi
  python3 "${REPO_ROOT}/scripts/generate_harris_lundquist_report.py" \
    --summary-csv "${SUMMARY_CSV}" \
    --plots-dir "${PLOT_OUTPUT_DIR}" \
    --report-path "${REPORT_PATH}" \
    --psi0-enhancement-min "${PSI0_ENHANCEMENT_MIN}" \
    --max-abs-dpsi0-enhancement-min "${MAX_ABS_DPSI0_ENHANCEMENT_MIN}" \
    --controlled-max-jw-ew-abs "${CONTROLLED_MAX_JW_EW_ABS}" \
    --controlled-max-s-leak-abs "${CONTROLLED_MAX_S_LEAK_ABS}" \
    --controlled-max-jw-mode-activity-proxy "${CONTROLLED_MAX_JW_MODE_ACTIVITY_PROXY}" \
    --full-min-jw-ew-abs "${FULL_MIN_JW_EW_ABS}" \
    --full-min-s-leak-abs "${FULL_MIN_S_LEAK_ABS}" \
    --full-min-jw-mode-activity-proxy "${FULL_MIN_JW_MODE_ACTIVITY_PROXY}"
fi
