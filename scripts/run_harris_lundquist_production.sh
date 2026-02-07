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
  --scan-min-abs-corr corr_logS_full_final_psi0_span=7.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_s_leak_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_jw_ew_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_mixed_ew2=7.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_em_leak_w_abs=2.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_helicity_sub_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_edotb_sub_abs=5.0e-1 \
  --scan-min-abs-corr corr_logS_full_final_jw_mode_activity_proxy=5.0e-1 \
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
    "S,psi0_span,jw_ew,s_leak_abs,mixed_ew2,mixed_c2,jw_mode_l2_1,jw_mode_activity_proxy,"
    "closure,transport,correlation,activity"
)
for row in rows:
    print(
        ",".join(
            [
                row.get("S", ""),
                row.get("full_final_psi0_span", ""),
                row.get("full_final_jw_ew", ""),
                row.get("full_final_s_leak_abs", ""),
                row.get("full_final_mixed_ew2", ""),
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
    --report-path "${REPORT_PATH}"
fi
