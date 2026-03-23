#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SOURCE_ROOT="/projects/fluid-engine/out/step218_postpatch_proxy_calibration_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step219_expanded_proxy_search_local"
GATE_T_ON="3.0e-1"
GATE_T_OFF="9.0e-1"
SIGMA_STAR="2.07614e-01"
EPS="1.0e-30"

usage() {
  cat <<'USAGE'
Run Step-219 expanded proxy search on existing Step-218 outputs.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-root) SOURCE_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --gate-on) GATE_T_ON="$2"; shift 2 ;;
    --gate-off) GATE_T_OFF="$2"; shift 2 ;;
    --sigma-star) SIGMA_STAR="$2"; shift 2 ;;
    --eps) EPS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for path in \
  "${SOURCE_ROOT}" \
  "${SOURCE_ROOT}_summary.csv" \
  "${SOURCE_ROOT}_phase2.csv"; do
  if [[ ! -e "${path}" ]]; then
    echo "Missing required Step-218 artifact: ${path}" >&2
    exit 1
  fi
done

mkdir -p "$(dirname "${OUTPUT_ROOT}")"

python3 "${REPO_ROOT}/scripts/summarize_step219_expanded_proxy_search.py" \
  --glob "${SOURCE_ROOT}/*" \
  --summary-csv "${SOURCE_ROOT}_summary.csv" \
  --phase2-csv "${SOURCE_ROOT}_phase2.csv" \
  --gate-on "${GATE_T_ON}" \
  --gate-off "${GATE_T_OFF}" \
  --sigma-star "${SIGMA_STAR}" \
  --eps "${EPS}" \
  --out-csv "${OUTPUT_ROOT}.csv" \
  --out-md "${OUTPUT_ROOT}.md"
