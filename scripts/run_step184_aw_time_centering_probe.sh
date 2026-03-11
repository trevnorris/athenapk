#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step183_aw_update_stage_probe_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step184_aw_time_centering_probe_local"
BASELINE_CASE="baseline_quadrature"
FOCUS_WINDOW="late"

usage() {
  cat <<'USAGE'
Run Step-184 A_w time-centering probe on existing Step-183 histories.

Goal:
  Test whether the surviving quadrature-baseline bulk residual is best explained by
  a time-centered or stage-shifted version of the full discrete A_w update-power sum.

Usage:
  ./scripts/run_step184_aw_time_centering_probe.sh [options]

Options:
  --input-root DIR       Existing output root containing baseline_quadrature
                         (default: /projects/fluid-engine/out/step183_aw_update_stage_probe_local)
  --output-root DIR      Output prefix root
                         (default: /projects/fluid-engine/out/step184_aw_time_centering_probe_local)
  --baseline-case NAME   Case name to analyze (default: baseline_quadrature)
  --focus-window NAME    Window to summarize: early, mid, late, all
                         (default: late)
  -h, --help             Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --baseline-case) BASELINE_CASE="$2"; shift 2 ;;
    --focus-window) FOCUS_WINDOW="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ ! -d "${INPUT_ROOT}" ]]; then
  echo "Input root missing: ${INPUT_ROOT}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUTPUT_ROOT}")"

python3 "${REPO_ROOT}/scripts/summarize_step184_aw_time_centering_probe.py" \
  --glob "${INPUT_ROOT}/*" \
  --baseline-case "${BASELINE_CASE}" \
  --focus-window "${FOCUS_WINDOW}" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "summary_csv,${OUTPUT_ROOT}_summary.csv"
