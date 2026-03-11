#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step183_aw_update_stage_probe_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step186_z_normalization_probe_local"
BASELINE_CASE="baseline_quadrature"
FOCUS_WINDOW="late"
LAMBDA_VALUE="1.0"

usage() {
  cat <<'USAGE'
Run Step-186 z-channel normalization probe on existing Step-183 histories.

Goal:
  Test whether the surviving bulk residual is primarily a sqrt(pi)*lambda scaling
  mismatch on the z-directed full discrete A_w update-power channel.

Usage:
  ./scripts/run_step186_z_normalization_probe.sh [options]

Options:
  --input-root DIR       Existing output root containing baseline_quadrature
                         (default: /projects/fluid-engine/out/step183_aw_update_stage_probe_local)
  --output-root DIR      Output prefix root
                         (default: /projects/fluid-engine/out/step186_z_normalization_probe_local)
  --baseline-case NAME   Case name to analyze (default: baseline_quadrature)
  --focus-window NAME    Window to summarize: early, mid, late, all
                         (default: late)
  --lambda VALUE         Lambda used in the baseline (default: 1.0)
  -h, --help             Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --baseline-case) BASELINE_CASE="$2"; shift 2 ;;
    --focus-window) FOCUS_WINDOW="$2"; shift 2 ;;
    --lambda) LAMBDA_VALUE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ ! -d "${INPUT_ROOT}" ]]; then
  echo "Input root missing: ${INPUT_ROOT}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUTPUT_ROOT}")"

python3 "${REPO_ROOT}/scripts/summarize_step186_z_normalization_probe.py" \
  --glob "${INPUT_ROOT}/*" \
  --baseline-case "${BASELINE_CASE}" \
  --focus-window "${FOCUS_WINDOW}" \
  --lambda "${LAMBDA_VALUE}" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "summary_csv,${OUTPUT_ROOT}_summary.csv"
