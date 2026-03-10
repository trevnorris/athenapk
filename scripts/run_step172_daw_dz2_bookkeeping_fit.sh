#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step170_local_gradient_bounded_signed_ab_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step172_daw_dz2_bookkeeping_fit_local"
BASELINE_CASE="baseline_quadrature"
FOCUS_WINDOW="late"

usage() {
  cat <<'USAGE'
Run Step-172 diagonal (d_z a_w)^2 bookkeeping coefficient fits on existing quadrature-baseline histories.

Goal:
  Fit the surviving ledger residual against the diagonal d_z a_w energy channels to determine
  whether the remaining mismatch is consistent with a unit-coefficient bookkeeping omission,
  a scaled diagonal mismatch, or a multi-mode combination.

Usage:
  ./scripts/run_step172_daw_dz2_bookkeeping_fit.sh [options]

Options:
  --input-root DIR       Existing output root containing baseline_quadrature
                         (default: /projects/fluid-engine/out/step170_local_gradient_bounded_signed_ab_local)
  --output-root DIR      Output prefix root
                         (default: /projects/fluid-engine/out/step172_daw_dz2_bookkeeping_fit_local)
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

python3 "${REPO_ROOT}/scripts/summarize_step172_daw_dz2_bookkeeping_fit.py" \
  --glob "${INPUT_ROOT}/*" \
  --baseline-case "${BASELINE_CASE}" \
  --focus-window "${FOCUS_WINDOW}" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "summary_csv,${OUTPUT_ROOT}_summary.csv"
