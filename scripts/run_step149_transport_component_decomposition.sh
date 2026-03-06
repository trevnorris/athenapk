#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

CASE_GLOB="/projects/fluid-engine/out/step148_ledger_decomposition_local/*"
OUTPUT_ROOT="/projects/fluid-engine/out/step149_transport_component_decomposition_local"
BASELINE_CASE="baseline_decomp"

usage() {
  cat <<'USAGE'
Analyze Step-148 transport closure by conserved component.

Goal:
  Recompute the mode-0 transport-balance residual for each conserved variable
  directly from the Step-148 full-case HST outputs, so we can see which
  component is actually dominating the remaining transport failure.

Usage:
  ./scripts/run_step149_transport_component_decomposition.sh [options]

Options:
  --glob PATTERN         Case directory glob (default: Step-148 output dirs)
  --output-root PATH     Output root prefix
  --baseline-case NAME   Baseline case name (default: baseline_decomp)
  -h, --help             Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --glob) CASE_GLOB="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --baseline-case) BASELINE_CASE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

python3 "${REPO_ROOT}/scripts/summarize_step149_transport_components.py" \
  --glob "${CASE_GLOB}" \
  --baseline-case "${BASELINE_CASE}" \
  --out-csv "${OUTPUT_ROOT}_components.csv" \
  --summary-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "components_csv,${OUTPUT_ROOT}_components.csv"
echo "summary_csv,${OUTPUT_ROOT}_summary.csv"
