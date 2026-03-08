#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step157_quadrature_ledger_decomposition_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step158_quadrature_ledger_terms_local"
FOCUS_WINDOW="late"
TOP_N=10

usage() {
  cat <<'USAGE'
Run Step-158 conservative-EM ledger term decomposition on existing Step-157 histories.

Goal:
  Rank the signed EM history-rate terms that best explain the remaining
  quadrature-baseline EM bulk ledger residual after transport has been fixed.

Usage:
  ./scripts/run_step158_quadrature_ledger_terms.sh [options]

Options:
  --input-root DIR       Step-157 output root (default:
                         /projects/fluid-engine/out/step157_quadrature_ledger_decomposition_local)
  --output-root DIR      Output prefix root (default:
                         /projects/fluid-engine/out/step158_quadrature_ledger_terms_local)
  --focus-window NAME    Window to summarize: early, mid, late, all
                         (default: late)
  --top-n N              Number of ranked terms per case/window (default: 10)
  -h, --help             Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --focus-window) FOCUS_WINDOW="$2"; shift 2 ;;
    --top-n) TOP_N="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ ! -d "${INPUT_ROOT}" ]]; then
  echo "Input root missing: ${INPUT_ROOT}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUTPUT_ROOT}")"

python3 "${REPO_ROOT}/scripts/summarize_step158_quadrature_ledger_terms.py" \
  --glob "${INPUT_ROOT}/*" \
  --focus-window "${FOCUS_WINDOW}" \
  --top-n "${TOP_N}" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "summary_csv,${OUTPUT_ROOT}_summary.csv"
