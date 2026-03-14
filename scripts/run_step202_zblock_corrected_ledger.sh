#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step200_zblock_exact_node_work_probe_local"
OUTPUT_PREFIX="/projects/fluid-engine/out/step202_zblock_corrected_ledger_local"

usage() {
  cat <<'USAGE'
Run Step-202 corrected-ledger A/B analysis on an existing post-fix baseline output.

Goal:
  Compare the legacy EM bulk ledger against fixed-coefficient z-block and xz-block
  corrected ledgers.

Usage:
  ./scripts/run_step202_zblock_corrected_ledger.sh [options]

Options:
  --input-root PATH
  --output-prefix PATH
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

python3 "${REPO_ROOT}/scripts/summarize_step202_zblock_corrected_ledger.py" \
  --glob "${INPUT_ROOT}/*" \
  --out-csv "${OUTPUT_PREFIX}_summary.csv" \
  --out-md "${OUTPUT_PREFIX}_summary.md"

echo "summary_csv,${OUTPUT_PREFIX}_summary.csv"
echo "summary_md,${OUTPUT_PREFIX}_summary.md"
