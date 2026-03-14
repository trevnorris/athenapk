#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step200_zblock_exact_node_work_probe_local"
OUTPUT_PREFIX="/projects/fluid-engine/out/step201_postfix_zblock_fixedrule_local"

usage() {
  cat <<'USAGE'
Run Step-201 fixed-coefficient z-block validation on an existing post-fix baseline output.

Goal:
  Test whether the stable post-fix z-block coefficients from Steps 198/199 generalize
  across early/mid/late windows without refitting.

Usage:
  ./scripts/run_step201_postfix_zblock_fixedrule.sh [options]

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

python3 "${REPO_ROOT}/scripts/summarize_step201_postfix_zblock_fixedrule.py" \
  --glob "${INPUT_ROOT}/*" \
  --out-csv "${OUTPUT_PREFIX}_summary.csv" \
  --out-md "${OUTPUT_PREFIX}_summary.md"

echo "summary_csv,${OUTPUT_PREFIX}_summary.csv"
echo "summary_md,${OUTPUT_PREFIX}_summary.md"
