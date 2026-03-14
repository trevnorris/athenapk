#!/usr/bin/env bash
set -euo pipefail

INPUT_ROOT="/projects/fluid-engine/out/step204_zblock_ledger_patch_local"
OUTPUT_PREFIX="/projects/fluid-engine/out/step206_postpatch_physics_reassessment_local"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
    -h|--help)
      cat <<'USAGE'
Run Step-206 post-patch physics reassessment on the existing Step-204 patched baseline.
USAGE
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

CASE_DIR="${INPUT_ROOT}/zblock_patch_candidate"

python3 "$SCRIPT_DIR/summarize_step206_postpatch_physics_reassessment.py" \
  --case-dir "$CASE_DIR" \
  --summary-csv "${INPUT_ROOT}_summary.csv" \
  --phase2-csv "${INPUT_ROOT}_phase2.csv" \
  --zblock-csv "${INPUT_ROOT}_zblock_patch_ab.csv" \
  --out-csv "${OUTPUT_PREFIX}_summary.csv" \
  --out-md "${OUTPUT_PREFIX}_summary.md"

echo "summary_csv,${OUTPUT_PREFIX}_summary.csv"
echo "summary_md,${OUTPUT_PREFIX}_summary.md"
