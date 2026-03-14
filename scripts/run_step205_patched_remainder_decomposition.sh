#!/usr/bin/env bash
set -euo pipefail

ROOT="/projects/fluid-engine/out/step204_zblock_ledger_patch_local"
OUT_PREFIX="/projects/fluid-engine/out/step205_patched_remainder_decomposition_local"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$SCRIPT_DIR/summarize_step205_patched_remainder_decomposition.py" \
  --glob "$ROOT/*" \
  --out-csv "${OUT_PREFIX}_summary.csv" \
  --out-md "${OUT_PREFIX}_summary.md"

echo "summary_csv,${OUT_PREFIX}_summary.csv"
echo "summary_md,${OUT_PREFIX}_summary.md"
