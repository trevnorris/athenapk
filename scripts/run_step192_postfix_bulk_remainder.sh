#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_ROOT="/projects/fluid-engine/out/step190_cx_cz_exact_delta_probe_local_rerun"
OUTPUT_ROOT="/projects/fluid-engine/out/step192_postfix_bulk_remainder_local"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--input-root DIR] [--output-root PREFIX]"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done
mkdir -p "$(dirname "${OUTPUT_ROOT}")"
python3 "${REPO_ROOT}/scripts/summarize_step192_postfix_bulk_remainder.py" \
  --glob "${INPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
