#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_ROOT="/projects/fluid-engine/out/step196_mixed_block_exact_probe_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step199_zblock_window_stability_local"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    -h|--help)
      cat <<'USAGE'
Post-fix analysis-only stability probe for the z-block energy/work partition across early/mid/late windows.
USAGE
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

mkdir -p "${OUTPUT_ROOT}"
python3 "${REPO_ROOT}/scripts/summarize_step199_zblock_window_stability.py" \
  --glob "${INPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

printf 'output_root,%s\n' "${OUTPUT_ROOT}"
printf 'summary_md,%s\n' "${OUTPUT_ROOT}_summary.md"
