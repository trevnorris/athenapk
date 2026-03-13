#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_ROOT="/projects/fluid-engine/out/step196_mixed_block_exact_probe_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step197_postfix_energy_work_split_local"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    -h|--help)
      cat <<'USAGE'
Post-fix analysis-only split of the remaining ledger residual into exact Cx/Cz energy and full discrete A_w update-work channels.
USAGE
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

mkdir -p "${OUTPUT_ROOT}"
python3 "${REPO_ROOT}/scripts/summarize_step197_postfix_energy_work_split.py" \
  --glob "${INPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

printf 'output_root,%s\n' "${OUTPUT_ROOT}"
printf 'summary_md,%s\n' "${OUTPUT_ROOT}_summary.md"
