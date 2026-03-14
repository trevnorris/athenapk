#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
INPUT_ROOT="/projects/fluid-engine/out/step200_zblock_exact_node_work_probe_local"
OUTPUT_PREFIX="/projects/fluid-engine/out/step203_postz_remainder_decomposition_local"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/run_step203_postz_remainder_decomposition.sh [options]

Options:
  --input-root PATH     Input step root containing case directories.
  --output-prefix PATH  Output prefix for summary artifacts.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root)
      INPUT_ROOT="$2"
      shift 2
      ;;
    --output-prefix)
      OUTPUT_PREFIX="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

python3 "${REPO_ROOT}/scripts/summarize_step203_postz_remainder_decomposition.py" \
  --glob "${INPUT_ROOT}/*" \
  --out-csv "${OUTPUT_PREFIX}_summary.csv" \
  --out-md "${OUTPUT_PREFIX}_summary.md"

echo "summary_csv,${OUTPUT_PREFIX}_summary.csv"
echo "summary_md,${OUTPUT_PREFIX}_summary.md"
