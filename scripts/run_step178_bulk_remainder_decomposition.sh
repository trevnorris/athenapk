#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_ROOT="/projects/fluid-engine/out/step177_aw_exact_delta_probe_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step178_bulk_remainder_decomposition_local"

usage() {
  cat <<'USAGE'
Run Step-178 bulk remainder decomposition on existing Step-177 outputs.

Goal:
  Subtract the exact all-mode A_w,z discrete energy increment from the EM bulk-ledger residual,
  then rank the surviving remainder against the other mixed/resolved bulk-energy channels.

Usage:
  ./scripts/run_step178_bulk_remainder_decomposition.sh [options]
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

mkdir -p "${OUTPUT_ROOT}"

python3 "${REPO_ROOT}/scripts/summarize_step178_bulk_remainder_decomposition.py" \
  --glob "${INPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_bulk_remainder.csv" \
  --out-md "${OUTPUT_ROOT}_bulk_remainder.md"

echo "output_root,${OUTPUT_ROOT}"
echo "bulk_remainder_md,${OUTPUT_ROOT}_bulk_remainder.md"
