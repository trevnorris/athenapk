#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step157_quadrature_ledger_decomposition_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step159_em_u_bulk_partition_local"
FOCUS_WINDOW="late"
TOP_N=8
MU0=1.0

usage() {
  cat <<'USAGE'
Run Step-159 EM bulk-energy partition analysis on existing Step-157 histories.

Goal:
  Decompose m4d_em_u_bulk into resolved, mixed, and leftover components under the
  corrected quadrature baseline to identify which energy partition still drives
  the conservative EM ledger residual.

Usage:
  ./scripts/run_step159_em_u_bulk_partition.sh [options]

Options:
  --input-root DIR       Step-157 output root (default:
                         /projects/fluid-engine/out/step157_quadrature_ledger_decomposition_local)
  --output-root DIR      Output prefix root (default:
                         /projects/fluid-engine/out/step159_em_u_bulk_partition_local)
  --focus-window NAME    Window to summarize: early, mid, late, all
                         (default: late)
  --top-n N              Number of ranked components per case/window (default: 8)
  --mu0 VALUE            EM permeability scaling for mixed/brane channels (default: 1.0)
  -h, --help             Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --focus-window) FOCUS_WINDOW="$2"; shift 2 ;;
    --top-n) TOP_N="$2"; shift 2 ;;
    --mu0) MU0="$2"; shift 2 ;;
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

python3 "${REPO_ROOT}/scripts/summarize_step159_em_u_bulk_partition.py" \
  --glob "${INPUT_ROOT}/*" \
  --focus-window "${FOCUS_WINDOW}" \
  --top-n "${TOP_N}" \
  --mu0 "${MU0}" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "summary_csv,${OUTPUT_ROOT}_summary.csv"
