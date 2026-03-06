#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_ROOT="/projects/fluid-engine/out/step150_em_pi_transport_probe_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step152_em_pi_source_decomposition_local"
BASELINE_CASE="baseline_pi_transport"

usage() {
  cat <<'USAGE'
Run Step-152 EM pi source decomposition on existing Step-150 histories.

Goal:
  Compare the failing EM pi source integrals against the existing signed
  mode-0 RHS diagnostics (`a0`, `ay`, `aw`) and expose any unaccounted gap.

Usage:
  ./scripts/run_step152_em_pi_source_decomposition.sh [options]

Options:
  --input-root DIR       Step-150 output root (default:
                         /projects/fluid-engine/out/step150_em_pi_transport_probe_local)
  --output-root DIR      Output prefix root (default:
                         /projects/fluid-engine/out/step152_em_pi_source_decomposition_local)
  --baseline-case NAME   Baseline case name (default: baseline_pi_transport)
  -h, --help             Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --baseline-case) BASELINE_CASE="$2"; shift 2 ;;
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

mkdir -p "${OUTPUT_ROOT}"

python3 "${REPO_ROOT}/scripts/summarize_step152_em_pi_source_decomposition.py" \
  --glob "${INPUT_ROOT}/*" \
  --baseline-case "${BASELINE_CASE}" \
  --out-csv "${OUTPUT_ROOT}_channels.csv" \
  --summary-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${OUTPUT_ROOT}_summary.md"
echo "summary_csv,${OUTPUT_ROOT}_summary.csv"
echo "channels_csv,${OUTPUT_ROOT}_channels.csv"
