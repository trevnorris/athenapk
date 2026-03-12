#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_ROOT="/projects/fluid-engine/out/step193_exact_current_cxz_equivalence_probe_local"
OUTPUT_ROOT="/projects/fluid-engine/out/step194_postfix_aw_update_reaudit_local"

usage() {
  cat <<'USAGE'
Run Step-194 post-fix A_w update re-audit.

Goal:
  Re-evaluate the A_w update-power decomposition on the post-fix baseline so we
  can separate obsolete pre-fix conclusions from the current authoritative branch.
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
python3 "${REPO_ROOT}/scripts/summarize_step194_postfix_aw_update_reaudit.py" \
  --input-root "${INPUT_ROOT}" \
  --out-csv "${OUTPUT_ROOT}_summary.csv" \
  --out-md "${OUTPUT_ROOT}_summary.md"

printf 'output_root,%s\n' "${OUTPUT_ROOT}"
printf 'summary_md,%s\n' "${OUTPUT_ROOT}_summary.md"
