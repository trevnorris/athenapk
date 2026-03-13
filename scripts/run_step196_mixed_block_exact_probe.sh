#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR="${REPO_ROOT}"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full_symbreak.in"
OUTPUT_ROOT="/projects/fluid-engine/out/step196_mixed_block_exact_probe_local"
NX1="128"
NX2="64"
NX3="64"
MB1="32"
MB2="16"
MB3="16"
HISTORY_DT="0.002"
TLIM="1.00"
NLIM="100000"
RESUME="true"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --controlled-input) CONTROLLED_INPUT="$2"; shift 2 ;;
    --full-input) FULL_INPUT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --nx1) NX1="$2"; shift 2 ;;
    --nx2) NX2="$2"; shift 2 ;;
    --nx3) NX3="$2"; shift 2 ;;
    --mb1) MB1="$2"; shift 2 ;;
    --mb2) MB2="$2"; shift 2 ;;
    --mb3) MB3="$2"; shift 2 ;;
    --history-dt) HISTORY_DT="$2"; shift 2 ;;
    --tlim) TLIM="$2"; shift 2 ;;
    --nlim) NLIM="$2"; shift 2 ;;
    --resume) RESUME="$2"; shift 2 ;;
    -h|--help)
      cat <<'USAGE'
Run Step-196 mixed-block exact probe on the corrected quadrature baseline.
USAGE
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

"${REPO_ROOT}/scripts/run_step190_cx_cz_exact_delta_probe.sh" \
  --binary "${BINARY}" \
  --workdir "${WORKDIR}" \
  --controlled-input "${CONTROLLED_INPUT}" \
  --full-input "${FULL_INPUT}" \
  --output-root "${OUTPUT_ROOT}" \
  --nx1 "${NX1}" --nx2 "${NX2}" --nx3 "${NX3}" \
  --mb1 "${MB1}" --mb2 "${MB2}" --mb3 "${MB3}" \
  --history-dt "${HISTORY_DT}" --tlim "${TLIM}" --nlim "${NLIM}" --resume "${RESUME}"

python3 "${REPO_ROOT}/scripts/summarize_step196_mixed_block_exact_probe.py" \
  --glob "${OUTPUT_ROOT}/*" \
  --out-csv "${OUTPUT_ROOT}_mixed_block_exact.csv" \
  --out-md "${OUTPUT_ROOT}_mixed_block_exact.md"

printf 'mixed_block_exact_md,%s\n' "${OUTPUT_ROOT}_mixed_block_exact.md"
