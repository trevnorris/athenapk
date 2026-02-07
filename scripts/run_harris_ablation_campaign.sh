#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-baseline/bin/athenaPK"
WORKDIR="${REPO_ROOT}/build-baseline/modes4d_harris_ablation_campaign"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full.in"
OUTPUT_DIR="${WORKDIR}/outputs"
TLIM="0.10"
NLIM="1000"
OUTPUT_DT="0.01"
ABLATION_STRONG_RATIO_THRESHOLD="0.5"
FAIL_ON_ABLATION_CHECK=0

usage() {
  cat <<'EOF'
Run the Harris 4D ablation campaign and generate CSV/plot/report artifacts.

Usage:
  ./scripts/run_harris_ablation_campaign.sh [options]

Options:
  --binary PATH            athenaPK binary (default: ./build-baseline/bin/athenaPK)
  --workdir DIR            Working directory (default: ./build-baseline/modes4d_harris_ablation_campaign)
  --controlled-input PATH  Controlled input deck
  --full-input PATH        Full input deck
  --output-dir DIR         Output directory (default: <workdir>/outputs)
  --tlim FLOAT             Runtime limit per run (default: 0.10)
  --nlim INT               Cycle cap per run (default: 1000)
  --output-dt FLOAT        History dt override (default: 0.01)
  --ablation-strong-ratio-threshold FLOAT
                          Expected strong-reduction ratio threshold (default: 0.5)
  --fail-on-ablation-check Exit nonzero when ablation checks fail
  -h, --help               Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary)
      BINARY="$2"
      shift 2
      ;;
    --workdir)
      WORKDIR="$2"
      shift 2
      ;;
    --controlled-input)
      CONTROLLED_INPUT="$2"
      shift 2
      ;;
    --full-input)
      FULL_INPUT="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --tlim)
      TLIM="$2"
      shift 2
      ;;
    --nlim)
      NLIM="$2"
      shift 2
      ;;
    --output-dt)
      OUTPUT_DT="$2"
      shift 2
      ;;
    --ablation-strong-ratio-threshold)
      ABLATION_STRONG_RATIO_THRESHOLD="$2"
      shift 2
      ;;
    --fail-on-ablation-check)
      FAIL_ON_ABLATION_CHECK=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

mkdir -p "${WORKDIR}" "${OUTPUT_DIR}"

CMD=(
  python3 "${REPO_ROOT}/scripts/harris_ablation_campaign.py"
  --binary "${BINARY}"
  --workdir "${WORKDIR}"
  --controlled-input "${CONTROLLED_INPUT}"
  --full-input "${FULL_INPUT}"
  --output-dir "${OUTPUT_DIR}"
  --tlim "${TLIM}"
  --nlim "${NLIM}"
  --output-dt "${OUTPUT_DT}"
  --ablation-strong-ratio-threshold "${ABLATION_STRONG_RATIO_THRESHOLD}"
  --check-transport-closure
  --check-em-bulk-ledger
)

if [[ "${FAIL_ON_ABLATION_CHECK}" -eq 1 ]]; then
  CMD+=(--fail-on-ablation-check)
fi

"${CMD[@]}"
