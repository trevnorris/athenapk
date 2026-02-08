#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-baseline/bin/athenaPK"
WORKDIR="${REPO_ROOT}/build-baseline/modes4d_harris_lundquist_envelope_campaign"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full.in"
OUTPUT_DIR=""
S_GROUPS="250,500,1000,2000;250,500,1000,2000,4000"
TLIMS="0.15,0.20,0.25,0.30"
BASE_TLIM="0.20"
BASE_NLIM="2000"
BASE_OUTPUT_DT="0.02"
MAX_CASES="0"
FAIL_ON_NO_PASS=0
SKIP_PLOTS=1
SKIP_TOPOLOGY=1
OUTPUT_DIR_EXPLICIT=0

usage() {
  cat <<'EOF'
Run a Lundquist-envelope campaign over S/tlim grids.

Usage:
  ./scripts/run_harris_lundquist_envelope_campaign.sh [options]

Options:
  --binary PATH            athenaPK binary (default: ./build-baseline/bin/athenaPK)
  --workdir DIR            Campaign working directory
                           (default: ./build-baseline/modes4d_harris_lundquist_envelope_campaign)
  --controlled-input PATH  Controlled input deck
  --full-input PATH        Full input deck
  --output-dir DIR         Campaign output directory (default: <workdir>/outputs)
  --s-groups TEXT          Semicolon-separated S groups
                           (default: 250,500,1000,2000;250,500,1000,2000,4000)
  --tlims CSV              Comma-separated tlim values (default: 0.15,0.20,0.25,0.30)
  --base-tlim FLOAT        Baseline tlim used for nlim/output_dt scaling (default: 0.20)
  --base-nlim INT          Baseline nlim at base-tlim (default: 2000)
  --base-output-dt FLOAT   Baseline output dt at base-tlim (default: 0.02)
  --max-cases INT          Limit number of campaign cases (default: 0 = no limit)
  --fail-on-no-pass        Exit nonzero if no passing envelope is found
  --with-plots             Enable per-case production plots (default: off)
  --with-topology          Enable per-case topology exports (default: off)
  --extra-arg ARG          Extra arg forwarded to production wrapper (repeatable)
  -h, --help               Show this help
EOF
}

EXTRA_ARGS=()

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
      OUTPUT_DIR_EXPLICIT=1
      shift 2
      ;;
    --s-groups)
      S_GROUPS="$2"
      shift 2
      ;;
    --tlims)
      TLIMS="$2"
      shift 2
      ;;
    --base-tlim)
      BASE_TLIM="$2"
      shift 2
      ;;
    --base-nlim)
      BASE_NLIM="$2"
      shift 2
      ;;
    --base-output-dt)
      BASE_OUTPUT_DT="$2"
      shift 2
      ;;
    --max-cases)
      MAX_CASES="$2"
      shift 2
      ;;
    --fail-on-no-pass)
      FAIL_ON_NO_PASS=1
      shift
      ;;
    --with-plots)
      SKIP_PLOTS=0
      shift
      ;;
    --with-topology)
      SKIP_TOPOLOGY=0
      shift
      ;;
    --extra-arg)
      EXTRA_ARGS+=("$2")
      shift 2
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

if [[ "${OUTPUT_DIR_EXPLICIT}" -eq 0 ]]; then
  OUTPUT_DIR="outputs"
fi

if [[ "${OUTPUT_DIR_EXPLICIT}" -eq 1 ]]; then
  mkdir -p "${WORKDIR}" "${OUTPUT_DIR}"
else
  mkdir -p "${WORKDIR}" "${WORKDIR}/outputs"
fi

CMD=(
  python3 "${REPO_ROOT}/scripts/harris_lundquist_envelope_campaign.py"
  --binary "${BINARY}"
  --workdir "${WORKDIR}"
  --controlled-input "${CONTROLLED_INPUT}"
  --full-input "${FULL_INPUT}"
  --output-dir "${OUTPUT_DIR}"
  --s-groups "${S_GROUPS}"
  --tlims "${TLIMS}"
  --base-tlim "${BASE_TLIM}"
  --base-nlim "${BASE_NLIM}"
  --base-output-dt "${BASE_OUTPUT_DT}"
  --max-cases "${MAX_CASES}"
)

if [[ "${FAIL_ON_NO_PASS}" -eq 1 ]]; then
  CMD+=(--fail-on-no-pass)
fi
if [[ "${SKIP_PLOTS}" -eq 1 ]]; then
  CMD+=(--skip-plots)
fi
if [[ "${SKIP_TOPOLOGY}" -eq 1 ]]; then
  CMD+=(--skip-topology)
fi
for arg in "${EXTRA_ARGS[@]}"; do
  # Keep values that begin with '--' from being parsed as top-level options.
  CMD+=("--extra-arg=${arg}")
done

"${CMD[@]}"
