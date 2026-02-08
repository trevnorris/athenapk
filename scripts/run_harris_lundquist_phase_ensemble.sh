#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-baseline/bin/athenaPK"
WORKDIR="${REPO_ROOT}/build-baseline/modes4d_harris_lundquist_phase_ensemble"
CONTROLLED_INPUT="${REPO_ROOT}/inputs/harris_4d_controlled.in"
FULL_INPUT="${REPO_ROOT}/inputs/harris_4d_full.in"
OUTPUT_DIR="outputs"
S_VALUES="250,500,1000,2000"
TLIM="0.20"
NLIM="2000"
OUTPUT_DT="0.02"
SEEDS="0,1,2,3,4"
SKIP_PLOTS=1
SKIP_TOPOLOGY=1
FAIL_ON_ANY_FAIL=0
EDOTB_GATE_MODE="strict"

usage() {
  cat <<'EOF'
Run deterministic phase-offset ensemble for Harris Lundquist production scans.

Usage:
  ./scripts/run_harris_lundquist_phase_ensemble.sh [options]

Options:
  --binary PATH            athenaPK binary (default: ./build-baseline/bin/athenaPK)
  --workdir DIR            Working directory
                           (default: ./build-baseline/modes4d_harris_lundquist_phase_ensemble)
  --controlled-input PATH  Controlled input deck
  --full-input PATH        Full input deck
  --output-dir DIR         Output directory (default: <workdir>/outputs)
  --s-values CSV           Lundquist scan values (default: 250,500,1000,2000)
  --tlim FLOAT             Runtime limit per scan run (default: 0.20)
  --nlim INT               Cycle cap per scan run (default: 2000)
  --output-dt FLOAT        History output dt (default: 0.02)
  --seeds CSV              Integer seeds mapped to deterministic phases (default: 0,1,2,3,4)
  --edotb-gate-mode MODE   `strict` (default) or `informational` for edotb_sub
                           scan-level correlation gates during ensemble runs
  --with-plots             Enable production plots for each seed (default: off)
  --with-topology          Enable topology exports for each seed (default: off)
  --fail-on-any-fail       Exit nonzero if any seed case fails
  --scan-arg ARG           Extra argument forwarded to production scan layer (repeatable)
  -h, --help               Show this help
EOF
}

SCAN_ARGS=()
OUTPUT_DIR_EXPLICIT=0

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
    --s-values)
      S_VALUES="$2"
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
    --seeds)
      SEEDS="$2"
      shift 2
      ;;
    --edotb-gate-mode)
      EDOTB_GATE_MODE="$2"
      shift 2
      ;;
    --with-plots)
      SKIP_PLOTS=0
      shift
      ;;
    --with-topology)
      SKIP_TOPOLOGY=0
      shift
      ;;
    --fail-on-any-fail)
      FAIL_ON_ANY_FAIL=1
      shift
      ;;
    --scan-arg)
      SCAN_ARGS+=("$2")
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

mkdir -p "${WORKDIR}"
WORKDIR_ABS="$(cd "${WORKDIR}" && pwd)"
if [[ "${OUTPUT_DIR_EXPLICIT}" -eq 1 ]]; then
  if [[ "${OUTPUT_DIR}" = /* ]]; then
    OUTPUT_DIR_RESOLVED="${OUTPUT_DIR}"
  else
    OUTPUT_DIR_RESOLVED="${WORKDIR_ABS}/${OUTPUT_DIR}"
  fi
else
  OUTPUT_DIR_RESOLVED="${WORKDIR_ABS}/outputs"
fi
mkdir -p "${OUTPUT_DIR_RESOLVED}"

CMD=(
  python3 "${REPO_ROOT}/scripts/harris_lundquist_phase_ensemble.py"
  --binary "${BINARY}"
  --workdir "${WORKDIR_ABS}"
  --controlled-input "${CONTROLLED_INPUT}"
  --full-input "${FULL_INPUT}"
  --output-dir "${OUTPUT_DIR_RESOLVED}"
  --s-values "${S_VALUES}"
  --tlim "${TLIM}"
  --nlim "${NLIM}"
  --output-dt "${OUTPUT_DT}"
  --seeds "${SEEDS}"
  --edotb-gate-mode "${EDOTB_GATE_MODE}"
)

if [[ "${SKIP_PLOTS}" -eq 1 ]]; then
  CMD+=(--skip-plots)
fi
if [[ "${SKIP_TOPOLOGY}" -eq 1 ]]; then
  CMD+=(--skip-topology)
fi
if [[ "${FAIL_ON_ANY_FAIL}" -eq 1 ]]; then
  CMD+=(--fail-on-any-fail)
fi
for arg in "${SCAN_ARGS[@]}"; do
  CMD+=(--scan-arg "$arg")
done

"${CMD[@]}"
