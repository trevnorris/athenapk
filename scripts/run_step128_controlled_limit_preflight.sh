#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY_DEFAULT="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
WORKDIR_DEFAULT="${REPO_ROOT}/build-gpu-cuda124-ada89"
OUTPUT_ROOT_DEFAULT="/projects/fluid-engine/out/step128_controlled_limit_preflight_local"

BINARY="${BINARY_DEFAULT}"
WORKDIR="${WORKDIR_DEFAULT}"
OUTPUT_ROOT="${OUTPUT_ROOT_DEFAULT}"

usage() {
  cat <<'EOF'
Run Step-128 controlled-limit preflight checks (fast correctness baseline).

Checks:
1) modes4d_unit_tests
2) projection-operator kernel contracts
3) controlled-limit channel inactivity regression (includes smoke deck runs)

Usage:
  ./scripts/run_step128_controlled_limit_preflight.sh [options]

Options:
  --binary PATH      Path to athenaPK binary
  --workdir PATH     Work directory for athenaPK runs
  --output-root DIR  Output directory root
  -h, --help         Show this help
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
    --output-root)
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

BINARY="$(realpath "${BINARY}")"
WORKDIR="$(realpath "${WORKDIR}")"
OUTPUT_ROOT="$(realpath -m "${OUTPUT_ROOT}")"

mkdir -p "${OUTPUT_ROOT}"

SUMMARY_CSV="${OUTPUT_ROOT}/preflight_summary.csv"
SUMMARY_MD="${OUTPUT_ROOT}/preflight_summary.md"

echo "check,status,log" > "${SUMMARY_CSV}"

run_check() {
  local check_name="$1"
  shift
  local log_path="${OUTPUT_ROOT}/${check_name}.log"
  local status="PASS"
  if ! "$@" >"${log_path}" 2>&1; then
    status="FAIL"
  fi
  echo "${check_name},${status},${log_path}" >> "${SUMMARY_CSV}"
  echo "${check_name},${status},${log_path}"
}

find_unit_test_bin() {
  local candidate_gpu="$(dirname "${BINARY}")/modes4d_unit_tests"
  local candidate_baseline="${REPO_ROOT}/build-baseline/bin/modes4d_unit_tests"
  local candidate_workdir="${WORKDIR}/bin/modes4d_unit_tests"

  if [[ -x "${candidate_gpu}" ]]; then
    echo "${candidate_gpu}"
    return 0
  fi
  if [[ -x "${candidate_workdir}" ]]; then
    echo "${candidate_workdir}"
    return 0
  fi
  if [[ -x "${candidate_baseline}" ]]; then
    echo "${candidate_baseline}"
    return 0
  fi
  return 1
}

run_modes4d_unit_tests() {
  local unit_bin
  if unit_bin="$(find_unit_test_bin)"; then
    "${unit_bin}"
    return 0
  fi

  # If no unit-test binary is present, try building the dedicated target.
  cmake --build "${WORKDIR}" --target modes4d_unit_tests -j"${JOBS:-$(nproc)}"

  unit_bin="$(find_unit_test_bin)"
  "${unit_bin}"
}

run_check "modes4d_unit_tests" \
  run_modes4d_unit_tests

run_check "projection_operator_contracts" \
  python3 "${REPO_ROOT}/scripts/modes4d_projection_operator_tests.py" \
    --source-root "${REPO_ROOT}"

run_check "controlled_limit_regression" \
  python3 "${REPO_ROOT}/scripts/modes4d_controlled_regression.py" \
    --binary "${BINARY}" \
    --workdir "${WORKDIR}" \
    --output-dir "${OUTPUT_ROOT}/controlled_regression_outputs"

python3 - "${SUMMARY_CSV}" "${SUMMARY_MD}" <<'PY'
import csv
import sys
from pathlib import Path

summary_csv = Path(sys.argv[1])
summary_md = Path(sys.argv[2])

rows = list(csv.DictReader(summary_csv.open()))
passed = sum(1 for r in rows if r["status"] == "PASS")
failed = len(rows) - passed
overall = "PASS" if failed == 0 else "FAIL"

lines = [
    "# Step-128 Controlled-Limit Preflight",
    "",
    f"- overall_status: `{overall}`",
    f"- checks_passed: `{passed}/{len(rows)}`",
    "",
    "| check | status | log |",
    "|---|---|---|",
]
for row in rows:
    lines.append(f"| `{row['check']}` | `{row['status']}` | `{row['log']}` |")

summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"summary_csv,{summary_csv}")
print(f"summary_md,{summary_md}")
print(f"overall_status,{overall}")
PY

echo ""
echo "Step-128 controlled-limit preflight complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${SUMMARY_MD}"
