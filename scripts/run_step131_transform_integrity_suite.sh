#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-gpu-cuda124-ada89/bin/athenaPK"
TEST_BUILD_DIR="${REPO_ROOT}/build-baseline"
WORKDIR="${REPO_ROOT}"
OUTPUT_ROOT="/projects/fluid-engine/out/step131_transform_integrity_suite_local"
JOBS="${JOBS:-2}"

AUDIT_STD_SAMPLES="256"
AUDIT_STD_ROUNDTRIPS="64"
AUDIT_DEEP_SAMPLES="1024"
AUDIT_DEEP_ROUNDTRIPS="128"

usage() {
  cat <<'USAGE'
Run Step-131 transform-integrity diagnostic suite.

This is a focused correctness gate before more long-horizon physics sweeps:
  1) C++ modes4d_unit_tests
  2) C++ mode-table Gram/Parseval/roundtrip audit (standard + deep)
  3) projection-operator contract checks
  4) controlled-limit regression smoke

Usage:
  ./scripts/run_step131_transform_integrity_suite.sh [options]

Options:
  --binary PATH            Path to athenaPK binary
  --test-build-dir DIR     Build directory for test/audit targets
  --workdir DIR            Runtime working directory
  --output-root DIR        Output directory root
  --jobs N                 Build jobs (default: 2)
  --std-samples N          Standard audit samples (default: 256)
  --std-roundtrips N       Standard audit repeated roundtrips (default: 64)
  --deep-samples N         Deep audit samples (default: 1024)
  --deep-roundtrips N      Deep audit repeated roundtrips (default: 128)
  -h, --help               Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="$2"; shift 2 ;;
    --test-build-dir) TEST_BUILD_DIR="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --std-samples) AUDIT_STD_SAMPLES="$2"; shift 2 ;;
    --std-roundtrips) AUDIT_STD_ROUNDTRIPS="$2"; shift 2 ;;
    --deep-samples) AUDIT_DEEP_SAMPLES="$2"; shift 2 ;;
    --deep-roundtrips) AUDIT_DEEP_ROUNDTRIPS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

BINARY="$(realpath "${BINARY}")"
TEST_BUILD_DIR="$(realpath -m "${TEST_BUILD_DIR}")"
WORKDIR="$(realpath "${WORKDIR}")"
OUTPUT_ROOT="$(realpath -m "${OUTPUT_ROOT}")"

if [[ ! -x "${BINARY}" ]]; then
  echo "Binary not found/executable: ${BINARY}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_ROOT}"

SUMMARY_CSV="${OUTPUT_ROOT}/transform_integrity_summary.csv"
SUMMARY_MD="${OUTPUT_ROOT}/transform_integrity_summary.md"
echo "check,status,log" > "${SUMMARY_CSV}"
TEST_BUILD_CONFIGURED=0

run_check() {
  local name="$1"
  shift
  local log_path="${OUTPUT_ROOT}/${name}.log"
  local status="PASS"
  if ! "$@" >"${log_path}" 2>&1; then
    status="FAIL"
  fi
  echo "${name},${status},${log_path}" >> "${SUMMARY_CSV}"
  echo "${name},${status},${log_path}"
}

ensure_test_build_configured() {
  if [[ "${TEST_BUILD_CONFIGURED}" -eq 1 ]]; then
    return 0
  fi
  cmake -S "${REPO_ROOT}" -B "${TEST_BUILD_DIR}" \
    -DPARTHENON_DISABLE_MPI=ON \
    -DPARTHENON_DISABLE_HDF5=ON \
    -DPARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF \
    -DPARTHENON_ENABLE_TESTING=OFF \
    -DCMAKE_BUILD_TYPE=Release >&2
  TEST_BUILD_CONFIGURED=1
}

ensure_target_bin() {
  local target="$1"
  local path="${TEST_BUILD_DIR}/bin/${target}"
  ensure_test_build_configured
  if [[ -x "${path}" ]]; then
    echo "${path}"
    return 0
  fi
  cmake --build "${TEST_BUILD_DIR}" --target "${target}" -j"${JOBS}" >&2
  if [[ ! -x "${path}" ]]; then
    echo "Expected target binary missing after build: ${path}" >&2
    return 1
  fi
  echo "${path}"
}

UNIT_BIN="$(ensure_target_bin modes4d_unit_tests)"
AUDIT_BIN="$(ensure_target_bin modes4d_transform_audit)"

run_check "modes4d_unit_tests" \
  "${UNIT_BIN}"

run_check "modes4d_transform_audit_standard" \
  "${AUDIT_BIN}" \
  --samples "${AUDIT_STD_SAMPLES}" \
  --roundtrips "${AUDIT_STD_ROUNDTRIPS}" \
  --gram-tol 5e-10 \
  --parseval-tol 5e-10 \
  --roundtrip-tol 1e-9 \
  --out-csv "${OUTPUT_ROOT}/transform_audit_standard.csv"

run_check "modes4d_transform_audit_deep" \
  "${AUDIT_BIN}" \
  --samples "${AUDIT_DEEP_SAMPLES}" \
  --roundtrips "${AUDIT_DEEP_ROUNDTRIPS}" \
  --gram-tol 5e-10 \
  --parseval-tol 5e-10 \
  --roundtrip-tol 1e-9 \
  --out-csv "${OUTPUT_ROOT}/transform_audit_deep.csv"

run_check "projection_operator_contracts" \
  python3 "${REPO_ROOT}/scripts/modes4d_projection_operator_tests.py" \
    --source-root "${REPO_ROOT}"

run_check "controlled_limit_regression_smoke" \
  python3 "${REPO_ROOT}/scripts/modes4d_controlled_regression.py" \
    --binary "${BINARY}" \
    --workdir "${WORKDIR}" \
    --output-dir "${OUTPUT_ROOT}/controlled_regression_smoke"

python3 - "${SUMMARY_CSV}" "${SUMMARY_MD}" <<'PY'
import csv
import sys
from pathlib import Path

summary_csv = Path(sys.argv[1])
summary_md = Path(sys.argv[2])
rows = list(csv.DictReader(summary_csv.open(encoding="utf-8")))
passed = sum(1 for r in rows if r["status"] == "PASS")
failed = len(rows) - passed
overall = "PASS" if failed == 0 else "FAIL"

lines = [
    "# Step-131 Transform Integrity Suite",
    "",
    f"- overall_status: `{overall}`",
    f"- checks_passed: `{passed}/{len(rows)}`",
    "",
    "| check | status | log |",
    "|---|---|---|",
]
for row in rows:
    lines.append(f"| `{row['check']}` | `{row['status']}` | `{row['log']}` |")

lines.extend(
    [
        "",
        "## Key Artifacts",
        f"- standard audit csv: `{summary_csv.parent / 'transform_audit_standard.csv'}`",
        f"- deep audit csv: `{summary_csv.parent / 'transform_audit_deep.csv'}`",
        f"- controlled smoke outputs: `{summary_csv.parent / 'controlled_regression_smoke'}`",
    ]
)

summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"summary_csv,{summary_csv}")
print(f"summary_md,{summary_md}")
print(f"overall_status,{overall}")
PY

echo ""
echo "Step-131 transform integrity suite complete."
echo "output_root,${OUTPUT_ROOT}"
echo "summary_md,${SUMMARY_MD}"
