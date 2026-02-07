#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BASELINE_BUILD_DIR="${REPO_ROOT}/build-baseline"
MPI_BUILD_DIR="${REPO_ROOT}/build-mpi-hdf5"
LOG_ROOT_DIR="${REPO_ROOT}/build-validation-logs"

SKIP_MPI=0
SKIP_CONFIGURE=0
VERBOSE=0
JOBS="${JOBS:-$(nproc)}"

usage() {
  cat <<'EOF'
Run full modes4d validation with concise PASS/FAIL output.

Usage:
  ./scripts/run_modes4d_validation.sh [options]

Options:
  --baseline-build DIR    Baseline build directory (default: ./build-baseline)
  --mpi-build DIR         MPI+HDF5 build directory (default: ./build-mpi-hdf5)
  --log-dir DIR           Log root directory (default: ./build-validation-logs)
  --jobs N                Build jobs for cmake --build (default: nproc)
  --skip-configure        Skip cmake configure steps
  --skip-mpi              Run only baseline + phase10 checks
  --verbose               Stream full command output in addition to logs
  -h, --help              Show this help
EOF
}

run_step() {
  local name="$1"
  shift
  local slug
  slug="$(echo "${name}" | tr '[:upper:]' '[:lower:]' | tr ' /+' '___')"
  local logfile="${RUN_LOG_DIR}/${slug}.log"
  local start_ts
  local elapsed
  local rc

  start_ts="$(date +%s)"
  printf '[RUN ] %s\n' "${name}"

  if [[ "${VERBOSE}" -eq 1 ]]; then
    set +e
    "$@" 2>&1 | tee "${logfile}"
    rc="${PIPESTATUS[0]}"
    set -e
  else
    set +e
    "$@" >"${logfile}" 2>&1
    rc="$?"
    set -e
  fi

  elapsed="$(( $(date +%s) - start_ts ))"
  if [[ "${rc}" -eq 0 ]]; then
    printf '[PASS] %s (%ss)\n' "${name}" "${elapsed}"
    return 0
  fi

  printf '[FAIL] %s (%ss)\n' "${name}" "${elapsed}"
  printf '       log: %s\n' "${logfile}"
  if [[ "${VERBOSE}" -eq 0 ]]; then
    tail -n 40 "${logfile}" || true
  fi
  return "${rc}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseline-build)
      BASELINE_BUILD_DIR="$2"
      shift 2
      ;;
    --mpi-build)
      MPI_BUILD_DIR="$2"
      shift 2
      ;;
    --log-dir)
      LOG_ROOT_DIR="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    --skip-configure)
      SKIP_CONFIGURE=1
      shift
      ;;
    --skip-mpi)
      SKIP_MPI=1
      shift
      ;;
    --verbose)
      VERBOSE=1
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

mkdir -p "${LOG_ROOT_DIR}"
RUN_LOG_DIR="${LOG_ROOT_DIR}/modes4d_validation_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${RUN_LOG_DIR}"

printf 'Logs: %s\n' "${RUN_LOG_DIR}"

if [[ "${SKIP_CONFIGURE}" -eq 0 ]]; then
  run_step "Configure baseline" \
    cmake -S "${REPO_ROOT}" -B "${BASELINE_BUILD_DIR}" \
      -DPARTHENON_DISABLE_MPI=ON \
      -DPARTHENON_DISABLE_HDF5=ON \
      -DPARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF \
      -DPARTHENON_ENABLE_TESTING=OFF \
      -DCMAKE_BUILD_TYPE=Release
fi

run_step "Build baseline" \
  cmake --build "${BASELINE_BUILD_DIR}" -j"${JOBS}"

run_step "Modes4D unit tests" \
  "${BASELINE_BUILD_DIR}/bin/modes4d_unit_tests"

run_step "Phase-10 baseline regression bundle" \
  cmake --build "${BASELINE_BUILD_DIR}" --target modes4d_phase10_regression

if [[ "${SKIP_MPI}" -eq 0 ]]; then
  if [[ "${SKIP_CONFIGURE}" -eq 0 ]]; then
    run_step "Configure MPI+HDF5" \
      cmake -S "${REPO_ROOT}" -B "${MPI_BUILD_DIR}" \
        -DCMAKE_C_COMPILER=mpicc \
        -DCMAKE_CXX_COMPILER=mpicxx \
        -DPARTHENON_DISABLE_MPI=OFF \
        -DPARTHENON_DISABLE_HDF5=OFF \
        -DPARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF \
        -DPARTHENON_ENABLE_TESTING=OFF \
        -DCMAKE_BUILD_TYPE=Release
  fi

  run_step "Build MPI+HDF5" \
    cmake --build "${MPI_BUILD_DIR}" -j"${JOBS}"

  run_step "MPI+HDF5 regression bundle" \
    cmake --build "${MPI_BUILD_DIR}" --target modes4d_mpi_hdf5_regression
fi

printf '\nAll requested validation stages passed.\n'
printf 'Log directory: %s\n' "${RUN_LOG_DIR}"
