#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${BINARY:-${REPO_ROOT}/build-mpi-hdf5/bin/athenaPK}"
INPUT="${INPUT:-${REPO_ROOT}/inputs/harris_4d_full.in}"
WORKDIR="${WORKDIR:-${REPO_ROOT}/build-mpi-hdf5/runs/harris_full_mpi}"
MPIRUN_BIN="${MPIRUN_BIN:-mpirun}"
RANKS=2
TLIM="0.1"
NLIM="1000"
OUTPUT_DT="0.01"
EXTRA_ARGS=()

usage() {
  cat <<'EOF'
Usage: ./scripts/run_harris_full_mpi_hdf5.sh [options] [-- <athenaPK overrides...>]

Run a full-channel Harris 4D case with MPI and HDF5 outputs.

Options:
  --binary <path>      athenaPK binary (default: build-mpi-hdf5/bin/athenaPK)
  --input <path>       Input deck (default: inputs/harris_4d_full.in)
  --workdir <path>     Run directory (default: build-mpi-hdf5/runs/harris_full_mpi)
  --ranks <int>        MPI ranks (default: 2)
  --tlim <real>        Simulation end time override (default: 0.1)
  --nlim <int>         Cycle limit override (default: 1000)
  --output-dt <real>   HDF5 output cadence override (default: 0.01)
  --mpirun <path>      mpirun executable (default: mpirun)
  -h, --help           Show this help

Examples:
  ./scripts/run_harris_full_mpi_hdf5.sh --ranks 4 --tlim 0.05
  ./scripts/run_harris_full_mpi_hdf5.sh --workdir /tmp/harris_mpi -- --modes4d/n_modes=6
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary)
      shift
      BINARY="$1"
      ;;
    --input)
      shift
      INPUT="$1"
      ;;
    --workdir)
      shift
      WORKDIR="$1"
      ;;
    --ranks)
      shift
      RANKS="$1"
      ;;
    --tlim)
      shift
      TLIM="$1"
      ;;
    --nlim)
      shift
      NLIM="$1"
      ;;
    --output-dt)
      shift
      OUTPUT_DT="$1"
      ;;
    --mpirun)
      shift
      MPIRUN_BIN="$1"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      EXTRA_ARGS+=("$@")
      break
      ;;
    *)
      echo "error: unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
  shift
done

if [[ ! -x "${BINARY}" ]]; then
  echo "error: athenaPK binary not executable: ${BINARY}" >&2
  exit 1
fi
if [[ ! -f "${INPUT}" ]]; then
  echo "error: input file not found: ${INPUT}" >&2
  exit 1
fi

mkdir -p "${WORKDIR}"

CMD=(
  "${MPIRUN_BIN}" -n "${RANKS}" "${BINARY}" -i "${INPUT}"
  "parthenon/time/tlim=${TLIM}"
  "parthenon/time/nlim=${NLIM}"
  "parthenon/output0/dt=${OUTPUT_DT}"
)
CMD+=("${EXTRA_ARGS[@]}")

echo "[run] workdir=${WORKDIR}"
echo "[run] ${CMD[*]}"
(
  cd "${WORKDIR}"
  "${CMD[@]}"
)

echo "[run] outputs:"
ls -1 "${WORKDIR}"/parthenon.out0*.phdf 2>/dev/null || true
ls -1 "${WORKDIR}"/parthenon.out1.hst 2>/dev/null || true
