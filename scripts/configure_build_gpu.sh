#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ -d "/usr/local/cuda-12.4" ]]; then
  CUDA_HOME_DEFAULT="/usr/local/cuda-12.4"
else
  CUDA_HOME_DEFAULT="/usr/local/cuda"
fi

CUDA_HOME="${CUDA_HOME:-${CUDA_HOME_DEFAULT}}"
BUILD_DIR="${REPO_ROOT}/build-gpu-cuda124"
BUILD_TYPE="Release"
KOKKOS_ARCH="AMPERE86"
TARGET="athenaPK"
JOBS="2"
PARTHENON_DISABLE_MPI="ON"
RUN_CONFIGURE=1
RUN_BUILD=1

usage() {
  cat <<'EOF'
Usage: ./scripts/configure_build_gpu.sh [options]

Configure and build AthenaPK with CUDA via Kokkos nvcc_wrapper.

Options:
  --cuda-home DIR       CUDA toolkit root (default: /usr/local/cuda-12.4 if present, else /usr/local/cuda)
  --build-dir DIR       Build directory (default: ./build-gpu-cuda124)
  --build-type TYPE     CMake build type (default: Release)
  --kokkos-arch ARCH    Kokkos GPU arch flag suffix (default: AMPERE86)
  --target NAME         Build target (default: athenaPK)
  --jobs N              Parallel build jobs (default: 2)
  --with-mpi            Enable MPI in this build (default: disabled)
  --configure-only      Run only CMake configure
  --build-only          Run only CMake build
  -h, --help            Show this help

Examples:
  ./scripts/configure_build_gpu.sh
  ./scripts/configure_build_gpu.sh --jobs 4
  ./scripts/configure_build_gpu.sh --build-dir ./build-gpu-custom --kokkos-arch AMPERE86
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cuda-home)
      CUDA_HOME="$2"
      shift 2
      ;;
    --build-dir)
      BUILD_DIR="$2"
      shift 2
      ;;
    --build-type)
      BUILD_TYPE="$2"
      shift 2
      ;;
    --kokkos-arch)
      KOKKOS_ARCH="$2"
      shift 2
      ;;
    --target)
      TARGET="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    --with-mpi)
      PARTHENON_DISABLE_MPI="OFF"
      shift
      ;;
    --configure-only)
      RUN_CONFIGURE=1
      RUN_BUILD=0
      shift
      ;;
    --build-only)
      RUN_CONFIGURE=0
      RUN_BUILD=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "error: unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

NVCC="${CUDA_HOME}/bin/nvcc"
NVCC_WRAPPER="${REPO_ROOT}/external/Kokkos/bin/nvcc_wrapper"

if [[ ! -x "${NVCC}" ]]; then
  echo "error: nvcc not executable at ${NVCC}" >&2
  echo "hint: set --cuda-home or CUDA_HOME to your CUDA 12.x install." >&2
  exit 1
fi
if [[ ! -x "${NVCC_WRAPPER}" ]]; then
  echo "error: nvcc_wrapper not found at ${NVCC_WRAPPER}" >&2
  exit 1
fi

mkdir -p "${BUILD_DIR}"

if [[ "${RUN_CONFIGURE}" -eq 1 ]]; then
  env -u CXXFLAGS -u CFLAGS -u CPPFLAGS -u LDFLAGS \
    PATH="${CUDA_HOME}/bin:${PATH}" \
    CUDACXX="${NVCC}" \
    cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
      -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
      -DPARTHENON_DISABLE_MPI="${PARTHENON_DISABLE_MPI}" \
      -DCMAKE_CXX_COMPILER="${NVCC_WRAPPER}" \
      -DKokkos_ENABLE_CUDA=ON \
      -DKokkos_ARCH_AMPERE86=OFF \
      -DKokkos_ARCH_ADA89=OFF \
      -DKokkos_ARCH_${KOKKOS_ARCH}=ON
fi

if [[ "${RUN_BUILD}" -eq 1 ]]; then
  env -u CXXFLAGS -u CFLAGS -u CPPFLAGS -u LDFLAGS \
    PATH="${CUDA_HOME}/bin:${PATH}" \
    CUDACXX="${NVCC}" \
    cmake --build "${BUILD_DIR}" --target "${TARGET}" -j"${JOBS}"
fi

echo "gpu_build_dir,${BUILD_DIR}"
echo "gpu_binary,${BUILD_DIR}/bin/athenaPK"
