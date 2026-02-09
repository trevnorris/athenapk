#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'EOF'
Check local CUDA runtime visibility before launching GPU scans.

Usage:
  ./scripts/check_cuda_runtime.sh [options]

Options:
  --cuda-home PATH   CUDA toolkit root (default: /usr/local/cuda-12.4 if present,
                     else /usr/local/cuda)
  --binary PATH      Optional athenaPK GPU binary to run a 1-step smoke test
  --help             Show this help
EOF
}

CUDA_HOME_DEFAULT="/usr/local/cuda"
if [[ -d "/usr/local/cuda-12.4" ]]; then
  CUDA_HOME_DEFAULT="/usr/local/cuda-12.4"
fi

CUDA_HOME="${CUDA_HOME_DEFAULT}"
BINARY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cuda-home)
      CUDA_HOME="$2"
      shift 2
      ;;
    --binary)
      BINARY="$2"
      shift 2
      ;;
    --help|-h)
      show_help
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      show_help
      exit 1
      ;;
  esac
done

if [[ ! -x "${CUDA_HOME}/bin/nvcc" ]]; then
  echo "Missing nvcc at ${CUDA_HOME}/bin/nvcc" >&2
  exit 1
fi

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "nvidia-smi not found in PATH." >&2
  exit 1
fi

echo "cuda_home,${CUDA_HOME}"
echo "nvcc,$(${CUDA_HOME}/bin/nvcc --version | tail -n 1)"
echo "nvidia_smi,$(nvidia-smi --query-gpu=name,driver_version,compute_mode --format=csv,noheader | head -n 1)"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

cat > "${tmpdir}/cuda_probe.cu" <<'EOF'
#include <cstdio>
#include <cuda_runtime.h>

int main() {
  int count = 0;
  cudaError_t err = cudaGetDeviceCount(&count);
  if (err != cudaSuccess) {
    std::fprintf(stderr, "cudaGetDeviceCount failed: %s\n", cudaGetErrorString(err));
    return 2;
  }
  std::printf("cuda_device_count,%d\n", count);
  return (count > 0) ? 0 : 3;
}
EOF

"${CUDA_HOME}/bin/nvcc" -O2 "${tmpdir}/cuda_probe.cu" -o "${tmpdir}/cuda_probe"
"${tmpdir}/cuda_probe"

if [[ -n "${BINARY}" ]]; then
  if [[ ! -x "${BINARY}" ]]; then
    echo "Binary is not executable: ${BINARY}" >&2
    exit 1
  fi
  repo_root="$(cd "$(dirname "${BINARY}")/.." && pwd -P)"
  input_guess="${repo_root}/../inputs/harris_4d_controlled.in"
  if [[ ! -f "${input_guess}" ]]; then
    input_guess="$(cd "$(dirname "$0")/.." && pwd -P)/inputs/harris_4d_controlled.in"
  fi
  if [[ ! -f "${input_guess}" ]]; then
    echo "Unable to locate inputs/harris_4d_controlled.in for binary smoke test." >&2
    exit 1
  fi
  "${BINARY}" -i "${input_guess}" \
    parthenon/time/nlim=1 \
    parthenon/time/tlim=1.0e-6 \
    parthenon/output0/dt=-1 \
    parthenon/output1/dt=-1 \
    problem/harris_4d/perturb_phase_x=0.0 \
    problem/harris_4d/perturb_phase_z=0.0 >/dev/null
  echo "athenapk_smoke,PASS"
fi
