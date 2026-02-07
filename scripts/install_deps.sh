#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REQUIREMENTS_FILE="${REPO_ROOT}/requirements.txt"
PYTHON_BIN="${PYTHON_BIN:-python3}"

INSTALL_SYSTEM=true
INSTALL_PYTHON=true
ASSUME_YES=false

usage() {
  cat <<'EOF'
Usage: ./scripts/install_deps.sh [options]

Install AthenaPK simulation dependencies (system + Python).

Options:
  --python-only      Install only Python dependencies
  --system-only      Install only system dependencies
  --python <bin>     Python executable to use (default: PYTHON_BIN env or python3)
  -y, --yes          Pass -y to apt-get install
  -h, --help         Show this help message

Examples:
  ./scripts/install_deps.sh --yes
  ./scripts/install_deps.sh --python-only
  ./scripts/install_deps.sh --system-only --yes
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python-only)
      INSTALL_SYSTEM=false
      ;;
    --system-only)
      INSTALL_PYTHON=false
      ;;
    --python)
      shift
      if [[ $# -eq 0 ]]; then
        echo "error: --python requires an argument" >&2
        exit 2
      fi
      PYTHON_BIN="$1"
      ;;
    -y|--yes)
      ASSUME_YES=true
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
  shift
done

if ! ${INSTALL_SYSTEM} && ! ${INSTALL_PYTHON}; then
  echo "error: nothing to do (both system and python installs disabled)" >&2
  exit 2
fi

if [[ ! -f "${REQUIREMENTS_FILE}" ]]; then
  echo "error: requirements file not found: ${REQUIREMENTS_FILE}" >&2
  exit 1
fi

if ${INSTALL_SYSTEM}; then
  echo "[deps] Installing system packages..."
  if ! command -v apt-get >/dev/null 2>&1; then
    echo "error: apt-get not found; install system dependencies manually." >&2
    exit 1
  fi

  SUDO_CMD=()
  if [[ ${EUID} -ne 0 ]]; then
    if command -v sudo >/dev/null 2>&1; then
      SUDO_CMD=(sudo)
    else
      echo "error: root or sudo required for apt-get installs." >&2
      exit 1
    fi
  fi

  APT_PACKAGES=(
    build-essential
    cmake
    ninja-build
    python3
    python3-pip
    python3-venv
    python3-dev
    openmpi-bin
    libopenmpi-dev
    libhdf5-dev
    libhdf5-openmpi-dev
  )

  "${SUDO_CMD[@]}" apt-get update
  APT_INSTALL_CMD=(apt-get install)
  if ${ASSUME_YES}; then
    APT_INSTALL_CMD+=(-y)
  fi
  APT_INSTALL_CMD+=("${APT_PACKAGES[@]}")
  "${SUDO_CMD[@]}" "${APT_INSTALL_CMD[@]}"
fi

if ${INSTALL_PYTHON}; then
  echo "[deps] Installing Python packages from ${REQUIREMENTS_FILE}..."
  "${PYTHON_BIN}" -m pip install --upgrade pip
  "${PYTHON_BIN}" -m pip install -r "${REQUIREMENTS_FILE}"
  echo "[deps] Verifying Python imports..."
  "${PYTHON_BIN}" -c "import h5py, unyt, numpy, scipy, matplotlib; print('python deps ok')"
fi

echo "[deps] Done."
