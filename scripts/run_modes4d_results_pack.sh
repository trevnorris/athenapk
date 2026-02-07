#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BUILD_DIR="${REPO_ROOT}/build-baseline"
OUTPUT_DIR="${BUILD_DIR}/modes4d_results_pack"
SKIP_PREREQS=0

usage() {
  cat <<'EOF'
Generate a final results pack combining production and ablation artifacts.

Usage:
  ./scripts/run_modes4d_results_pack.sh [options]

Options:
  --build-dir DIR      Baseline build directory (default: ./build-baseline)
  --output-dir DIR     Results-pack output directory (default: <build-dir>/modes4d_results_pack)
  --skip-prereqs       Skip running production/ablation prerequisite targets
  -h, --help           Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      BUILD_DIR="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --skip-prereqs)
      SKIP_PREREQS=1
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

mkdir -p "${OUTPUT_DIR}"

if [[ "${SKIP_PREREQS}" -eq 0 ]]; then
  cmake --build "${BUILD_DIR}" --target modes4d_harris_lundquist_production
  cmake --build "${BUILD_DIR}" --target modes4d_harris_ablation_campaign
fi

python3 "${REPO_ROOT}/scripts/generate_modes4d_results_pack.py" \
  --lundquist-summary-csv "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/lundquist_scan_summary.csv" \
  --lundquist-report "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/production_report.md" \
  --lundquist-primary-plot "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/plots/lundquist_primary_channels.png" \
  --lundquist-subscale-plot "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/plots/lundquist_subscale_channels.png" \
  --topology-report "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/topology/topology_comparison_report.md" \
  --topology-plot "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/topology/topology_timeseries.png" \
  --topology-status-csv "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/topology/topology_status.csv" \
  --topology-gates-csv "${BUILD_DIR}/modes4d_harris_lundquist_production/outputs/topology/topology_gate_details.csv" \
  --ablation-summary-csv "${BUILD_DIR}/modes4d_harris_ablation_campaign/outputs/ablation_summary.csv" \
  --ablation-report "${BUILD_DIR}/modes4d_harris_ablation_campaign/outputs/ablation_report.md" \
  --ablation-plot "${BUILD_DIR}/modes4d_harris_ablation_campaign/outputs/plots/ablation_channels.png" \
  --output-dir "${OUTPUT_DIR}"
