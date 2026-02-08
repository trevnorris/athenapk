#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BINARY="${REPO_ROOT}/build-baseline/bin/athenaPK"
WORKROOT="${REPO_ROOT}/build-baseline"
PROFILE_TAG="highs_pub_v2"
S_VALUES="250,500,1000,2000,4000"
TLIM="0.20"
NLIM="2000"
OUTPUT_DT="0.02"
CLOSURE_LOCAL_MODE0_ABS_RATE_TOL="2.0e-6"
SEEDS="0,1,2,3,4,5,6,7,8,9"
RUN_MPI_SMOKE=1
MPI_BINARY="${REPO_ROOT}/build-mpi-hdf5/bin/athenaPK"
MPI_WORKDIR="${REPO_ROOT}/build-mpi-hdf5/runs/harris_full_mpi_highs_pub_v2"
MPI_RANKS="2"
MPI_TLIM="0.02"
MPI_NLIM="200"
MPI_OUTPUT_DT="0.01"
MPIRUN_BIN="mpirun"

usage() {
  cat <<'EOF'
Run the frozen highs_pub_v2 publication bundle end-to-end.

Usage:
  ./scripts/run_harris_highs_pub_v2_bundle.sh [options]

Options:
  --binary PATH            athenaPK binary for production/ablation/ensemble
                           (default: ./build-baseline/bin/athenaPK)
  --workroot DIR           Build/output root for artifacts
                           (default: ./build-baseline)
  --profile-tag TAG        Artifact suffix tag (default: highs_pub_v2)
  --s-values CSV           Lundquist values (default: 250,500,1000,2000,4000)
  --tlim FLOAT             Runtime limit per scan run (default: 0.20)
  --nlim INT               Cycle cap per scan run (default: 2000)
  --output-dt FLOAT        History dt per scan run (default: 0.02)
  --closure-local-mode0-abs-rate-tol FLOAT
                           Full-case local closure abs-rate tolerance (default: 2.0e-6)
  --seeds CSV              Phase-ensemble seeds (default: 0..9)
  --skip-mpi-smoke         Skip MPI+HDF5 smoke run
  --mpi-binary PATH        athenaPK MPI binary (default: ./build-mpi-hdf5/bin/athenaPK)
  --mpi-workdir DIR        MPI smoke workdir (default: ./build-mpi-hdf5/runs/harris_full_mpi_highs_pub_v2)
  --mpi-ranks INT          MPI smoke ranks (default: 2)
  --mpi-tlim FLOAT         MPI smoke tlim (default: 0.02)
  --mpi-nlim INT           MPI smoke nlim (default: 200)
  --mpi-output-dt FLOAT    MPI smoke output dt (default: 0.01)
  --mpirun PATH            MPI launcher (default: mpirun)
  -h, --help               Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary)
      BINARY="$2"
      shift 2
      ;;
    --workroot)
      WORKROOT="$2"
      shift 2
      ;;
    --profile-tag)
      PROFILE_TAG="$2"
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
    --closure-local-mode0-abs-rate-tol)
      CLOSURE_LOCAL_MODE0_ABS_RATE_TOL="$2"
      shift 2
      ;;
    --seeds)
      SEEDS="$2"
      shift 2
      ;;
    --skip-mpi-smoke)
      RUN_MPI_SMOKE=0
      shift
      ;;
    --mpi-binary)
      MPI_BINARY="$2"
      shift 2
      ;;
    --mpi-workdir)
      MPI_WORKDIR="$2"
      shift 2
      ;;
    --mpi-ranks)
      MPI_RANKS="$2"
      shift 2
      ;;
    --mpi-tlim)
      MPI_TLIM="$2"
      shift 2
      ;;
    --mpi-nlim)
      MPI_NLIM="$2"
      shift 2
      ;;
    --mpi-output-dt)
      MPI_OUTPUT_DT="$2"
      shift 2
      ;;
    --mpirun)
      MPIRUN_BIN="$2"
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

PROD_WORKDIR="${WORKROOT}/modes4d_harris_lundquist_production_${PROFILE_TAG}"
ABL_WORKDIR="${WORKROOT}/modes4d_harris_ablation_campaign_${PROFILE_TAG}"
PHASE_WORKDIR="${WORKROOT}/modes4d_harris_lundquist_phase_ensemble_${PROFILE_TAG}"
RESULTS_DIR="${WORKROOT}/modes4d_results_pack_${PROFILE_TAG}"
UNCERTAINTY_DIR="${RESULTS_DIR}/uncertainty_bands"

mkdir -p "${WORKROOT}" "${RESULTS_DIR}" "${UNCERTAINTY_DIR}"

echo "[run] production (${PROFILE_TAG})"
"${REPO_ROOT}/scripts/run_harris_lundquist_production.sh" \
  --binary "${BINARY}" \
  --workdir "${PROD_WORKDIR}" \
  --s-values "${S_VALUES}" \
  --tlim "${TLIM}" \
  --nlim "${NLIM}" \
  --output-dt "${OUTPUT_DT}" \
  --closure-local-mode0-abs-rate-tol "${CLOSURE_LOCAL_MODE0_ABS_RATE_TOL}"

echo "[run] ablation (${PROFILE_TAG})"
"${REPO_ROOT}/scripts/run_harris_ablation_campaign.sh" \
  --binary "${BINARY}" \
  --workdir "${ABL_WORKDIR}" \
  --tlim 0.10 \
  --nlim 1000 \
  --output-dt 0.01

echo "[run] phase ensemble (${PROFILE_TAG})"
"${REPO_ROOT}/scripts/run_harris_lundquist_phase_ensemble.sh" \
  --binary "${BINARY}" \
  --workdir "${PHASE_WORKDIR}" \
  --s-values "${S_VALUES}" \
  --tlim "${TLIM}" \
  --nlim "${NLIM}" \
  --output-dt "${OUTPUT_DT}" \
  --seeds "${SEEDS}" \
  --edotb-gate-mode informational \
  --fail-on-any-fail \
  --scan-arg "--closure-local-mode0-abs-rate-tol=${CLOSURE_LOCAL_MODE0_ABS_RATE_TOL}"

echo "[run] uncertainty bands (${PROFILE_TAG})"
python3 "${REPO_ROOT}/scripts/plot_harris_lundquist_uncertainty_bands.py" \
  --phase-ensemble-summary-csv "${PHASE_WORKDIR}/outputs/phase_ensemble_summary.csv" \
  --output-dir "${UNCERTAINTY_DIR}"

echo "[run] results pack (${PROFILE_TAG})"
python3 "${REPO_ROOT}/scripts/generate_modes4d_results_pack.py" \
  --lundquist-summary-csv "${PROD_WORKDIR}/outputs/lundquist_scan_summary.csv" \
  --lundquist-report "${PROD_WORKDIR}/outputs/production_report.md" \
  --lundquist-primary-plot "${PROD_WORKDIR}/outputs/plots/lundquist_primary_channels.png" \
  --lundquist-subscale-plot "${PROD_WORKDIR}/outputs/plots/lundquist_subscale_channels.png" \
  --topology-report "${PROD_WORKDIR}/outputs/topology/topology_comparison_report.md" \
  --topology-plot "${PROD_WORKDIR}/outputs/topology/topology_timeseries.png" \
  --topology-status-csv "${PROD_WORKDIR}/outputs/topology/topology_status.csv" \
  --topology-gates-csv "${PROD_WORKDIR}/outputs/topology/topology_gate_details.csv" \
  --ablation-summary-csv "${ABL_WORKDIR}/outputs/ablation_summary.csv" \
  --ablation-report "${ABL_WORKDIR}/outputs/ablation_report.md" \
  --ablation-plot "${ABL_WORKDIR}/outputs/plots/ablation_channels.png" \
  --phase-ensemble-summary-csv "${PHASE_WORKDIR}/outputs/phase_ensemble_summary.csv" \
  --phase-ensemble-report "${PHASE_WORKDIR}/outputs/phase_ensemble_report.md" \
  --output-dir "${RESULTS_DIR}"

if [[ "${RUN_MPI_SMOKE}" -eq 1 ]]; then
  echo "[run] mpi+hdf5 smoke (${PROFILE_TAG})"
  "${REPO_ROOT}/scripts/run_harris_full_mpi_hdf5.sh" \
    --binary "${MPI_BINARY}" \
    --workdir "${MPI_WORKDIR}" \
    --ranks "${MPI_RANKS}" \
    --tlim "${MPI_TLIM}" \
    --nlim "${MPI_NLIM}" \
    --output-dt "${MPI_OUTPUT_DT}" \
    --mpirun "${MPIRUN_BIN}" \
    -- \
    "diffusion/resistivity=ohmic" \
    "diffusion/resistivity_coeff=fixed" \
    "diffusion/integrator=rkl2" \
    "diffusion/rkl2_max_dt_ratio=100.0" \
    "diffusion/ohm_diff_coeff_code=1.118033988750e-03"
fi

echo "profile_tag,${PROFILE_TAG}"
echo "production_report,${PROD_WORKDIR}/outputs/production_report.md"
echo "phase_ensemble_report,${PHASE_WORKDIR}/outputs/phase_ensemble_report.md"
echo "ablation_report,${ABL_WORKDIR}/outputs/ablation_report.md"
echo "uncertainty_bands_report,${UNCERTAINTY_DIR}/uncertainty_bands_report.md"
echo "uncertainty_bands_plot,${UNCERTAINTY_DIR}/uncertainty_bands.png"
echo "results_pack_report,${RESULTS_DIR}/modes4d_results_pack.md"
echo "results_pack_status_csv,${RESULTS_DIR}/modes4d_results_status.csv"
