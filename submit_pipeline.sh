#!/bin/bash
#SBATCH --job-name=scDNA_pipeline_plate17_4 # EDIT THIS
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --ntasks=2
#SBATCH --partition=ncpu
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lincoln.xu@crick.ac.uk # EDIT THIS

# Single-Cell DNA Pipeline SLURM Submission Script
# Edit the USER SETTINGS section below, then submit with:
#   sbatch submit_pipeline.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Prefer the real script location when it contains the pipeline. Fall back to
# SLURM_SUBMIT_DIR only when the script is executed from a spool copy.
if [[ -f "${SCRIPT_DIR}/Snakefile" && -f "${SCRIPT_DIR}/config.yaml" ]]; then
  PIPELINE_DIR="${SCRIPT_DIR}"
else
  PIPELINE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
fi
cd "${PIPELINE_DIR}"

# ============================================================================
# USER SETTINGS - EDIT THESE
# ============================================================================

# Mode to run (choose one):
# - "full"       : Run complete pipeline
# - "preprocessing"  : Demux, dedup, dimer removal, align
# - "qc"         : Generate QC reports (preprocessing + FastQC + MultiQC)
# - "aneufinder" : Run blacklist, GC template generation, and AneuFinder
MODE="aneufinder" # EDIT THIS

# Plate directory EDIT THIS
PLATE_DIR="/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/384_well/plate17/plate17_4"

# Snakemake launcher environment
CONDA_ENV="/camp/project/tracerX/working/CRENAL/kl_scripts/software/anaconda3/envs/snakemake_scDNA"
CONDA_ROOT="/camp/project/tracerX/working/CRENAL/kl_scripts/software/anaconda3"

# Number of jobs to submit
JOBS=50

# Keep all rule environments and conda package downloads shared at the pipeline root.
SNAKEMAKE_CONDA_PREFIX="${PIPELINE_DIR}/.snakemake/conda"
CONDA_PKGS_DIRS="${PIPELINE_DIR}/.conda/pkgs"
export SNAKEMAKE_CONDA_PREFIX
export CONDA_PKGS_DIRS
mkdir -p "${SNAKEMAKE_CONDA_PREFIX}" "${CONDA_PKGS_DIRS}" "${PIPELINE_DIR}/logs"

# Load snakemake launcher environment
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

EXTRA_FLAGS=""
case "${MODE}" in
  blacklist)
    TARGET="generate_blacklist_plots"
    echo "Generating blacklist diagnosis plots..."
    ;;
  full)
    TARGET="all"
    echo "Running complete pipeline..."
    ;;
  qc)
    TARGET="all_qc"
    echo "Generating QC reports..."
    ;;
  aneufinder)
    TARGET="all_aneufinder"
    EXTRA_FLAGS=""
    echo "Running AneuFinder analysis..."
    ;;
  *)
    echo "ERROR: Unknown mode: ${MODE}"
    exit 1
    ;;
esac

echo "============================================================"
echo "Single-Cell DNA Pipeline - SLURM Execution"
echo "============================================================"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Pipeline root: ${PIPELINE_DIR}"
echo "Mode: ${MODE}"
echo "Plate: ${PLATE_DIR}"
echo "Shared conda env prefix: ${SNAKEMAKE_CONDA_PREFIX}"
echo "Shared conda package cache: ${CONDA_PKGS_DIRS}"
echo "============================================================"

snakemake \
  --snakefile "${PIPELINE_DIR}/Snakefile" \
  --configfile "${PIPELINE_DIR}/config.yaml" \
  --config plate_dir="${PLATE_DIR}" \
  --executor slurm \
  --jobs "${JOBS}" \
  --use-conda \
  --profile /nemo/project/proj-tracerX/working/PIPELINES/Snakemake \
  --conda-prefix "${SNAKEMAKE_CONDA_PREFIX}" \
  ${TARGET} ${EXTRA_FLAGS}

EXIT_CODE=$?

echo ""
echo "============================================================"
if [ ${EXIT_CODE} -eq 0 ]; then
  echo "Pipeline completed successfully!"
  echo "Results in: ${PLATE_DIR}/"
  case "${MODE}" in
    blacklist)
      echo "Review plots at: ${PLATE_DIR}/mappability/mappability_plot.pdf"
      ;;
    full|qc)
      echo "QC report at: ${PLATE_DIR}/multiqc/multiqc_report.html"
      ;;
  esac
else
  echo "Pipeline failed with exit code: ${EXIT_CODE}"
  echo "Check logs in: ${PIPELINE_DIR}/logs/ and ${PLATE_DIR}/logs/"
fi
echo "Completed at: $(date)"
echo "============================================================"

exit ${EXIT_CODE}
