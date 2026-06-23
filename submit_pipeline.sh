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

# Mode to run (choose one). The pipeline is a two-pass workflow with a human review
# step in the middle. The NORMAL operator sequence is just two runs:
#
#   1. MODE=pre_review    -> preprocessing + first-pass AneuFinder + review HTML, then STOP
#   2. (human) open qc_review/review.html, mark decisions, save <PLATE_DIR>/qc_decisions.csv
#      (the report has a "Copy terminal save command" button that writes it for you)
#   3. MODE=post_review   -> validate -> second-pass AneuFinder (PASS wells) -> cn_review.html
#
# Primary modes:
# - "pre_review"         : everything up to and INCLUDING the first-pass review HTML,
#                          then stops for human review (this is the old "full")
# - "post_review"        : ALL post-review steps in one go -> qc_review/cn_review.html
#                          (requires a saved <PLATE_DIR>/qc_decisions.csv)
#
# Other / lower-level modes (useful for debugging or partial runs):
# - "preprocessing"      : Demux, dimer filtering, alignment, dedup, QC summaries
# - "qc"                 : preprocessing + FastQC + MultiQC
# - "aneufinder_first"   : blacklist + first AneuFinder pass (all wells)
# - "review"             : (re)build the first-pass review HTML only
# - "validate_review"    : validate a saved <PLATE_DIR>/qc_decisions.csv
# - "aneufinder_reviewed": second AneuFinder pass on PASS wells only (needs the CSV)
# - "cn_review"          : final copy-number review viewer over the second pass
# - "blacklist"          : blacklist diagnosis plots only
# Defaults can be overridden at submit time, e.g.
#   MODE=pre_review PLATE_DIR=/path/to/plate sbatch submit_pipeline.sh
MODE="${MODE:-pre_review}" # EDIT THIS

# Plate directory EDIT THIS
PLATE_DIR="${PLATE_DIR:-/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/384_well/plate17_umi/plate17_4}"

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
  preprocessing)
    TARGET="all_preprocessing"
    echo "Running preprocessing (demux -> dedup)..."
    ;;
  pre_review|full)
    TARGET="all"
    echo "Running pipeline through the first-pass review report (stops for human review)..."
    ;;
  post_review)
    TARGET="post_review"
    echo "Running all post-review steps (validate -> second-pass AneuFinder -> CN viewer)..."
    ;;
  qc)
    TARGET="all_qc"
    echo "Generating QC reports..."
    ;;
  aneufinder_first|aneufinder)
    TARGET="all_aneufinder_first"
    EXTRA_FLAGS=""
    echo "Running first-pass AneuFinder analysis..."
    ;;
  review)
    TARGET="qc_review_report"
    echo "Building first-pass static HTML QC review report..."
    ;;
  validate_review)
    TARGET="validate_qc_decisions"
    echo "Validating ${PLATE_DIR}/qc_decisions.csv..."
    ;;
  aneufinder_reviewed)
    TARGET="all_aneufinder_reviewed"
    echo "Running second (post-review) AneuFinder pass on PASS wells..."
    ;;
  cn_review)
    TARGET="cn_review_report"
    echo "Building final copy-number review viewer..."
    ;;
  *)
    echo "ERROR: Unknown mode: ${MODE}"
    exit 1
    ;;
esac

# Post-review modes consume the human-saved decisions file. Fail early with a clear
# message rather than a deep Snakemake MissingInputException if it is not there yet.
case "${MODE}" in
  post_review|validate_review|aneufinder_reviewed|cn_review)
    DECISIONS_CSV="${PLATE_DIR}/qc_decisions.csv"
    if [ ! -f "${DECISIONS_CSV}" ]; then
      echo "ERROR: MODE='${MODE}' needs the human decisions file, which was not found:"
      echo "       ${DECISIONS_CSV}"
      echo "Run MODE=pre_review first, open ${PLATE_DIR}/qc_review/review.html,"
      echo "make decisions, then use the report's 'Copy terminal save command' button"
      echo "(or save the downloaded CSV) to write ${DECISIONS_CSV}."
      exit 1
    fi
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
    pre_review|full|qc)
      echo "QC report at: ${PLATE_DIR}/multiqc/multiqc_report.html"
      ;;
  esac
  case "${MODE}" in
    pre_review|full|review)
      echo "First-pass review report at: ${PLATE_DIR}/qc_review/review.html"
      echo "  Open it, make decisions, then use its 'Copy terminal save command' button"
      echo "  (or save the downloaded CSV) to write: ${PLATE_DIR}/qc_decisions.csv"
      echo "  Next: MODE=post_review (runs validate -> second pass -> final CN viewer)"
      ;;
    post_review|cn_review)
      echo "Final CN review viewer at: ${PLATE_DIR}/qc_review/cn_review.html"
      echo "Second-pass AneuFinder at: ${PLATE_DIR}/aneufinder_reviewed/"
      ;;
    aneufinder_reviewed)
      echo "Second-pass AneuFinder at: ${PLATE_DIR}/aneufinder_reviewed/"
      echo "  Build the final viewer with: MODE=cn_review (or MODE=post_review)"
      ;;
  esac
else
  echo "Pipeline failed with exit code: ${EXIT_CODE}"
  echo "Check logs in: ${PIPELINE_DIR}/logs/ and ${PLATE_DIR}/logs/"
fi
echo "Completed at: $(date)"
echo "============================================================"

exit ${EXIT_CODE}
