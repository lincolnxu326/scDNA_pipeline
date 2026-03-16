#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}"
cd "${PIPELINE_DIR}"

PLATE_DIR="/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/plate9"
SNAKEMAKE_CONDA_PREFIX="${PIPELINE_DIR}/.snakemake/conda"
CONDA_PKGS_DIRS="${PIPELINE_DIR}/.conda/pkgs"
mkdir -p "${SNAKEMAKE_CONDA_PREFIX}" "${CONDA_PKGS_DIRS}"
export CONDA_PKGS_DIRS

snakemake \
  --snakefile "${PIPELINE_DIR}/Snakefile" \
  --configfile "${PIPELINE_DIR}/config.yaml" \
  --config plate_dir="${PLATE_DIR}" \
  --executor slurm   --jobs 50   --use-conda   --profile /nemo/project/proj-tracerX/working/PIPELINES/Snakemake   --conda-prefix "${SNAKEMAKE_CONDA_PREFIX}" \
  --forcerun align   all_preprocessing
