#!/bin/bash
set -euo pipefail

module purge
ml Graphviz/8.1.0-GCCcore-12.3.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}"
cd "${PIPELINE_DIR}"

CONDA_ENV="/camp/project/tracerX/working/CRENAL/kl_scripts/software/anaconda3/envs/snakemake_scDNA"
CONDA_ROOT="/camp/project/tracerX/working/CRENAL/kl_scripts/software/anaconda3"
PLATE_DIR="/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/plate10"
SNAKEMAKE_CONDA_PREFIX="${PIPELINE_DIR}/.snakemake/conda"
CONDA_PKGS_DIRS="${PIPELINE_DIR}/.conda/pkgs"

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
mkdir -p "${SNAKEMAKE_CONDA_PREFIX}" "${CONDA_PKGS_DIRS}"
export CONDA_PKGS_DIRS

snakemake \
  --snakefile "${PIPELINE_DIR}/Snakefile" \
  --configfile "${PIPELINE_DIR}/config.yaml" \
  --config plate_dir="${PLATE_DIR}" \
  --use-conda   --conda-prefix "${SNAKEMAKE_CONDA_PREFIX}" \
  --profile /nemo/project/proj-tracerX/working/PIPELINES/Snakemake   --dag > "${PIPELINE_DIR}/rulegraph_raw.dot"

# dot -Tpdf "${PIPELINE_DIR}/rulegraph_raw.dot" -o "${PIPELINE_DIR}/rulegraph_plate10.pdf"
