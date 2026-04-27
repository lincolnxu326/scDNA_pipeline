#!/bin/bash
# Pre-flight checks for Single-Cell DNA Pipeline
# Run from anywhere with: bash pre_check.sh

set -euo pipefail

RED='[0;31m'
GREEN='[0;32m'
YELLOW='[1;33m'
NC='[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}"
SUBMIT_SCRIPT="${PIPELINE_DIR}/submit_pipeline.sh"
cd "${PIPELINE_DIR}"

echo ""
echo "========================================="
echo "Single-Cell DNA Pipeline - Pre-flight Check"
echo "========================================="
echo ""

if [ ! -f "${PIPELINE_DIR}/Snakefile" ]; then
    echo -e "${RED}x Snakefile not found${NC}"
    exit 1
fi
echo -e "${GREEN}OK pipeline directory found${NC}"

if [ ! -f "${PIPELINE_DIR}/config.yaml" ]; then
    echo -e "${RED}x config.yaml not found${NC}"
    exit 1
fi
echo -e "${GREEN}OK config file found${NC}"

PLATE_DIR=$(grep '^PLATE_DIR=' "${SUBMIT_SCRIPT}" | cut -d'"' -f2)
MODE=$(grep '^MODE=' "${SUBMIT_SCRIPT}" | cut -d'"' -f2)
CONDA_PREFIX=$(grep '^SNAKEMAKE_CONDA_PREFIX=' "${SUBMIT_SCRIPT}" | cut -d'"' -f2)
CONDA_PKGS=$(grep '^CONDA_PKGS_DIRS=' "${SUBMIT_SCRIPT}" | cut -d'"' -f2)

if [ -z "${PLATE_DIR}" ]; then
    echo -e "${YELLOW}! Could not read PLATE_DIR from submit_pipeline.sh${NC}"
else
    echo -e "${GREEN}OK plate directory set to: ${PLATE_DIR}${NC}"
    if [ -d "${PLATE_DIR}" ]; then
        echo -e "${GREEN}OK plate directory exists${NC}"
        PLATE_NAME=$(basename "${PLATE_DIR}")
        if [ -f "${PLATE_DIR}/${PLATE_NAME}_R1.fastq.gz" ]; then
            echo -e "${GREEN}OK R1 FASTQ found${NC}"
        else
            echo -e "${RED}x R1 FASTQ not found: ${PLATE_DIR}/${PLATE_NAME}_R1.fastq.gz${NC}"
        fi
        if [ -f "${PLATE_DIR}/${PLATE_NAME}_R2.fastq.gz" ]; then
            echo -e "${GREEN}OK R2 FASTQ found${NC}"
        else
            echo -e "${RED}x R2 FASTQ not found: ${PLATE_DIR}/${PLATE_NAME}_R2.fastq.gz${NC}"
        fi
    else
        echo -e "${RED}x Plate directory not found: ${PLATE_DIR}${NC}"
    fi
fi

if [ -f "${PIPELINE_DIR}/resources/barcodes.tsv" ] || \
   [ -f "${PIPELINE_DIR}/resources/barcodes/barcodes.tsv" ] || \
   [ -n "${PLATE_DIR:-}" ] && [ -f "${PLATE_DIR}/barcodes.tsv" ]; then
    echo -e "${GREEN}OK barcodes file found${NC}"
else
    echo -e "${YELLOW}! barcodes.tsv not found in shared resources or plate directory${NC}"
fi

REF_GENOME=$(grep 'fasta:' "${PIPELINE_DIR}/config.yaml" | head -1 | awk '{print $2}' | tr -d '"')
if [ -f "${REF_GENOME}" ]; then
    echo -e "${GREEN}OK reference genome found${NC}"
else
    echo -e "${RED}x Reference genome not found: ${REF_GENOME}${NC}"
fi

echo ""
echo "Mode set to: ${YELLOW}${MODE:-unknown}${NC}"
case "${MODE:-}" in
    blacklist)
        echo "  -> Will generate blacklist diagnosis plots"
        ;;
    full)
        if [ -f "${PLATE_DIR}/mappability/blacklist.bed.gz" ]; then
            echo -e "${GREEN}OK blacklist exists${NC}"
        else
            echo -e "${YELLOW}! blacklist not found yet; consider running MODE="blacklist" first${NC}"
        fi
        ;;
    qc)
        echo "  -> Will generate QC reports only"
        ;;
    aneufinder)
        echo "  -> Will run AneuFinder only"
        ;;
    *)
        echo -e "${YELLOW}! MODE is unset or unrecognized${NC}"
        ;;
esac

EMAIL=$(grep 'mail-user=' "${SUBMIT_SCRIPT}" | head -1 | cut -d'=' -f2)
echo -e "${GREEN}OK email notifications to: ${EMAIL}${NC}"
echo -e "${GREEN}OK shared conda prefix: ${CONDA_PREFIX}${NC}"
echo -e "${GREEN}OK shared conda package cache: ${CONDA_PKGS}${NC}"

echo ""
echo "========================================="
echo "Summary"
echo "========================================="
echo "Pipeline directory: ${PIPELINE_DIR}"
echo "Plate directory: ${PLATE_DIR}"
echo "Mode: ${MODE}"
echo ""
echo "Submit with: sbatch submit_pipeline.sh"
echo "Monitor with: tail -f logs/pipeline_*.out"
echo ""
