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

# submit_pipeline.sh writes these as `VAR="${VAR:-default}"`, so unwrap the default
# out of the parameter expansion rather than reporting it literally.
read_default() {
    grep "^$1=" "${SUBMIT_SCRIPT}" | head -1 | cut -d'"' -f2 \
        | sed -E "s/^\\\$\{$1:-(.*)\}$/\1/"
}

PLATE_DIR="${PLATE_DIR:-$(read_default PLATE_DIR)}"
MODE="${MODE:-$(read_default MODE)}"
PLATE_FORMAT="${PLATE_FORMAT:-$(read_default PLATE_FORMAT)}"
PLATE_FORMAT="${PLATE_FORMAT:-96}"
CONDA_PREFIX=$(grep '^SNAKEMAKE_CONDA_PREFIX=' "${SUBMIT_SCRIPT}" | cut -d'"' -f2)
CONDA_PKGS=$(grep '^CONDA_PKGS_DIRS=' "${SUBMIT_SCRIPT}" | cut -d'"' -f2)

# Every place a (sub)plate's own inputs are checked. In 384 mode this is each
# subplate directory; in 96 mode it is just the plate directory itself.
UNITS=()

if [ -z "${PLATE_DIR}" ]; then
    echo -e "${YELLOW}! Could not read PLATE_DIR from submit_pipeline.sh${NC}"
else
    echo -e "${GREEN}OK plate directory set to: ${PLATE_DIR}${NC}"
    if [ -d "${PLATE_DIR}" ]; then
        echo -e "${GREEN}OK plate directory exists${NC}"
        PLATE_NAME=$(basename "${PLATE_DIR}")

        if [ "${PLATE_FORMAT}" = "384" ]; then
            while IFS= read -r d; do
                [ -n "$d" ] && UNITS+=("$d")
            done < <(ls -d "${PLATE_DIR}/${PLATE_NAME}"_[0-9]* 2>/dev/null | sort -V)
            if [ ${#UNITS[@]} -eq 0 ]; then
                echo -e "${RED}x PLATE_FORMAT=384 but no ${PLATE_NAME}_<n>/ subplate dirs found${NC}"
            else
                echo -e "${GREEN}OK ${#UNITS[@]} subplate director(ies): $(for u in "${UNITS[@]}"; do basename "$u"; done | tr '\n' ' ')${NC}"
            fi
        else
            # The likeliest operator mistake: a 384 plate dir left at PLATE_FORMAT=96.
            if [ ! -f "${PLATE_DIR}/${PLATE_NAME}_R1.fastq.gz" ] && \
               ls -d "${PLATE_DIR}/${PLATE_NAME}"_[0-9]* >/dev/null 2>&1; then
                echo -e "${RED}x ${PLATE_DIR} looks like a 384 plate directory (subplate dirs, no plate-level FASTQ).${NC}"
                echo -e "${RED}  Set PLATE_FORMAT=384, or point PLATE_DIR at a single subplate.${NC}"
            fi
            UNITS=("${PLATE_DIR}")
        fi

        for u in "${UNITS[@]}"; do
            un=$(basename "$u")
            for r in R1 R2; do
                if [ -f "${u}/${un}_${r}.fastq.gz" ]; then
                    echo -e "${GREEN}OK ${un} ${r} FASTQ found${NC}"
                else
                    echo -e "${RED}x ${un} ${r} FASTQ not found: ${u}/${un}_${r}.fastq.gz${NC}"
                fi
            done
        done
    else
        echo -e "${RED}x Plate directory not found: ${PLATE_DIR}${NC}"
    fi
fi

# Barcodes: shared resources win, otherwise each (sub)plate needs its own.
if [ -f "${PIPELINE_DIR}/resources/barcodes.tsv" ] || \
   [ -f "${PIPELINE_DIR}/resources/barcodes/barcodes.tsv" ]; then
    echo -e "${GREEN}OK barcodes file found (shared resources)${NC}"
elif [ ${#UNITS[@]} -gt 0 ]; then
    missing=0
    for u in "${UNITS[@]}"; do
        [ -f "${u}/barcodes.tsv" ] || { echo -e "${YELLOW}! barcodes.tsv missing: ${u}/barcodes.tsv${NC}"; missing=1; }
    done
    [ ${missing} -eq 0 ] && echo -e "${GREEN}OK barcodes.tsv present for all ${#UNITS[@]} (sub)plate(s)${NC}"
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
    pre_review|full)
        echo "  -> Preprocessing + first-pass AneuFinder + review HTML, then stops for human review"
        if [ -f "${PLATE_DIR}/mappability/blacklist.bed.gz" ]; then
            echo -e "${GREEN}OK blacklist exists${NC}"
        else
            echo -e "${YELLOW}! blacklist not found yet; consider running MODE=blacklist first${NC}"
        fi
        ;;
    preprocessing)
        echo "  -> Demux, dimer filtering, alignment, dedup, QC summaries"
        ;;
    qc)
        echo "  -> Will generate QC reports only"
        ;;
    aneufinder|aneufinder_first)
        echo "  -> Will run the first AneuFinder pass only"
        ;;
    review)
        echo "  -> Will (re)build the first-pass review HTML only"
        ;;
    post_review|validate_review|aneufinder_reviewed|cn_review)
        echo "  -> Post-review step; needs ${PLATE_DIR}/qc_decisions.csv"
        if [ -f "${PLATE_DIR}/qc_decisions.csv" ]; then
            echo -e "${GREEN}OK qc_decisions.csv found${NC}"
        else
            echo -e "${RED}x qc_decisions.csv not found: ${PLATE_DIR}/qc_decisions.csv${NC}"
        fi
        ;;
    cellenone)
        echo "  -> Will ingest + render CellenONE cell images only"
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
echo "Plate format: ${PLATE_FORMAT}-well"
echo "Mode: ${MODE}"
echo ""
echo "Submit with: sbatch submit_pipeline.sh"
echo "  (384 plate: PLATE_FORMAT=384 PLATE_DIR=/path/to/plate sbatch submit_pipeline.sh)"
echo "Monitor with: tail -f logs/pipeline_*.out"
echo ""
