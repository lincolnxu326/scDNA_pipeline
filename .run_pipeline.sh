#!/bin/bash
# SLURM Submission Script for Single-Cell DNA Sequencing Pipeline
# Directory structure: ./scDNA_pipeline/ for pipeline, ../output/plate*/ for data

set -e

# Default values
MODE="run"
CONFIG="config.yaml"
JOBS=50
EMAIL=""
CONDA_ENV="/camp/project/tracerX/working/CRENAL/kl_scripts/software/anaconda3/envs/snakemake_env"
CONDA_PREFIX="/camp/project/tracerX/working/CRENAL/kl_scripts/software/anaconda3/"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m'

# Function to print colored messages
print_message() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# Header
print_header() {
    echo ""
    print_message "$PURPLE" "╔═══════════════════════════════════════════════════════╗"
    print_message "$PURPLE" "║     Single-Cell DNA Sequencing Pipeline v2.0         ║"
    print_message "$PURPLE" "║           SLURM Cluster Execution Manager            ║"
    print_message "$PURPLE" "╚═══════════════════════════════════════════════════════╝"
    echo ""
}

# Help function
show_help() {
    print_header
    cat << EOF
${GREEN}USAGE:${NC} $0 [MODE] [OPTIONS]

${BLUE}PIPELINE MODES:${NC}
    ${GREEN}-r, --run${NC}           Run complete pipeline (default)
    ${GREEN}-b, --blacklist${NC}     Generate blacklist diagnosis plots only
    ${GREEN}-a, --aneufinder${NC}    Run AneuFinder analysis (after blacklist review)
    ${GREEN}-q, --qc${NC}            Generate QC reports only (MultiQC)
    ${GREEN}-d, --dry-run${NC}       Perform a dry run
    ${GREEN}-u, --unlock${NC}        Unlock working directory
    ${GREEN}-c, --clean${NC}         Clean intermediate files

${BLUE}OPTIONS:${NC}
    ${GREEN}--config FILE${NC}       Config file (default: config.yaml)
    ${GREEN}--jobs NUM${NC}          Number of SLURM jobs (default: 50)
    ${GREEN}--email EMAIL${NC}       Email for notifications
    ${GREEN}--conda-env PATH${NC}    Conda environment path
    ${GREEN}--conda-prefix PATH${NC} Conda prefix for environments

${BLUE}SLURM OPTIONS:${NC}
    ${GREEN}--partition NAME${NC}    SLURM partition (default: ncpu)
    ${GREEN}--account NAME${NC}      SLURM account
    ${GREEN}--time HOURS${NC}        Runtime in hours (default: 24)
    ${GREEN}--mem GB${NC}            Memory in GB (default: 128)

${BLUE}EXAMPLES:${NC}
    # Generate blacklist (FIRST STEP)
    $0 --blacklist --email user@institute.com

    # Review blacklist plots
    evince ../output/mappability/mappability_plot.pdf

    # Run full pipeline
    $0 --run --email user@institute.com

    # Generate QC reports only
    $0 --qc

    # Dry run
    $0 --dry-run

${BLUE}DIRECTORY STRUCTURE:${NC}
    scDNA_pipeline/         # Pipeline files (current directory)
    ../output/              # Data directory
        ├── plate1/         # Plate directory
        │   ├── plate1_R1.fastq.gz
        │   ├── plate1_R2.fastq.gz
        │   └── [outputs]
        └── plate2/

EOF
}

# Parse arguments
PARTITION="ncpu"
ACCOUNT=""
TIME="24"
MEM="128"

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -r|--run)
            MODE="run"
            shift
            ;;
        -b|--blacklist)
            MODE="blacklist"
            shift
            ;;
        -a|--aneufinder)
            MODE="aneufinder"
            shift
            ;;
        -q|--qc)
            MODE="qc"
            shift
            ;;
        -d|--dry-run)
            MODE="dry-run"
            shift
            ;;
        -u|--unlock)
            MODE="unlock"
            shift
            ;;
        -c|--clean)
            MODE="clean"
            shift
            ;;
        --config)
            CONFIG="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
            ;;
        --email)
            EMAIL="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --conda-prefix)
            CONDA_PREFIX="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --account)
            ACCOUNT="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        *)
            print_message "$RED" "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check directory structure
check_directory_structure() {
    print_message "$BLUE" "Checking directory structure..."
    
    # Check if we're in scDNA_pipeline directory
    if [ ! -f "Snakefile" ]; then
        print_message "$RED" "Error: Not in scDNA_pipeline directory!"
        print_message "$YELLOW" "Please run this script from the scDNA_pipeline directory."
        exit 1
    fi
    
    # Check for output directory
    if [ ! -d "../output" ]; then
        print_message "$RED" "Error: ../output directory not found!"
        print_message "$YELLOW" "Please create the output directory and add plate folders."
        exit 1
    fi
    
    # Check for plate directories
    PLATES=$(ls -d ../output/plate* 2>/dev/null | wc -l)
    if [ "$PLATES" -eq 0 ]; then
        print_message "$RED" "Error: No plate directories found in ../output/"
        print_message "$YELLOW" "Please create plate directories (e.g., ../output/plate1/)"
        exit 1
    fi
    
    print_message "$GREEN" "✓ Found $PLATES plate directories"
    
    # Check for FASTQ files
    for plate_dir in ../output/plate*/; do
        plate=$(basename $plate_dir)
        if [ ! -f "$plate_dir/${plate}_R1.fastq.gz" ] || [ ! -f "$plate_dir/${plate}_R2.fastq.gz" ]; then
            print_message "$YELLOW" "Warning: FASTQ files not found in $plate_dir"
            print_message "$YELLOW" "Expected: ${plate}_R1.fastq.gz and ${plate}_R2.fastq.gz"
        fi
    done
}

# Create logs directory
mkdir -p logs

# Function to create SLURM script
create_slurm_script() {
    local mode=$1
    local script_name="submit_${mode}_$$.sh"
    
    cat > $script_name << EOF
#!/bin/bash
#SBATCH --job-name=scDNA_${mode}
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --partition=${PARTITION}
#SBATCH --mem=${MEM}G
#SBATCH --time=${TIME}:00:00
EOF

    if [ -n "$ACCOUNT" ]; then
        echo "#SBATCH --account=${ACCOUNT}" >> $script_name
    fi
    
    if [ -n "$EMAIL" ]; then
        cat >> $script_name << EOF
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${EMAIL}
EOF
    fi
    
    cat >> $script_name << EOF

source ~/.bashrc

# Load modules if needed
# module load Bowtie2/2.5.1-GCC-12.3.0
# module load SAMtools/1.18-GCC-12.3.0

# Activate conda environment
source activate ${CONDA_ENV}

echo "═══════════════════════════════════════════════════════"
echo "Single-Cell DNA Pipeline - Mode: ${mode}"
echo "═══════════════════════════════════════════════════════"
echo "Time: \$(date)"
echo "Node: \$(hostname)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "Working directory: \$(pwd)"
echo "═══════════════════════════════════════════════════════"
echo ""

EOF
    
    echo $script_name
}

# Add snakemake command
add_snakemake_command() {
    local script=$1
    local target=$2
    local extra_params=$3
    
    cat >> $script << EOF
snakemake \\
  --configfile ${CONFIG} \\
  --executor slurm \\
  --jobs ${JOBS} \\
  --use-conda \\
  --conda-prefix "${CONDA_PREFIX}" \\
  --default-resources \\
    mem_mb=12000 \\
    cpus=4 \\
    slurm_partition=${PARTITION} \\
    runtime=240 \\
    slurm_output="logs/{rule}.{jobid}.out" \\
    slurm_error="logs/{rule}.{jobid}.err" \\
EOF

    if [ -n "$ACCOUNT" ]; then
        echo "    slurm_account=${ACCOUNT} \\" >> $script
    fi
    
    cat >> $script << EOF
  --latency-wait 60 \\
  --rerun-incomplete \\
  --keep-going \\
  --printshellcmds \\
  --show-failed-logs \\
  ${target} ${extra_params}

echo ""
echo "═══════════════════════════════════════════════════════"
echo "Pipeline completed at: \$(date)"
echo "═══════════════════════════════════════════════════════"
EOF
}

# Main execution
print_header

case $MODE in
    run)
        check_directory_structure
        print_message "$GREEN" "🚀 Preparing to run complete pipeline..."
        
        script=$(create_slurm_script "full")
        add_snakemake_command $script "--until all" ""
        
        print_message "$YELLOW" "📤 Submitting job to SLURM..."
        job_id=$(sbatch --parsable $script)
        print_message "$GREEN" "✅ Job submitted with ID: $job_id"
        print_message "$BLUE" "📊 Monitor progress: squeue -j $job_id"
        print_message "$BLUE" "📄 Check logs: tail -f logs/scDNA_full_${job_id}.out"
        ;;
        
    blacklist)
        check_directory_structure
        print_message "$GREEN" "🔬 Generating blacklist diagnosis plots..."
        
        script=$(create_slurm_script "blacklist")
        add_snakemake_command $script "--until generate_blacklist_plots" ""
        
        print_message "$YELLOW" "📤 Submitting blacklist job..."
        job_id=$(sbatch --parsable $script)
        print_message "$GREEN" "✅ Job submitted with ID: $job_id"
        print_message "$YELLOW" "⚠️  Review plots before running AneuFinder:"
        print_message "$BLUE" "   📊 Plots: ../output/mappability/mappability_plot.pdf"
        print_message "$BLUE" "   📝 Stats: ../output/mappability/blacklist_diagnosis.txt"
        ;;
        
    aneufinder)
        check_directory_structure
        print_message "$GREEN" "🧬 Running AneuFinder analysis..."
        
        # Check blacklist exists
        if [ ! -f "../output/mappability/blacklist.bed.gz" ]; then
            print_message "$RED" "❌ Error: Blacklist not found!"
            print_message "$YELLOW" "Run --blacklist first to generate blacklist."
            exit 1
        fi
        
        script=$(create_slurm_script "aneufinder")
        add_snakemake_command $script "--until all" "--forcerun run_aneufinder"
        
        print_message "$YELLOW" "📤 Submitting AneuFinder job..."
        job_id=$(sbatch --parsable $script)
        print_message "$GREEN" "✅ Job submitted with ID: $job_id"
        ;;
        
    qc)
        check_directory_structure
        print_message "$GREEN" "📊 Generating QC reports..."
        
        script=$(create_slurm_script "qc")
        add_snakemake_command $script "--until all_qc" ""
        
        print_message "$YELLOW" "📤 Submitting QC job..."
        job_id=$(sbatch --parsable $script)
        print_message "$GREEN" "✅ Job submitted with ID: $job_id"
        print_message "$BLUE" "📊 Reports will be in:"
        print_message "$BLUE" "   - Plate reports: ../output/plate*/multiqc/"
        print_message "$BLUE" "   - Combined: ../output/combined_qc_summary.html"
        ;;
        
    dry-run)
        check_directory_structure
        print_message "$GREEN" "🔍 Performing dry run..."
        source activate ${CONDA_ENV}
        
        snakemake \
            --configfile ${CONFIG} \
            --dry-run \
            --printshellcmds \
            --reason \
            --summary
        
        if [ $? -eq 0 ]; then
            print_message "$GREEN" "✅ Dry run successful!"
            print_message "$YELLOW" "Pipeline is ready to run."
        else
            print_message "$RED" "❌ Dry run failed. Check configuration."
        fi
        ;;
        
    unlock)
        print_message "$YELLOW" "🔓 Unlocking working directory..."
        source activate ${CONDA_ENV}
        snakemake --unlock --configfile ${CONFIG}
        print_message "$GREEN" "✅ Directory unlocked!"
        ;;
        
    clean)
        print_message "$YELLOW" "🧹 Cleaning intermediate files..."
        read -p "Remove all intermediate files? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            source activate ${CONDA_ENV}
            snakemake --configfile ${CONFIG} clean
            print_message "$GREEN" "✅ Cleanup complete!"
        else
            print_message "$YELLOW" "Cleanup cancelled."
        fi
        ;;
        
    *)
        print_message "$RED" "Unknown mode: $MODE"
        show_help
        exit 1
        ;;
esac

# Clean up submission scripts
if [ "$DEBUG" != "1" ] && [ "$MODE" != "dry-run" ]; then
    rm -f submit_*.sh
fi

print_message "$PURPLE" "═══════════════════════════════════════════════════════"
