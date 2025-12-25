#!/bin/bash

# ==============================================================================
#                 METATRIMX: NEURAL ADAPTIVE ENGINE v3.0
#      Powered by Random Forest Optimization & Deep Neural Networks
# ==============================================================================
# Author: SUBRAMANIAM VIJAYAKUMAR
# ==============================================================================
#
# SYSTEM LOGIC:
# This pipeline uses a Hybrid AI Architecture:
#   1. ADAPTIVE CORE (Mode 1): A Random Forest Regressor that predicts optimal
#      trimming parameters based on raw data entropy.
#   2. NEURAL SENTINEL (Mode 2): A Deep Learning Classifier (TensorFlow) that
#      detects and removes PCR artifacts from final sequences.
#
# ==============================================================================

# --- COLOR DEFINITIONS FOR UI ---
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# --- FUNCTION: PROFESSIONAL BANNER ---
print_banner() {
    clear
    echo -e "${CYAN}${BOLD}"
    echo "  __  __      _         _______   _            __  __ "
    echo " |  \/  |    | |       |__   __| (_)          \ \/ / "
    echo " | \  / | ___| |_ __ _    | |_ __ _ _ __ ___   \  /  "
    echo " | |\/| |/ _ \ __/ _\` |   | | '__| | '_ \` _ \  /  \  "
    echo " | |  | |  __/ || (_| |   | | |  | | | | | | |/    \ "
    echo " |_|  |_|\___|\__\__,_|   |_|_|  |_|_| |_| |_/_/\_\ "
    echo -e "${NC}"
    echo -e "${BLUE}  :: ADAPTIVE BIO-CORTEX ENGINE v3.0                ::${NC}"
    echo -e "${BLUE}  :: Author: SUBRAMANIAM VIJAYAKUMAR                ::${NC}"
    echo ""
}

# --- FUNCTION: LOADING BAR ---
show_loading() {
    local pid=$1
    local delay=0.1
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}

# ==============================================================================
# SECTION 1: SYSTEM & I/O CONFIGURATION
# ==============================================================================

# --- 1.1 Output Location ---
# All results will be stored here. A timestamp prevents overwriting previous runs.
export OUTPUT_BASE_DIR="./MetaTrimX_Results_$(date +%Y%m%d_%H%M%S)"

# --- 1.2 Raw Data Inputs ---
# [USER ATTENTION]: Update these paths to point to your specific data files.
export RAW_R1="" # e.g., "./raw_data/Sample_F.fq"
export RAW_R2="" # e.g., "./raw_data/Sample_R.fq"

# --- 1.3 System Resources ---
# Adjust based on your server's capacity.
export MAX_PARALLEL_JOBS=""      # e.g., 4
export CORES_PER_JOB=""          # e.g., 12

# ==============================================================================
# SECTION 2: USER PREFERENCE SWITCHES
# ==============================================================================

# --- 2.1 Demultiplexing & Trimming Switches ---

# [PRIMER ANCHOR]
# TRUE:  Strict. Primer MUST be at the very start of the read.
# FALSE: Loose. Primer can be found anywhere in the read.
export PRIMER_ANCHOR="" # TRUE or FALSE

# [DISCARD UNTRIMMED]
# This switch ONLY controls Primer Trimming (Step 2).
# TRUE:  Strict. If primer is not found, DELETE the read.
# FALSE: Loose. If primer is not found, KEEP the read.
export DISCARD_UNTRIMMED="" # TRUE or FALSE

# --- 2.2 Clustering Switches ---

# [REMOVE SINGLETONS]
# TRUE:  Removes sequences that appear only once.
# FALSE: Keeps everything.
export REMOVE_SINGLETONS="" # TRUE or FALSE

# --- 2.3 Analysis Mode ---
# Choose "OTU" (Standard) or "ASV" (High-Res Denoising)
export ANALYSIS_MODE="" # OTU or ASV

# ==============================================================================
# SECTION 3: ADVANCED BIOINFORMATICS CONFIGURATION
# ==============================================================================

# --- 3.1 Demultiplexing ---
export DEMUX_ERROR_RATE=""    # e.g., 0.15 (15%)

# --- 3.2 Adapter Trimming ---
export ADAPTER_F=""           # e.g., AGATCGGAAGAGC
export ADAPTER_R=""           # e.g., AGATCGGAAGAGC

# --- 3.3 Primer Trimming Strategy ---
export TRIM_ERROR_RATE=""     # e.g., 0.3 (30%)

# --- 3.4 Quality Control ---
export QUALITY_CUTOFF=""      # e.g., 20
export MIN_PREPROCESS_LEN=""  # e.g., 50

# --- 3.5 Merging (VSEARCH) ---
export MERGE_MIN_OVLEN=""     # e.g., 16
export MERGE_MAX_DIFFS=""     # e.g., 2

# --- 3.6 Final Filtering ---
export MIN_FINAL_LEN=""       # e.g., 100
# Stricter filtering for high-quality OTUs
export MAX_EXPECTED_ERRORS="" # e.g., 1.0

# --- 3.7 Clustering Settings ---
export CLUSTER_IDENTITY=""    # e.g., 0.97
export ASV_IDENTITY=""        # e.g., 0.99

# >>>>> MACHINE LEARNING SETTINGS <<<<<
# TRUE: Performs deep scan of raw data to detect spacers/orientation.
export AUTO_DIAGNOSE=""       # TRUE or FALSE

# [UPDATE]: Master Switch for TensorFlow Engine
export USE_ML_FILTER=""       # TRUE or FALSE

# >>>>> PRIMER CONFIGURATION <<<<<

# SET 1 (Mandatory)
export FWD_PRIMER_1=""
export REV_PRIMER_1=""

# SET 2 (Optional)
# If you are NOT using a second set, leave these blank: ""
export FWD_PRIMER_2=""
export REV_PRIMER_2=""

export THREADS="" # e.g., 12

# ==============================================================================
# SECTION 4: SAMPLE DEFINITIONS
# ==============================================================================
# Define your samples below using Associative Arrays.
# Format: TAGS["SampleID"]="TagSequence"

declare -A TAGS

# >>>>> USER: ADD YOUR SAMPLES HERE <<<<<

# Example:
# TAGS["Sample_ID_1"]="Tag_Sequence_1"

# >>>>> END OF CONFIGURATION <<<<<

# ==============================================================================
#  ENGINE INITIALIZATION
# ==============================================================================

print_banner
echo -e "${YELLOW}[SYSTEM]${NC} Initializing MetaTrimX Engine..."
sleep 1

# Define Core Directory
CORE_DIR="./core"

# Check if Core Script exists
if [ ! -f "$CORE_DIR/metatrimx_core.py" ]; then
    echo ""
    echo -e "${RED}[CRITICAL ERROR] Core script not found!${NC}"
    echo "Missing file: $CORE_DIR/metatrimx_core.py"
    echo "Please ensure the 'core' folder is in the same directory as this script."
    exit 1
fi

# # Pass sample–tag pairs to the downstream script
SAMPLE_DATA=""
for key in "${!TAGS[@]}"; do
    # Format: ID|TAG
    ITEM="${key}|${TAGS[$key]}"
    SAMPLE_DATA+="$ITEM"$'\n'
done
export SAMPLE_DATA

# Check inputs
if [ ! -f "$RAW_R1" ]; then
    echo -e "${RED}[ERROR] Input file R1 not found: $RAW_R1${NC}"
    echo -e "Please edit this script and set the correct RAW_R1 path."
    exit 1
fi

# ------------------------------------------------------------------------------
# LOGIC INSERTION
# ------------------------------------------------------------------------------

# Set minimum cluster size based on singleton removal flag
if [ "$REMOVE_SINGLETONS" == "TRUE" ]; then
    export MIN_SIZE=2
    echo -e "${CYAN}[CONFIG]${NC} Singleton Removal: ON (Min Size = 2)"
else
    export MIN_SIZE=1
    echo -e "${YELLOW}[CONFIG]${NC} Singleton Removal: OFF (Min Size = 1)"
fi

# ==============================================================================
# SECTION 5: INTERACTIVE MODE SELECTION
# ==============================================================================

# Sync environment variables for the Configuration Mode (Mode 1 requirements)
export PRIMER_F_1="$FWD_PRIMER_1"
export PRIMER_R_1="$REV_PRIMER_1"
export PRIMER_F_2="$FWD_PRIMER_2"
export PRIMER_R_2="$REV_PRIMER_2"
export ADAPTER_SEQ="$ADAPTER_F"

echo -e "\n${BOLD}${BLUE}>>> SELECT EXECUTION MODE <<<${NC}"
echo -e "${YELLOW}1) Random Forest Optimizer:${NC} (Scans data, predicts parameters, and generates custom scripts)"
echo -e "${GREEN}2) Deep Neural Network Pipeline:${NC} (Full Pipeline - Runs Core, Merging, Clustering, and TensorFlow Filtration)"

read -r -p "Enter choice [1 or 2]: " choice

case $choice in
    1)
        export RUN_MODE="1"
        echo -e "${CYAN}Running in Mode 1: RF OPTIMIZER. Generating custom scripts...${NC}"
        ;;
    2)
        export RUN_MODE="2"
        echo -e "${CYAN}Running in Mode 2: DEEP LEARNING PIPELINE. Full pipeline execution.${NC}"
        ;;
    *)
        echo -e "${RED}[ERROR] Invalid choice. Exiting.${NC}"
        exit 1
        ;;
esac

# ==============================================================================
# SECTION 6: CONDITIONAL EXECUTION LOGIC
# ==============================================================================

if [ "$RUN_MODE" == "1" ]; then
    # Run configuration mode
    echo -e "${YELLOW}[SYSTEM]${NC} Launching Random Forest Optimizer..."
    echo -e "${CYAN}-------------------------------------------------------------${NC}"
    
    # Execute the compiler script from the core directory
    # The compiler handles the scanning, prediction, and script generation.
    python3 "$CORE_DIR/metatrimx_compiler.py"
    
    echo -e "\n=============================================================================="
    echo -e "  ${GREEN}[SUCCESS] Configuration Finished.${NC}"
    echo -e "  Please check the generated scripts: step1_...py and step2_...py"
    echo -e "=============================================================================="
    exit 0 # Exit cleanly after Mode 1 is finished

fi

# Execute Python Core (Mode 2)
echo -e "${YELLOW}[SYSTEM]${NC} Handing over control to Python Core..."
echo -e "${CYAN}-------------------------------------------------------------${NC}"

# Run the Python Core
python3 "$CORE_DIR/metatrimx_core.py"
CORE_STATUS=$?

# Run the HTML Dashboard Generator
echo -e "${YELLOW}[SYSTEM]${NC} Generating Interactive Dashboard..."
python3 "$CORE_DIR/metatrimx_report.py"

echo -e "${CYAN}-------------------------------------------------------------${NC}"

# Check status
if [ $CORE_STATUS -eq 0 ]; then
    echo ""
    echo -e "=============================================================================="
    echo -e "  ${GREEN}[SUCCESS] Pipeline completed successfully.${NC}"
    echo -e "  Results: $OUTPUT_BASE_DIR"
    echo -e "  Report:  $OUTPUT_BASE_DIR/MetaTrimX_Interactive_Report.html"
    echo -e "=============================================================================="
else
    echo ""
    echo -e "=============================================================================="
    echo -e "  ${RED}[FAILURE] Pipeline stopped due to errors.${NC}"
    echo -e "  Check the logs in: $OUTPUT_BASE_DIR/Logs"
    echo -e "=============================================================================="
    exit 1
fi