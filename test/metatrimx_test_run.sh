#!/bin/bash

# ==============================================================================
#           METATRIMX TEST ENV: NEURAL ADAPTIVE ENGINE v3.0
#      Powered by Random Forest Optimization & Deep Neural Networks
# ==============================================================================
#                      Author: SUBRAMANIAM VIJAYAKUMAR
# ==============================================================================
#
# SYSTEM LOGIC (TEST ENVIRONMENT):
# This script executes the Hybrid Machine Learning Architecture in a sandbox mode:
#   1. ADAPTIVE CORE (Mode 1): A Random Forest Regressor that predicts optimal
#      trimming parameters based on raw data entropy.
#   2. NEURAL SENTINEL (Mode 2): A Deep Learning Classifier (TensorFlow) that
#      detects and removes PCR artifacts from final sequences.
#
# NOTE: Results are saved to a local test folder to protect production data.
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
    echo -e "${BLUE}  :: METATRIMX: ADAPTIVE BIO-CORTEX ENGINE TEST RUN :: v3.0   ::${NC}"
    echo -e "${BLUE}  :: Author: SUBRAMANIAM VIJAYAKUMAR                          ::${NC}"
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
# Results will be generated INSIDE the 'test' folder to keep your main folder clean.
export OUTPUT_BASE_DIR="./MetaTrimX_Test_Results_$(date +%Y%m%d_%H%M%S)"

# --- 1.2 Raw Data Inputs ---
# [USER ATTENTION]: We use ".." to go up one level to find the raw_data folder.
export RAW_R1="../raw_data/Fish_diet_R.FASTQ"
export RAW_R2="../raw_data/Fish_diet_F.FASTQ"

# --- 1.3 System Resources ---
# Adjust based on your server's capacity.
export MAX_PARALLEL_JOBS=4      # How many samples to process simultaneously
export CORES_PER_JOB=12          # Threads per sample

# ==============================================================================
# SECTION 2: USER PREFERENCE SWITCHES (EASY MODE)
# ==============================================================================

# --- 2.1 Demultiplexing & Trimming Switches ---

# [PRIMER ANCHOR]
# TRUE:  Strict. Primer MUST be at the very start of the read (e.g. ^Primer).
# FALSE: Loose. Primer can be found anywhere in the read (Good for variable spacers).
export PRIMER_ANCHOR="FALSE"

# [DISCARD UNTRIMMED]
# NOTE: Demultiplexing (Step 1) ALWAYS discards untrimmed reads (Strict).
# This switch ONLY controls Primer Trimming (Step 2).
# TRUE:  Strict. If primer is not found, DELETE the read. (Cleanest data).
# FALSE: Loose. If primer is not found, KEEP the read. (Maximum data rescue).
export DISCARD_UNTRIMMED="FALSE"

# --- 2.2 Clustering Switches ---

# [REMOVE SINGLETONS]
# TRUE:  Removes sequences that appear only once (Recommended to remove sequencing error).
# FALSE: Keeps everything. (Use only for extremely rare species detection).
export REMOVE_SINGLETONS="TRUE"

# --- 2.3 Analysis Mode ---
# Choose "OTU" (Standard 97%) or "ASV" (High-Res Denoising)
export ANALYSIS_MODE="OTU"

# ==============================================================================
# SECTION 3: ADVANCED BIOINFORMATICS CONFIGURATION
# ==============================================================================

# --- 3.1 Demultiplexing ---
export DEMUX_ERROR_RATE=0.15    # Error rate for tag detection (0.15 = 15%)

# --- 3.2 Adapter Trimming ---
export ADAPTER_F="AGATCGGAAGAGC"
export ADAPTER_R="AGATCGGAAGAGC"

# --- 3.3 Primer Trimming Strategy ---
export TRIM_ERROR_RATE=0.3     # Primer error tolerance (0.15 = 15%)

# --- 3.4 Quality Control ---
export QUALITY_CUTOFF=20          # Phred Score
export MIN_PREPROCESS_LEN=50      # Discard reads <50bp after adapter trimming

# --- 3.5 Merging (VSEARCH) ---
export MERGE_MIN_OVLEN=16         # Minimum overlap required (bp)
export MERGE_MAX_DIFFS=2          # Maximum mismatches allowed in the overlap region

# --- 3.6 Final Filtering ---
export MIN_FINAL_LEN=100          # Minimum merged length to keep
# Stricter filtering for high-quality OTUs
export MAX_EXPECTED_ERRORS=1.0    

# --- 3.7 Clustering Settings ---
export CLUSTER_IDENTITY=0.97      # 97% Similarity threshold (Ignored in ASV mode)
export ASV_IDENTITY=0.99          # For ASV mode (Standard is 0.99 or 1.00)

# >>>>> INTELLIGENCE SETTINGS <<<<<
# TRUE: Performs deep scan of raw data to detect spacers/orientation.
export AUTO_DIAGNOSE="TRUE" 

# Master Switch for TensorFlow Engine
export USE_ML_FILTER="TRUE"

# >>>>> PRIMER CONFIGURATION (FLEXIBLE) <<<<<

# SET 1 (Mandatory): e.g., MiFish-E
export FWD_PRIMER_1="GTCGGTAAAACTCGTGCCAGC"
export REV_PRIMER_1="CATAGTGGGGTATCTAATCCCAGTTTG"

# SET 2 (Optional): e.g., MiFish-U
# If you are NOT using a second set, leave these blank: ""
export FWD_PRIMER_2="GTTGGTAAATCTCGTGCCAGC"
export REV_PRIMER_2="CATAGTGGGGTATCTAATCCTAGTTTG"

export THREADS=12

# ==============================================================================
# SECTION 4: SAMPLE DEFINITIONS
# ==============================================================================
# Define your samples below using Associative Arrays.
# Format: TAGS["SampleID"]="TagSequence"

declare -A TAGS

# >>>>> USER: ADD YOUR SAMPLES HERE <<<<<

# Sample 1
TAGS["Des08a"]="ggtaag"

# Sample 2
TAGS["Des10a"]="cactct"

# >>>>> END OF CONFIGURATION <<<<<

# ==============================================================================
#  ENGINE INITIALIZATION
# ==============================================================================

print_banner
echo -e "${YELLOW}[SYSTEM]${NC} Initializing MetaTrimX Test Engine..."
sleep 1

# # Define Core Directory
CORE_DIR="../core"

# Check if Core Script exists
if [ ! -f "$CORE_DIR/metatrimx_core.py" ]; then
    echo ""
    echo -e "${RED}[CRITICAL ERROR] Core script not found!${NC}"
    echo "Expected location: $CORE_DIR/metatrimx_core.py"
    echo "Please ensure the 'core' folder is in the parent directory."
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
    echo -e "Please ensure your test script is inside the 'test' folder and raw_data is in the parent."
    exit 1
fi

# ------------------------------------------------------------------------------
# >>>>>> DUAL-MODE LOGIC INSERTION <<<<<<
# ------------------------------------------------------------------------------

# LOGIC: Convert User Switch "REMOVE_SINGLETONS" to "MIN_SIZE"
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

# Set compiler input variables (Mode 1 requirements)
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
        echo -e "${CYAN}Running in Mode 1: ML COMPILER. Generating custom scripts...${NC}"
        ;;
    2)
        export RUN_MODE="2"
        echo -e "${CYAN}Running in Mode 2: PRODUCTION. Full pipeline execution.${NC}"
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
    # Mode 1: Run compiler mode
    echo -e "${YELLOW}[SYSTEM]${NC} Launching Random Forest Optimizer..."
    echo -e "${CYAN}-------------------------------------------------------------${NC}"
    
    # Run compiler from core
    python3 "$CORE_DIR/metatrimx_compiler.py"
    
    echo -e "\n=============================================================================="
    echo -e "  ${GREEN}[SUCCESS] Compiler Finished.${NC}"
    echo -e "  Please check the generated scripts in the test folder."
    echo -e "=============================================================================="
    exit 0 # Exit after Mode 1 is finished

fi

# ------------------------------------------------------------------------------
# >>>>>> EXECUTION PHASE <<<<<<
# ------------------------------------------------------------------------------

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
    echo -e "  ${GREEN}[SUCCESS] Test Run completed successfully.${NC}"
    echo -e "  Results: $OUTPUT_BASE_DIR"
    echo -e "  Report:  $OUTPUT_BASE_DIR/MetaTrimX_Interactive_Report.html"
    echo -e "=============================================================================="
else
    echo ""
    echo -e "=============================================================================="
    echo -e "  ${RED}[FAILURE] Test Run stopped due to errors.${NC}"
    echo -e "  Check the logs in: $OUTPUT_BASE_DIR/Logs"
    echo -e "=============================================================================="
    exit 1
fi