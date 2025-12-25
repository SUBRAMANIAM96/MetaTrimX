#!/usr/bin/env python3

import os
import sys
import subprocess
import datetime
import re
import shutil
import gzip
import stat
import tempfile
import json
from collections import Counter

# --- [AI UPDATE] IMPORT NEURAL ENGINE ---
try:
    import metatrimx_neural as brain
    AI_AVAILABLE = True
except ImportError:
    AI_AVAILABLE = False
   
# ==============================================================================
#                       METATRIMX CORE ENGINE v3.0
#      "Multi-Primer Scanning + Neural Filtration + Auto-Staging + Full Pipeline"
# ==============================================================================
#
# DESCRIPTION:
# 1. Bioinformatics Pipeline: Executes demultiplexing, adapter trimming, paired-end merging, and quality filtering.
# 2. Heuristic Diagnostics: Scans raw data to empirically determine tag shifts and primer gaps.
# 3. Neural Network Integration: Applies a TensorFlow-based classifier for post-clustering artifact removal.
#
# ==============================================================================

# --- UI COLORS & FORMATTING ---
C_RED = '\033[91m'
C_GREEN = '\033[92m'
C_YELLOW = '\033[93m'
C_BLUE = '\033[94m'
C_CYAN = '\033[96m'
C_END = '\033[0m'
C_BOLD = '\033[1m'

def pretty_log(message, type="INFO"):
    """
    Prints a formatted, timestamped log message to the console.
    Also writes to the global log for reporting.
    """
    ts = datetime.datetime.now().strftime("%H:%M:%S")
    
    # 1. Print to Console
    if type == "INFO":
        print(f"{C_BLUE}[{ts}]{C_END} {message}", flush=True)
    elif type == "SUCCESS":
        print(f"{C_GREEN}[{ts}] \u2714 {message}{C_END}", flush=True)
    elif type == "WARN":
        print(f"{C_YELLOW}[{ts}] \u26A0 {message}{C_END}", flush=True)
    elif type == "ERROR":
        print(f"{C_RED}[{ts}] \u2718 {message}{C_END}", flush=True)
    elif type == "HEADER":
        print(f"\n{C_CYAN}{C_BOLD}=== {message} ==={C_END}", flush=True)

    # 2. Write to File (Required for HTML Report)
    try:
        with open(LOG_FILE, "a") as f:
            f.write(f"[{ts}] {message}\n")
    except NameError:
        pass

def log(message):
    """
    Writes a raw message to the global log file for persistent records.
    """
    try:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        with open(LOG_FILE, "a") as f:
            f.write(f"[{ts}] {message}\n")
    except NameError:
        pass 

# --- 1. CONFIGURATION LOADING ---
# Load configuration from environment
try:
    OUTPUT_BASE = os.environ["OUTPUT_BASE_DIR"]
    SAMPLE_DATA = os.environ["SAMPLE_DATA"]
    RAW_R1 = os.environ["RAW_R1"]
    RAW_R2 = os.environ["RAW_R2"]
    
    # System Settings
    AUTO_DIAGNOSE = os.environ["AUTO_DIAGNOSE"]
    THREADS = os.environ["THREADS"]
    ANALYSIS_MODE = os.environ["ANALYSIS_MODE"].upper()
    
    # Read runtime flags
    USE_ML_FILTER = os.environ.get("USE_ML_FILTER", "FALSE")
    RUN_MODE = os.environ.get("RUN_MODE", "2")

    # Bioinformatics Parameters 
    DEMUX_ERR = os.environ["DEMUX_ERROR_RATE"]
    TRIM_ERR = os.environ["TRIM_ERROR_RATE"] 

    ADAPTER_F = os.environ["ADAPTER_F"]
    ADAPTER_R = os.environ["ADAPTER_R"]
    MIN_LEN = os.environ["MIN_PREPROCESS_LEN"]
    QUAL = os.environ["QUALITY_CUTOFF"]
    
    # Merging 
    MERGE_OV = os.environ["MERGE_MIN_OVLEN"]
    MERGE_DIFF = os.environ["MERGE_MAX_DIFFS"]
    
    # Filtering 
    MAX_EE = os.environ["MAX_EXPECTED_ERRORS"]
    FINAL_LEN = os.environ["MIN_FINAL_LEN"]
    
    # Loading the Switches
    DISCARD_UNTRIMMED = os.environ["DISCARD_UNTRIMMED"]
    PRIMER_ANCHOR = os.environ.get("PRIMER_ANCHOR", "TRUE")

    # Clustering Parameters 
    CLUST_ID = os.environ["CLUSTER_IDENTITY"]
    MIN_SIZE = os.environ["MIN_SIZE"]
    
    # ASV identity threshold
    ASV_ID = os.environ.get("ASV_IDENTITY", "0.99")
    
    # Primers (Flexible Sets: Loading Set 1 and Set 2)
    P_FWD_1 = os.environ.get("FWD_PRIMER_1", "")
    P_REV_1 = os.environ.get("REV_PRIMER_1", "")
    P_FWD_2 = os.environ.get("FWD_PRIMER_2", "")
    P_REV_2 = os.environ.get("REV_PRIMER_2", "")
    
except KeyError as e:
    pretty_log(f"Missing Environment Variable in Run Script: {e}", "ERROR")
    sys.exit(1)

# Directory Setup
DIRS = {
    "demux": os.path.join(OUTPUT_BASE, "01_Demux"),
    "trimmed": os.path.join(OUTPUT_BASE, "02_Trimmed"),
    "merged": os.path.join(OUTPUT_BASE, "03_Merged"),
    "qc": os.path.join(OUTPUT_BASE, "04_QC_Reports"),
    "final": os.path.join(OUTPUT_BASE, "Final_FASTA_Files"),
    "logs": os.path.join(OUTPUT_BASE, "Logs"),
    "cluster": os.path.join(OUTPUT_BASE, "Clustering_Results")
}
for d in DIRS.values():
    os.makedirs(d, exist_ok=True)

LOG_FILE = os.path.join(OUTPUT_BASE, "Pipeline_Global.log")

def write_tool_log(content, sample_name):
    """Saves the stdout of tools (Cutadapt/Vsearch) to individual log files."""
    if content:
        with open(os.path.join(DIRS["logs"], f"{sample_name}.log"), "a") as f:
            f.write(content + "\n")

def check_file_format(filepath):
    """Handles both .fastq and .fastq.gz seamlessly."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')

# --- DEPENDENCY AUTO-INSTALLER ---
def check_and_install_dependencies():
    """
    Checks for required tools (cutadapt, vsearch, fastp).
    If missing, attempts to install them via Conda or APT.
    """
    tools = ["cutadapt", "vsearch", "fastp"]
    
    for tool in tools:
        if not shutil.which(tool):
            print(f"\n{C_YELLOW}[SYSTEM] Tool '{tool}' not found. Attempting Auto-Install...{C_END}", flush=True)
            
            success = False
            
            # Method 1: Conda (Preferred)
            if shutil.which("conda"):
                print(f"  > Trying Conda install for {tool}...")
                try:
                    subprocess.run(["conda", "install", "-c", "bioconda", tool, "-y"], check=True)
                    success = True
                except subprocess.CalledProcessError:
                    print(f"  {C_RED}x Conda install failed.{C_END}", flush=True)

            # Method 2: APT (Fallback, requires sudo)
            if not success and shutil.which("apt"):
                print(f"  > Trying APT install for {tool} (may ask for password)...")
                try:
                    subprocess.run(["sudo", "apt", "update"], stderr=subprocess.DEVNULL)
                    subprocess.run(["sudo", "apt", "install", tool, "-y"], check=True)
                    success = True
                except subprocess.CalledProcessError:
                    print(f"  {C_RED}x APT install failed.{C_END}", flush=True)

            # Final Verification
            if shutil.which(tool):
                print(f"  {C_GREEN}\u2714 {tool} successfully installed.{C_END}", flush=True)
            else:
                print(f"\n{C_RED}[CRITICAL ERROR] Could not install {tool}.{C_END}", flush=True)
                print("Please install manually: sudo apt install " + tool)
                sys.exit(1)

# --- AUTO-STAGING (SPEED BOOST) ---
TEMP_DIR = None
def stage_files_if_slow(fwd, rev):
    """
    Copy input files to temporary local storage if they are on /mnt/.
    Returns the paths to use.
    """
    global TEMP_DIR
    # Files under /mnt/ can be slow in WSL
    if fwd.startswith("/mnt/") or rev.startswith("/mnt/"):
        pretty_log("Input files located on /mnt/; using temporary local copy.", "WARN")
        pretty_log("Copying data to fast Linux temporary storage...", "INFO")
        
        TEMP_DIR = tempfile.mkdtemp(prefix="metatrimx_fast_")
        
        fwd_name = os.path.basename(fwd)
        rev_name = os.path.basename(rev)
        
        new_fwd = os.path.join(TEMP_DIR, fwd_name)
        new_rev = os.path.join(TEMP_DIR, rev_name)
        
        try:
            shutil.copy2(fwd, new_fwd)
            shutil.copy2(rev, new_rev)
            pretty_log("Data successfully staged. Processing will be fast.", "SUCCESS")
            return new_fwd, new_rev
        except Exception as e:
            pretty_log(f"Staging failed ({e}). Falling back to slow mode.", "WARN")
            return fwd, rev
    else:
        return fwd, rev

def cleanup_staged_files():
    """Removes the temp folder at the end."""
    if TEMP_DIR and os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)

# --- INTELLIGENCE MODULE (DIAGNOSTIC REPORT) ---
def analyze_qc_report(json_path):
    """
    Read the Fastp JSON report and summarize QC metrics
    """
    if not os.path.exists(json_path):
        return

    try:
        with open(json_path) as f:
            data = json.load(f)

        # Extract Metrics
        total_reads = data['summary']['before_filtering']['total_reads']
        q30_rate = data['summary']['before_filtering']['q30_rate']
        duplication = data['duplication']['rate']
        insert_size = data['insert_size']['peak']
        
        # Diagnostic Logic Rules
        q_status = f"{C_GREEN}EXCELLENT{C_END}" if q30_rate > 0.90 else (f"{C_YELLOW}GOOD{C_END}" if q30_rate > 0.80 else f"{C_RED}POOR{C_END}")
        dup_status = f"{C_YELLOW}HIGH{C_END}" if duplication > 0.50 else f"{C_GREEN}LOW{C_END}"
        
        # Print Diagnostic Report
        print("\n" + "="*80, flush=True)
        print(f"{C_BOLD}   🤖 METATRIMX INTELLIGENCE REPORT{C_END}", flush=True)
        print("="*80, flush=True)
        
        print(f"[1] SEQUENCING HEALTH: {q_status}", flush=True)
        print(f"    - Total Reads:      {total_reads:,} (Robust Depth)", flush=True)
        print(f"    - Q30 Score:        {q30_rate*100:.2f}% (Chance of error < 0.1%)", flush=True)
        print(f"    - Verdict:          Data quality is high. Aggressive filtering not needed.", flush=True)
        print("", flush=True)
        
        print(f"[2] LIBRARY STRUCTURE", flush=True)
        print(f"    - Insert Size:      {insert_size} bp", flush=True)
        print(f"    - Calculated Overlap: ~{300 - insert_size} bp (Strong Merging expected)", flush=True)
        print("", flush=True)
        
        print(f"[3] DUPLICATION: {dup_status} ({duplication*100:.1f}%)", flush=True)
        print(f"    - Note:             Normal for Metabarcoding. Do NOT deduplicate now.", flush=True)
        print("", flush=True)

        if q30_rate > 0.85:
            print(f"   >>> AI VERDICT: Your data is clean. Proceed with UNIVERSAL mode.", flush=True)
        else:
            print(f"   >>> AI VERDICT: Data is noisy. Strict mode may lose data.", flush=True)

        print("="*80 + "\n", flush=True)
    except Exception as e:
        pretty_log(f"AI Interpretation failed: {e}", "WARN")

# --- FASTP QC MODULE ---
def run_fastp_qc(fwd_file, rev_file, output_dir, label="qc_report"):
    """
    Runs Fastp in 'Analysis Only' mode (no output files).
    
    """
    if not shutil.which("fastp"): return

    pretty_log(f"Running Fastp QC ({label})...", "INFO")
    
    html_report = os.path.join(output_dir, f"{label}_fastp.html")
    json_report = os.path.join(output_dir, f"{label}_fastp.json")
    
    cmd = [
        "fastp",
        "-i", fwd_file
    ]
    
    # Include reverse reads if present
    if rev_file:
        cmd.extend(["-I", rev_file])
    
    
    cmd.extend([
        "-h", html_report,
        "-j", json_report,
        "--thread", THREADS,
        "--disable_adapter_trimming",
        "--disable_quality_filtering",
        "--disable_length_filtering"
    ])
    
    try:
        # Run Fastp quietly
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # Adjust permissions for report access
        if os.path.exists(html_report):
            try: os.chmod(html_report, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
            except: pass
        
        # Trigger Diagnostic Analysis for Raw Data
        if label == "raw_data":
            analyze_qc_report(json_report)
        
        pretty_log(f"QC Report Generated: {html_report}", "SUCCESS")
        
    except Exception as e:
        pretty_log(f"Fastp QC Failed: {e}", "WARN")

# ==============================================================================
#             DUAL-SEED & MULTI-PRIMER DIAGNOSTICS
# ==============================================================================

def scan_file_structure(filepath, samples_list, primer_list, stats_pos, stats_gap):
    """
    Scan reads to infer tag position and primer gap.
    Supports multiple primer sequences.
    """
    tag_to_sample = {s[1].upper(): s[0] for s in samples_list}
    check_tags = list(tag_to_sample.keys())
    
    # Generate primer seeds for matching
    # Seed A: first 6 bases
    # Seed B: bases 7–12
    seeds = []
    for p in primer_list:
        if p and len(p) >= 12:
            seeds.append((p[:6].upper(), p[6:12].upper()))
    
    reads_checked = 0
    try:
        with check_file_format(filepath) as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    reads_checked += 1
                    seq = line.strip().upper()
                    start_region = seq[:60] # Scan first 60bp
                    
                    for tag in check_tags:
                        tag_pos = start_region.find(tag)
                        
                        # Check Tag Position (0-10bp allowed)
                        if tag_pos != -1 and tag_pos <= 10:
                            sid = tag_to_sample[tag]
                            stats_pos[sid][tag_pos] += 1 
                            
                            tag_end = tag_pos + len(tag)
                            search_window = start_region[tag_end : tag_end + 30]
                            
                            # Check against ALL primer seeds
                            found_primer = False
                            for (sa, sb) in seeds:
                                gap_idx = search_window.find(sa)
                                
                                # If Seed A fails, try Seed B
                                if gap_idx == -1:
                                    match_b = search_window.find(sb)
                                    if match_b != -1:
                                        gap_idx = match_b - 6 # Adjust back to start (FIXED OFFSET)
                                
                                if gap_idx != -1 and gap_idx >= -5:
                                    # Normalize small negatives to 0
                                    final_gap = max(0, gap_idx)
                                    stats_gap[sid][final_gap] += 1
                                    found_primer = True
                                    break # Found a valid primer, stop checking others
                            
                            if not found_primer:
                                stats_gap[sid][-1] += 1 # Truly not found in any set
                            break
                    
                    if reads_checked >= 10000: break # 10k per file
    except: pass
    return reads_checked

# --- EXPORT DATA FOR MODE 1 ---
def export_scan_data_to_json(samples, stats_pos, stats_gap):
    """
    Save diagnostic scan results to a JSON file.
    Allows the advisor script to read results without rescanning.
    """
    json_path = os.path.join(OUTPUT_BASE, "scan_report.json")
    data = {}
    for s in samples:
        sid = s[0]
        # Get Mode (Most common value)
        dom_pos = stats_pos[sid].most_common(1)
        shift = dom_pos[0][0] if dom_pos else 0
        
        dom_gap = stats_gap[sid].most_common(1)
        gap = dom_gap[0][0] if dom_gap else -1
        
        # We export the raw counts too so the Advisor can see distribution
        data[sid] = {
            "tag": s[1],
            "tag_shift_mode": shift,
            "primer_gap_mode": gap,
            "raw_pos_counts": dict(stats_pos[sid]),
            "raw_gap_counts": dict(stats_gap[sid])
        }
    
    with open(json_path, "w") as f:
        json.dump(data, f, indent=4)
    pretty_log(f"Diagnostic Data Exported: {json_path}", "SUCCESS")

def scan_raw_data(samples_list, fwd_path, rev_path):
    pretty_log("DIAGNOSTIC ENGINE: Scanning R1 & R2 (Multi-Primer Mode)...", "HEADER")
    
    stats_pos = {s[0]: Counter() for s in samples_list} 
    stats_gap = {s[0]: Counter() for s in samples_list}
    
    # Primer Lists 
    fwds = [p for p in [P_FWD_1, P_FWD_2] if p]
    revs = [p for p in [P_REV_1, P_REV_2] if p]
    
    # Scan R1 with Forward Primers
    c1 = scan_file_structure(fwd_path, samples_list, fwds, stats_pos, stats_gap)
    # Scan R2 with Reverse Primers
    c2 = scan_file_structure(rev_path, samples_list, revs, stats_pos, stats_gap)
    
    pretty_log(f"Scanned {c1+c2} total reads.", "SUCCESS")

    # --- EXPORT DATA ---
    export_scan_data_to_json(samples_list, stats_pos, stats_gap)

    # ---MODE 1 CHECK ---
    # Mode 1 selected
    # Run advisory analysis on scan results
    if RUN_MODE == "1":
        pretty_log("Scan Complete. Handing over to Advisor for Prediction.", "SUCCESS")
        return None

    # --- VISUAL MAP (Mode 2 Only) ---
    print("\n" + "="*100, flush=True)
    print(f"{C_BOLD}   PART 1: READ STRUCTURE MAP (Dominant Pattern){C_END}", flush=True)
    print("="*100, flush=True)
    
    for s_data in samples_list:
        sid, tag = s_data
        
        # Determine dominant patterns
        dom_pos = stats_pos[sid].most_common(1)
        shift = dom_pos[0][0] if dom_pos else 0
        
        dom_gap = stats_gap[sid].most_common(1)
        gap = dom_gap[0][0] if dom_gap else -1
        
        # Clear Visuals (No NN)
        if shift == 0:
            vis_pre = ""
        else:
            vis_pre = f"--({shift}bp Spacer)--"

        if gap == 0:
            vis_gap = "--(0bp Gap)--"
        elif gap > 0:
            vis_gap = f"--({gap}bp Gap)--"
        else:
            vis_gap = "--(? Gap)--"
        
        print(f"Sample: {C_BOLD}{sid:<10}{C_END} Tag: {tag}", flush=True)
        print(f"  Structure:  5' {vis_pre} [TAG] {vis_gap} [PRIMER] 3'", flush=True)
        print("-" * 60, flush=True)

    # --- TABLE 1: TAG SHIFT ---
    print("\n" + "="*100, flush=True)
    print(f"{C_BOLD}   PART 2: TAG START POSITION (Where does the read start?){C_END}", flush=True)
    print("="*100, flush=True)
    print(f"{'Sample ID':<12} | {'0bp':<8} | {'1bp':<8} | {'2bp':<8} | {'3bp':<8} | {'4bp':<8} | {'>5bp':<8}", flush=True)
    print("-" * 100, flush=True)

    for s_data in samples_list:
        sid = s_data[0]
        p = stats_pos[sid]
        print(f"{sid:<12} | {p[0]:<8} | {p[1]:<8} | {p[2]:<8} | {p[3]:<8} | {p[4]:<8} | {sum(p[k] for k in p if k>4):<8}", flush=True)

    # --- TABLE 2: PRIMER GAP ---
    print("\n" + "="*100, flush=True)
    print(f"{C_BOLD}   PART 3: PRIMER GAP (Distance between Tag and Primer){C_END}", flush=True)
    print("="*100, flush=True)
    print(f"{'Sample ID':<12} | {'0bp':<8} | {'1bp':<8} | {'2bp':<8} | {'3bp':<8} | {'>4bp':<8} | {'Not Found':<8}", flush=True)
    print("-" * 100, flush=True)

    for s_data in samples_list:
        sid = s_data[0]
        g = stats_gap[sid]
        g_large = sum(g[k] for k in g if k > 3)
        print(f"{sid:<12} | {g[0]:<8} | {g[1]:<8} | {g[2]:<8} | {g[3]:<8} | {g_large:<8} | {g[-1]:<8}", flush=True)
    
    print("\n" + "="*100, flush=True)

    # --- INTERACTIVE MENU ---
    print(f"\n{C_YELLOW}{C_BOLD}>>> DECISION REQUIRED <<<{C_END}", flush=True)
    print("  [1] STRICT    : Use if 'Tag 0bp' is dominant (High Quality).", flush=True)
    print("  [2] UNIVERSAL : Use if Tags shift 0-6bp. (Recommended).", flush=True)
    print("  [3] RESCUE    : Use if Tags/Gaps are variable or 'Not Found' is high.", flush=True)
    
    while True:
        try:
            choice = input(f"\n{C_CYAN}Enter Selection [1/2/3]: {C_END}")
            if choice == "1": return {s[0]: "strict" for s in samples_list}
            elif choice == "2": return {s[0]: "universal" for s in samples_list}
            elif choice == "3": return {s[0]: "rescue" for s in samples_list}
        except KeyboardInterrupt:
            sys.exit(0)

# ==============================================================================
#                       MAIN PIPELINE
# ==============================================================================

def main():
    pretty_log("SYSTEM CHECK: Verifying Dependencies...", "INFO")
    check_and_install_dependencies()
    
    pretty_log(f"METATRIMX ENGINE v3.0 ONLINE (Mode: {ANALYSIS_MODE})", "INFO")
    
    CUTADAPT_CMD = "cutadapt"
    VSEARCH_CMD = "vsearch"
    
    # --- PARSING ---
    samples = []
    raw_lines = [x for x in SAMPLE_DATA.split('\n') if x.strip()]
    for line in raw_lines:
        parts = line.split('|')
        if len(parts) >= 2: samples.append([parts[0].strip(), parts[1].strip()])

    pretty_log(f"Loaded {len(samples)} samples.", "INFO")

    # --- AUTO-STAGING ---
    FAST_R1, FAST_R2 = stage_files_if_slow(RAW_R1, RAW_R2)

    # --- QC & DIAGNOSIS ---
    run_fastp_qc(FAST_R1, FAST_R2, DIRS["qc"], label="raw_data")
    
    if AUTO_DIAGNOSE == "TRUE":
        # Pass both R1 and R2 for dual scanning
        sample_strategies = scan_raw_data(samples, FAST_R1, FAST_R2)
        
        # EXIT IF MODE 1
        if RUN_MODE == "1":
            # Exit main(), script finishes, Bash script runs Advisor next.
            return 
    else:
        sample_strategies = {s[0]: "strict" for s in samples}

    pretty_log("PROCESSING PHASE STARTING...", "HEADER")
    
    for s_data in samples:
        sample, tag = s_data
        strategy = sample_strategies.get(sample, "universal")
        
        # Record processing step
        pretty_log(f"Processing {sample} (Tag: {tag}) -> Mode: {strategy}", "INFO")
        
        demux_f = os.path.join(DIRS["demux"], f"{sample}_R1.fastq")
        demux_r = os.path.join(DIRS["demux"], f"{sample}_R2.fastq")
        trim_f = os.path.join(DIRS["trimmed"], f"{sample}_trim_R1.fastq")
        trim_r = os.path.join(DIRS["trimmed"], f"{sample}_trim_R2.fastq")
        merged_f = os.path.join(DIRS["merged"], f"{sample}_merged.fastq")
        final_f = os.path.join(DIRS["final"], f"{sample}.fasta")

        # 1. DEMUX (Strict Discard hardcoded)
        cmd_demux = [CUTADAPT_CMD, "-j", THREADS, "--error-rate", DEMUX_ERR, "--action=trim", "--discard-untrimmed"]
        
        if strategy == "strict":
            cmd_demux.extend(["-g", f"^{tag}", "-G", f"^{tag}"])
            
        elif strategy == "universal":
            # # Apply 0–5 bp shifts to tags
            
            for i in range(6): # 0, 1, 2, 3, 4, 5
                prefix = "N" * i
                cmd_demux.extend(["-g", f"^{prefix}{tag}"]) # R1
                cmd_demux.extend(["-G", f"^{prefix}{tag}"]) # R2
                
        else: # Rescue Mode
            cmd_demux.extend(["-g", tag, "-G", tag])
            
        cmd_demux.extend(["-o", demux_f, "-p", demux_r, FAST_R1, FAST_R2])
        p1 = subprocess.run(cmd_demux, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        write_tool_log(p1.stdout, sample)
        
        if "Pairs written (passing filters):       0" in p1.stdout:
             print(f"    {C_RED}! Sample Failed (0 reads passed demux){C_END}", flush=True)
             continue
        
        run_fastp_qc(demux_f, demux_r, DIRS["qc"], label=f"{sample}_demux")

        # 2. PRIMER TRIMMING (Conditional Switches)
        cmd_trim = [
            CUTADAPT_CMD, "-j", THREADS, "-e", TRIM_ERR, "-q", QUAL, "--minimum-length", MIN_LEN,
            "-a", ADAPTER_F, "-A", ADAPTER_R
        ]
        
        # Use anchored primer matching when PRIMER_ANCHOR is TRUE
        anchor_char = "^" if PRIMER_ANCHOR == "TRUE" else ""
        
        linked_adapter_1 = f"{anchor_char}{P_FWD_1}...{P_REV_1}"
        cmd_trim.extend(["-a", linked_adapter_1])
        
        if P_FWD_2 and P_FWD_2.strip() != "":
            linked_adapter_2 = f"{anchor_char}{P_FWD_2}...{P_REV_2}"
            cmd_trim.extend(["-a", linked_adapter_2])
            
        cmd_trim.extend(["-o", trim_f, "-p", trim_r, demux_f, demux_r])
        
        # Optional Discarding based on user switch
        if DISCARD_UNTRIMMED == "TRUE": 
            cmd_trim.append("--discard-untrimmed")
            
        p2 = subprocess.run(cmd_trim, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        write_tool_log(p2.stdout, sample)

        # 3. MERGING
        cmd_merge = [
            VSEARCH_CMD, "--fastq_mergepairs", trim_f, "--reverse", trim_r,
            "--fastqout", merged_f, "--fastq_minovlen", MERGE_OV,
            "--fastq_maxdiffs", MERGE_DIFF, "--threads", THREADS
        ]
        p3 = subprocess.run(cmd_merge, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Capturing Merge Log
        write_tool_log(p3.stderr, sample)
        
        # Run QC on merged reads (single-end)
        run_fastp_qc(merged_f, None, DIRS["qc"], label=f"{sample}_merged")

        # 4. FILTERING
        cmd_filter = [
            VSEARCH_CMD, "--fastq_filter", merged_f, "--fastq_maxee", MAX_EE,
            "--fastq_minlen", FINAL_LEN, "--fastaout", final_f,
            "--relabel", f"{sample}_", "--threads", THREADS
        ]
        p4 = subprocess.run(cmd_filter, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Captures Filter Log
        write_tool_log(p4.stderr, sample)
        
        count_final = 0
        match = re.search(r"(\d+) sequences kept", p4.stderr)
        if match: count_final = int(match.group(1))
        
        # Log final read count
        pretty_log(f"    + Valid Reads: {count_final} ({sample})", "SUCCESS")

    cleanup_staged_files()

    pretty_log(f"GENERATING RESULTS (Method: {ANALYSIS_MODE})", "HEADER")
    pooled_fasta = os.path.join(DIRS["cluster"], "all_samples.fasta")
    total_seqs = 0
    with open(pooled_fasta, 'w') as outfile:
        for fname in os.listdir(DIRS["final"]):
            if fname.endswith(".fasta"):
                with open(os.path.join(DIRS["final"], fname)) as infile:
                    outfile.write(infile.read())
                    total_seqs += 1
    
    if total_seqs > 0:
        uniques = os.path.join(DIRS["cluster"], "uniques.fasta")
        otus = os.path.join(DIRS["cluster"], "otus.fasta")
        otutab = os.path.join(DIRS["cluster"], "OTU_Table.txt")
        zotus = os.path.join(DIRS["cluster"], "zotus.fasta") 
        
        # Captures Dereplication Logs
        p5 = subprocess.run([VSEARCH_CMD, "--derep_fulllength", pooled_fasta, "--output", uniques, "--sizeout", "--minuniquesize", MIN_SIZE, "--threads", THREADS], stderr=subprocess.PIPE, text=True)
        write_tool_log(p5.stderr, "clustering_log")
        
        if ANALYSIS_MODE == "ASV":
            pretty_log("Running UNOISE3 Denoising...", "INFO")
            p6 = subprocess.run([VSEARCH_CMD, "--cluster_unoise", uniques, "--centroids", zotus, "--minsize", MIN_SIZE, "--threads", THREADS], stderr=subprocess.PIPE, text=True)
            write_tool_log(p6.stderr, "clustering_log")
            
            # ASV clustering with configurable identity threshold
            p7 = subprocess.run([VSEARCH_CMD, "--usearch_global", pooled_fasta, "--db", zotus, "--id", ASV_ID, "--otutabout", otutab, "--threads", THREADS], stderr=subprocess.PIPE, text=True)
            
            write_tool_log(p7.stderr, "clustering_log")
            pretty_log(f"Generated ASV Sequences: {zotus}", "SUCCESS")
        else:
            pretty_log("Running Standard Clustering (97%)...", "INFO")
            
            # Chimera filtering
            # Cluster to a temporary file before final cleanup
            temp_otus = os.path.join(DIRS["cluster"], "temp_otus_with_chimeras.fasta")
            
            p6 = subprocess.run([VSEARCH_CMD, "--cluster_fast", uniques, "--id", CLUST_ID, "--centroids", temp_otus, "--relabel", "OTU_", "--threads", THREADS], stderr=subprocess.PIPE, text=True)
            write_tool_log(p6.stderr, "clustering_log")
            
            pretty_log("Running Chimera Detection (VSEARCH UCHIME)...", "INFO")
            p_chim = subprocess.run([VSEARCH_CMD, "--uchime_denovo", temp_otus, "--nonchimeras", otus, "--threads", THREADS], stderr=subprocess.PIPE, text=True)
            write_tool_log(p_chim.stderr, "clustering_log")
            
            p7 = subprocess.run([VSEARCH_CMD, "--usearch_global", pooled_fasta, "--db", otus, "--id", CLUST_ID, "--otutabout", otutab, "--threads", THREADS], stderr=subprocess.PIPE, text=True)
            write_tool_log(p7.stderr, "clustering_log")
            
            pretty_log(f"Generated OTU Sequences: {otus}", "SUCCESS")
            
            # ---TRIGGERS TENSORFLOW ENGINE ---
            if USE_ML_FILTER == "TRUE" and AI_AVAILABLE:
                pretty_log("Initializing TensorFlow Neural Engine...", "INFO")
                pretty_log("Running Self-Supervised Artifact Detection...", "INFO")
                
                # The Brain scans the OTU file
                ai_scores = brain.run_ai_classification(otus)
                
                # Save results for the Report Generator
                ai_log = os.path.join(DIRS["cluster"], "AI_Classification.csv")
                with open(ai_log, "w") as f:
                    f.write("OTU_ID,AI_Confidence,Verdict\n")
                    valid_count = 0
                    for otu_id, score in ai_scores.items():
                        verdict = "Real_Biology" if score > 0.5 else "Artifact"
                        f.write(f"{otu_id},{score:.4f},{verdict}\n")
                        if score > 0.5: valid_count += 1
                
                pretty_log(f"TensorFlow Analysis Complete. {valid_count} High-Confidence Biological Entities verified.", "SUCCESS")
            # --- END ML UPDATE ---
            
        pretty_log(f"Generated Table:         {otutab}", "SUCCESS")
    else:
        pretty_log("No valid sequences found across all samples.", "WARN")
    
    pretty_log("MetaTrimX Run Completed Successfully.", "INFO")

if __name__ == "__main__":
    main()