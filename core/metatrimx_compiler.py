#!/usr/bin/env python3

import os
import sys
import glob
import time
import joblib
import statistics
import subprocess
import datetime
import json
import random
import numpy as np

# ==============================================================================
#                 🧬 METATRIMX COMPILER: ULTIMATE ML EDITION 🧬
#          "The Orchestrator: Scans, Predicts, Writes, and Guides"
# ==============================================================================

# --- 1. SYSTEM CONFIGURATION & IMPORTS ---
# Dynamic Path for the AI Brain (looks in the same folder as this script)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BRAIN_FILE = os.path.join(SCRIPT_DIR, "metatrimx_brain.pkl")  # Your Trained AI Model

# Point to the Vizier script (Mode 1 Reporter)
REPORT_SCRIPT = os.path.join(SCRIPT_DIR, "metatrimx_vizier.py")

OUTPUT_SCRIPT_1 = "step1_cleaning_pipeline.py"
OUTPUT_SCRIPT_2 = "step2_clustering_pipeline.py"

# ANSI Colors for Professional UI
C_RED   = "\033[91m"
C_GREEN = "\033[92m"
C_YELLOW= "\033[93m"
C_CYAN  = "\033[96m"
C_BLUE  = "\033[94m"
C_RESET = "\033[0m"

def log(msg, level="INFO"):
    """
    Prints timestamped, colored logs to the console to track progress.
    """
    timestamp = datetime.datetime.now().strftime("%H:%M:%S")
    if level == "INFO":
        print(f"{C_BLUE}[{timestamp}]{C_RESET} {msg}")
    elif level == "SUCCESS":
        print(f"{C_GREEN}[{timestamp}] ✅ {msg}{C_RESET}")
    elif level == "WARN":
        print(f"{C_YELLOW}[{timestamp}] ⚠️  {msg}{C_RESET}")
    elif level == "ERROR":
        print(f"{C_RED}[{timestamp}] ❌ {msg}{C_RESET}")
    elif level == "AI":
        print(f"{C_CYAN}[{timestamp}] 🧠 {msg}{C_RESET}")

# --- 2. INPUT HARVESTING (FROM SHELL ENVIRONMENT) ---
# This ensures NO Hardcoding. It grabs what you typed in metatrimx_run.sh
PRIMER_F_1 = os.environ.get("PRIMER_F_1", "")
PRIMER_R_1 = os.environ.get("PRIMER_R_1", "")
PRIMER_F_2 = os.environ.get("PRIMER_F_2", "")
PRIMER_R_2 = os.environ.get("PRIMER_R_2", "")
ADAPTER    = os.environ.get("ADAPTER_SEQ", "AGATCGGAAGAGC")

# Import Raw Files and Sample Data from Bash
RAW_R1_PATH = os.environ.get("RAW_R1", "")
RAW_R2_PATH = os.environ.get("RAW_R2", "")
SAMPLE_DATA_RAW = os.environ.get("SAMPLE_DATA", "")

# ==============================================================================
#                     PHASE 1: THE INTELLIGENT SCANNER (THE EYES)
# ==============================================================================

def run_fastp_diagnostics(r1, r2):
    """
    Runs Fastp to accurately measure Q30 Score and Read Length.
    SAVES the HTML report for the user to see later.
    """
    log("Initializing Fastp Diagnostic Scan (Quality Check)...", "INFO")
    
    # Ensure output dir exists for the report
    report_dir = "MetaTrimX_Output/Reports"
    os.makedirs(report_dir, exist_ok=True)
    
    html_out = os.path.join(report_dir, "Fastp_Diagnostic_Report.html")
    json_out = "scan_temp.json"
    
    # We scan 20,000 reads to get a statistically significant sample
    cmd = [
        "fastp", "-i", r1, "-I", r2,
        "-j", json_out, "-h", html_out,
        "--reads_to_process", "20000", 
        "--disable_adapter_trimming", "--disable_trim_poly_g"
    ]
    
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        with open(json_out) as f: 
            data = json.load(f)
        
        # Extract Key Metrics for AI Input
        q30 = data['summary']['before_filtering']['q30_rate']
        avg_len = data['summary']['before_filtering']['read1_mean_length']
        
        # Cleanup JSON (keep HTML for user)
        if os.path.exists(json_out): os.remove(json_out)
        
        log(f"Fastp Analysis: Q30={q30:.2f} | AvgLen={int(avg_len)}bp", "SUCCESS")
        log(f"Diagnostic Report Saved: {html_out}", "INFO")
        
        return q30, int(avg_len)
        
    except Exception as e:
        log(f"Fastp Scan Failed ({e}). Using safety fallbacks.", "WARN")
        return 0.85, 250 # Fallback values

def scan_structural_features(r1, known_tags, p1_seq, p2_seq):
    """
    Physically scans raw reads to find Structural Features for the Brain:
    1. Tag Shift Mode (Exact Integer Mode)
    2. Primer Gap (Exact Integer Mode)
    
    [UPGRADE]: Checks BOTH Primer Set 1 and Primer Set 2.
    """
    log("Scanning Read Structure (Hunter Mode - Dual Primer Check)...", "INFO")
    
    scan_limit = 20000 
    shifts = []
    gaps = []
    
    # Prepare Stubs (First 10bp of primers for quick matching)
    p1_stub = p1_seq[:10].upper() if p1_seq else ""
    p2_stub = p2_seq[:10].upper() if p2_seq else ""
    
    with open(r1, 'r') as f:
        for i, line in enumerate(f):
            if i >= scan_limit * 4: break
            
            # Line 2: Sequence
            if i % 4 == 1:
                seq = line.strip().upper()
                
                # 1. FIND TAG
                found_pos = -1
                found_len = 0
                for t in known_tags:
                    pos = seq.find(t.upper())
                    # Look at start (0-25bp window allows for 2bp shifts easily)
                    if pos != -1 and pos < 25: 
                        found_pos = pos
                        found_len = len(t)
                        break
                
                if found_pos != -1:
                    shifts.append(found_pos)
                    
                    # 2. VIRTUAL TRIM & PRIMER SEARCH
                    # We define the window starting IMMEDIATELY after the tag
                    search_start = found_pos + found_len
                    # Look in the next 60bp for the primer (Allows for spacers + primers)
                    search_window = seq[search_start : search_start+60]
                    
                    p_idx = -1
                    
                    # Check Primer Set 1
                    if p1_stub:
                        p_idx = search_window.find(p1_stub)
                    
                    # Check Primer Set 2 (if Set 1 not found)
                    if p_idx == -1 and p2_stub:
                         p_idx = search_window.find(p2_stub)
                    
                    if p_idx != -1:
                        gaps.append(p_idx)
    
    # [VERIFICATION]: Log exact hit counts so the user knows it worked
    log(f"Scan Stats: Found {len(shifts)} tags and {len(gaps)} primer matches in {scan_limit} reads.", "INFO")

    # --- DETERMINE FEATURES (USING MODE/INTEGER) ---
    
    # Feature 1: Shift Mode
    if not shifts: 
        final_shift = 0 
        log("Warning: No tags found in raw data. Assuming 0bp shift.", "WARN")
    else:
        try: final_shift = int(statistics.mode(shifts))
        except: final_shift = int(shifts[0]) # Fallback if multiple modes

    # Feature 2: Primer Gap
    if not gaps: 
        final_gap = 0
    else: 
        try: final_gap = int(statistics.mode(gaps))
        except: final_gap = int(gaps[0])

    return final_shift, final_gap

# ==============================================================================
#                 PHASE 2: AI PREDICTION (THE BRAIN)
# ==============================================================================

def consult_brain(q30, shift, gap, length):
    """
    Loads trained .pkl model and predicts ALL 7 PARAMETERS.
    Input Vector: [Q30, Shift, Gap, Length]
    """
    if not os.path.exists(BRAIN_FILE):
        log("CRITICAL: Brain file missing! Using default safe parameters.", "ERROR")
        return 0.15, 0.2, 10, 2.0, 50, 20, 10 # Defaults

    try:
        model = joblib.load(BRAIN_FILE)
        
        # Construct Input Vector
        input_features = np.array([[q30, shift, gap, length]])
        
        log(f"Brain Input: Q30={q30:.2f} | Shift={shift} | Gap={gap} | Len={length}", "AI")
        
        # Run Inference
        prediction = model.predict(input_features)[0]
        
        # Unpack the 7 Targets defined in the Trainer
        p_demux   = round(float(prediction[0]), 2)
        p_trim    = round(float(prediction[1]), 2)
        p_merge   = int(prediction[2])
        p_maxee   = round(float(prediction[3]), 1)
        p_min_len = int(prediction[4])
        p_qual    = int(prediction[5])
        p_min_ov  = int(prediction[6])
        
        log(f"AI Prescription Generated Successfully:", "SUCCESS")
        print(f"   > Demux Err: {p_demux} | Trim Err: {p_trim}")
        print(f"   > Merge Tol: {p_merge} | MaxEE: {p_maxee}")
        print(f"   > Safety: Len>{p_min_len}, Q>{p_qual}, Overlap>{p_min_ov}")
        
        return p_demux, p_trim, p_merge, p_maxee, p_min_len, p_qual, p_min_ov

    except Exception as e:
        log(f"Brain Inference Failed ({e}). Using safe defaults.", "WARN")
        return 0.15, 0.2, 10, 2.0, 50, 20, 10

# ==============================================================================
#                 PHASE 3: SCRIPT GENERATION (THE HANDS)
# ==============================================================================

def generate_step1_script(file_data, use_anchor):
    """
    Writes the Python script for Demultiplexing, Trimming, Merging, and Filtering.
    [NEW] INCLUDES FASTP QC CALLS AFTER EVERY STEP.
    """
    
    anchor_str = "^" if use_anchor else ""
    has_dual_primers = "True" if (PRIMER_F_2 and PRIMER_R_2) else "False"

    script_content = f'''#!/usr/bin/env python3
import os
import subprocess
import sys
import time

# ==============================================================================
#           STEP 1: CLEANING, MERGING & RESCUE (AUTO-GENERATED)
# ==============================================================================
# AUTO-GENERATED BY METATRIMX COMPILER (ML EDITION)
# BASED ON AI DIAGNOSIS OF RAW DATA

# --- BIOLOGICAL CONSTANTS (INJECTED FROM SHELL) ---
PRIMER_F_1 = "{PRIMER_F_1}"
PRIMER_R_1 = "{PRIMER_R_1}"
PRIMER_F_2 = "{PRIMER_F_2}"
PRIMER_R_2 = "{PRIMER_R_2}"
ADAPTER    = "{ADAPTER}"
HAS_DUAL   = {has_dual_primers}

# --- SYSTEM SETTINGS ---
BASE_OUT = "MetaTrimX_Output"
LOG_DIR  = os.path.join(BASE_OUT, "Logs")
QC_DIR   = os.path.join(BASE_OUT, "06_QC_Reports")
DIRS = ["01_Demux", "02_Primers", "03_Adapters", "04_Merged", "05_Filtered", "06_QC_Reports"]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#                         USER CONFIGURATION ZONE
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# [AI DIAGNOSIS]: Tags detected as "{ "Anchored (^)" if use_anchor else "Floating (No ^)" }"
# [INSTRUCTION]: Please VERIFY the tags below.

SAMPLES = [
'''
    
    # DYNAMICALLY ADD SAMPLES WITH AI PREDICTIONS
    for sample in file_data:
        script_content += f'''    {{
        "id": "{sample['id']}",
        "r1": "{sample['r1']}", 
        "r2": "{sample['r2']}",
        "f_tag": "{anchor_str}{sample['tag']}",
        "r_tag": "{anchor_str}{sample['tag']}",
        # AI PREDICTIONS (Optimized by Neural Network)
        "p_demux": {sample['p_demux']},
        "p_trim":  {sample['p_trim']},
        "p_merge": {sample['p_merge']},
        "p_maxee": {sample['p_maxee']},
        "p_minlen": {sample['p_min_len']},
        "p_qual": {sample['p_qual']},
        "p_minov": {sample['p_min_ov']}
    }},
'''

    script_content += f'''
]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_cmd(cmd, log_file):
    with open(log_file, "a") as f:
        f.write(f"\\nCMD: {{' '.join(cmd)}}\\n")
        f.flush()
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

def run_fastp_qc(r1, r2, step_name, log_file):
    """Generates QC Report for this specific pipeline step"""
    html_report = os.path.join(QC_DIR, f"{{step_name}}_fastp.html")
    json_report = os.path.join(QC_DIR, f"{{step_name}}_fastp.json")
    cmd = [
        "fastp", "-i", r1,
        "-h", html_report, "-j", json_report,
        "--disable_adapter_trimming", "--disable_quality_filtering", "--disable_length_filtering",
        "--reads_to_process", "10000"
    ]
    if r2: cmd.extend(["-I", r2])
    
    with open(log_file, "a") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

def main():
    print("\\n==================================================")
    print("      RUNNING STEP 1: AI-OPTIMIZED CLEANING")
    print("==================================================")

    # 1. Setup Directories
    for d in DIRS:
        os.makedirs(os.path.join(BASE_OUT, d), exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    # 2. Process Loop
    for s in SAMPLES:
        print(f"Processing {{s['id']}}... ", end="", flush=True)
        log_f = os.path.join(LOG_DIR, f"{{s['id']}}_pipeline.log")
        
        # --- A. DEMULTIPLEXING (AI Error Rate) ---
        out_r1 = os.path.join(BASE_OUT, "01_Demux", f"{{s['id']}}_R1.fq")
        out_r2 = os.path.join(BASE_OUT, "01_Demux", f"{{s['id']}}_R2.fq")
        cmd_demux = [
            "cutadapt", "-j", "8",
            "-e", str(s['p_demux']), "--discard-untrimmed",
            "-g", s['f_tag'], "-G", s['r_tag'],
            "-o", out_r1, "-p", out_r2,
            s['r1'], s['r2']
        ]
        run_cmd(cmd_demux, log_f)
        run_fastp_qc(out_r1, out_r2, f"{{s['id']}}_01_demux", log_f) # QC

        # --- B. PRIMER REMOVAL (AI Trim Rate + Dual Support) ---
        nop_r1 = os.path.join(BASE_OUT, "02_Primers", f"{{s['id']}}_noP_R1.fq")
        nop_r2 = os.path.join(BASE_OUT, "02_Primers", f"{{s['id']}}_noP_R2.fq")
        
        cmd_p = [
            "cutadapt", "-j", "8",
            "-e", str(s['p_trim']), "--discard-untrimmed",
            "-g", PRIMER_F_1, "-G", PRIMER_R_1
        ]
        if HAS_DUAL:
            cmd_p.extend(["-g", PRIMER_F_2, "-G", PRIMER_R_2])
            
        cmd_p.extend(["-o", nop_r1, "-p", nop_r2, out_r1, out_r2])
        run_cmd(cmd_p, log_f)
        run_fastp_qc(nop_r1, nop_r2, f"{{s['id']}}_02_primers", log_f) # QC

        # --- C. ADAPTER REMOVAL (AI Quality & Length) ---
        clean_r1 = os.path.join(BASE_OUT, "03_Adapters", f"{{s['id']}}_cl_R1.fq")
        clean_r2 = os.path.join(BASE_OUT, "03_Adapters", f"{{s['id']}}_cl_R2.fq")
        cmd_adapt = [
            "cutadapt", "-j", "8",
            "-a", ADAPTER, "-A", ADAPTER,
            "-q", str(s['p_qual']),              # <-- AI DECISION
            "--minimum-length", str(s['p_minlen']), # <-- AI DECISION
            "-o", clean_r1, "-p", clean_r2,
            nop_r1, nop_r2
        ]
        run_cmd(cmd_adapt, log_f)

        # --- D. MERGING (AI Tolerance & Overlap) ---
        merged = os.path.join(BASE_OUT, "04_Merged", f"{{s['id']}}_merged.fq")
        cmd_merge = [
            "vsearch", "--fastq_mergepairs", clean_r1, "--reverse", clean_r2,
            "--fastq_maxdiffs", str(s['p_merge']),
            "--fastq_minovlen", str(s['p_minov']),  # <-- AI DECISION
            "--fastqout", merged
        ]
        run_cmd(cmd_merge, log_f)
        run_fastp_qc(merged, None, f"{{s['id']}}_04_merged", log_f) # QC

        # --- E. FILTERING (AI MaxEE) ---
        final = os.path.join(BASE_OUT, "05_Filtered", f"{{s['id']}}_final.fasta")
        cmd_filt = [
            "vsearch", "--fastq_filter", merged,
            "--fastq_maxee", str(s['p_maxee']),
            "--fastq_minlen", str(s['p_minlen']),   # <-- AI DECISION
            "--fastaout", final
        ]
        run_cmd(cmd_filt, log_f)
        
        print(" Done.")

if __name__ == "__main__":
    main()
'''
    with open(OUTPUT_SCRIPT_1, "w") as f:
        f.write(script_content)
    log(f"Generated Cleaning Script: {OUTPUT_SCRIPT_1}", "SUCCESS")


def generate_step2_script():
    """
    Writes the Python script for Replicate Merging, Dereplication, Chimera Removal, OTU Clustering,
    AND [IMPORTANT] OTU Table Creation (Mapping).
    """
    script_content = f'''#!/usr/bin/env python3
import os
import subprocess
import glob
import sys

# ==============================================================================
#           STEP 2: CLUSTERING & TAXONOMY (AUTO-GENERATED)
# ==============================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#                         USER CONFIGURATION ZONE
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# [INSTRUCTION]: Combine your technical replicates here.
# Format: "Biological_ID": ["File_A", "File_B"]
# Example: "Des01": ["Des01a", "Des01b"]

REPLICATE_MAP = {{
    # FILL THIS IN MANUALLY BELOW:
    
}}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

INPUT_DIR  = "MetaTrimX_Output/05_Filtered"
OUTPUT_DIR = "MetaTrimX_Output/06_OTUs"

def run(cmd, log):
    with open(log, "a") as f: 
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

def main():
    print("\\n==================================================")
    print("      RUNNING STEP 2: OTU CLUSTERING")
    print("==================================================")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Fallback if user leaves map empty
    if not REPLICATE_MAP:
        print("[WARN] No Replicate Map found. Treating files as individual samples.")
        files = glob.glob(os.path.join(INPUT_DIR, "*_final.fasta"))
        for f in files:
            name = os.path.basename(f).replace("_final.fasta", "")
            REPLICATE_MAP[name] = [name]

    for bio_name, reps in REPLICATE_MAP.items():
        print(f"  Clustering {{bio_name}} ({{len(reps)}} replicates)... ", end="", flush=True)
        log_f = os.path.join(OUTPUT_DIR, f"{{bio_name}}.log")
        
        # 1. Combine Technical Replicates
        combined_fasta = os.path.join(OUTPUT_DIR, f"{{bio_name}}_combined.fasta")
        with open(combined_fasta, "w") as outfile:
            for r in reps:
                fname = os.path.join(INPUT_DIR, f"{{r}}_final.fasta")
                if os.path.exists(fname):
                    with open(fname) as infile: outfile.write(infile.read())
        
        # 2. Dereplication (VSEARCH)
        derep = os.path.join(OUTPUT_DIR, f"{{bio_name}}_derep.fasta")
        run(["vsearch", "--derep_fulllength", combined_fasta, "--output", derep, "--sizeout", "--minuniquesize", "2"], log_f)

        # 3. Chimera Removal (UCHIME)
        nonchim = os.path.join(OUTPUT_DIR, f"{{bio_name}}_nonchimera.fasta")
        run(["vsearch", "--uchime_denovo", derep, "--nonchimeras", nonchim], log_f)

        # 4. OTU Clustering (97%)
        otu = os.path.join(OUTPUT_DIR, f"{{bio_name}}_OTUs.fasta")
        run(["vsearch", "--cluster_size", nonchim, "--id", "0.97", "--centroids", otu, "--sizein", "--sizeout"], log_f)

        # 5. Mapping (Create OTU Table) - [ADDED FOR REPORTING]
        otutab = os.path.join(OUTPUT_DIR, f"{{bio_name}}_OTU_Table.txt")
        run(["vsearch", "--usearch_global", combined_fasta, "--db", otu, "--id", "0.97", "--otutabout", otutab, "--threads", "8"], log_f)

        print("Done.")

    print(f"\\n[SUCCESS] Final OTU files are in {{OUTPUT_DIR}}")

if __name__ == "__main__":
    main()
'''
    with open(OUTPUT_SCRIPT_2, "w") as f:
        f.write(script_content)
    log(f"Generated Clustering Script: {OUTPUT_SCRIPT_2}", "SUCCESS")

# ==============================================================================
#                 MAIN EXECUTION
# ==============================================================================

def main():
    print(f"{C_CYAN}")
    print("===============================================================")
    print("        🧬 METATRIMX COMPILER: ML-INFERENCE EDITION 🧬")
    print("===============================================================")
    print(f"{C_RESET}")

    # 1. Check AI Brain
    if not os.path.exists(BRAIN_FILE):
        log("CRITICAL ERROR: 'metatrimx_brain.pkl' not found!", "ERROR")
        # Note: We continue cautiously to let defaults handle it if needed
    
    # 2. Check Primers
    if not PRIMER_F_1 or not PRIMER_R_1:
        log("CRITICAL ERROR: Primers not detected in environment.", "ERROR")
        sys.exit(1)
    
    # 3. Check for Imported Raw Files (from Bash)
    if not RAW_R1_PATH or not os.path.exists(RAW_R1_PATH):
        log("CRITICAL ERROR: Raw Input Files (R1) not found.", "ERROR")
        sys.exit(1)

    log(f"Raw Input Source: {os.path.basename(RAW_R1_PATH)}", "INFO")
    
    # 4. Load Sample Tags from Env
    tags = []
    file_data = []
    
    if SAMPLE_DATA_RAW:
        lines = SAMPLE_DATA_RAW.strip().split('\n')
        for line in lines:
            if "|" in line:
                s_id, s_tag = line.split('|')
                s_id = s_id.strip()
                s_tag = s_tag.strip()
                tags.append(s_tag)
                file_data.append({
                    "id": s_id, 
                    "r1": RAW_R1_PATH, 
                    "r2": RAW_R2_PATH,
                    "tag": s_tag
                })
    else:
        log("No Sample Data found.", "ERROR")
        sys.exit(1)

    # 5. GATHER FEATURES
    q30, avg_len = run_fastp_diagnostics(RAW_R1_PATH, RAW_R2_PATH)
    shift, gap = scan_structural_features(RAW_R1_PATH, tags, PRIMER_F_1, PRIMER_F_2)

    # 6. PREDICT
    p_demux, p_trim, p_merge, p_maxee, p_min_len, p_qual, p_min_ov = consult_brain(q30, shift, gap, avg_len)

    # 7. UPDATE DATA
    use_anchor = True if shift == 0 else False
    for s in file_data:
        s.update({
            "p_demux": p_demux, "p_trim": p_trim, 
            "p_merge": p_merge, "p_maxee": p_maxee,
            "p_min_len": p_min_len, "p_qual": p_qual, "p_min_ov": p_min_ov
        })

    # --- EXECUTE STEP 1 ---
    print("")
    generate_step1_script(file_data, use_anchor)
    
    print(f"\n{C_YELLOW}" + "="*60)
    print(f"✋ HUMAN-IN-THE-LOOP ACTION REQUIRED")
    print(f"="*60 + f"{C_RESET}")
    print(f"1. Open {C_CYAN}{OUTPUT_SCRIPT_1}{C_RESET} in your editor.")
    print(f"2. Verify Tags. 3. Save.")
    input(f"\n{C_GREEN}Press ENTER to run Step 1 (Cleaning)...{C_RESET}")

    log("Executing Step 1 Pipeline...", "INFO")
    subprocess.run(["python3", OUTPUT_SCRIPT_1])

    # --- EXECUTE STEP 2 ---
    print("")
    generate_step2_script()
    
    print(f"\n{C_YELLOW}" + "="*60)
    print(f"✋ HUMAN-IN-THE-LOOP ACTION REQUIRED")
    print(f"="*60 + f"{C_RESET}")
    print(f"1. Open {C_CYAN}{OUTPUT_SCRIPT_2}{C_RESET} in your editor.")
    print(f"2. Define REPLICATE_MAP. 3. Save.")
    input(f"\n{C_GREEN}Press ENTER to run Step 2 (Clustering)...{C_RESET}")

    log("Executing Step 2 Pipeline...", "INFO")
    subprocess.run(["python3", OUTPUT_SCRIPT_2])

    # --- FINAL REPORTING (THE UPGRADE) ---
    if os.path.exists(REPORT_SCRIPT):
        print(f"\n{C_CYAN}>>> 📊 Auto-Generating Dual Dashboards...{C_RESET}")
        subprocess.run(["python3", REPORT_SCRIPT])
    else:
        log("Vizier script (metatrimx_vizier.py) not found. Skipping HTML generation.", "WARN")
    
    log("METATRIMX PIPELINE COMPLETED SUCCESSFULLY.", "SUCCESS")

if __name__ == "__main__":
    main()