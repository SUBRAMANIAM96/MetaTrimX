#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import joblib  # For saving the brain

# Check for Machine Learning Libraries
try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.multioutput import MultiOutputRegressor
except ImportError:
    print("\n[CRITICAL ERROR] Missing AI libraries.")
    print("Please run this command first: pip install scikit-learn pandas numpy joblib")
    sys.exit(1)

# ------------------------------------------------------------------------------
# MetaTrimX - Brain Training Utility
# Author: SUBRAMANIAM VIJAYAKUMAR
# 
# This script was used to "teach" MetaTrimX how to pick the right settings. 
# It looks at the CSV file and builds a Random Forest model so the pipeline 
# can automate itself later on.
# ------------------------------------------------------------------------------

# --- CONFIGURATION ---
CSV_FILE = "metatrimx_training_data.csv"
MODEL_FILE = "metatrimx_brain.pkl"

# --- 1. INPUTS ---

FEATURES = [
    "q30_score",        # Quality (0.0 - 1.0)
    "tag_shift_mode",   # Shift (0, 1, 2)
    "primer_gap",       # Gap (0, 2, 4)
    "avg_read_len"      # Length (e.g. 250)
]

# --- 2. OUTPUTS  ---
TARGETS = [
    "best_demux_e",      # Demux Error Rate (0.1 - 0.3)
    "best_trim_e",       # Primer Trim Error (0.05 - 0.25)
    "best_merge_diffs",  # Merge Tolerance (2 - 20)
    "best_max_ee",       # Quality Filter (0.5 - 5.0)
    "best_min_len",      # Minimum Length (40 - 50)
    "best_qual_cutoff",  # Tail Trimming (15 - 25)
    "best_min_overlap"   # Merge Overlap (5 - 20)
]

def train_brain():
    print(f"\n" + "="*60)
    print(f"      INITIALIZING METATRIMX TRAINER")
    print(f"="*60)
    
    # 1. Check File Existence
    if not os.path.exists(CSV_FILE):
        print(f"[ERROR] Could not find '{CSV_FILE}'")
        print("Please ensure the CSV file is in this folder.")
        return None

    # 2. Load Data
    try:
        df = pd.read_csv(CSV_FILE)
        print(f"  -> Loaded Training Data: {len(df)} scenarios found.")
    except Exception as e:
        print(f"[ERROR] Failed to read CSV: {e}")
        return None

    # 3. Validation Check
    missing_features = [f for f in FEATURES if f not in df.columns]
    missing_targets = [t for t in TARGETS if t not in df.columns]
    
    if missing_features or missing_targets:
        print(f"[ERROR] CSV is missing columns!")
        print(f"  Missing Inputs: {missing_features}")
        print(f"  Missing Outputs: {missing_targets}")
        return None

    # 4. Prepare Training Sets
    X = df[FEATURES]
    y = df[TARGETS]

    # 5. Train Random Forest
    # The model is trained with 200 trees to ensure stable and reliable performance. A MultiOutputRegressor is used to predict all seven targets simultaneously.
    print("  -> Training Neural Decision Forest (200 Trees)...", end="")
    model = MultiOutputRegressor(RandomForestRegressor(n_estimators=200, random_state=42))
    model.fit(X, y)
    print(" Done.")

    # 6. Save the Brain
    joblib.dump(model, MODEL_FILE)
    print(f"\n[SUCCESS] ML Model compiled and saved to: '{MODEL_FILE}'")
    print("The pipeline is now ready to use this brain.")
    
    return model

def test_interaction(model):
    print("\n" + "="*60)
    print("      ML DIAGNOSTIC TERMINAL (TEST MODE)")
    print("      Type values to verify the AI's logic.")
    print("="*60)

    while True:
        try:
            print("\n--- Enter Simulated Diagnostics ---")
            q30_input = input("1. Q30 Score (e.g. 0.95, 0.85): ")
            if q30_input.lower() in ['q', 'quit', 'exit']: break
            q30 = float(q30_input)
            
            shift = int(input("2. Tag Shift (0, 1, 2): "))
            gap = int(input("3. Primer Gap (0, 2, 4): "))
            length = int(input("4. Read Length (e.g. 250): "))

            # Run Prediction
            # Note: We wrap inputs in [[]] because the model expects a 2D array
            prediction = model.predict([[q30, shift, gap, length]])[0]

            # Parse Results
            demux_e   = round(prediction[0], 2)
            trim_e    = round(prediction[1], 2)
            merge_d   = int(prediction[2])
            max_ee    = round(prediction[3], 1)
            min_len   = int(prediction[4])
            qual_cut  = int(prediction[5])
            min_over  = int(prediction[6])

            # Display "Prescription"
            print("\n>>> AI PRESCRIPTION:")
            print(f"    [Demux Error]:      {demux_e}   (Tag Strictness)")
            print(f"    [Trim Error]:       {trim_e}   (Primer Strictness)")
            print(f"    [Merge Tolerance]:  {merge_d}    (Mismatches allowed)")
            print(f"    [MaxEE Filter]:     {max_ee}    (Quality Threshold)")
            print("-" * 40)
            print(f"    [Safety Rails]:     Len>{min_len}, Qual>{qual_cut}, Overlap>{min_over}")

            # Logic Sanity Check
            if shift > 0 and demux_e < 0.20:
                print("\n    ⚠️  WARNING: AI suggests strict settings despite a shift. Check training data.")
            elif q30 > 0.95 and demux_e > 0.15:
                print("\n    ⚠️  WARNING: AI suggests loose settings for perfect data. Check training data.")
            else:
                print("\n    ✅  Logic Verified: Parameters match the scenario.")

        except ValueError:
            print("[!] Invalid input. Please enter numbers.")
        except KeyboardInterrupt:
            print("\nExiting.")
            break

if __name__ == "__main__":
    trained_model = train_brain()
    if trained_model:
        test_interaction(trained_model)