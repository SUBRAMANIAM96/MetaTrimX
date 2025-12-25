# MetaTrimX v3.0: Precision Metabarcoding, Powered by Machine Learning
### The Universal "Bio-Cortex" Engine for Automated Amplicon Processing
---

##  Overview

**MetaTrimX** is a revolutionary, fully automated bioinformatics pipeline for metabarcoding analysis. It transforms raw paired-end Illumina sequencing reads into publication-ready Operational Taxonomic Units (OTUs) or Amplicon Sequence Variants (ASVs) through two distinct execution pathways powered by machine learning.

### **What Makes MetaTrimX Unique:**

MetaTrimX employs a **Dual-Path Machine Learning Architecture**:

1. **Mode 1 (Diagnostic - The Predictive Optimizer)**: Uses **Multi-Output Random Forest Regressor** to intelligently diagnose your raw sequencing data and automatically predict optimal processing parameters before execution
2. **Mode 2 The Bio-Sentinel (Quality Control)**: Deploys a **Hybrid Ensemble (RF + SVM)** to validate the pipeline's output. It acts as an automated critic, scoring every OTU to distinguish true biology from subtle artifacts.
---

##  Machine Learning Engine

MetaTrimX ships with an opinionated but **transparent ML layer** that behaves like a smart co-pilot rather than a black box. No magic, just learned patterns from real bioinformatics work.

### **Mode 1 – The Predictive Optimizer**

In **Mode 1**, MetaTrimX uses a trained **Multi-Output Random Forest Regressor** (stored in `metatrimx_brain.pkl`) to turn raw data signatures into a complete parameter set in one shot.

- **Model Type**: `RandomForestRegressor` wrapped in `MultiOutputRegressor` to predict a **7-dimensional vector** of optimal parameters simultaneously (demux error, trim tolerance, merge diffs, maxEE, min length, quality cutoff, min overlap).
  
- **Input Features** (4-vector `[Q30, tag_shift, primer_gap, read_length]`):
  - **Q30 rate** from focused `fastp` scan of 20,000 reads
  - **Tag shift** and **primer gap** estimated directly from raw reads (dual-primer aware)
  - **Mean read length** from same diagnostic pass

- **How It Works**: The model treats **pipeline configuration as a supervised regression problem**, learning from historical runs where these 7 parameters already worked well for similar data profiles. It scans your data, computes 4 features, feeds them to the RF "brain", and writes **ready-to-run scripts** (`step1_cleaning_pipeline.py`, `step2_clustering_pipeline.py`) plus a JSON config explaining what the model decided.

- **Output**: 7 simultaneous predictions:
  ```
  [demux_error, trim_error, merge_diffs, maxEE, min_len, quality_cutoff, min_overlap]
  ```

---

### **Mode 2 – The Neural Sentinel**

In **Mode 2**, MetaTrimX promotes every OTU to a **candidate biological entity** and asks a hybrid ensemble, the *Neural Sentinel*, to vote on whether it looks real or suspicious.

**Path A – Statistical Brain (Random Forest Classifier)**

- **Features per sequence**: length, GC fraction, Shannon entropy
- **Training**: Top ~20% abundant sequences = "real biology" proxies; matching number of synthetic 200-bp random sequences = "noise" class
- **Model**: `RandomForestClassifier` (200 trees, `max_depth=10`)
- **Output**: P(Real_Biology) ∈ [0,1]
- **Why**: Catches statistically impossible sequences (wrong GC%, weird length distribution)

**Path B – Motif Brain (One-Class SVM + PCA)**

- **K-mer Decomposition**: Each sequence → overlapping 5-mers
- **Vectorization**: `CountVectorizer` creates feature matrix (~4,096 k-mer types)
- **Compression**: PCA reduces to 10 principal components (retains ~95.3% variance)
- **Anomaly Detection**: `OneClassSVM` (RBF kernel, `nu=0.05`) learns "normal" k-mer boundary
- **Scoring**: Decision scores normalized via sigmoid → probability-like score in [0,1]
- **Why**: Catches unusual motif combinations (chimeras, reading errors)

### ⚖️ The Ensemble Verdict
For each OTU, MetaTrimX calculates a consensus score to determine biological validity.

**1. Calculation:**

$$
\text{Final Confidence} = \frac{\text{RF}_{\text{score}} + \text{SVM}_{\text{score}}}{2}
$$

**2. Decision Rule:**

$$
\text{Verdict} = \begin{cases} 
\text{Real Biology} & \text{if } \text{Confidence} > 0.5 \\
\text{Artifact} & \text{if } \text{Confidence} \leq 0.5 
\end{cases}
$$

```

Results exported to `AI_Classification.csv` with per-OTU RF score, SVM score, ensemble confidence, and verdict—so you **inspect, filter, or override at will**.

---

## 🚀 Key Features

| Feature | Description |
|---------|-------------|
| **General-Purpose Tool** | Works with any metabarcoding marker gene (12S, 16S, ITS, COI, etc.) |
| **Dual Execution Modes** | Mode 1 (Predictive Optimizer) + Mode 2 (Neural Sentinel) |
| **ML Parameter Prediction** | Multi-Output Random Forest diagnoses data and prescribes optimal settings |
| **ML Confidence Scoring** | Every OTU receives ensemble confidence metrics (Random Forest + One-Class SVM) |
| **Parallel Processing** | Process multiple samples simultaneously |
| **Interactive Dashboards** | HTML reports with real-time QC visualizations |
| **Dual Primer Support** | Handle multiple primer sets in single run |
| **Chimera + Error Detection** | VSEARCH UCHIME + Machine Learning ensemble filtering |
| **Flexible Analysis Modes** | OTU (97% similarity) or ASV (exact matching) |
| **Comprehensive Logging** | Detailed logs for every processing step |

---

## 📂 Project Structure

The MetaTrimX architecture is organized into the **Root** (Execution & Training), **Core** (Internal Logic), and **Data** directories.

```
MetaTrimX/
├── README.md                    # Comprehensive documentation and manual
├── requirements.txt             # List of required Python libraries
├── workflow.png                 # Visual diagram of the AI architecture
│
├── metatrimx_run.sh             # COMMAND CENTER: Run this to start the pipeline
├── metatrimx_trainer.py         # AI TRAINER: Script to retrain the model (Author: S. Vijayakumar)
├── metatrimx_training_data.csv  # GROUND TRUTH: The dataset used to train the brain
│
├── core/                        # ⚙️ ENGINE ROOM (Internal Logic)
│   ├── metatrimx_brain.pkl         # The Serialized "Brain" (Random Forest Model)
│   ├── metatrimx_compiler.py       # Mode 1: Optimizer & Script Generator
│   ├── metatrimx_core.py           # Mode 2: Main Bioinformatics Processing
│   ├── metatrimx_neural.py         # Mode 2: Neural Sentinel (Artifact Detection)
│   ├── metatrimx_report.py         # Generates the Interactive HTML Dashboards
│   └── metatrimx_vizier.py         # Logic for Mode 1 Visualizations
│
├── raw_data/                    # 📂 INPUT: Place your raw FASTQ files here
└── test/                        #     SANDBOX: Contains test scripts and sample data

```

#### **🛠️ Mode 1 Workflow: The Step-by-Step Guide**

**Step 1: Initialization & Environment Harvesting**

- **Script**: `metatrimx_run.sh`
- **Action**: You launch the pipeline and select Option 1 (Random Forest Optimizer).
- **Process**: The engine retrieves your primers, adapter sequences, and sample metadata directly from the environment variables you configured in the shell script. This ensures that all subsequent steps are 100% customized to your biological targets with no hardcoded fallbacks.

**Step 2: Diagnostic Scanning (The "Eyes")**

- **Script**: `metatrimx_compiler.py`
- **Action**: The system performs a deep-scan of a subset (typically 20,000 read pairs) of your raw FASTQ files.
- **Extracted Metrics**:
  - **Quality Profiling**: Uses `fastp` to calculate the Q30 Quality Score and Mean Read Length.
  - **Structural Mapping**: Physically scans reads for barcodes and primers to calculate Tag Shift (ΔT) (detecting if barcodes are "floating") and the Primer Gap (Γ) (the exact distance between the barcode and the biological insert).

**Step 3: Neural Parameter Prediction (The "Brain")**

- **Script/File**: `metatrimx_compiler.py` consulting `metatrimx_brain.pkl`
- **Action**: The extracted metrics (Q30, Shift, Gap, and Length) are fed into the Predictive Engine.
- **Process**: The Random Forest Regressor (stored in `metatrimx_brain.pkl`) performs a multi-output regression to predict the Optimal 7 Parameters for your dataset:
  - **Demux Error**: Threshold for barcode identification.
  - **Trim Error**: Precision for primer removal.
  - **Merge Tolerance**: Maximum differences allowed in R1/R2 overlap.
  - **MaxEE**: The "Expected Error" quality filter threshold.
  - **Min Length**: The shortest acceptable biological sequence.
  - **Quality Cutoff**: Tail-trimming stringency.
  - **Min Overlap**: Minimum required overlap for contig assembly.

**Step 4: Script Synthesis & Mandatory Pauses**

- **Script**: `metatrimx_compiler.py`
- **Action**: The compiler dynamically writes two standalone Python scripts tailored specifically to your data.
- **Generated Files**:
  - **step1_cleaning_pipeline.py**: Contains the demultiplexing, trimming, merging, and filtering logic pre-loaded with your tags and ML-predicted parameters.
  - **step2_clustering_pipeline.py**: Contains the logic for dereplication, chimera removal, and OTU/ASV clustering.

**Step 5: Human-in-the-Loop (HITL) Verification**

- **Action**: The system pauses execution twice to ensure human oversight.
- **Checkpoint Alpha**: You must open `step1_cleaning_pipeline.py`, verify your sample IDs and barcodes, save, and press ENTER to authorize processing.
- **Checkpoint Beta**: After cleaning, you must open `step2_clustering_pipeline.py` to define your REPLICATE_MAP (grouping technical replicates into biological samples) and press ENTER to start clustering.

**Step 6: Preliminary Analytics (The Vizier)**

- **Script**: `metatrimx_vizier.py`
- **Action**: Once the cleaning and clustering phases are complete, the Vizier module activates.
- **Process**: It parses the terminal logs from every tool to calculate exact read retention rates. It then renders an Interactive HTML Dashboard (e.g., `MetaTrimX_Pipeline_Dashboard.html`), providing a visual audit of where data was lost and confirming the pipeline's health.

---

## 🧠 Mathematical Model Training & Calibration

The "brain" of MetaTrimX (`metatrimx_brain.pkl`) is a serialized Multi-Output Random Forest Regressor that acts as the pipeline's decision-making core. This model was calibrated to translate raw sequencing metrics into high-precision bioinformatics parameters, eliminating the need for manual trial-and-error.

**1. Training Methodology**
The model was trained using the `metatrimx_trainer.py` utility on a ground-truth dataset (`metatrimx_training_data.csv`) containing diverse sequencing scenarios.
- **Ensemble Logic**: It utilizes an ensemble of 200 Decision Trees to ensure statistical stability and prevent over-fitting to noisy sequencing data.
- **Optimization Criterion**: Each tree is constructed through recursive partitioning to minimize the Mean Squared Error (MSE) across the multi-dimensional output space:
  ```
  MSE = (1/n) * Σ(yᵢ - ŷᵢ)²
  ```

**2. Feature Engineering (The Inputs)**
During the diagnostic phase, the engine extracts four primary features from the raw FASTQ files to represent the dataset's complexity:
- **Q30 Score (X₁)**: Overall base-calling reliability, representing the probability of accuracy P = 1 - 10^(-Q/10).
- **Tag Shift Mode (ΔT)**: Positional variance of sample barcodes (detecting sequencer phasing errors).
- **Primer Gap (Γ)**: The physical nucleotide distance between the barcode and the primer sequence.
- **Average Read Length (L)**: The mean length of raw sequencing reads.

**3. Parameter Prediction (The Outputs)**
The model treats pipeline configuration as a multidimensional regression problem, simultaneously predicting a 7-dimensional vector (Ŷ) of optimized thresholds:
```
Ŷ = [Ed, Et, Md, EE, Lmin, Qcut, Omin]
```
- **Ed / Et**: Optimal stringency for barcode and primer identification (error rates).
- **Md / Omin**: Mathematical tolerances for paired-end overlap and mismatches.
- **EE**: Predicted MaxEE (Expected Error) for high-fidelity filtering.
- **Lmin / Qcut**: Safety rails for sequence length and quality trimming.

**4. The metatrimx_brain.pkl File**
Once training is complete, the model's logic is compressed into this binary file using the joblib serialization format. The Compiler loads this file in Mode 1 to perform real-time inference on new data, effectively functioning as a "digital consultant" that prescribes the best workflow for your specific sequencing run.

**5. Retraining the Engine**
Users can update the model with lab-specific data to maintain peak performance on their institution's unique sequencing platforms:
- **Modify Data**: Add your empirical "best" parameters to `metatrimx_training_data.csv`.
- **Execute Trainer**: Run `python3 metatrimx_trainer.py` to overwrite the existing brain file.
- **Validate Logic**: Use the interactive Diagnostic Terminal that launches after training to verify that the model's predictions align with your expected biological outcomes.

---

## 📂 Mode 1: Output Directory Glossary

The `MetaTrimX_Output` folder is the primary results hub for the Predictive Optimizer. It contains all intermediate processing stages, high-resolution quality reports, and the interactive analytics dashboards generated during the run.

| Directory / File | Description |
|------------------|-------------|
| **`01_Demux/`** | Contains FastQ files separated by their sample-specific barcode tags. |
| **`02_Primers/`** | Includes reads that have successfully undergone biological primer removal. |
| **`03_Adapters/`** | Stores sequences after sequencing adapters and low-quality tails have been trimmed. |
| **`04_Merged/`** | Contains paired-end R1 and R2 reads that were successfully assembled into single contigs. |
| **`05_Filtered/`** | High-fidelity sequences that passed the ML-predicted quality and MaxEE filters. |
| **`06_OTUs/`** | The final synthesis directory containing dereplication, chimera removal, and OTU clustering data. |
| **`06_QC_Reports/`** | Houses step-by-step Fastp quality reports (HTML/JSON) generated after every pipeline stage. |
| **`Logs/`** | A comprehensive archive of raw terminal outputs and audit logs from every tool used. |
| **`Reports/`** | Contains the `Fastp_Diagnostic_Report.html`, which profiles the quality of the raw input files. |
| **`MetaTrimX_Stats.json`** | A machine-readable file containing all parsed read counts and efficiency metrics. |
| **`MetaTrimX_Pipeline_Dashboard.html`** | The interactive dashboard for the "Cleaning" phase (Demux, Primers, Merging, and Filtering). |
| **`MetaTrimX_Clustering_Dashboard.html`** | The interactive dashboard for the "Clustering" phase (Dereplication, Chimeras, and OTUs). |

---

#### **🛡️ Mode 2: The Neural Sentinel Verification (Detailed Workflow)**

In Mode 2, MetaTrimX transitions from the "Architect" phase into the "Executor" phase. This mode performs high-throughput bioinformatics processing—guided by the AI-prescribed parameters from Mode 1—and culminates in the Neural Sentinel verification, where a hybrid machine learning ensemble purges PCR artifacts and biological noise.

**Step 1: Core Execution Engine (`metatrimx_core.py`)**

This script orchestrates the primary bioinformatics workflow using the optimized parameters generated in Mode 1.

- **Demultiplexing**: Splits raw FASTQ files into sample-specific reads based on barcode tags.
- **Primer & Adapter Trimming**: Removes biological primers and sequencing adapters using predicted precision tolerances.
- **Paired-End Merging**: Utilizes VSEARCH to assemble R1 and R2 reads into single high-fidelity contigs.
- **Quality Filtering**: Applies the MaxEE (Expected Error) filter to discard low-confidence sequences.

**Step 2: Clustering & Replicate Integration**

The pipeline moves toward biological synthesis by grouping sequences:

- **Dereplication**: Collapses identical sequences into unique entries while tracking abundance.
- **Chimera Removal**: Identifies and discards PCR-generated hybrid artifacts.
- **OTU/ASV Generation**: Groups sequences by identity (e.g., 97% for OTUs) to define biological entities.

**Step 3: The Neural Sentinel Verification (`metatrimx_neural.py`)**

This is the flagship "Neural Firewall" that validates every generated sequence using a dual-path ensemble.
**Mathematical Framework**

### 🧠 Mathematical Architecture

#### 1. Statistical Features (Random Forest)
The Sentinel Engine calculates complexity metrics for every sequence to detect anomalies.

- **Shannon Entropy ($H$):** Measures the randomness of the nucleotide sequence.

$$
H = - \sum p_i \log_2 p_i
$$

*(Where $p_i$ is the frequency of each base A, C, G, T)*

- **GC Content (GC%):** The proportion of Guanine and Cytosine bases, critical for distinguishing biological sequences from sequencing artifacts.

#### 2. Motif Analysis (One-Class SVM)
- **K-mer Decomposition:** Sequences are vectorized into overlapping substrings of length $k$.
- **Feature Vector Size:** For a read of length $L$, the number of k-mers generated is:

$$
N = L - k + 1
$$

#### 3. Ensemble Scoring Logic
The final classification is a consensus between the Statistical Brain and the Motif Brain.

- **Confidence Score Formula:**

$$
\text{Confidence Score} = \frac{P(\text{RF}) + P(\text{SVM})}{2}
$$

- **Decision Thresholds**:
  - **Real Biology**: Score > 0.5
  - **Artifact**: Score ≤ 0.5

**Step 4: Final Master Reporting (`metatrimx_report.py`)**
Action: Compiles the results into the MetaTrimX Interactive Report.

- **Final OTU Table:** Quantitative matrix of samples versus biological entities.
- **AI Classification Map:** Detailed CSV showing the confidence scores and biological verdicts for every sequence.
- **Interactive HTML Dashboard:** A visual master-report for the entire production run.

---

## 📂 Mode 2: Output Directory Glossary

The `MetaTrimX_Results` folder is the comprehensive results ecosystem for the Production Run. It contains all stages of the bioinformatics pipeline, the final Neural Sentinel classification data, and the high-fidelity interactive reports.

| Directory / File | Description |
|------------------|-------------|
| **`01_Demux/`** | Contains raw reads separated by their sample-specific barcode tags. |
| **`02_Trimmed/`** | Stores sequences after biological primer and sequencing adapter removal. |
| **`03_Merged/`** | Contains paired-end R1 and R2 reads successfully assembled into single contigs via VSEARCH. |
| **`04_QC_Reports/`** | Houses detailed Fastp quality reports generated after every major processing stage. |
| **`Clustering_Results/`** | The final analytical hub containing OTU tables, FASTA files, and AI validation data. |
| **`Final_FASTA_Files/`** | High-quality, filtered sequences ready for downstream ecological analysis. |
| **`Logs/`** | Individual audit logs for every sample and specific pipeline event. |
| **`MetaTrimX_Interactive_Report.html`** | The master visual dashboard for the production run, summarizing quality and neural verification. |
| **`MetaTrimX_Master_Log.csv`** | A comprehensive quantitative audit log of all read transitions across the entire pipeline. |
| **`Pipeline_Global.log`** | The main system log capturing the overall execution flow and engine status. |
| **`scan_report.json`** | Machine-readable metrics from the initial data diagnostic scan. |

---

## 🦠 Inside `Clustering_Results`

This sub-directory contains the definitive biological outputs of the pipeline, verified by the Neural Sentinel.

| File | Description |
|------|-------------|
| **`AI_Classification.csv`** | The output of the Neural Sentinel ensemble; contains confidence scores and biological verdicts for every sequence. |
| **`OTU_Table.txt`** | The final quantitative matrix representing Sample IDs versus OTU/ASV counts. |
| **`otus.fasta`** | The representative sequences for each verified biological cluster. |
| **`uniques.fasta`** | The results of the dereplication phase, containing all unique sequences before clustering. |
| **`all_samples.fasta`** | A pooled file containing all quality-filtered reads from all samples before dereplication. |
| **`temp_otus_with_chimeras.fasta`** | An intermediate file containing potential biological entities before chimera and artifact purging. |


---

## 🛠 System Requirements & Dependencies

### Minimum System Requirements
- **OS**: Linux (Ubuntu 18.04+) or macOS (10.12+), WSL2 on Windows
- **Python**: 3.6 or higher
- **RAM**: 8GB minimum (16GB recommended for >1M reads)
- **Disk Space**: 50GB+ (depends on raw dataset size)
- **CPU**: 8+ cores recommended

### Complete Dependency List

#### **Core System Tools**

| Tool | Version | Installation | Purpose |
|------|---------|--------------|---------|
| **Python3** | ≥3.6 | `apt install python3` / `brew install python3` | Pipeline runtime |
| **Bash** | ≥4.0 | Pre-installed on Linux/Mac | Shell orchestration |
| **Git** | ≥2.0 | `apt install git` / `brew install git` | Repository management |

#### **Bioinformatics Tools**

| Tool | Version | Installation | Purpose |
|------|---------|--------------|---------|
| **cutadapt** | ≥4.0 | `pip install cutadapt` or `conda install -c bioconda cutadapt` | Primer/adapter trimming |
| **vsearch** | ≥2.18 | `conda install -c bioconda vsearch` or `apt install vsearch` | Merging, dereplication, clustering, chimera detection |
| **fastp** | ≥0.23 | `conda install -c bioconda fastp` or `apt install fastp` | Real-time QC reporting per step |
| **perl** | ≥5.0 | Pre-installed on Linux/Mac | USEARCH script compatibility |

#### **Python Libraries**

| Library | Version | Installation | Purpose |
|---------|---------|--------------|---------|
| **numpy** | ≥1.19 | `pip install numpy` | Numerical operations |
| **scikit-learn** | ≥0.24.0 | `pip install scikit-learn` | ML models (Random Forest, SVM, PCA, MultiOutputRegressor) |
| **pandas** | ≥1.0 | `pip install pandas` | Data manipulation, CSV processing |
| **joblib** | ≥1.0 | `pip install joblib` | Parallel job processing |
| **matplotlib** | ≥3.3 | `pip install matplotlib` | Graph generation for reports |
| **jinja2** | ≥2.11 | `pip install jinja2` | HTML template rendering |

---

# MetaTrimX v3.0: Comprehensive User Operation Manual

This manual provides a detailed, step-by-step guide to configuring, launching, and managing the MetaTrimX Neural Adaptive Engine.

## 1. Master Configuration

All user inputs are centralized in the `metatrimx_run.sh` script. Open this file in your text editor. You must configure the following four key sections before execution.

### Section I: System & Data Paths
Located in "Section 1" of the script.

| Variable | Description | Input Requirement |
|----------|-------------|-------------------|
| **RAW_R1** | Full path to your raw Forward reads (fastq/fastq.gz). | **Mandatory** |
| **RAW_R2** | Full path to your raw Reverse reads (fastq/fastq.gz). | **Mandatory** |
| **OUTPUT_BASE_DIR** | Folder name for results. Defaults to timestamped folder. | **Optional** |
| **MAX_PARALLEL_JOBS** | Number of samples to process simultaneously. | Adjust based on RAM |
| **CORES_PER_JOB** | CPU threads assigned to each sample. | Adjust based on CPU |

### Section II: Pipeline Switches
Located in "Section 2" of the script. Use `TRUE` or `FALSE` for toggles.

| Switch | Function | Recommendation |
|--------|----------|----------------|
| **ANALYSIS_MODE** | Choose clustering method: OTU (97%) or ASV (Denoising). | ASV for high-res; OTU for standard |
| **REMOVE_SINGLETONS** | Discards reads that appear only once. | **TRUE** (Reduces sequencing noise) |
| **PRIMER_ANCHOR** | **TRUE**: Primer must be at the very start (0bp).<br>**FALSE**: Primer can be "floating" inside the read. | **TRUE** for standard Illumina; **FALSE** for ragged ends |
| **DISCARD_UNTRIMMED** | **TRUE**: Deletes read if primer is missing.<br>**FALSE**: Keeps read even if primer isn't found. | **TRUE** (Ensures data purity) |

### Section III: Biological Definitions
Located in "Section 3" of the script.

* **Primers**: Enter your specific sequences.
  * **FWD_PRIMER_1 / REV_PRIMER_1**: Mandatory. The main primer set.
  * **FWD_PRIMER_2 / REV_PRIMER_2**: Optional. Use only if you have a second primer set mixed in.

* **Adapters**:
  * **ADAPTER_F / ADAPTER_R**: Enter your sequencing adapters (default is Illumina Universal).

* **Machine Learning**:
  * **USE_ML_FILTER**: Set to `TRUE` to activate the TensorFlow Artifact Detector post-clustering.

### Section IV: Sample Registration
Located in "Section 4" of the script.

You must map every Sample ID to its unique barcode (Tag). This allows the engine to demultiplex your pooled data.

* **Syntax**: `TAGS["Your_Sample_Name"]="Barcode_Sequence"`
* **Example**:
  ```bash
  TAGS["River_Site_A"]="ATGCGT"
  TAGS["River_Site_B"]="GCATGA"
  ```

## 2. Launching the Engine

Once configured, save the file. Open your terminal in the script's directory and run:

```bash
# 1. Make the script executable (only needed once)
chmod +x metatrimx_run.sh

# 2. Start the pipeline
./metatrimx_run.sh
```

## 3. Execution Modes (Critical Decision)

Upon running the script, you will be prompted to choose an execution mode. This choice determines how parameters are handled.

### MODE 1: Random Forest Optimizer (The "Advisor")
* **Best For**: New datasets, complex environmental samples, or when you are unsure of the best settings.

**How it Works**:
1. **Scanning**: The engine scans your raw data to measure entropy, quality (Q30), and structural anomalies (Tag Shifts/Primer Gaps).
2. **Prediction**: It feeds this data into the Random Forest Brain (`metatrimx_brain.pkl`) to predict the optimal error rates and filtering thresholds for your specific run.
3. **Generation**: It creates two custom Python scripts for you:
   * `step1_cleaning_pipeline.py` (Demux/Trim/Filter)
   * `step2_clustering_pipeline.py` (OTU/ASV Clustering)

**User Responsibility (The "Human-in-the-Loop")**:
* The pipeline stops after generating the scripts.
* **Action**: You must open `step1_cleaning_pipeline.py` to verify the AI's choices.
* **Action**: You must open `step2_clustering_pipeline.py` to fill in the Replicate Map (if you have technical replicates to merge).
* **Run**: Execute these python scripts manually when ready.
* **Note**: In this mode, the numerical parameters (error rates) in `metatrimx_run.sh` are **IGNORED** in favor of AI predictions.

### MODE 2: Deep Neural Network Pipeline (The "Auto-Pilot")
* **Best For**: Routine runs, high-throughput processing, or when you have validated parameters you want to enforce.

**How it Works**:
1. **Direct Execution**: It skips the prediction phase and immediately launches the Core Engine (`metatrimx_core.py`).
2. **Parameter Enforcement**: It uses the **exact numeric values** you entered in `metatrimx_run.sh` (e.g., `DEMUX_ERROR_RATE`, `MAX_EXPECTED_ERRORS`, `TRIM_ERROR_RATE`).
3. **End-to-End**: It runs Demultiplexing → Trimming → Merging → Clustering → Reporting without stopping.

**User Responsibility**:
* Ensure all Error Rate and Quality Cutoff variables in `metatrimx_run.sh` are filled out correctly before starting.

## 4. Strategic Parameter Guide

Use this table to decide which parameters to edit based on your chosen mode.

| Variable in Script | Mode 1 (AI Advisor) | Mode 2 (Auto-Pilot) |
|--------------------|---------------------|---------------------|
| **DEMUX_ERROR_RATE** | Ignored (AI predicts this). | **Required**. (Standard: 0.15) |
| **TRIM_ERROR_RATE** | Ignored (AI predicts this). | **Required**. (Standard: 0.1 - 0.2) |
| **MAX_EXPECTED_ERRORS** | Ignored (AI predicts this). | **Required**. (Standard: 1.0 - 2.0) |
| **MIN_PREPROCESS_LEN** | Ignored (AI predicts this). | **Required**. (Standard: 50) |
| **PRIMER_SEQUENCES** | **Required**. (AI needs these to scan). | **Required**. (Engine needs these to cut). |
| **TAGS (Samples)** | **Required**. | **Required**. |

## 5. Output Interpretation

Results are saved in a directory named `MetaTrimX_Results_[Timestamp]`.

* **01_Demux/**: Raw FastQ files split by sample.
* **02_Trimmed/**: Data with primers and adapters removed.
* **03_Merged/**: Paired-end reads merged into single amplicons.
* **Clustering_Results/**:
  * `OTU_Table.txt`: The final count table.
  * `otus.fasta` / `zotus.fasta`: Representative sequences.
  * `AI_Classification.csv`: (If ML Active) "Real" vs "Artifact" scores for every sequence.
* **MetaTrimX_Interactive_Report.html**: A full visual dashboard of run quality.

---

# 🧬 MetaTrimX Workflow Diagram

![MetaTrimX Dual-Mode Workflow](workflow.png)

---

## 🤝 Social & Community

MetaTrimX is designed to be **shared, discussed, and tweaked** in real projects, not isolated in a black box.

- **Reproducible ML Configs**: Every ML decision (Mode 1 predictions and Mode 2 ensemble scores) exported as JSON/CSV—easy to share exact runs with collaborators, supervisors, or reviewers.

- **Human-in-the-Loop**: Mode 1 writes human-readable Python scripts, not hidden binaries. Labs can fork, comment, and version-control their favorite parameter profiles.

- **Lab-to-Lab Knowledge Transfer**: The RF "brain" (`metatrimx_brain.pkl`) can be retrained or swapped. Institutes maintain their own tuned brain and share it across teams.

- **Community-Driven Evolution**: Open an issue on GitHub with your dataset, ML predictions, and verdicts. The ensemble and brain get retrained on real edge cases (extreme GC, ultra-short markers, rare-species enrichment) so everyone benefits.

---

## 🔬 Use Cases

### **Marine Metabarcoding (12S/Fish)**
Identify fish species and diet composition from water samples or stomach contents using 12S mitochondrial DNA barcoding.

### **Bacterial/Archaeal Surveys (16S)**
Profile microbial communities with automatic parameter optimization for your specific quality characteristics.

### **Fungal Ecology (ITS)**
Analyze fungal diversity in soils, plants, or microbiomes with ITS1/ITS2 marker support.

### **Arthropod Detection (COI)**
Environmental DNA surveys for insects, crustaceans, and other arthropods in aquatic or terrestrial systems.

---

## 🐛 Troubleshooting

### Common Issues

**Issue**: "ModuleNotFoundError: No module named 'sklearn'"
- **Solution**: Install scikit-learn: `pip install scikit-learn`

**Issue**: "fastp command not found"
- **Solution**: Install via conda: `conda install -c bioconda fastp`

**Issue**: "FastQ parsing error" with your custom files
- **Solution**: Ensure FastQ files are properly formatted (4-line format). Check with: `head -4 your_file.fastq`

**Issue**: Out of memory on large datasets
- **Solution**: Increase available RAM or reduce read subset via `--max_reads` parameter

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 👨‍💻 Authors & Acknowledgments

**Original Development**: SUBRAMANIAM VIJAYAKUMAR 

**ML Model Training**: Trained on real metabarcoding datasets from marine ecosystems obtained from CIIMAR.

**Collaborating Institutions**: Acknowledging contributions from CIIMAR - Interdisciplinary Centre of Marine and Environmental Research of The University of Porto.

---

## 📞 Support & Contact

- **Email**: vijayakumar.subraman@mail.ucv.es

---

## 🔮 Roadmap

- [ ] GPU-accelerated clustering (Mode 2)
- [ ] Integration with QIIME 2 ecosystem
- [ ] Web-based dashboard (standalone or Docker)
- [ ] Multi-marker gene analysis in single pipeline
- [ ] Real-time parameter tuning based on run history
- [ ] Community brain re-training pipeline
- [ ] Extended documentation with video tutorials

---

**Last Updated**: December 2025 | **Version**: 3.0.0 (Stable)

