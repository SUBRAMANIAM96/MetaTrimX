# 📌 MetaTrimX  
*A first-of-its-kind automated pipeline for 12S rRNA metabarcoding analysis*  

---

## 🔥 Overview  
MetaTrimX is a **novel** and **fully automated** bioinformatics pipeline specifically designed for **12S rRNA metabarcoding**. It converts **raw paired-end sequencing reads** into high-quality **Operational Taxonomic Units (OTUs) with 99% accuracy**.  

Unlike conventional workflows that require multiple commands at different stages, MetaTrimX simplifies the entire process into a **one-time input command**—no additional user intervention is required!  

### ✅ Features  
✔️ **One-of-a-kind pipeline for 12S rRNA analysis**  
✔️ **Fully automated** from raw reads to OTU tables  
✔️ **Handles multiple samples & multiple primer pairs**  
✔️ **Supports high-throughput sequencing data**  
✔️ **Trimming, demultiplexing, merging, dereplication, deionization in one go**  
✔️ **Specially designed for paired-end reads**  
✔️ **No separate commands required—just input primer & tag information**  

### 🚀 What MetaTrimX Does:  
🔹 **Removes adapters & primers** (forward & reverse)  
🔹 **Demultiplexes** cleaned amplicons based on sample-assigned tags  
🔹 **Merges forward & reverse reads** for downstream analysis  
🔹 **Dereplicates & deionizes** reads  
🔹 **Generates high-quality OTUs** for final analysis  

---

## 🛠 Installation  

### **Prerequisites**  
Ensure the following dependencies are installed:  
- **Python** (>=3.6)  
- **Cutadapt** (>=4.0)  
- **Vsearch**  
- **Bash Shell**  

You can install Cutadapt and Vsearch using:  
```bash
pip install cutadapt
conda install -c bioconda vsearch


##Install MetaTrimX
git clone https://github.com/SUBRAMANIAM96/MetaTrimX.git
cd MetaTrimX
pip install .


##Project Structure
MetaTrimX/
│── metatrimx/
│   ├── cli.py                 # CLI interface for running the pipeline
│   ├── core/
│   │   ├── metatrimx.sh        # Interface script (DO NOT MODIFY)
│   ├── trim_pipeline.sh        # User-editable script (Modify as needed)
│   ├── __init__.py             # Python package initialization
│   ├── Test/                   # Test data directory
│── setup.py                    # Installation script
│── README.md                   # Project documentation



## 🚀 Running the Pipeline  

MetaTrimX requires only a **single command** to process your sequencing data:  

```bash
bash metatrimx/trim_pipeline.sh



📝 Editing trim_pipeline.sh

The trim_pipeline.sh script is the main user-editable file where you must define:

    Primer sequences (Forward & Reverse)
    Tag information (for demultiplexing)
    Cutadapt & Vsearch parameters
   
    To edit the script, use:

nano metatrimx/trim_pipeline.sh





⚙️ How It Works

1️⃣ The trim_pipeline.sh script contains all necessary parameters, including primers, tags, and settings.
2️⃣ It calls metatrimx.sh, which acts as a fixed interface to execute Cutadapt and Vsearch commands.
3️⃣ The pipeline runs automatically without requiring any further inputs.




📂 Output Structure

MetaTrimX generates sequential output files at different stages:

📁 Trimmed Reads:

    trimmed_F.fastq (Forward reads)
    trimmed_R.fastq (Reverse reads)

📁 Demultiplexed Reads:

    demultiplexed_F.fastq (Forward reads per sample)
    demultiplexed_R.fastq (Reverse reads per sample)

📁 Merged Reads:

    merged_<tag>.fastq (Merged reads with assigned sample tag)

📁 Dereplicated Reads:

    dereplicated_<tag>.fastq

📁 Deionized Reads:

    deionized_<tag>.fastq

📁 Final OTU Output:

    otu_<tag>.fastq



## 🧪 Running the Test

A **test dataset** is included in the `test/` folder, which contains:  

✅ **`forward_reads.fastq`** – Sample forward reads  
✅ **`reverse_reads.fastq`** – Sample reverse reads  
✅ **`test_trim_pipeline.sh`** – Pre-configured test pipeline script  

To run the test, navigate to the `test/` folder and execute the test pipeline:  

```bash
cd metatrimx/test
bash test_trim_pipeline.sh
## Extracting FASTQ Files

To ensure the extracted FASTQ files remain in the same folder, use:

```bash
gunzip -k forward_reads.fastq.gz reverse_reads.fastq.gz



📜 License

MetaTrimX is licensed under the MIT License. Feel free to use, modify, and contribute!
👥 Contributors



Developed by: Subramaniam Vijayakumar
📧 Email: subramanyamvkumar@gmail.com
🔗 GitHub: SUBRAMANIAM96
⭐ Support & Contributions


💡 Found a bug? Have a feature request? Open an issue on GitHub Issues
📢 If you like this project, consider starring ⭐ the repository!
🔥 Key Updates in This Version:


✅ No need for separate files for primers & tags—everything is inside trim_pipeline.sh
✅ Clarifies that metatrimx.sh acts as an interface
✅ Explains step-by-step how the pipeline processes data
✅ Detailed breakdown of output files at each stage
✅ Single-command execution without interruptions


