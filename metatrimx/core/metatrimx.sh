#!/bin/bash

# Check for the correct number of arguments
if [[ $# -ne 10 ]]; then
    echo "Usage: $0 <file_path> <forward_read> <reverse_read> <forward_primer_1> <reverse_primer_1> <forward_primer_2> <reverse_primer_2> <sample_tag> <quality_cutoff> <min_length>"
    exit 1
fi

# Assign input parameters
file_path="$1"
forward_read="$2"
reverse_read="$3"
forward_primer_1="$4"
reverse_primer_1="$5"
forward_primer_2="$6"
reverse_primer_2="$7"
sample_tag="$8"
quality_cutoff="$9"
min_length="${10}"

# Check if required tools are installed
for tool in cutadapt vsearch; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is not installed. Please install it."
        exit 1
    fi
done

# Define output files
trimmed_forward="${file_path}/trimmed_${sample_tag}_F.fastq"
trimmed_reverse="${file_path}/trimmed_${sample_tag}_R.fastq"
no_adapters_forward="${file_path}/no_adapters_${sample_tag}_F.fastq"
no_adapters_reverse="${file_path}/no_adapters_${sample_tag}_R.fastq"
demux_forward="${file_path}/demux_${sample_tag}_F.fastq"
demux_reverse="${file_path}/demux_${sample_tag}_R.fastq"
merged_output="${file_path}/merged_${sample_tag}.fastq"
unique_output="${file_path}/unique_${sample_tag}.fastq"
denovo_derep="${file_path}/derep_${sample_tag}.fastq"
chimera_free="${file_path}/no_chimera_${sample_tag}.fastq"
otu_output="${file_path}/otu_${sample_tag}.txt"

# Step 1: Trim reads
cutadapt -q "$quality_cutoff" --minimum-length "$min_length" \
    --pair-filter=any \
    -g "$forward_primer_1" -G "$reverse_primer_1" \
    -a "$forward_primer_1" -A "$reverse_primer_1" \
    -g "$forward_primer_2" -G "$reverse_primer_2" \
    -a "$forward_primer_2" -A "$reverse_primer_2" \
    -o "$trimmed_forward" -p "$trimmed_reverse" \
    "$forward_read" "$reverse_read"

# Step 2: Remove remaining adapters
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
    -o "$no_adapters_forward" -p "$no_adapters_reverse" \
    "$trimmed_forward" "$trimmed_reverse"

# Step 3: Demultiplex reads based on the tag
cutadapt -g "$sample_tag" -G "$sample_tag" \
    -o "$demux_forward" -p "$demux_reverse" \
    "$no_adapters_forward" "$no_adapters_reverse"

# Verify demultiplexing output before proceeding
if [[ ! -s "$demux_forward" || ! -s "$demux_reverse" ]]; then
    echo "Error: Demultiplexed files are missing or empty. Exiting."
    exit 1
else
    echo "Sample $sample_tag processed successfully."
fi

# Step 4: Merge paired-end reads using vsearch with 4 cores
vsearch --fastq_mergepairs "$demux_forward" --reverse "$demux_reverse" --fastqout "$merged_output" --threads 4

# Check if merging was successful
if [[ ! -s "$merged_output" ]]; then
    echo "Error: Merging failed. No merged reads found for $sample_tag. Exiting."
    exit 1
fi

# Step 5: Remove duplicate sequences using 4 cores
vsearch --fastx_uniques "$merged_output" --sizeout --fastaout "$unique_output" --threads 12

# Step 6: Denoise sequences using 4 cores
vsearch --cluster_size "$unique_output" --id 0.99 --centroids "$denovo_derep" --threads 4

# Step 7: Remove chimeras using 4 cores
vsearch --uchime_denovo "$denovo_derep" --nonchimeras "$chimera_free" --threads 4

# Step 8: Cluster into OTUs with 99% similarity using 4 cores
vsearch --cluster_fast "$chimera_free" --id 0.97 --centroids "$otu_output" --threads 4

# Final output message
echo "MetaTrimX Pipeline completed successfully for sample: $sample_tag"
echo "Final OTUs are saved in: $otu_output"
