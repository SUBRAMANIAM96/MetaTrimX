#!/bin/bash

# Check for correct number of arguments
if [[ $# -ne 10 ]]; then
    echo "Usage: $0 <file_path> <forward_read> <reverse_read> <forward_primer_1> <reverse_primer_1> <forward_primer_2> <reverse_primer_2> <sample_tag> <quality_cutoff> <min_length>"
    exit 1
fi

# Assign input parameters
file_path="$1"
forward_read="$2"
reverse_read="$3"
fwd1="$4"
rev1="$5"
fwd2="$6"
rev2="$7"
tag="$8"
quality="$9"
min_len="${10}"

# Output file paths
demux_R1="${file_path}/demux_${tag}_F.fastq"
demux_R2="${file_path}/demux_${tag}_R.fastq"
trimmed_R1="${file_path}/trimmed_${tag}_F.fastq"
trimmed_R2="${file_path}/trimmed_${tag}_R.fastq"
clean_R1="${file_path}/no_adapters_${tag}_F.fastq"
clean_R2="${file_path}/no_adapters_${tag}_R.fastq"
merged="${file_path}/merged_${tag}.fastq"
filtered="${file_path}/filtered_${tag}.fastq"
fasta="${file_path}/filtered_${tag}.fasta"
derep="${file_path}/dereplicated_${tag}.fasta"
nonchim="${file_path}/nonchim_${tag}.fasta"
centroids="${file_path}/otu_centroids_${tag}.fasta"
otu_map="${file_path}/otu_map_${tag}.uc"

# Check for required tools
for tool in cutadapt vsearch seqtk; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool not found in PATH"
        exit 1
    fi
done

echo "[1] Demultiplexing by tag: $tag"
cutadapt -j 12 -g "^${tag}" -G "^${tag}" --action=trim \
    -o "$demux_R1" -p "$demux_R2" \
    "$forward_read" "$reverse_read"

if [[ $? -ne 0 ]]; then
    echo "Error: Demultiplexing failed for $tag"
    exit 1
fi

echo "[2] Primer trimming (2 primer sets)"
cutadapt -j 12 -q "$quality" --minimum-length "$min_len" \
    --pair-filter=any \
    -g "$fwd1" -G "$rev1" \
    -g "$fwd2" -G "$rev2" \
    -o "$trimmed_R1" -p "$trimmed_R2" \
    "$demux_R1" "$demux_R2"

if [[ $? -ne 0 ]]; then
    echo "Error: Primer trimming failed for $tag"
    exit 1
fi

echo "[3] Adapter removal"
cutadapt -j 12 -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
    -o "$clean_R1" -p "$clean_R2" \
    "$trimmed_R1" "$trimmed_R2"

if [[ $? -ne 0 ]]; then
    echo "Error: Adapter removal failed for $tag"
    exit 1
fi

echo "[4] Merging paired reads..."
vsearch --fastq_mergepairs "$clean_R1" \
    --reverse "$clean_R2" \
    --fastqout "$merged" \
    --fastq_minovlen 20 \
    --fastq_maxdiffs 5 \
    --relabel all \
    --threads 4

if [[ ! -s "$merged" ]]; then
    echo "Error: Merging failed for $tag"
    exit 1
fi

read_count=$(grep -c "^@" "$merged")
echo "[INFO] Merged read count for $tag: $read_count"

echo "[5] Quality filtering..."
vsearch --fastq_filter "$merged" \
    --fastqout "$filtered" \
    --fastq_minlen 150 \
    --fastq_maxee 1.0 \
    --threads 4

echo "[6] Convert to FASTA..."
seqtk seq -A "$filtered" > "$fasta"

echo "[7] Dereplicating..."
vsearch --derep_fulllength "$fasta" \
    --output "$derep" \
    --sizeout \
    --threads 4

echo "[8] Chimera removal..."
vsearch --uchime_denovo "$derep" \
    --nonchimeras "$nonchim" \
    --threads 4

echo "[9] OTU clustering..."
vsearch --cluster_fast "$nonchim" \
    --id 0.97 \
    --centroids "$centroids" \
    --uc "$otu_map" \
    --relabel OTU_

echo " Processing completed for sample: $tag"
