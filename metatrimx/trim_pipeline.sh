#!/bin/bash

# Define the file path where the input files are located 
file_path="# User provides the folder path as the first argument"
forward_read="${file_path}/forward read .fastq"
reverse_read="${file_path}/reverse read .fastq"

# Define quality filtering options
quality_cutoff=25  # Minimum quality score for trimming
min_length=100      # Minimum read length after trimming

# Define primer sequences and tags for each sample
declare -A forward_primers_1
declare -A reverse_primers_1
declare -A forward_primers_2
declare -A reverse_primers_2
declare -A sample_tags

# Add the primer sequences and tags for each sample
forward_primers_1["sample_name"]="Forward primer"
reverse_primers_1["sample_name "]="Reverse primer "
forward_primers_2["sample_name"]="Forward primer"
reverse_primers_2["sample_name"]="Reverse primer"
sample_tags["sample_name"]="tags"

# Path to the MetaTrimX script
metatrimx_script="/mnt/c/Users/User/Downloads/MetaTrimX/metatrimx/core/metatrimx.sh"

# Check if the main processing script exists
if [[ ! -f "$metatrimx_script" ]]; then
    echo "Error: The main processing script '$metatrimx_script' is missing!"
    exit 1
fi

# Loop through the samples and process each one using the external script
for sample in "${!forward_primers_1[@]}"; do
    echo "Processing sample: $sample"

    # Run the metatrimx script with necessary parameters
    bash "$metatrimx_script" \
        "$file_path" \
        "$forward_read" "$reverse_read" \
        "${forward_primers_1[$sample]}" "${reverse_primers_1[$sample]}" \
        "${forward_primers_2[$sample]}" "${reverse_primers_2[$sample]}" \
        "${sample_tags[$sample]}" \
        "$quality_cutoff" "$min_length"

    # Check if the script executed successfully
    if [[ $? -eq 0 ]]; then
        echo "Sample $sample processed successfully."
    else
        echo "Error processing sample $sample."
    fi
done

echo "MetaTrimX execution completed."
