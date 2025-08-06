#!/bin/bash

# Define the file paths
file_path="/mnt/c/Users/User/Downloads/MetaTrimX/metatrimx/Test"
forward_read="${file_path}/forward_reads.fastq"
reverse_read="${file_path}/reverse_reads.fastq"

# Define quality filtering options
quality_cutoff=20  # Minimum quality score for trimming
min_length=100      # Minimum read length after trimming

# Define primer sequences and tags for each sample
declare -A forward_primers_1
declare -A reverse_primers_1
declare -A forward_primers_2
declare -A reverse_primers_2
declare -A sample_tags

# Add the primer sequences and full tag sequence for each sample
forward_primers_1["Seal_dietDes10a"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Seal_dietDes10a"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Seal_dietDes10a"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Seal_dietDes10a"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Seal_dietDes10a"]="cactct"

forward_primers_1["Seal_dietDes08a"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Seal_dietDes08a"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Seal_dietDes08a"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Seal_dietDes08a"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Seal_dietDes08a"]="ggtaag"

# Path to the MetaTrimX script
metatrimx_script="/mnt/c/Users/User/Downloads/MetaTrimX/metatrimx/core/metatrimx.sh"

# Check if the main processing script exists
if [[ ! -f "$metatrimx_script" ]]; then
    echo "Error: The main processing script '$metatrimx_script' is missing!"
    exit 1
fi

# Function to run MetaTrimX for a sample
run_metatrimx() {
    sample="$1"
    tag="${sample_tags[$sample]}"

    if [[ -z "$tag" ]]; then
        echo " Error: Tag for sample '$sample' is empty. Skipping..."
        return
    fi

    echo "\n--- Starting: $sample (tag: $tag) ---"

    bash "$metatrimx_script" \
        "$file_path" \
        "$forward_read" "$reverse_read" \
        "${forward_primers_1[$sample]}" "${reverse_primers_1[$sample]}" \
        "${forward_primers_2[$sample]}" "${reverse_primers_2[$sample]}" \
        "$tag" \
        "$quality_cutoff" "$min_length"

    if [[ $? -eq 0 ]]; then
        echo "--- Completed: $sample ---"
    else
        echo "*** Error in processing: $sample ***"
    fi
}

# Max number of parallel jobs (tweak if needed)
MAX_JOBS=3

# Counter for jobs
job_count=0

# Loop through samples
for sample in "${!forward_primers_1[@]}"; do
    run_metatrimx "$sample" &
    ((job_count++))

    # If we've hit the max parallel jobs, wait for all to finish
    if [[ $job_count -ge $MAX_JOBS ]]; then
        wait
        job_count=0
    fi
done

# Final wait in case remaining jobs are still running
wait

echo "All MetaTrimX sample jobs completed."
