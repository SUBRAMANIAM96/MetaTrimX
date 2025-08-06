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
metatrimx_script="${base_path}/metatrimx/core/metatrimx.sh"

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
=======
#!/bin/bash

# Define the file paths
file_path="/path/to/your/data"
forward_read="${file_path}/forward read .fastq"
reverse_read="${file_path}/reverse read .fastq"

# Define quality filtering options
quality_cutoff=20  # Minimum quality score for trimming
min_length=50      # Minimum read length after trimming

# Define primer sequences and tags for each sample (2 sets of primers per sample)
declare -A forward_primers_1
declare -A reverse_primers_1
declare -A forward_primers_2
declare -A reverse_primers_2
declare -A sample_tags

# Add the primer sequences and full tag sequence for each sample
forward_primers_1["Sample_name"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Sample_name"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Sample_name"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Sample_name"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Sample_name"]="cactct" 

# Path to the MetaTrimX script
metatrimx_script="{File_path}/metatrimx/core/metatrimx.sh" # User provides the path 

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
