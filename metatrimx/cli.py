import argparse
import os
import subprocess
import sys
import shutil  # For checking external dependencies

def check_dependencies():
    """
    Check if required external tools are installed.
    """
    if not shutil.which("cutadapt"):
        print("❌ Error: Cutadapt is not installed. Please install it manually (e.g., `pip install cutadapt`).")
        sys.exit(1)
    
    if not shutil.which("vsearch"):
        print("❌ Error: Vsearch is not installed. Please install it manually (e.g., `sudo apt install vsearch`).")
        sys.exit(1)

def run_pipeline(input_file, output_dir, min_length, quality_cutoff):
    """
    Run the MetaTrimX pipeline using the trim_pipeline.sh script.
    The user edits trim_pipeline.sh, which then calls metatrimx.sh.
    """
    check_dependencies()  # Ensure Cutadapt and Vsearch are installed

    script_path = os.path.join(os.path.dirname(__file__), "trim_pipeline.sh")

    # Ensure trim_pipeline.sh exists
    if not os.path.exists(script_path):
        print(f"❌ Error: Script not found at {script_path}")
        sys.exit(1)

    # Ensure input file exists
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file '{input_file}' not found.")
        sys.exit(1)

    # Ensure output directory exists, if not, create it
    if not os.path.exists(output_dir):
        print(f"📂 Output directory '{output_dir}' does not exist. Creating it...")
        os.makedirs(output_dir)

    # Ensure trim_pipeline.sh is executable
    os.chmod(script_path, 0o755)

    # Construct command
    cmd = [
        "bash", script_path, input_file, output_dir, str(min_length), str(quality_cutoff)
    ]

    print(f"🚀 Running MetaTrimX pipeline with the following command:\n{' '.join(cmd)}\n")

    try:
        subprocess.run(cmd, check=True)
        print("\n🎉 MetaTrimX has successfully completed! 🎉\n")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error while executing the pipeline: {e}")
        sys.exit(1)

def main():
    """
    Entry point for the CLI tool.
    """
    parser = argparse.ArgumentParser(
        description="MetaTrimX - A CLI tool for automated metabarcoding data processing"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input sequencing file"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output directory"
    )

    parser.add_argument(
        "-l", "--min_length",
        type=int,
        default=100,
        help="Minimum read length to retain (default: 100)"
    )

    parser.add_argument(
        "-q", "--quality",
        type=int,
        default=25,
        help="Quality cutoff for filtering (default: 25)"
    )

    args = parser.parse_args()

    run_pipeline(args.input, args.output, args.min_length, args.quality)

if __name__ == "__main__":
    main()
