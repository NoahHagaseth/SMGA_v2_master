#!/bin/bash
#SBATCH --job-name=NanoPlotQC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=84G
#SBATCH --time=08:00:00
#SBATCH --output=nanoplot_%j.out
#SBATCH --error=nanoplot_%j.err

echo "Starting NanoPlot QC"
echo "Start time: $(date)"

# === Absolute path to NanoPlot ===
NANOPLOT="/path/to/NanoPlot"

# === Arguments ===
INPUT_FASTQ=$1
OUTPUT_DIR=$2

if [ -z "$INPUT_FASTQ" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: sbatch nanoplot_qc.sh <input_fastq.gz> <output_dir>"
    exit 1
fi

echo "Input file: $INPUT_FASTQ"
echo "Output dir:  $OUTPUT_DIR"

# === Create output directory ===
mkdir -p "$OUTPUT_DIR"

# === Run NanoPlot ===
/mnt/orion/conda/CIGENE/ONPTools/bin/NanoPlot \
    --fastq "$INPUT_FASTQ" \
    --outdir "$OUTPUT_DIR" \
    --threads $SLURM_CPUS_PER_TASK \
    --loglength \
    --N50 \
    --plots dot \
    --verbose

echo "NanoPlot finished"
echo "End time: $(date)"