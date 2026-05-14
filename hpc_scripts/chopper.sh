#!/bin/bash
#SBATCH --job-name=Chopper_merged
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=orion
#SBATCH --output=chopper_%j.out
##SBATCH --out slurm-chopper-%A_%a.out

#################################
# Usage:
# sbatch chopper_single_merged.sh <input.fastq.gz>
#################################

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: sbatch $0 <merged.fastq.gz>"
    exit 1
fi

INPUT_FILE="$1"

if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: File not found: $INPUT_FILE"
    exit 1
fi

echo "Starting Chopper QC"
date
echo "Input file: $INPUT_FILE"

#################################
# Load Chopper from ONPTools
#################################

export PATH=/path/to/chopper

echo "Chopper version:"
chopper --version || true

#################################
# Scratch setup
#################################

SCRATCH_DIR="/mnt/SCRATCH/$USER/chopper_single"
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

#################################
# File handling
#################################

BASENAME=$(basename "$INPUT_FILE" .fastq.gz)
OUTDIR=$(dirname "$INPUT_FILE")

echo "Base name: $BASENAME"
echo "Output dir: $OUTDIR"

cp "$INPUT_FILE" "$BASENAME.fq.gz"
pigz -d "$BASENAME.fq.gz"

#################################
# Run Chopper
#################################

echo "Running Chopper..."
date

chopper \
  --threads $SLURM_CPUS_ON_NODE \
  -q 10 \
  -l 1000 \
  < "$BASENAME.fq" \
  > "$BASENAME.cleaned.fq"

#################################
# Compress & copy back
#################################

pigz -p $SLURM_CPUS_ON_NODE "$BASENAME.cleaned.fq"

rsync -av "$BASENAME.cleaned.fq.gz" "$OUTDIR/"

echo "Chopper finished successfully"
date

#################################
# Cleanup
#################################

rm -rf "$SCRATCH_DIR"