#!/bin/bash
#SBATCH --job-name=coverm
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x_%A.out

# =============================================================================
# CoverM genome coverage script
# Usage:
#   sbatch coverm.sh <genome_dir> <reads.fq.gz> <output.tsv> [fasta_extension]
#
# Arguments:
#   $1  Path to directory containing genome FASTA files
#   $2  Path to single-end reads (fastq.gz)
#   $3  Output TSV file name/path
#   $4  (Optional) FASTA file extension, default: fasta
#
# Example:
#   sbatch coverm.sh \
#     /path/to/dereplicated_genomes \
#     /path/to/reads.fq.gz \
#     freshwater_coverm.tsv \
#     fasta
# =============================================================================

# --- Argument parsing ---------------------------------------------------------
GENOME_DIR="${1}"
READS_DIR="${2}"
OUTPUT="${3}"
FASTA_EXT="${4:-fasta}"   # Default to 'fasta' if not provided

# --- Validate inputs ----------------------------------------------------------
if [[ -z "$GENOME_DIR" || -z "$READS_DIR" || -z "$OUTPUT" ]]; then
    echo "ERROR: Missing required arguments."
    echo "Usage: sbatch coverm.sh <genome_dir> <reads_dir> <output.tsv> [fasta_extension]"
    exit 1
fi

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "ERROR: Genome directory not found: $GENOME_DIR"
    exit 1
fi

if [[ ! -d "$READS_DIR" ]]; then
    echo "ERROR: Reads directory not found: $READS_DIR"
    exit 1
fi

# --- Log run info -------------------------------------------------------------
echo "============================================"
echo "CoverM run started: $(date)"
echo "Job ID:             $SLURM_JOB_ID"
echo "Genome directory:   $GENOME_DIR"
echo "FASTA extension:    $FASTA_EXT"
echo "Reads:              $READS_DIR"
echo "Reads found:        $(ls "$READS_DIR"/*.gz | wc -l) files"
echo "Output:             $OUTPUT"
echo "Threads:            $SLURM_CPUS_PER_TASK"
echo "============================================"

# --- Activate conda environment -----------------------------------------------
source ~/miniforge3/etc/profile.d/conda.sh
conda activate coverm_env

# --- Run CoverM ---------------------------------------------------------------
coverm genome \
    --genome-fasta-directory "$GENOME_DIR" \
    --genome-fasta-extension "$FASTA_EXT" \
    --single "$READS_DIR"/*.gz \
    --threads "$SLURM_CPUS_PER_TASK" \
    -p minimap2-ont \
    --methods relative_abundance covered_fraction \
    --min-covered-fraction 0 \
    --output-file "$OUTPUT"

# --- Exit status check --------------------------------------------------------
EXIT_CODE=$?
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "============================================"
    echo "CoverM finished successfully: $(date)"
    echo "Output written to: $OUTPUT"
    echo "============================================"
else
    echo "============================================"
    echo "ERROR: CoverM failed with exit code $EXIT_CODE"
    echo "============================================"
    exit $EXIT_CODE
fi