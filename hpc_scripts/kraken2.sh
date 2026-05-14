#!/bin/bash
############## SLURM SCRIPT ###################################
#SBATCH --job-name=kraken2
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH -p orion
#SBATCH --output=slurm-%x_%A.out
###############################################################

set -euo pipefail

#####################
# Usage check
#####################
if [[ $# -ne 3 ]]; then
  echo "Usage:"
  echo "  sbatch kraken2.sh SAMPLE READS OUTDIR"
  echo
  echo "Example:"
  echo "  sbatch kraken2.sh saltwater /path/to/saltwater.chopper.fq.gz /path/to/kraken2_out"
  exit 1
fi

#####################
# Inputs
#####################
SAMPLE="$1"     # saltwater | freshwater | mucus
READS="$2"      # cleaned chopper reads (.fq.gz)
OUTDIR="$3"     # output directory

#Kraken database path
DB="/path/to/kraken2/database"

#####################
# Conda environment
#####################
module --quiet purge
module load Anaconda3/2021.11

source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null

conda activate /path/to/kraken2

echo "🧬 Using conda env: $CONDA_PREFIX"
command -v kraken2 >/dev/null || { echo "ERROR: kraken2 not found"; exit 2; }

#####################
# Debug info
#####################
echo "👋 Hello $USER"
echo "🧾 Job ID: $SLURM_JOB_ID"
echo "🖥️ Node: $SLURM_NODELIST"
echo "🧵 CPUs: $SLURM_CPUS_ON_NODE"
echo "💾 Memory: 200G"
echo "📅 Date:"; date
echo "🧬 Sample: $SAMPLE"
echo "📥 Reads: $READS"
echo "🧬 Database: $DB"
echo "📤 Output dir: $OUTDIR"

#####################
# Safety checks
#####################
if [[ ! -s "$READS" ]]; then
  echo "ERROR: Reads file not found: $READS" >&2
  exit 3
fi

if [[ ! -d "$DB" ]]; then
  echo "ERROR: Kraken2 DB not found: $DB" >&2
  exit 4
fi

mkdir -p "$OUTDIR"

#####################
# Run Kraken2
#####################
echo "🚀 Starting Kraken2 for $SAMPLE"
date +%d\ %b\ %T

kraken2 \
  --db "$DB" \
  --threads "$SLURM_CPUS_ON_NODE" \
  --memory-mapping \
  --output "$OUTDIR/${SAMPLE}.kraken2.out" \
  --report "$OUTDIR/${SAMPLE}.kraken2.report.tsv" \
  "$READS"

echo "✅ Kraken2 finished for $SAMPLE"
date +%d\ %b\ %T

echo "📍 Results:"
echo " - $OUTDIR/${SAMPLE}.kraken2.out"
echo " - $OUTDIR/${SAMPLE}.kraken2.report.tsv"