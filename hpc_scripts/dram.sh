#!/bin/bash
############## SLURM SCRIPT ###################################
#SBATCH --job-name=DRAM
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH -p orion,RStudio
#SBATCH --output=slurm-%x_%A.out
###############################################################

set -euo pipefail

#####################
# Inputs
#####################
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 INDIR EXT OUTDIR" >&2
  echo "Example: sbatch $0 /path/to/mags fasta /path/to/results" >&2
  exit 2
fi

indir=$1     # directory containing MAGs
ext=$2       # extension (e.g. fasta, fa, fna)
outdir=$3    # final output directory

# Your DRAM conda env on ORION
CONDA_ENV="/path/to/dram"

RSYNC='rsync -aLhv --no-perms --no-owner --no-group'

#####################
# Conda (ORION)
#####################

module load Miniconda3
eval "$(conda shell.bash hook)"

conda activate "$CONDA_ENV"

echo "I am working with this conda env: $CONDA_PREFIX"
command -v DRAM.py >/dev/null 2>&1 || { echo "ERROR: DRAM.py not found in env PATH" >&2; exit 10; }

#####################
# Debug info
#####################
echo "👋 Hello $USER"
echo "📌 submit dir: $SLURM_SUBMIT_DIR"
echo "🧾 job: ${SLURM_JOB_ID}"
echo "🖥️ nodelist: $SLURM_NODELIST"
echo "🧵 cpus: $SLURM_CPUS_ON_NODE"
echo "📅 date:"; date
echo "📥 INDIR:  $indir"
echo "🧬 ext:    $ext"
echo "📤 OUTDIR: $outdir"

#####################
# Local scratch setup
#####################
TMPROOT="/work/users"
mkdir -p "$TMPROOT/$USER"
LOCALSCRATCH="$TMPROOT/$USER/dram.${SLURM_JOB_ID}"
mkdir -p "$LOCALSCRATCH"

echo "🧪 Working dir (LOCALSCRATCH): $LOCALSCRATCH"
cd "$LOCALSCRATCH"

# Output folder name on scratch (and synced to outdir)
RESULTS_DIR="DRAM.Results.dir"

cleanup() {
  set +e
  echo
  echo "📦 Syncing results to $outdir ..."
  mkdir -p "$outdir"
  if [[ -d "$LOCALSCRATCH/$RESULTS_DIR" ]]; then
    $RSYNC "$LOCALSCRATCH/$RESULTS_DIR/" "$outdir/$RESULTS_DIR/" || true
  else
    echo "⚠️ No results dir found at $LOCALSCRATCH/$RESULTS_DIR"
  fi
  echo "🧹 Cleaning scratch: $LOCALSCRATCH"
  rm -rf "$LOCALSCRATCH" || true
}
trap cleanup EXIT

#####################
# Stage MAGs to node
#####################
echo "📥 Copying MAGs → node (.$ext)"
shopt -s nullglob
src_files=( "$indir"/*."$ext" )
if (( ${#src_files[@]} == 0 )); then
  echo "ERROR: No files found matching: $indir/*.$ext" >&2
  exit 3
fi

mkdir -p MAGs
time $RSYNC "$indir"/*."$ext" ./MAGs/

echo "Number of genomes to annotate:"
ls -1 MAGs/*."$ext" | wc -l

#####################
# Run DRAM annotate
#####################

echo "🚀 DRAM annotate start:"
date +%d\ %b\ %T

time DRAM.py annotate \
  -i "MAGs/*.$ext" \
  -o dram.annotation.dir \
  --min_contig_size 500 \
  --threads "$SLURM_CPUS_ON_NODE"

#####################
# Run DRAM distill
#####################
echo "🧪 DRAM distill start:"
date +%d\ %b\ %T

time DRAM.py distill \
  -i dram.annotation.dir/annotations.tsv \
  -o dram.genome_summaries.dir \
  --trna_path dram.annotation.dir/trnas.tsv \
  --rrna_path dram.annotation.dir/rrnas.tsv

echo "✅ DRAM finished at:"
date +%d\ %b\ %T

#####################
# Pack results
#####################
mkdir -p "$RESULTS_DIR"
mv dram.annotation.dir "$RESULTS_DIR/"
mv dram.genome_summaries.dir "$RESULTS_DIR/"

echo "📍 Results will be in: $outdir/$RESULTS_DIR"
echo "✅ I've done"; date
