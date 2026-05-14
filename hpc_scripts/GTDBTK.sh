#!/bin/bash
#########################################################################
# SLURM script for running GTDB-Tk classify_wf on Orion cluster
#
# - Stages MAGs to local scratch for speed
# - Runs gtdbtk classify_wf locally
# - Syncs results back to an output directory
#
# Updated: Feb 2026
# Conda env: /mnt/orion/conda/CIGENE/GTDBTK_2.6.1
#
# Usage:
#   sbatch gtdbtk.sh INDIR EXT OUTDIR
#########################################################################

############## SLURM SCRIPT ###################################
#SBATCH --job-name=gtdbtk_classifywf
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH -p orion,RStudio
#SBATCH --output=slurm-%x_%A.out
###############################################################

set -euo pipefail

#####################
# Inputs
#####################
print_usage() {
  echo "Usage: sbatch $0 INDIR EXT OUTDIR" >&2
  echo "Example: sbatch $0 /path/to/mags fna /path/to/results" >&2
}

if [[ $# -ne 3 ]]; then
  print_usage
  exit 2
fi

indir=$1     # directory containing MAGs
ext=$2       # extension (e.g. fasta, fa, fna)
outdir=$3    # final output directory to copy results into

#####################
# Conda (Orion Miniforge template style)
#####################
module --quiet purge

# If Orion provides a standard conda init, keep the miniforge style from your template.
# If your site uses a different base, change CONDA_BASE accordingly.
CONDA_BASE="/mnt/users/auve/miniforge3"
source "$CONDA_BASE/etc/profile.d/conda.sh"
export PATH="$CONDA_BASE/bin:$PATH"

CONDA_ENV="/path/to/gtdbtk"

conda deactivate >/dev/null 2>&1 || true
conda activate "$CONDA_ENV"

echo "I am working with this conda env: $CONDA_PREFIX"
command -v gtdbtk >/dev/null 2>&1 || { echo "ERROR: gtdbtk not found in env PATH" >&2; exit 10; }
gtdbtk --version || true

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
echo "🧬 EXT:    $ext"
echo "📤 OUTDIR: $outdir"

#####################
# Local scratch setup
#####################
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'

TMPROOT="/work/users"
mkdir -p "$TMPROOT/$USER"
LOCALSCRATCH="$TMPROOT/$USER/gtdbtk.${SLURM_JOB_ID}"
mkdir -p "$LOCALSCRATCH"

echo "🧪 Working dir (LOCALSCRATCH): $LOCALSCRATCH"
cd "$LOCALSCRATCH"

# Names
MAGDIR="MAGs"
OUTNAME="MAGs_gtdbtk.dir"

cleanup() {
  set +e
  echo
  echo "📦 Syncing results to $outdir ..."
  mkdir -p "$outdir"
  if [[ -d "$LOCALSCRATCH/$OUTNAME" ]]; then
    $RSYNC "$LOCALSCRATCH/$OUTNAME/" "$outdir/$OUTNAME/" || true
  fi
  echo "🧹 Cleaning scratch: $LOCALSCRATCH"
  rm -rf "$LOCALSCRATCH" || true
}
trap cleanup EXIT

#####################
# Stage MAGs to node
#####################
echo "📥 Copying MAGs files → node (*.$ext)"

shopt -s nullglob
src_files=( "$indir"/*."$ext" )
if (( ${#src_files[@]} == 0 )); then
  echo "ERROR: No files found matching: $indir/*.$ext" >&2
  exit 3
fi

time $RSYNC "$indir"/*."$ext" ./

mkdir -p "$MAGDIR"
mv *."$ext" "$MAGDIR/"

echo "Number of genomes to GTDB-Tk:"
ls -1 "$MAGDIR"/*."$ext" | wc -l

#####################
# Run GTDB-Tk classify_wf
#####################
echo "🚀 start gtdbtk classify_wf at:"
date +%d\ %b\ %T

time gtdbtk classify_wf \
  --genome_dir "$MAGDIR" \
  --out_dir "$OUTNAME" \
  -x "$ext" \
  --cpus "$SLURM_CPUS_ON_NODE"

echo "✅ GTDB-Tk finished at:"
date +%d\ %b\ %T

echo "📍 Results will be in: $outdir/$OUTNAME"
