#!/bin/bash
###################################
#SBATCH --job-name=dRep
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=200G
#SBATCH -p orion
#SBATCH --output=slurm-%x_%A.out
###################################

set -euo pipefail

#####################
# Inputs
#####################
if [[ $# -ne 6 ]]; then
  echo "Usage: $0 INPUT_NAME INDIR EXT COMP CON OUTDIR" >&2
  exit 2
fi

input=$1     # e.g. saltwater
indir=$2     # directory with MAGs
ext=$3       # fasta extension (e.g. fasta)
comp=$4      # completeness threshold
con=$5       # contamination threshold
outdir=$6    # output directory

RSYNC='rsync -aLhv --no-perms --no-owner --no-group'

#####################
# Modules + Conda
#####################
module --quiet purge
module load Anaconda3/2021.11

source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null

CONDA_ENV="/path/to/drep"

conda activate "$CONDA_ENV"
echo "🧬 Using conda env: $CONDA_PREFIX"

command -v dRep || { echo "ERROR: dRep not found in PATH"; exit 1; }

#####################
# Debug info
#####################
echo "👋 Hello $USER"
echo "📌 submit dir: $SLURM_SUBMIT_DIR"
echo "🧾 job ID: $SLURM_JOB_ID"
echo "🖥️ node: $SLURM_NODELIST"
echo "🧵 CPUs: $SLURM_CPUS_ON_NODE"
echo "📅 date:"; date
echo "📥 INDIR:  $indir"
echo "📤 OUTDIR: $outdir"

#####################
# Local scratch
#####################
TMPROOT="/work/users"
mkdir -p "$TMPROOT/$USER"
cd "$TMPROOT/$USER"

WORKDIR="drep.${input}.${SLURM_JOB_ID}"
mkdir "$WORKDIR"
cd "$WORKDIR"

LOCALSCRATCH=$(pwd)
echo "🧪 Working dir: $LOCALSCRATCH"

# Always sync results back
cleanup() {
  set +e
  echo "📦 Syncing results to $outdir"
  mkdir -p "$outdir"
  $RSYNC "$input.DREP.$comp.$con.out" "$outdir/" || true
}
trap cleanup EXIT

#####################
# Stage MAGs
#####################
echo "📥 Copying MAGs"
time $RSYNC "$indir"/*."$ext" .

#####################
# Remove unwanted files
#####################
echo "🧹 Removing unbinned / junk MAGs"

rm -f *unbin*.fasta \
      *contigs.fasta \
      *lowDepth*.fasta \
      *tooShort*.fasta || true

echo "🧮 Number of genomes to dRep:"
ls -1 *."$ext" 2>/dev/null | wc -l

#####################
# Organize MAGs
#####################
mkdir "$input.MAGs"
mv *."$ext" "$input.MAGs/"

#####################
# Run dRep
#####################
echo "🚀 Starting dRep"
date +%d\ %b\ %T

time dRep dereplicate \
  "$input.DREP.$comp.$con.out" \
  -g "$input.MAGs"/*."$ext" \
  -p "$SLURM_CPUS_ON_NODE" \
  -comp "$comp" \
  -con "$con"

#####################
# Done
#####################
echo "✅ dRep finished"
echo "📁 Results: $outdir/$input.DREP.$comp.$con.out"
date