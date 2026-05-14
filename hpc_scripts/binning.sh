#!/bin/bash
############## SLURM SCRIPT ###################################
#SBATCH --job-name=Binning
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=300G
#SBATCH -p orion,RStudio
#SBATCH --output=slurm-%x_%A.out
###############################################################
#
# Run MetaBAT2 and MaxBin2 in parallel on an ONT assembly
#
# USAGE
#   sbatch slurm_binning_metabat2_maxbin2.sh \
#          READSDIR MEDAKADIR OUTDIR [CONDA_ENV]
#
# ARGUMENTS
#   READSDIR   : directory containing *.fq.gz (e.g., *.chopper.fq.gz)
#   MEDAKADIR  : directory that contains the assembly file: consensus.fasta
#   OUTDIR     : output directory for final results
#   CONDA_ENV  : optional; path to conda/mamba env with tools
#                 default=/net/fs-2/scale/OrionStore/Orion03/conda/CIGENE/BINNING_TOOLS
#
# TOOLS (expected in env)
#   minimap2, samtools, jgi_summarize_bam_contig_depths (MetaBAT2), metabat2,
#   MaxBin2 (run_MaxBin.pl), parallel, pigz
#
set -euo pipefail

#####################
# Inputs
#####################
if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 READSDIR MEDAKADIR OUTDIR [CONDA_ENV]" >&2
  exit 2
fi
READSDIR=$1          # directory with reads (*.fq.gz)
MEDAKADIR=$2         # directory that contains consensus.fasta
OUTDIR=$3            # where final Binning.dir is synced
CONDA_ENV=${4:-path/to/binning_tools}}

RSYNC='rsync -aLhv --no-perms --no-owner --no-group'

#####################
# Modules + Conda/Mamba
#####################
module --quiet purge

#####################
# Conda (your Miniforge)
#####################
module --quiet purge   # optional: keeps the environment clean

# Point to your Miniforge install
CONDA_BASE="/mnt/users/auve/miniforge3"

# Initialize conda in this non-interactive SLURM shell
source "$CONDA_BASE/etc/profile.d/conda.sh"

# (Optional but nice) ensure conda is the one from Miniforge
export PATH="$CONDA_BASE/bin:$PATH"

# Activate your env
conda deactivate >/dev/null 2>&1 || true
conda activate $CONDA_ENV

echo "I am working with this conda env: $CONDA_PREFIX"

echo "I am working with this env: $CONDA_PREFIX"
command -v minimap2   && minimap2 --version   || true
command -v samtools   && samtools --version   || true
command -v metabat2   && metabat2 2>/dev/null || true
command -v run_MaxBin.pl >/dev/null 2>&1 || echo "WARNING: MaxBin2 not found in PATH"
command -v jgi_summarize_bam_contig_depths >/dev/null 2>&1 || echo "WARNING: jgi_summarize_bam_contig_depths not in PATH"
command -v parallel >/dev/null 2>&1 || echo "WARNING: GNU parallel not in PATH"

#####################
# Debug info
#####################
echo "👋 Hello $USER"
echo "📌 submit dir: $SLURM_SUBMIT_DIR"
echo "🧾 job: ${SLURM_JOB_ID}"
echo "🖥️ nodelist: $SLURM_NODELIST"
echo "🧵 cpus: $SLURM_CPUS_ON_NODE"
echo "📅 date:"; date
echo "📥 READSDIR: $READSDIR"
echo "🧬 MEDAKADIR: $MEDAKADIR"
echo "📤 OUTDIR:   $OUTDIR"

#####################
# Local scratch setup
#####################
TMPROOT='/work/users'
available_space=$(df -BG "$TMPROOT" | awk 'NR==2 {print $4}' | sed 's/G//')
required_space=300
if [[ "${available_space:-0}" -lt "$required_space" ]]; then
  echo "ERROR: Insufficient space in $TMPROOT. Available: ${available_space}GB, Required: ${required_space}GB." >&2
  exit 1
fi

mkdir -p "$TMPROOT/$USER" && cd "$TMPROOT/$USER"
mkdir -p "binning.${SLURM_JOB_ID}" && cd "binning.${SLURM_JOB_ID}"
LOCALSCRATCH=$(pwd)
echo "🧪 Working dir (LOCALSCRATCH): $LOCALSCRATCH"

# Always try to ship results back on exit
cleanup() {
  set +e
  echo "\n📦 Syncing results to $OUTDIR ..."
  mkdir -p "$OUTDIR"
  $RSYNC "$LOCALSCRATCH/$SAMPLE.Binning.dir/" "$OUTDIR/" || true
}
trap cleanup EXIT

#####################
# Validate inputs and stage files
#####################
shopt -s nullglob
mapfile -t READS < <(printf '%s\n' "$READSDIR"/*.fq.gz "$READSDIR"/*.fastq.gz 2>/dev/null | sort -u)
if (( ${#READS[@]} == 0 )); then
  echo "ERROR: No *.fq.gz or *.fastq.gz files in $READSDIR" >&2
  exit 3
fi

ASSEMBLY_SRC=$(ls "$MEDAKADIR"/*consensus*.fasta 2>/dev/null | head -n 1)

if [[ -z "$ASSEMBLY_SRC" || ! -s "$ASSEMBLY_SRC" ]]; then
  echo "ERROR: No consensus fasta found in $MEDAKADIR" >&2
  exit 4
fi

# Sample name derived from MEDAKADIR basename unless SAMPLE pre-set
SAMPLE=${SAMPLE:-$(basename "$MEDAKADIR")}

# Working directories
WORKDIR="$LOCALSCRATCH/$SAMPLE.Binning.dir"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

echo "📥 Copy reads → node"
mkdir -p reads
for f in "${READS[@]}"; do
  echo " - $f"
  time $RSYNC "$f" reads/
done

echo "🔗 Merge reads → $SAMPLE.fq.gz"
zcat reads/*.f*q.gz > "$SAMPLE.fq"  # supports fq.gz or fastq.gz
#pigz -p "$SLURM_CPUS_ON_NODE" "$SAMPLE.fq"
rm -rf reads

echo "🧬 Copy assembly → $SAMPLE.assembly.fasta"
$RSYNC "$ASSEMBLY_SRC" "$SAMPLE.assembly.fasta"

#####################
# Alignment (minimap2) + BAM
#####################
echo "🚀 minimap2"
date +%d\ %b\ %T

time minimap2 \
  -ax map-ont \
  -t "$SLURM_CPUS_ON_NODE" \
  "$SAMPLE.assembly.fasta" \
  "$SAMPLE.fq" > "$SAMPLE.assembly.sam"

# BAM conversion

echo "🧰 samtools view → BAM"
time samtools view -@ "$SLURM_CPUS_ON_NODE" -bS "$SAMPLE.assembly.sam" > "$SAMPLE.assembly.bam"
rm -f "$SAMPLE.assembly.sam" "$SAMPLE.fq.gz"

echo "🧰 samtools sort"
time samtools sort -@ "$SLURM_CPUS_ON_NODE" "$SAMPLE.assembly.bam" > "$SAMPLE.assembly.sorted.bam"
rm -f "$SAMPLE.assembly.bam"

#####################
# Depth table
#####################
echo "📊 jgi_summarize_bam_contig_depths"
time jgi_summarize_bam_contig_depths \
  --outputDepth "$SAMPLE.depth.txt" \
  "$SAMPLE.assembly.sorted.bam"

echo "🔧 Depth for MaxBin2"
# MaxBin wants contig and average depth columns, no header
cut -f1,3 "$SAMPLE.depth.txt" | tail -n+2 > "$SAMPLE.depth_maxbin.txt"

#####################
# Binning: MetaBAT2 + MaxBin2 (in parallel)
#####################
cpu_half=$(( SLURM_CPUS_ON_NODE / 2 ))
(( cpu_half < 1 )) && cpu_half=1

echo "🏃 Running MetaBAT2 and MaxBin2 with $cpu_half threads each"

time parallel -j2 ::: \
  "metabat2 -i $SAMPLE.assembly.fasta -a $SAMPLE.depth.txt -m 1500 --seed 100 -t $cpu_half --unbinned -o $SAMPLE.Metabat2" \
  "run_MaxBin.pl -contig $SAMPLE.assembly.fasta -out $SAMPLE.MaxBin.out -abund $SAMPLE.depth_maxbin.txt -thread $cpu_half"

# Normalize extensions for downstream tools (.fa -> .fasta)
for i in *.fa; do
  [[ -e "$i" ]] || continue
  mv "$i" "${i%.fa}.fasta"
done

#####################
# Tidy up + report
#####################
rm -f "$SAMPLE.assembly.fasta"

# Create a simple summary
{
  echo "# Binning summary for $SAMPLE"
  echo "date: $(date)"
  echo "node: $SLURM_NODELIST"
  echo "cpus: $SLURM_CPUS_ON_NODE (each binning tool: $cpu_half)"
  echo "contigs: $(grep -c '^>' $SAMPLE.assembly.sorted.bam 2>/dev/null || echo 'n/a')"
  echo "MetaBAT2 bins: $(ls -1 $SAMPLE.Metabat2.*.fasta 2>/dev/null | wc -l)"
  echo "MaxBin2 bins:  $(ls -1 $SAMPLE.MaxBin.out.*.fasta 2>/dev/null | wc -l)"
} > "$SAMPLE.binning.summary.txt"

echo "✅ Done. Results in $WORKDIR and synced to $OUTDIR/"

#####################
# Final cleanup
#####################
cd "$TMPDIR/$USER"
rm -rf "tmpDir_of.${SLURM_JOB_ID}"
echo "✅ I've done"; date