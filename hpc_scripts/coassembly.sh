#!/bin/bash
#SBATCH --job-name=MetaFlye
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --out slurm-coassembly-%A_%a.out

############################

############################
# CHOOSE SAMPLE HERE
############################
# Options: freshwater | mucus | saltwater

SAMPLE="freshwater"

############################
# Paths
############################

PROJECTS=/mnt/project
BASE=$PROJECTS/path/to/folder/with/inputfiles
OUTDIR=$BASE/assemblies

############################
# Input files
############################

case "$SAMPLE" in
  freshwater)
    INPUT="$BASE/Freshwater/freshwater_combined2.cleaned.fq.gz"
    PREFIX="freshwater"
    ;;
  mucus)
    INPUT="$BASE/Mucus/mucus_combined.cleaned.fq.gz"
    PREFIX="mucus"
    ;;
  saltwater)
    INPUT="$BASE/Saltwater/saltwater_combined.cleaned.fq.gz"
    PREFIX="saltwater"
    ;;
  *)
    echo "ERROR: Unknown SAMPLE: $SAMPLE"
    exit 1
    ;;
esac

############################
# Sanity check
############################

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

############################
# Activate Flye environment
############################

module --quiet purge
module load Anaconda3/2021.11

source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null
conda activate /path/to/FLYE_2.9.3.0

echo "Using Flye from:"
which flye
flye --version

############################
# Scratch setup
############################

SCRATCH=/mnt/SCRATCH/$USER
WORKDIR=$SCRATCH/${PREFIX}_flye

mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

echo "Working directory: $WORKDIR"
echo "Sample: $SAMPLE"
echo "Input file: $INPUT"
echo "CPUs: $SLURM_CPUS_ON_NODE"
date

############################
# Copy input to scratch
############################

cp "$INPUT" "${PREFIX}.fq.gz"

############################
# Run MetaFlye
############################

echo "Starting MetaFlye assembly..."
date

time flye \
  --nano-raw "${PREFIX}.fq.gz" \
  --meta \
  --out-dir "${PREFIX}.flye.outdir" \
  -t $SLURM_CPUS_ON_NODE

############################
# Copy results back
############################

mkdir -p "$OUTDIR"
rsync -a "${PREFIX}.flye.outdir" "$OUTDIR/"

echo "Assembly finished for $SAMPLE"
date

############################
# Clean scratch
############################

rm -rf "$WORKDIR"