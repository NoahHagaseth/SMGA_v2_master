#!/bin/bash
#SBATCH --job-name=medaka_polish
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=36:00:00
#SBATCH --out slurm-polish-%A_%a.out
################

#Variables
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'
input=$1          # e.g. saltwater
READFILE=$2       # full path to reads .fq.gz
ASSDIR=$3         # flye output dir (contains assembly.fasta)
OUTDIR=$4         # output directory

##Activate conda environments

module --quiet purge  # Reset the modules to the system default
echo "SLURM partition assigned: $SLURM_JOB_PARTITION"

# Auto-detect the partition and load the correct module environment
if [[ "$SLURM_JOB_PARTITION" == "a100" ]]; then
    echo "Running on A100 partition - Swapping to Zen2Env"
    module --force swap StdEnv Zen2Env
else
    echo "Running on Accel (P100) partition - Using default Intel environment"
fi


##Activate conda environments
module load Anaconda3/2021.11
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null

conda activate /path/to/medaka_env
echo "I am working with this" $CONDA_PREFIX

###Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date

## Copying data to local node for faster computation
SCRATCH=/mnt/SCRATCH/$USER/medaka_${input}
mkdir -p "$SCRATCH"
cd "$SCRATCH"


echo "copying files to" "$SCRATCH"

echo "Copy fq file"

cp "$READFILE" "$input.fq.gz"

echo "Copy assembly"

cp "$ASSDIR/assembly.fasta" "$input.assembly.fasta"

##MEdaking

echo "Starting Medaka..."
date +%d\ %b\ %T

time medaka_consensus \
-i $input.fq.gz \
-d $input.assembly.fasta \
-o "${input}.medaka.dir" \
-t $SLURM_CPUS_ON_NODE

echo "Cleaning and changing names..."

if [[ ! -f "${input}.medaka.dir/consensus.fasta" ]]; then
    echo "ERROR: consensus.fasta not produced!"
    exit 1
fi


###Cleaning

cd "${input}.medaka.dir"

mv consensus.fasta "${input}.medaka.consensus.fasta"

# Remove everything else
find . -mindepth 1 ! -name "${input}.medaka.consensus.fasta" -exec rm -rf {} +

# =========================
# Copy results back
# =========================
mkdir -p "$OUTDIR"
cd "$SCRATCH"

echo "Copying results to $OUTDIR"
$RSYNC "${input}.medaka.dir" "$OUTDIR/"

echo "Medaka polishing completed successfully"
date