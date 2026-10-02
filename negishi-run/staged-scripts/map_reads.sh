#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --account=rcac-rnaseq
#SBATCH --qos=standby
#SBATCH --partition=cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=read_mapping
#SBATCH --array=1-8
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load star

FASTQ_DIR="$SCRATCH/rnaseq-workshop/data"
GENOME_INDEX="$SCRATCH/rnaseq-workshop/data/star_index"
OUTDIR="$SCRATCH/rnaseq-workshop/results/mapping"

mkdir -p $OUTDIR

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

R1=${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz
R2=${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz

STAR \
    --runThreadN ${SLURM_CPUS_ON_NODE} \
    --genomeDir $GENOME_INDEX \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --outFileNamePrefix $OUTDIR/$SAMPLE \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate
