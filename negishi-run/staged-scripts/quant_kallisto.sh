#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=rcac-rnaseq
#SBATCH --qos=standby
#SBATCH --partition=cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=kallisto_quant
#SBATCH --array=1-8
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load kallisto

DATA="$SCRATCH/rnaseq-workshop/data"
INDEX="$SCRATCH/rnaseq-workshop/data/kallisto_index/transcripts.idx"
OUT="$SCRATCH/rnaseq-workshop/results/kallisto_quant"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

mkdir -p ${OUT}/${SAMPLE}

R1=${DATA}/${SAMPLE}_R1.fastq.gz
R2=${DATA}/${SAMPLE}_R2.fastq.gz

kallisto quant \
    -i $INDEX \
    -o $OUT/$SAMPLE \
    -b 0 \
    -t ${SLURM_CPUS_ON_NODE} \
    $R1 $R2 &> $OUT/$SAMPLE/${SAMPLE}.log
