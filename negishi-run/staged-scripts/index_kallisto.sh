#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --account=rcac-rnaseq
#SBATCH --qos=standby
#SBATCH --partition=cpu
#SBATCH --time=2:00:00
#SBATCH --job-name=kallisto_index
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load kallisto

cd $SCRATCH/rnaseq-workshop
mkdir -p data/kallisto_index

kallisto index \
    -i data/kallisto_index/transcripts.idx \
    data/gencode.vM38.transcripts-clean.fa
