#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --account=rcac-rnaseq
#SBATCH --qos=standby
#SBATCH --partition=cpu
#SBATCH --time=2:00:00
#SBATCH --job-name=star_index
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load star

data_dir="${SCRATCH}/rnaseq-workshop/data"
genome="${data_dir}/GRCm39.primary_assembly.genome.fa"
gtf="${data_dir}/gencode.vM38.primary_assembly.basic.annotation.gtf"
indexdir="${data_dir}/star_index"

STAR \
    --runMode genomeGenerate \
    --runThreadN ${SLURM_CPUS_ON_NODE} \
    --genomeDir ${indexdir} \
    --genomeFastaFiles ${genome} \
    --sjdbGTFfile ${gtf}

