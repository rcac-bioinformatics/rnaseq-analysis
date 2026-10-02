#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=rcac-rnaseq
#SBATCH --qos=standby
#SBATCH --partition=cpu
#SBATCH --time=0:30:00
#SBATCH --job-name=featurecounts
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load subread

GTF="$SCRATCH/rnaseq-workshop/data/gencode.vM38.primary_assembly.basic.annotation.gtf"
BAMS="$SCRATCH/rnaseq-workshop/results/mapping/*.bam"
output_dir="$SCRATCH/rnaseq-workshop/results/counts"

featureCounts \
  -T ${SLURM_CPUS_PER_TASK} \
  -p \
  -s 0 \
  -t exon \
  -g gene_id \
  -a ${GTF} \
  -o ${output_dir}/gene_counts.txt \
  ${BAMS}
