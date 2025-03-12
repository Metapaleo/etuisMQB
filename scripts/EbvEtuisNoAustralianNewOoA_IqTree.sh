#!/bin/bash

#SBATCH --job-name IQTreeHHV4Genome
#SBATCH -p common
#SBATCH -o IQTreeHHV4GenomeEtuis.out -e IQTreeHHV4GenomeEtuis.err
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

echo "**** Job starts ****"
date

echo "**** info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Task id: ${SLURM_TASK_PID}"

module load IQ-TREE/2.4.0

cd ${WORKINGDIR}

# Run phylogeny after estimating the best model
iqtree2 -s ${WORKINGDIR}/msaEbvEtuis.fa --seqtype DNA -T 24 -m TEST -alrt 1000 -B 10000

echo "**** Job ends ****"
date
