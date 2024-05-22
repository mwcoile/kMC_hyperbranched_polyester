#!/bin/bash
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-1
#SBATCH -t 48:00:00
#SBATCH -A p31595
#SBATCH --job-name="500 only"
#SBATCH --output=p.%A_%a.out

cd ${SLURM_SUBMIT_DIR}
NPROCS=`echo ${SLURM_JOB_NODELIST} | wc -l`

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

./program

