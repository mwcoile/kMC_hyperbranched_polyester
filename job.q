#!/bin/bash
#SBATCH -N 1
#SBATCH -p buyin
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-20
#SBATCH -t 480:00:00
#SBATCH -A b1039
#SBATCH --job-name="memAccess"
#SBATCH --output=p.%A_%a.out

cd ${SLURM_SUBMIT_DIR}
NPROCS=`echo ${SLURM_JOB_NODELIST} | wc -l`

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

./program
