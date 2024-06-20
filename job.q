#!/bin/bash
#SBATCH -N 1
#SBATCH -p shared
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-3
#SBATCH -t 48:00:00
#SBATCH -A TG-CTS120055
#SBATCH --job-name="unnamed"
#SBATCH --output=p.%A_%a.out

cd ${SLURM_SUBMIT_DIR}
NPROCS=`echo ${SLURM_JOB_NODELIST} | wc -l`

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

./program

