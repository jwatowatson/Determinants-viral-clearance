#!/bin/bash

#SBATCH -A moru-batty.prj
#SBATCH -D /well/moru-batty/users/lcd199/PLATCOV-SAP/Determinants-viral-clearance
#SBATCH -J res
#SBATCH -n 4
#SBATCH -o /well/moru-batty/users/lcd199/PLATCOV-SAP/Determinants-viral-clearance/o_and_e_files/output.o%A_%a.out
#SBATCH -e /well/moru-batty/users/lcd199/PLATCOV-SAP/Determinants-viral-clearance/o_and_e_files/output.e%A_%a.out
#SBATCH -p short
#SBATCH --array 1-1800


echo started=`date`
module purge
module load R/4.1.2-foss-2021b



echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`
Rscript /well/moru-batty/users/lcd199/PLATCOV-SAP/Determinants-viral-clearance/03_run_cluster_bootstraps_analysis_sampling_fq.R ${SLURM_ARRAY_TASK_ID} --no-save --no-restore
echo "finished="`date`
exit 0