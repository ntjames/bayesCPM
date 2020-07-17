#!/bin/bash
#SBATCH --job-name=cpm_ppr1_a2 
#SBATCH --mail-user=nathan.t.james@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:30:00
#SBATCH --mem=2G
#SBATCH --constraint=sandybridge|haswell
#SBATCH --partition=production ## change from debug to production for real run
#SBATCH --array=1-25
#SBATCH --output=out/sim_a2_%a.out

module restore R_360
R --version

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## run this first to compile Stan models
## Rscript --no-save 1_compile_mod_paper_a.R

## run array job there are 25 rows
Rscript --no-save 3_sim_a2.R $SLURM_ARRAY_TASK_ID
