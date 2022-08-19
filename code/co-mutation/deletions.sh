#!/bin/bash
#$ -cwd
#$ -j y
#$ -t 1-8
#$ -l mem_free=24G
#$ -l h_vmem=25G
module load conda_R
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript deletions.R $SGE_TASK_ID
