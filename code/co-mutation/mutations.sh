#!/bin/bash
#$ -cwd
#$ -j y
#$ -t 1-215
#$ -l mem_free=29G
#$ -l h_vmem=30G
unalias R
module load conda_R
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript mutations.R $SGE_TASK_ID
