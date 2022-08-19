#!/bin/bash
#$ -cwd
#$ -j y
#$ -t 1-118
#$ -l mem_free=23G
#$ -l h_vmem=24G
module load conda_R
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript mutations_mutsig.R $SGE_TASK_ID
