#!/bin/bash
#$ -cwd
#$ -j y
#$ -t 1-24
#$ -l mem_free=23G
#$ -l h_vmem=24G
unalias R
module load conda_R
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript rearrangements.R $SGE_TASK_ID
