#!/bin/bash
#$ -cwd
#$ -j y
#$ -t 4
#$ -l mem_free=24G
#$ -l h_vmem=26G
# #$ -pe local 16
#$ -l cancergen
module load conda_R/4.1
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript amplifications.R $SGE_TASK_ID
