#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=4G
#$ -l h_vmem=5G
unalias R
module load conda_R
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript summarize_models.R
