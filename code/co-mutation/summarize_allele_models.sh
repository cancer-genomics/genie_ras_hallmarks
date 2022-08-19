#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=7G
#$ -l h_vmem=8G
unalias R
module load conda_R
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript summarize_allele_models.R
