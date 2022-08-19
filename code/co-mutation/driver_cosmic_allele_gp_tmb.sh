#!/bin/bash
#$ -cwd
#$ -j y
#$ -t 1-7
#$ -l mem_free=25G
#$ -l h_vmem=26G
#$ -o /dcs04/scharpf/data/abalan/kras_cd/trial
# #$ -l cancergen
unalias R
module load conda_R/4.1
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib Rscript driver_cosmic_allele_gp_tmb.R  $SGE_TASK_ID
