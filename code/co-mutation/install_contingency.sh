#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=9G
#$ -l h_vmem=10G
unalias R
module load conda_R
rlib="${HOME}/Library/R/3.12-bioc-release-conda"
R_LIBS_USER=$rlib R CMD build ../contingency.table
R_LIBS_USER=$rlib R CMD INSTALL contingency.table_1.1.0.tar.gz
rm contingency.table_1.1.0.tar.gz
