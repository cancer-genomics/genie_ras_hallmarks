#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=1G
#$ -l h_vmem=2G
qsub install_contingency.sh
qsub -hold_jid install_contingency.sh mutations.sh
#qsub -hold_jid install_contingency.sh inactivating.sh
#qsub -hold_jid install_contingency.sh mutations_tmb.sh
#qsub -hold_jid install_contingency.sh inactivating_tmb.sh
#qsub -hold_jid install_contingency.sh mutations_mutsig.sh
qsub -hold_jid install_contingency.sh deletions.sh
qsub -hold_jid install_contingency.sh amplifications.sh
qsub -hold_jid install_contingency.sh any_deletion.sh
qsub -hold_jid install_contingency.sh any_amplification.sh
qsub -hold_jid install_contingency.sh rearrangements.sh
#qsub -hold_jid rearrangements.sh summarize_models.sh
