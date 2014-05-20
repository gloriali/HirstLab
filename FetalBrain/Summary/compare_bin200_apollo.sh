#!/bin/sh
#$ -S /bin/sh
#$ -m e
#$ -M lli@bcgsc.ca
#$ -N compare3methods_bin200
#$ -V
/gsc/software/linux-x86_64/R-2.14.2/bin/R CMD BATCH /home/lli/FetalBrain/MeDIP/compare_bin200.R