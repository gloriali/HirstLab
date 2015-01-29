#!/bin/sh
#$ -S /bin/sh
#$ -l mf=4G
#$ -m e
#$ -M lli@bcgsc.ca
#$ -N RunW2BW
#$ -V


/home/mbilenky/UCSCtools/wigToBigWig $1.wig.gz /home/mbilenky/UCSC_chr/$4.chrom.sizes $3/$2.bw
