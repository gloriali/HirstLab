#!/bin/sh
#$ -S /bin/sh
#$ -l mf=4G
#$ -m e
#$ -M lli@bcgsc.ca
#$ -N RunB2W
#$ -V

JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
SAMTOOLS=/gsc/software/linux-x86_64/samtools-0.1.13/samtools

# -bamFile <bam file name (BAM file has to be coordinate sorted!)> -out <output folder> [parameters]

par=$(echo $3 | sed -e 's/:/ /g' -e 's/,/ /g')
name=$(basename $1)
mkdir -p $2
$JAVA -jar -Xmx4G /home/mbilenky/bin/Solexa_Java/BAM2WIG.jar -samtools $SAMTOOLS -bamFile $1 -out $2 $par > $2/$name".log"
