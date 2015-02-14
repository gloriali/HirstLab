#!/bin/sh
#$ -S /bin/sh
#$ -m e
#$ -M lli@bcgsc.ca
#$ -N junction_FetalBrain
#$ -V

out=/home/lli/FetalBrain/RNAseq/junction

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A03484/A03484.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A07825/A07825.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out S

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A03473/A03473.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A03475/A03475.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A04599/A04599.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A15298/A15298.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out S

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A03474/A03474.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A03476/A03476.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A15295/A15295.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out S

bam=/projects/epigenomics/ep50/internal/jqc.1.7.6/A15299/A15299.bam
/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh $bam $out S

