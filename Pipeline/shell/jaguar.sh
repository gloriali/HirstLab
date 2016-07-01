#!/bin/sh

# Calculate RPKM from fastq files using Jaguar alignment
if [ "$1" == "-h" ] ; then
    echo -e "Calculate RPKM from fastq files using Jaguar alignment 
Usage: `basename $0` -i <dirIn> -o <dirOut> -f1 <fastq1> -f2 <fastq2> -n <name> -r <ref> -v <ens>
    <dirIn>: input directory
    <dirOut>: output directory
    <fastq1>: input fastq1 file     
    <fastq2>: input fastq2 file
    <name>: name of output files
    <ref>: Ensembl reference genome
    <ens>: Ensembl version, e.g. hg19_ens69"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f1) f1="$2"; shift;;
        -f2) f2="$2"; shift;;
        -n) name="$2"; shift;;
        -r) ref="$2"; shift;;
        -v) ens="$2"; shift;; 
    esac
    shift
done

####################
# jaguar alignments
####################

mkdir -p $dirOut/run;
AUTOMATION_HOME=/projects/wtsspipeline/programs/external_programs/Automation/
/projects/wtsspipeline/programs/external_programs/Automation/bin/runBWA.sh -1 $dirIn/$f1 -2 $dirIn/$f2 -n $name -o $dirOut -r $ref &> $dirOut/run/$name.bwa.log

####################
# sorting by name
####################

/gsc/software/linux-x86_64-centos5/samtools-0.1.13/samtools sort -n $dirOut/$name.sorted.bam $dirOut/$name".sortedByName";

#######################
# repositioning step using:
#######################

mkdir -p $dirOut/j/;
rm -rf $out/j/
/home/mbilenky/bin/Solexa_Shell/RunJR.sh $dirOut/$name".sortedByName.bam" $dirOut/j $ens &> $dirOut/run/$name.j.log
