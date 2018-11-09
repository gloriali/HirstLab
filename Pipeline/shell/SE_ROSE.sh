#!/bin/sh

# Run ROSE rank algorithm for superEnhnacers
if [ "$1" == "-h" ] ; then
    echo -e "Compute SuperEnhancres
Usage: `basename $0` -o <dirOut> -r <region> -f <file> -c <control> -g <genome>
    <dirOut>: output directory
    <region>: full path to enrich region bed file, format bedFile
    <file>: full path to input source file, ChIP bamFile
    <control>: full path to input control file, Input bamFile
    <genome>:Genome build: hg19 or mm10"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -o) dirOut="$2"; shift;;
        -r) region="$2"; shift;;
        -f) file="$2"; shift;;
        -c) control="$2"; shift;;
        -g) genome="$2"; shift;;
    esac
    shift
done

samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
mkdir -p $dirOut

# generate gff file
less $region | sed 's/chr//g' | awk '{OFS="\t"; print "chr"$1, "chr"$1":"$2"-"$3,".",$2,$3, ".",".",".", "chr"$1":"$2"-"$3}' > $dirOut/$region.gff

# Generating sorted bamFiles with Chr
for bam in $file $control; do
    name=$(basename $bam | sed 's/.bam//')
    echo "Generating sorted bamFiles with Chr for " $name
    samtools view -H $bam | sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | sed -e 's/SN:MT/SN:chrM/' | samtools reheader - $bam >  $dirOut/$name"_chr.bam"
    $sambamba index $dirOut/$name"_chr.bam" -t 8 
done

# running ROSE algorithm
name1=$(basename $file | sed 's/.bam//'); name2=$(basename $control | sed 's/.bam//');
cd /projects/epigenomics3/epigenomics3_results/users/alorzadeh/young_computation-rose-1a9bb86b5464/
/home/jyzhu/anaconda2/bin/python2.7 /projects/epigenomics3/epigenomics3_results/users/alorzadeh/young_computation-rose-1a9bb86b5464/ROSE_main.py -g $genome -i $dirOut/$region.gff -r $dirOut/$name1"_chr.bam" -o $dirOut -t 2500 -c $dirOut/$name2"_chr.bam" &> $dirOut/ROSE.log

rm $dirOut/$region.gff
rm $dirOut/*chr.bam
rm $dirOut/*bam.bai
