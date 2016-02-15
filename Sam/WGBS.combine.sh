#!/bin/sh

# Combine CpG coverage on both strand
if [ "$1" == "-h" ] ; then
    echo -e "Combine CpG coverage on both strand
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <fractional> -c <coverage> -s <CG.strand>
    <dirIn>: input directory
    <dirOut>: output directory
    <fractional>: WGBS .sam.bedGraph.gz input file, format chr\tstart\tend\tfractional_methylation
    <coverage>: WGBS .sam.bedGraph.gz input file, format chr\tstart\tend\tcoverage
    <CG.strand>: mapping of stranded CpG IDs and combined CpG IDs, default genome hg19. "
    exit 0
fi

#Initializing variables
strand=/home/lli/hg19/CG.strand 
while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) fractional="$2"; shift;;
        -c) coverage="$2"; shift;;
        -s) strand="$2"; shift;;
    esac
    shift
done



# Combine fractional and coverage and reformat to:  chr:start  fractional_methylation  coverage
echo "Combining strands for "$file" output to "$dirOut
less $dirIn/$fractional | awk '{if($1~/chr/){print $1":"$2"\t"$4} else{print "chr"$1":"$2"\t"$4}}' > $dirOut/$fractional.tmp
less $dirIn/$coverage | awk '{if($1~/chr/){print $1":"$2"\t"$4} else{print "chr"$1":"$2"\t"$4}}' > $dirOut/$coverage.tmp
awk 'NR==FNR {b[$1]=$2; next} NR!=FNR {if($1 in b){print $0"\t"b[$1]}}' $dirOut/$coverage.tmp $dirOut/$fractional.tmp > $dirOut/$fractional.coverage.tmp
rm $dirOut/$fractional.tmp $dirOut/$coverage.tmp

# Find matches to strand and append the appropriate combined CpG ID
awk 'NR==FNR {h[$1]=$2; next} NR!=FNR {if($1 in h){print $0"\t"h[$1]}}' $strand $dirOut/$fractional.coverage.tmp > $dirOut/$fractional.coverage.tmp.join

# Send data related to CpG IDs that appear more than twice to an error file
# Output format:   Combined_CpG_ID  Fractional_Methylation  Coverage    m
less $dirOut/$fractional.coverage.tmp.join | awk '{c[$4]=c[$4]+$2; t[$4]=t[$4]+$3; count[$4]=count[$4]+1} END{for(i in count){if(count[i]>2){print i"\t"c[i]"\t"t[i]"\t"count[i] >> "'$dirOut'""/ERROR.""'$fractional'"".error"}}
for(i in c){m=c[i]/(c[i]+t[i]); print i"\t"c[i]"\t"t[i]"\t"m}}' > $dirOut/$fractional.combine
rm $dirOut/$fractional.coverage.tmp $dirOut/$fractional.coverage.tmp.join
