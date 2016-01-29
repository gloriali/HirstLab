#!/bin/sh

# Combine CpG coverage on both strand
if [ "$1" == "-h" ] ; then
    echo -e "Combine CpG coverage on both strand
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <file> -s <CG.strand>
    <dirIn>: input directory
    <dirOut>: output directory
    <file>: WGBS .sam.bedGraph input file, format chr\tstart\tend\tfractional_methylation\tcoverage      
    <CG.strand>: mapping of stranded CpG IDs and combined CpG IDs, default genome hg19. "
    exit 0
fi

strand=/home/lli/hg19/CG.strand 
while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; shift;;
        -s) strand="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

echo "Combining strands for "$file" output to "$dirOut
less $dirIn/$file | awk '{if($1~/chr/){print $1":"$2"\t"$4"\t"$5} else{print "chr"$1":"$2"\t"$4"\t"$5}}' > $dirOut/$file.tmp
awk 'NR==FNR {h[$1]=$2; next} {if($1 in h){print $0"\t"h[$1]}}' $strand $dirOut/$file.tmp > $dirOut/$file.tmp.join
less $dirOut/$file.tmp.join | awk '{c[$4]=c[$4]+$2; t[$4]=t[$4]+$3; count[$4]=count[$4]+1} END{for(i in count){if(count[i]>2){print i"\t"c[i]"\t"t[i]"\t"count[i] >> "'$dirOut'""/ERROR.""'$file'"".error"}}; for(i in c){m=c[i]/(c[i]+t[i]); print i"\t"c[i]"\t"t[i]"\t"m}}' > $dirOut/$file.combine
rm $dirOut/$file.tmp $dirOut/$file.tmp.join

