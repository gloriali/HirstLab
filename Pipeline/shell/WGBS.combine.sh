#!/bin/sh

# Combine CpG coverage on both strand
if [ "$1" == "-h" ] ; then
    echo -e "Combine CpG coverage on both strand
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <file> -s <CG.strand> -c <cov>
    <dirIn>: input directory
    <dirOut>: output directory
    <file>: input file, format chr\tposition\tstrand\tconverted\tunconverted (novoalign output .5mC.CpG file)      
    <CG.strand>: mapping of stranded CpG IDs and combined CpG IDs, default genome hg19
    <cov>: minimal coverage required, default to 3. "
    exit 0
fi

strand=/home/lli/hg19/CG.strand 
cov=3
while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; shift;;
        -s) strand="$2"; shift;;
        -c) cov="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

echo "Combining strands for "$file" output to "$dirOut
less $dirIn/$file | awk '{gsub(/chr/, ""); print $1":"$2"\t"$4"\t"$5}' > $dirOut/$file.tmp
awk 'NR==FNR {h[$1]=$2; next} {if($1 in h){print $0"\t"h[$1]}}' $strand $dirOut/$file.tmp > $dirOut/$file.tmp.join
less $dirOut/$file.tmp.join | awk '{c[$4]=c[$4]+$3; t[$4]=t[$4]+$2; count[$4]=count[$4]+1} END{for(i in count){if(count[i]>2){print i"\t"c[i]"\t"t[i]"\t"count[i] >> "'$dirOut'""/ERROR.""'$file'"".error"}}; for(i in c){chr=gensub(":.+", "", "g", i); end=gensub(".+-", "", "g", i); if(t[i]+c[i] >= "'$cov'"+0){print chr"\t"end-2"\t"end"\t"t[i]"\t"c[i]"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 -k2,2n -T /projects/epigenomics/temp/ > $dirOut/$file.combine.5mC.CpG
rm $dirOut/$file.tmp $dirOut/$file.tmp.join

