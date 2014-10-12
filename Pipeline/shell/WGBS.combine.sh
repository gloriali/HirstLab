#!/bin/sh

# Combine CpG coverage on both strand
while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

echo "Combining strands for "$file" output to "$dirOut
less $dirIn/$file | awk '{print "chr"$1":"$2"\t"$5"\t"$6}' > $dirOut/$file.tmp
awk 'NR==FNR {h[$1]=$2; next} {if($1 in h){print $0"\t"h[$1]}}' /home/lli/hg19/CG.strand $dirOut/$file.tmp > $dirOut/$file.tmp.join
less $dirOut/$file.tmp.join | awk '{c[$4]=c[$4]+$2; t[$4]=t[$4]+$3; count[$4]=count[$4]+1} END{for(i in count){if(count[i]>2){print i"\t"c[i]"\t"t[i]"\t"count[i] >> "'$dirOut'""/ERROR.""'$file'"".error"}}; for(i in c){m=c[i]/(c[i]+t[i]); print i"\t"c[i]"\t"t[i]"\t"m}}' > $dirOut/$file.combine
rm $dirOut/$file.tmp $dirOut/$file.tmp.join

