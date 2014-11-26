#!/bin/sh

# Dynamic growth appraoch to join DM CpGs
if [ "$1" == "-h" ] ; then
    echo -e "Dynamic growth appraoch to join DM CpGs
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <file> -n <name> -s <size> -c <cut>
    <dirIn>: input directory
    <dirOut>: output directory
    <file>: input DM CpG bed file
    <name>: output name
    <size>: max distance to join adjacent CpGs
    <cut>: min No. of CpGs in each DMR"
    exit 0
fi

size=500  # max distance between two consecutive CpGs
cut=3     # minimum number of CpGs
while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; shift;;
        -n) name="$2"; shift;;
        -s) size="$2"; shift;;
        -c) cut="$2"; shift;;
    esac
    shift
done

echo -e "Processing "$file", output format chr\tstart\tend\tID\t(hyper/hypo)\tcount\tlength"

less $dirIn/$file | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirOut/DMR.$name.s$size.c$cut
dmr=($(less $dirOut/DMR.$name.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirOut'""/DMR.""'$name'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirOut'""/DMR.""'$name'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
hyper=($(wc -l $dirOut/DMR.$name.s$size.c$cut.hyper))
hypo=($(wc -l $dirOut/DMR.$name.s$size.c$cut.hypo))
length=($(less $dirOut/DMR.$name.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
count=($(less $dirOut/DMR.$name.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
echo -e $name"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirOut/DMR.summary.stats

