#!/bin/sh

# ChIP-seq check if regions are enriched from signal files
if [ "$1" == "-h" ] ; then
    echo -e "ChIP-seq unique enriched regions in pairwise comparisons 
Usage: `basename $0` -r <region> -w <wig> -o <dirOut> -n <name> -excl <excl> -q <quant>
    <region>: enriched regions in sample1, format: chr\tstart\tend\t<additional columns>
    <wig>: wig signalling file for sample2
    <dirOut>: output directory
    <name>: output name
    <excl>: enriched regions in sample2, format: chr\tstart\tend\t<additional columns>
    <quant>: quantile of background signal to use as cutoff, default to 0.9
Output: $dirOut/$name.unique, format: chr\tstart\tend\tID\taverage coverage\tmax coverage"
    exit 0
fi

quant=0.9
while [ $# -gt 0 ]
do
    case "$1" in
        -r) region="$2"; shift;;
        -w) wig="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -n) name="$2"; shift;;
        -excl) excl="$2"; shift;;
        -q) quant="$2"; shift;;
    esac
    shift
done

echo "Test whether $region are enriched in $wig, output to $dirOut/$name"
JAVA=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
RegCov=/home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/R
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes

## calculate region signal from wig
less $region | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/$name.signal.bed
$JAVA -jar -Xmx15G $RegCov -w $wig -r $dirOut/$name.signal.bed -o $dirOut -c $chr -n $name.signal > $dirOut/$name.log

## generate background regions
$BEDTOOLS/shuffleBed -i $region -g $chr -excl $excl | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/$name.background.bed

## calculate coverage for background regions
$JAVA -jar -Xmx15G $RegCov -w $wig -r $dirOut/$name.background.bed -o $dirOut -c $chr -n $name.background >> $dirOut/$name.log

## plot and cutoff: 90% quantile of coverage for background regions
less $excl | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/$name.ER.bed
$JAVA -jar -Xmx15G $RegCov -w $wig -r $dirOut/$name.ER.bed -o $dirOut -c $chr -n $name.ER > $dirOut/$name.log
echo "signal <- read.delim(\"$(ls $dirOut/*.$name.signal.coverage)\", as.is = T, head = F, col.names = c(\"chr\", \"start\", \"end\", \"ID\", \"cov\", \"max\"))
background <- read.delim(\"$(ls $dirOut/*.$name.background.coverage)\", as.is = T, head = F, col.names = c(\"chr\", \"start\", \"end\", \"ID\", \"cov\", \"max\"))
ER <- read.delim(\"$(ls $dirOut/*.$name.ER.coverage)\", as.is = T, head = F, col.names = c(\"chr\", \"start\", \"end\", \"ID\", \"cov\", \"max\"))
pdf(\"$dirOut/$name.pdf\")
plot(c(0, 20), c(0, 1), type = \"n\", main = \"$name\", xlab = \"average coverage\", ylab = \"ecdf\")
lines(ecdf(ER\$cov), col = \"black\")
lines(ecdf(signal\$cov), col = \"red\")
lines(ecdf(background\$cov), col = \"blue\")
abline(v = quantile(background\$cov, 0.9), col = \"blue\")
legend(\"bottomright\", c(\"Sample1 ER\", \"Sample2 ER\", \"Sample2 background\"), col = c(\"red\", \"black\", \"blue\"), lwd = 5, lty = 1, cex = 0.8)
dev.off()
" > $dirOut/$name.R
$R CMD BATCH $dirOut/$name.R
cutoff=$(less $dirOut/*.$name.background.coverage | awk '{print $5}' | sort -k1,1n | awk '{c[NR]=$1} END{print c[int(NR*'"$quant"')]}')
echo "Coverage cutoff: $cutoff"

## output
less $dirOut/*.$name.signal.coverage | awk '{if($5 <= '"$cutoff"'){print $0}}' > $dirOut/$name.unique
rm $dirOut/*$name.ER* $dirOut/*$name.signal* $dirOut/*$name.background* $dirOut/$name.R

