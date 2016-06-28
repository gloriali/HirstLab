#!/bin/sh

samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
java=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java
resource=/projects/epigenomics/resources/ERCC/mix1.class.concentration.length
out=/projects/edcc_prj2/RNAseq/PX0348/SpikeIn/
cd $out
for file in *.sorted.bam; do
    name=$(echo $file | sed -e 's/.sorted.bam//g')
    echo $name
    $java -jar -Xmx8G /home/pubseq/BioSw/picard/picard-tools-1.52/MarkDuplicates.jar I=$out/$name.sorted.bam O=$out/$name.sorted.dups_marked.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $out/$name.sorted.dups_marked.bam > $out/$name.flagstat
    rm -rf $out/*.sa* $out/$name.bam $out/$name.sorted.bam $out/$name.f.bam
done

for index in AACCCC ACCCAG AGCGCT CAAAAG; do
    echo $index
    $samtools merge $out/$index.sorted.dups_marked.bam $(ls $out/s*$index.sorted.dups_marked.bam)
    $samtools flagstat $out/$index.sorted.dups_marked.bam > $out/$index.flagstat
    
    $samtools view -F4 $out/$index.sorted.dups_marked.bam | cut -f3 | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1 > $out/x
    tot=$(less $resource | wc -l)
    num=$(less $out/x | wc -l)
    echo "Spike-in species recovered (at least one read): $num out of $tot"
    cat $out/x $resource | cut -f1 | sort | uniq -c | awk '{if($1==1) print $2"\t"0.1}' > $out/y # use 0.1 for unrecovered species
    cat $out/x $out/y | sort -k1,1 > $out/z
    join $out/z $resource | sed 's/ /\t/g' > $out/$index.count.class.concentration.length
    rm $out/x $out/y $out/z    
done

for file in *.count.class.concentration.length; do
    index=$(echo $file | sed -e 's/.count.class.concentration.length//g')
    depth=$(less $index.flagstat | awk 'NR==1 {print $1}')
    echo $index $depth
    less $file | awk '{print $0"\t"$2/$5/"'$depth'"*10^9"\t""'$index'"}' > $file.RPKM
done
cat *.RPKM > SpikeIn.concentration.RPKM

