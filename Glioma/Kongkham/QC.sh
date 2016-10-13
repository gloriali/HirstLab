#!/bin/sh

# QC alignment
## QC original alignment - unsorted bam
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirOut
for bam in *_trimmed.bam; do
    name=$(echo $bam | sed -e 's/.bam//g')
    echo $name
    $samtools sort $dirOut/$bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam  > $dirOut/$name.bamstats
done
## re-align a couple of libraries with Fetal Brain pipeline and our current pipeline
### current pipeline
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
genome=/home/lli/hg19/GRCh37-lite.fa
bwa=/home/pubseq/BioSw/bwa/bwa-0.7.5a/bwa
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/fq/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirIn
for file in *.fq; do
    name=$(echo $file | sed -e 's/.fastq.trimmed.fq/_current/g')
    echo $name
    $bwa mem -M -t 12 $genome $dirIn/$file  > $dirOut/$name.sam
    $samtools view -Sb $dirOut/$name.sam > $dirOut/$name.bam
    $samtools sort $dirOut/$name.bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.bamstats
done 
### Fetal Brain pipeline
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
genome=/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa
bwa=/home/pubseq/BioSw/bwa/bwa-0.5.7/bwa-0.5.7/bwa
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/fq/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirIn
for file in *.fq; do
    name=$(echo $file | sed -e 's/.fastq.trimmed.fq/_FB/g')
    echo $name
    $bwa aln $genome $file > $dirOut/$name.sai
    $bwa samse $genome $dirOut/$name.sai $file > $dirOut/$name.sam
    $samtools view -Sb $dirOut/$name.sam > $dirOut/$name.bam
    $samtools sort $dirOut/$name.bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat   
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.bamstats
done
## QC reports
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirOut
for file in *.bamstats; do
    lib=$(echo $file | sed -e 's/.bamstats//g')
    echo $lib
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $lib $file
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut

# sort mark dups and QC 
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirOut
for bam in *_trimmed.bam; do
    name=$(echo $bam | sed -e 's/.bam//g')
    echo $name
    $samtools sort $dirOut/$bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam  > $dirOut/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirOut/$name.bamstats
	if [ -s $name.sorted.dupsFlagged.bam ]; then
		rm $dirOut/$bam $dirOut/$name.sorted.bam
	fi
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
dirIn=/projects/epigenomics2/users/mmingay/all_bam/hmedip/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/Vitc_hMC/
mkdir -p $dirOut
cd $dirOut
for bam in $dirIn/*.bam; do
    name=$(basename $bam | sed -e 's/.bam//g' | sed 's/mutant/Vitc/g')
    echo $name
    $samtools sort $bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam  > $dirOut/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirOut/$name.bamstats
    rm $dirOut/$name.sorted.bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut

# mapping to repetitive regions
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/
echo -e "Library\tNumber_Repetitive_Mapping" > $dirOut/QC_repetitive
cd $dirIn
for bam in *.bam /projects/epigenomics2/users/lli/glioma/Kongkham/Vitc_hMC/*.sorted.dupsFlagged.bam /projects/epigenomics/users/lli/FetalBrain/MeDIP/bam/*.bam; do
    lib=$(basename $bam | sed 's/.bam//g'| sed 's/.sorted.dupsFlagged//g')
    echo $lib
    rep=$($samtools view $bam | awk '{if($5==0){n++}}END{print n}')
    echo -e $lib"\t"$rep >> $dirOut/QC_repetitive
done

# CpG coverage
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
dirNPC=/projects/epigenomics/users/lli/FetalBrain/MeDIP/bam/
CG=/home/lli/hg19/CG.BED # 28217448
cd $dirIn
## MC
echo "Coverage" > $dirIn/CG.coverage.MC
for i in {0..50}; do
    echo $i >> $dirIn/CG.coverage.MC
done
echo ">50" >> $dirIn/CG.coverage.MC
for bam in *_MC_*_trimmed.sorted.dupsFlagged.bam; do
    sample=$(echo $bam | sed 's/_MC_.*_trimmed.sorted.dupsFlagged.bam//g')
    echo "Processing" $sample
    echo $sample > x
    less $CG | awk '{print "chr"$0}' | $BEDTOOLS/coverageBed -a stdin -b <($samtools view -q 5 -F 1028 -b $bam) -counts | awk 'BEGIN{for(i=0;i<=51;i++){s[i]=0}} {if($5>50){s[51]++} else {s[$5]++}} END{for(i=0;i<=51;i++){print s[i]}}' >>x
    paste $dirIn/CG.coverage.MC x > y
    mv y $dirIn/CG.coverage.MC
done
for bam in $dirNPC/*.bam; do
    sample=$(basename $bam | sed 's/.bam//g')
    echo "Processing" $sample
    echo $sample > x
    $BEDTOOLS/coverageBed -a $CG -b <($samtools view -q 5 -F 1028 -b $bam) -counts | awk 'BEGIN{for(i=0;i<=51;i++){s[i]=0}} {if($5>50){s[51]++} else {s[$5]++}} END{for(i=0;i<=51;i++){print s[i]}}' >>x
    paste $dirIn/CG.coverage.MC x > y
    mv y $dirIn/CG.coverage.MC
done
rm x
## hMC
echo "Coverage" > $dirIn/CG.coverage.hMC
for i in {0..50}; do
    echo $i >> $dirIn/CG.coverage.hMC
done
echo ">50" >> $dirIn/CG.coverage.hMC
for bam in *_hMC_*_trimmed.sorted.dupsFlagged.bam; do
    sample=$(echo $bam | sed 's/_hMC_.*_trimmed.sorted.dupsFlagged.bam//g')
    echo "Processing" $sample
    echo $sample > a
    less $CG | awk '{print "chr"$0}' | $BEDTOOLS/coverageBed -a stdin -b <($samtools view -q 5 -F 1028 -b $bam) -counts | awk 'BEGIN{for(i=0;i<=51;i++){s[i]=0}} {if($5>50){s[51]++} else {s[$5]++}} END{for(i=0;i<=51;i++){print s[i]}}' >>a
    paste $dirIn/CG.coverage.hMC a > b
    mv b $dirIn/CG.coverage.hMC
done
for bam in /projects/epigenomics2/users/lli/glioma/Kongkham/Vitc_hMC/*.sorted.dupsFlagged.bam; do
    sample=$(basename $bam | sed 's/.sorted.dupsFlagged.bam//g')
    echo "Processing" $sample
    echo $sample > a
    $BEDTOOLS/coverageBed -a $CG -b <($samtools view -q 5 -F 1028 -b $bam) -counts | awk 'BEGIN{for(i=0;i<=51;i++){s[i]=0}} {if($5>50){s[51]++} else {s[$5]++}} END{for(i=0;i<=51;i++){print s[i]}}' >>a
    paste $dirIn/CG.coverage.hMC a > b
    mv b $dirIn/CG.coverage.hMC
done
rm a

# Merge bam files (IDH mut vs IDH wt)
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
CG=/home/lli/hg19/CG.BED # 28217448
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirIn
for file in *.bam; do
	sample=$(echo $file | sed 's/_.*//g');
	name=$(echo $file | sed "s/$sample//g");
	IDH=$(less ../Sample_info.txt | awk -F"\t" '{if($1=="'$sample'"){if($2=="positive"){print "_IDHmut"}else{print "_IDHwt"}}}');
	new=$sample$IDH$name;
	mv $file $new
done
$samtools merge IDHmut_MC_merge.trimmed.sorted.dupsFlagged.bam $(ls *IDHmut_MC*trimmed.sorted.dupsFlagged.bam)
$samtools merge IDHwt_MC_merge.trimmed.sorted.dupsFlagged.bam $(ls *IDHwt_MC*trimmed.sorted.dupsFlagged.bam)
$samtools merge IDHmut_hMC_merge.trimmed.sorted.dupsFlagged.bam $(ls *IDHmut_hMC*trimmed.sorted.dupsFlagged.bam)
$samtools merge IDHwt_hMC_merge.trimmed.sorted.dupsFlagged.bam $(ls *IDHwt_hMC*trimmed.sorted.dupsFlagged.bam)
## QC
for bam in *merge*.bam; do
    name=$(echo $bam | sed -e 's/.bam//g')
    echo $name
    $samtools flagstat $dirIn/$name.bam > $dirIn/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirIn/$name.bam  > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
echo "Coverage" > $dirIn/CG.coverage.merge
for i in {0..50}; do
    echo $i >> $dirIn/CG.coverage.merge
done
echo ">50" >> $dirIn/CG.coverage.merge
for bam in *merge*.bam; do
    sample=$(echo $bam | sed 's/.trimmed.sorted.dupsFlagged.bam//g')
    echo "Processing" $sample
    echo $sample > x
    less $CG | awk '{print "chr"$0}' | $BEDTOOLS/coverageBed -a stdin -b <($samtools view -q 5 -F 1028 -b $bam) -counts | awk 'BEGIN{for(i=0;i<=51;i++){s[i]=0}} {if($5>50){s[51]++} else {s[$5]++}} END{for(i=0;i<=51;i++){print s[i]}}' >>x
    paste $dirIn/CG.coverage.merge x > y
    mv y $dirIn/CG.coverage.merge
done
## wig
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
chrsize=/home/lli/hg19/hg19.chrom.sizes
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/wig/
mkdir -p $dirOut
fl=175
for bam in $dirIn/*merge*.bam; do
	name=$(basename $bam | sed -e 's/.trimmed.sorted.dupsFlagged.bam//g')
	echo $fl $name
	/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirOut -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
done
for file in $dirOut/*.wig.gz; do
    name=$(basename $file | sed -e 's/.wig.gz//g')
    echo "Processing" $name
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $file $chrsize /gsc/www/bcgsc.ca/downloads/mb/Kongkham/hg19/$name.bw
done
## generate 5mC calls
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
LIB=/home/mbilenky/bin/Solexa_Java
csizes=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes
dirr=/projects/epigenomics/users/mbilenky/CpG/hg19/CG_25_around_chr/
dirw=/projects/epigenomics2/users/lli/glioma/Kongkham/wig/
for file in $dirw/*.wig.gz; do
    name=$(basename $file | sed -e 's/.wig.gz//g')
    echo "$name"
    out="/projects/epigenomics2/users/lli/glioma/Kongkham/CG_25_around_chr/"$name
    mkdir -p $out
    for chr in {1..22} "X" "Y"; do
        chr="chr"$chr
        echo "$chr"
        mkdir -p $out/$chr
        $JAVA -jar -Xmx15G $LIB/RegionsCoverageFromWigCalculator.jar -w $file -r $dirr/$chr.gz -o $out/$chr -c $csizes -n $name
        less $out/$chr/*.coverage | awk '{gsub("chr", "", $1); print $1"_"$2"\t"$4}' > $out/$chr/$chr"."$name.cov
    done
done

