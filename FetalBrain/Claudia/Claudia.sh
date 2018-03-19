#!/bin/sh

# ================ RNAseq ===================
# align to hg19v65
## jaguar alignment
ref=/home/pubseq/genomes/Homo_sapiens/hg19a/jaguar/1.7.5/ens65only/bwa_ind/transcriptome/100/ref.fa
ens=hg19_ens65
dirIn=/projects/sftp/ckleinman/incoming/fetal_brain_RNAseq/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/RNAseq/bam/
mkdir -p $dirOut
cd $dirIn
for f1 in *_R1.fastq.gz; do
    f2=$(echo $f1 | sed 's/_R1./_R2./')
    id=$(echo $f1 | cut -d'.' -f5)
    age=$(less /projects/epigenomics3/users/lli/Claudia/SampleInfo.txt | awk '{if($1=="'$id'"){print $2}}')
    name=$id.$age
    echo $name
    /home/lli/HirstLab/Pipeline/shell/jaguar.sh -i $dirIn -o $dirOut -f1 $f1 -f2 $f2 -n $name -r $ref -v $ens
done
function JRalign {
    ens=hg19_ens65
    dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/RNAseq/bam/
    name=$1
    rm -rf $dirOut/$name/
    /home/mbilenky/bin/Solexa_Shell/RunJR.sh $dirOut/$name".sortedByName.bam" $dirOut/$name $ens &> $dirOut/run/$name.j.log
}
export -f JRalign
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/RNAseq/bam/
cd $dirOut
ls *.sortedByName.bam | awk '{gsub(".sortedByName.bam", ""); print $0}' > $dirOut/List.txt
cat List.txt | parallel --gnu JRalign 
## QC and RPKM
samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/RNAseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/RNAseq/RPKM/
ens=hg19v65
mkdir -p $dirOut
for bam in $dirIn/*_withJunctionsOnGenome_dupsFlagged.bam; do
    name=$(basename $bam | sed 's/_withJunctionsOnGenome_dupsFlagged.bam//g')
    echo $name
    rm -rf $dirOut/$name/
    /home/lli/bin/Solexa_Shell/src/RNAseqMaster.sh $bam $name $dirOut $ens S 0 "1,1,1,1,1" /projects/epigenomics/resources/ $JAVA
done
cd /projects/epigenomics3/epigenomics3_results/users/lli/Claudia/RNAseq/RPKM/
echo "ENSG" > Claudia.RPKM
echo -e "Sample\tENSG\tRPKM" > RPKM.long
less M291.77d_11w/coverage/M291.77d_11w.G.A.rpkm.pc | awk '{print $1}' >> Claudia.RPKM
for file in */coverage/*.G.A.rpkm.pc; do
    lib=$(basename $file | sed 's/.G.A.rpkm.pc//g')
    echo -e "ENSG\t$lib" > x
    less $file | awk '{print $1"\t"$3}' >> x
    join Claudia.RPKM x | sed 's/ /\t/g' > y
    mv y Claudia.RPKM
    less $file | awk '{print "'$lib'""\t"$1"\t"$3}' >> RPKM.long
done
rm x

# ================ 450K ===================
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirDMR=/projects/epigenomics/users/lli/FetalBrain/GW/DMR/
cd /projects/epigenomics3/epigenomics3_results/users/lli/Claudia/DNAme/
less DNAme450K.txt | awk '$4 ~ /cg/ {print "chr"$0}' | sort -k1,1 -k2,2n > DNAme450K.id
$BEDTOOLS/intersectBed -a DNAme450K.id -b $dirDMR/DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -u > DMR.HuFNSC02_HuFNSC04.450K.hyper.bed
$BEDTOOLS/intersectBed -a DNAme450K.id -b $dirDMR/DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -u > DMR.HuFNSC02_HuFNSC04.450K.hypo.bed
$BEDTOOLS/intersectBed -a DNAme450K.id -b <(echo -e "chr21\t34396154\t34400153\tENSG00000205927") -u | $BEDTOOLS/intersectBed -a stdin -b $dirDMR/DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -u > Olig2.promoter.450K.bed
