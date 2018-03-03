#!/bin/sh

export PATH=/home/rislam/anaconda2/bin/:$PATH
export PYTHONPATH=/home/rislam/anaconda2/lib/python2.7/site-packages
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirDE=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/DEfine/
dirPromoter=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/K4me3_K27me3/
dirEnhancer=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/K27ac_K4me1/
dirK27ac=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/
dir5mC=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/DMR/
dir5hmC=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique3/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/wt_mut_vitc/
mkdir -p $dirOut

## DE
mkdir -p $dirOut/DE/
cat $dirDE/DN.NHA* $dirDE/UP.NHA* | awk '{print $1}' | sort | uniq | join - <(less /home/lli/hg19/hg19v69_genes.TSS.pc.bed | sort -k4,4) -2 4 | awk -F' ' '{gsub("chr", ""); print $2"\t"$3"\t"$4"\t"$1}' > $dirOut/DE/DE.union.TSS.bed
cat <(less $dirDE/UP.NHAR_control_NHA_control.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $1"\tUP"}') <(less $dirDE/DN.NHAR_control_NHA_control.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $1"\tDN"}') | sort -k1,1 > $dirOut/DE/DE.wt-mut.id
awk 'NR==FNR {de[$1]=$2; next} {if($4 in de){print $0"\t"de[$4]}else{print $0"\tST"}}' $dirOut/DE/DE.wt-mut.id $dirOut/DE/DE.union.TSS.bed > $dirOut/DE/DE.wt-mut.TSS.bed
cat <(less $dirDE/UP.NHAR_vitc_NHAR_control.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $1"\tUP"}') <(less $dirDE/DN.NHAR_vitc_NHAR_control.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $1"\tDN"}') | sort -k1,1 > $dirOut/DE/DE.mut-vitc.id
awk 'NR==FNR {de[$1]=$2; next} {if($4 in de){print $0"\t"de[$4]}else{print $0"\tST"}}' $dirOut/DE/DE.mut-vitc.id $dirOut/DE/DE.wt-mut.TSS.bed | awk '{print $0"\t"$5"_"$6}' > $dirOut/DE/DE.wt-mut-vitc.TSS.bed
less $dirOut/DE/DE.wt-mut-vitc.TSS.bed | awk '{print $7}' | sort | uniq -c | sort -k1,1rn > $dirOut/DE/DE.wt-mut-vitc.summary
awk '{print $0 > "'$dirOut'""DE/DE.wt-mut-vitc.TSS."$7".bed"}' $dirOut/DE/DE.wt-mut-vitc.TSS.bed

## promoters: H3K4me3 & H3K27me3
mkdir -p $dirOut/promoter/
paste <(less $dirPromoter/NHA_control.promoter.K4me3_K27me3.bed | sort -k4,4) <(less $dirPromoter/NHAR_control.promoter.K4me3_K27me3.bed | sort -k4,4) <(less $dirPromoter/NHAR_vitc.promoter.K4me3_K27me3.bed | sort -k4,4) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10"\t"$15"\t"$5"_"$10"_"$15}' > $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.bed
less $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.bed | awk '{print $8}' | sort | uniq -c | sort -k1,1rn > $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.summary
less $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.bed | awk '{if($5==$7&&$5!=$6){print $0}}' > $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.reverse.bed
less $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.reverse.bed | awk '{print $8}' | sort | uniq -c | sort -k1,1rn > $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.reverse.summary

## enhancers: H3K27ac & H3K4me1
mkdir -p $dirOut/enhancer/
cat $dirEnhancer/NHA_control.*.bed $dirEnhancer/NHAR_control.*.bed $dirEnhancer/NHAR_vitc.*.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/enhancer/enhancer.uion.bed
echo -e "Sample\tActive\tWeak\tUnmarked" > $dirOut/enhancer/Enhancer.K27ac_K4me1.summary
for file in $dirEnhancer/NHA_control.active*.bed $dirEnhancer/NHAR_control.active*.bed $dirEnhancer/NHAR_vitc.active*.bed; do
    sample=$(basename $file | cut -d'.' -f1); file2=$(echo $file | sed 's/active/weak/')
    echo $sample
    $BEDTOOLS/intersectBed -a $dirOut/enhancer/enhancer.uion.bed -b $file -u | awk '{print $0"\tActive"}' > $dirOut/enhancer/$sample.enhancer.bed
    active=$(less $dirOut/enhancer/$sample.enhancer.bed | wc -l)
    $BEDTOOLS/intersectBed -a $dirOut/enhancer/enhancer.uion.bed -b $file2 -u | $BEDTOOLS/intersectBed -a stdin -b $dirOut/enhancer/$sample.enhancer.bed -v -f 1 | awk '{print $0"\tWeak"}' >> $dirOut/enhancer/$sample.enhancer.bed
    weak=$(expr $(less $dirOut/enhancer/$sample.enhancer.bed | wc -l) - $active)
    $BEDTOOLS/intersectBed -a $dirOut/enhancer/enhancer.uion.bed -b $dirOut/enhancer/$sample.enhancer.bed -v -f 1 | awk '{print $0"\tUnmarked"}' >> $dirOut/enhancer/$sample.enhancer.bed
    unmarked=$(expr $(less $dirOut/enhancer/$sample.enhancer.bed | wc -l) - $active - $weak)
    echo -e $sample"\t"$active"\t"$weak"\t"$unmarked >> $dirOut/enhancer/Enhancer.K27ac_K4me1.summary
done
paste <(less $dirOut/enhancer/NHA_control.enhancer.bed | sort -k1,1 -k2,2n) <(less $dirOut/enhancer/NHAR_control.enhancer.bed | sort -k1,1 -k2,2n) <(less $dirOut/enhancer/NHAR_vitc.enhancer.bed | sort -k1,1 -k2,2n) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10"\t"$15"\t"$5"_"$10"_"$15}' > $dirOut/enhancer/Enhancer.wt-mut-vitc.K27ac_K4me1.bed
less $dirOut/enhancer/Enhancer.wt-mut-vitc.K27ac_K4me1.bed | awk '{print $8}' | sort | uniq -c | sort -k1,1rn > $dirOut/enhancer/Enhancer.wt-mut-vitc.K27ac_K4me1.summary
less $dirOut/enhancer/Enhancer.wt-mut-vitc.K27ac_K4me1.bed | awk '{if($5==$7&&$5!=$6){print $0}}' > $dirOut/enhancer/Enhancer.wt-mut-vitc.K27ac_K4me1.reverse.bed
less $dirOut/enhancer/Enhancer.wt-mut-vitc.K27ac_K4me1.reverse.bed | awk '{print $8}' | sort | uniq -c | sort -k1,1rn > $dirOut/enhancer/Enhancer.wt-mut-vitc.K27ac_K4me1.reverse.summary

## H3K27ac: FC
mkdir -p $dirOut/H3K27ac/
cat <(less $dirK27ac/H3K27ac_NHA_control.vs.input_NHA_control.FDR_0.05.FindER.bed.gz) <(less $dirK27ac/H3K27ac_NHAR_control.vs.input_NHAR_control.FDR_0.05.FindER.bed.gz) <(less $dirK27ac/H3K27ac_NHAR_vitc.vs.input_NHAR_vitc.FDR_0.05.FindER.bed.gz) | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/H3K27ac/H3K27ac.ER.union.bed
dirBW=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/BW/H3K27ac/
multiBigwigSummary BED-file -b $dirBW/H3K27ac_NHA_control.bw $dirBW/H3K27ac_NHAR_control.bw $dirBW/H3K27ac_NHAR_vitc.bw --BED $dirOut/H3K27ac/H3K27ac.ER.union.bed --labels NHA_control NHAR_control NHAR_vitc -out $dirOut/H3K27ac/H3K27ac.ER.union.npz --outRawCounts $dirOut/H3K27ac/H3K27ac.ER.union.RPKM
cd $dirOut/H3K27ac/
less $dirOut/H3K27ac/H3K27ac.ER.union.RPKM | awk 'NR>1{fc=($5+0.0001)/($4+0.0001); if($5>=5&&fc>=2){print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"fc}}' > $dirOut/H3K27ac/H3K27ac.wt_mut.UP.bed
less $dirOut/H3K27ac/H3K27ac.ER.union.RPKM  | awk 'NR>1{fc=($4+0.0001)/($5+0.0001); if($4>=5&&fc>=2){print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"fc}}' > $dirOut/H3K27ac/H3K27ac.wt_mut.DN.bed
less $dirOut/H3K27ac/H3K27ac.ER.union.RPKM | awk 'NR>1{fc=($6+0.0001)/($5+0.0001); if($6>=5&&fc>=2){print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"fc}}' > $dirOut/H3K27ac/H3K27ac.mut_vitc.UP.bed
less $dirOut/H3K27ac/H3K27ac.ER.union.RPKM | awk 'NR>1{fc=($5+0.0001)/($6+0.0001); if($5>=5&&fc>=2){print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"fc}}' > $dirOut/H3K27ac/H3K27ac.mut_vitc.DN.bed
cat <(less $dirOut/H3K27ac/H3K27ac.wt_mut.UP.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\tUP"}') <(less $dirOut/H3K27ac/H3K27ac.wt_mut.DN.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\tDN"}') > $dirOut/H3K27ac/H3K27ac.wt_mut.bed
$BEDTOOLS/intersectBed -a $dirOut/H3K27ac/H3K27ac.ER.union.bed -b $dirOut/H3K27ac/H3K27ac.wt_mut.bed -v | awk '{print $0"\tST"}' >> $dirOut/H3K27ac/H3K27ac.wt_mut.bed
cat <(less $dirOut/H3K27ac/H3K27ac.mut_vitc.UP.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\tUP"}') <(less $dirOut/H3K27ac/H3K27ac.mut_vitc.DN.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\tDN"}') > $dirOut/H3K27ac/H3K27ac.mut_vitc.bed
$BEDTOOLS/intersectBed -a $dirOut/H3K27ac/H3K27ac.ER.union.bed -b $dirOut/H3K27ac/H3K27ac.mut_vitc.bed -v | awk '{print $0"\tST"}' >> $dirOut/H3K27ac/H3K27ac.mut_vitc.bed
paste <(less $dirOut/H3K27ac/H3K27ac.wt_mut.bed | sort -k1,1 -k2,2n) <(less $dirOut/H3K27ac/H3K27ac.mut_vitc.bed | sort -k1,1 -k2,2n) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10"\t"$5"_"$10}' > $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.bed
less $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.bed | awk '{print $7}' | sort | uniq -c | sort -k1,1rn > $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.summary
less $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.bed | awk '{if($7=="UP_DN"||$7=="DN_UP"){print $0}}' > $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.reverse.bed 
less $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.bed | awk '$7 !~ /ST_ST/ {print > "'$dirOut'""/H3K27ac/group."$7".bed"}' 
computeMatrix scale-regions -R $dirOut/H3K27ac/group.*.bed -S $dirBW/H3K27ac_NHA_control.bw $dirBW/H3K27ac_NHAR_control.bw $dirBW/H3K27ac_NHAR_vitc.bw -out $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.gz --startLabel start --endLabel end -bs 20
plotHeatmap -m $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.gz -out $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.png --colorMap coolwarm --xAxisLabel "H3K27ac ER (bp)" --startLabel start --endLabel end -max 50 --samplesLabel wt mut VitC --regionsLabel DN_DN DN_ST DN_UP ST_DN ST_UP UP_DN UP_ST UP_UP
computeMatrix scale-regions -R $dirOut/H3K27ac/group.UP_DN.bed $dirOut/H3K27ac/group.DN_UP.bed -S $dirBW/H3K27ac_NHA_control.bw $dirBW/H3K27ac_NHAR_control.bw $dirBW/H3K27ac_NHAR_vitc.bw -out $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.reverse.gz --startLabel start --endLabel end -bs 20
plotHeatmap -m $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.reverse.gz -out $dirOut/H3K27ac/H3K27ac.wt-mut-vitc.reverse.png --colorMap coolwarm --xAxisLabel "H3K27ac ER (bp)" --startLabel start --endLabel end --samplesLabel wt mut VitC --regionsLabel UP_DN DN_UP
