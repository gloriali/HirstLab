#!/bin/sh

dirDE=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/DEfine/
dirPromoter=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/K4me3_K27me3/
dirEnhancer=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/K27ac_K4me1/
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

## promoter: H3K4me3 & H3K27me3
mkdir -p $dirOut/promoter/
paste <(less $dirPromoter/NHA_control.promoter.K4me3_K27me3.bed | sort -k4,4) <(less $dirPromoter/NHAR_control.promoter.K4me3_K27me3.bed | sort -k4,4) <(less $dirPromoter/NHAR_vitc.promoter.K4me3_K27me3.bed | sort -k4,4) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10"\t"$15"\t"$5"_"$10"_"$15}' > $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.bed
less $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.bed | awk '{print $8}' | sort | uniq -c | sort -k1,1rn > $dirOut/promoter/Promoter.wt-mut-vitc.K4me3_K27me3.summary


