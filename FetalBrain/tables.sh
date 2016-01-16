#!/bin/sh

dirIn='/projects/epigenomics/users/lli/FetalBrain/'
dirOut='/projects/epigenomics/users/lli/FetalBrain/Tables/'

# Table S1. Summary of molecular libraries and accessions.
cp /home/lli/HirstLab/FetalBrain/FetalBrainLibrariesDetail.tsv $dirOut/TableS1.txt

# Table S2. List of UMRs between MZ twins, NPCs and across developmental stages.
echo -e "chr\tstart\tend\tID\tUMR\tNo.of.CpGs\tlength\tSample" > $dirOut/TableS2_MZ.txt
less $dirIn/MeDIP/DMR/DMR.Brain-HuFNSC01_Brain-HuFNSC02.m0.75.d0.6.s300.c4 | awk '{if($5==1){DM="UMR.in.Subject2"}else{DM="UMR.in.Subject1"} print $1"\t"$2"\t"$3"\t"$4"\t"DM"\t"$6"\t"$7"\tBrain"}' >> $dirOut/TableS2_MZ.txt
less $dirIn/MeDIP/DMR/DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.m0.75.d0.6.s300.c4 | awk '{if($5==1){DM="UMR.in.Subject2"}else{DM="UMR.in.Subject1"} print $1"\t"$2"\t"$3"\t"$4"\t"DM"\t"$6"\t"$7"\tCortex"}' >> $dirOut/TableS2_MZ.txt
less $dirIn/MeDIP/DMR/DMR.GE-HuFNSC01_GE-HuFNSC02.m0.75.d0.6.s300.c4 | awk '{if($5==1){DM="UMR.in.Subject2"}else{DM="UMR.in.Subject1"} print $1"\t"$2"\t"$3"\t"$4"\t"DM"\t"$6"\t"$7"\tGE"}' >> $dirOut/TableS2_MZ.txt
echo -e "chr\tstart\tend\tID\tUMR\tNo.of.CpGs\tlength\tGW" > $dirOut/TableS2_NPC.txt
less $dirIn/WGBS/DMR/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3 | awk '{if($5==1){DM="UMR.in.GE"}else{DM="UMR.in.Cortex"} print $1"\t"$2"\t"$3"\t"$4"\t"DM"\t"$6"\t"$7"\tGW13"}' >> $dirOut/TableS2_NPC.txt
less $dirIn/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3 | awk '{if($5==1){DM="UMR.in.GE"}else{DM="UMR.in.Cortex"} print $1"\t"$2"\t"$3"\t"$4"\t"DM"\t"$6"\t"$7"\tGW17"}' >> $dirOut/TableS2_NPC.txt
echo -e "chr\tstart\tend\tID\tUMR\tNo.of.CpGs\tlength\tSample" > $dirOut/TableS2_GW.txt
less $dirIn/GW/DMR/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3 | awk '{if($5==1){DM="UMR.in.GW13"}else{DM="UMR.in.GW17"} print $1"\t"$2"\t"$3"\t"$4"\t"DM"\t"$6"\t"$7"\tCortex"}' >> $dirOut/TableS2_GW.txt
less $dirIn/GW/DMR/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3 | awk '{if($5==1){DM="UMR.in.GW13"}else{DM="UMR.in.GW17"} print $1"\t"$2"\t"$3"\t"$4"\t"DM"\t"$6"\t"$7"\tGE"}' >> $dirOut/TableS2_GW.txt

# Table S3. Histone modification enriched regions from FindER.
for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3; do
    echo -e "chr\tstart\tend\tID\tSample\tSubject" > $dirOut'/TableS3_'$mark'.txt'
    cd $dirIn/ChIPseq/ER/$mark
    for file in FindER_scan.*.bed; do
        sample=$(echo $file | cut -d'.' -f4)
        cell=$(echo $sample | sed -e 's/[0-9]\+//g')
        donor='Subject'$(echo $sample | sed -e 's/[A-Za-z]\+0//g')
        echo $cell $donor $file
        less $file | awk '{print $0"\t"$1":"$2"-"$3"\t""'$cell'""\t""'$donor'"}' >> $dirOut'/TableS3_'$mark'.txt'
    done
done
echo -e "chr\tstart\tend\tID" > $dirOut/TableS3_coreEnhancer.txt
less $dirIn/ChIPseq/ER/H3K4me1/core/core_enhancers.bed >> $dirOut/TableS3_coreEnhancer.txt
echo -e "chr\tstart\tend\tID\tComparison\tunique" > $dirOut/TableS3_uniqueEnhancer.txt
less $dirIn/ChIPseq/ER/H3K4me1/unique/unique_enhancer.MZ.HuFNSC01.bed | awk '{print $0"\tMZ\tSubject1"}' >> $dirOut/TableS3_uniqueEnhancer.txt
less $dirIn/ChIPseq/ER/H3K4me1/unique/unique_enhancer.MZ.HuFNSC02.bed | awk '{print $0"\tMZ\tSubject2"}' >> $dirOut/TableS3_uniqueEnhancer.txt
less $dirIn/ChIPseq/ER/H3K4me1/unique/unique_enhancer.Neurospheres.Cortex.bed | awk '{print $0"\tNPCs\tCortex"}' >> $dirOut/TableS3_uniqueEnhancer.txt
less $dirIn/ChIPseq/ER/H3K4me1/unique/unique_enhancer.Neurospheres.GE.bed | awk '{print $0"\tNPCs\tGE"}' >> $dirOut/TableS3_uniqueEnhancer.txt
less $dirIn/ChIPseq/ER/H3K4me1/unique/unique_enhancer.GW.GW13.bed | awk '{print $0"\tGW\tGW13"}' >> $dirOut/TableS3_uniqueEnhancer.txt
less $dirIn/ChIPseq/ER/H3K4me1/unique/unique_enhancer.GW.GW17.bed | awk '{print $0"\tGW\tGW17"}' >> $dirOut/TableS3_uniqueEnhancer.txt

# Table S4. List of gene promoters differentially marked by H3K4me3 and H3K27me3 between MZ twins, NPCs and across developmental stages.
## see FindER.R 

# Table S5. Differentially expressed genes between MZ twins, NPCs and across developmental stages.
echo -e "ID\tRPKM1\tRPKM2\tp-value\tFDR\tDE\tSample" > $dirOut/TableS5_MZ.txt
cd $dirIn/RNAseq/DEfine/gene/MZ/
for file in *Brain* *Cortex* *GE*; do
    DE=$(echo $file | cut -d'.' -f1)
    cell=$(echo $file | cut -d'.' -f2 | cut -d'-' -f1)
    echo $DE $cell $file
    less $file | awk '{if("'$DE'"=="UP"){de="UP.in.Subject1"}else{de="UP.in.Subject2"}; print $0"\t"de"\t""'$cell'"}' >> $dirOut/TableS5_MZ.txt
done
echo -e "ID\tRPKM1\tRPKM2\tp-value\tFDR\tDE\tSubject" > $dirOut/TableS5_NPC.txt
cd $dirIn/RNAseq/DEfine/gene/NPC/
for file in *HuFNSC01* *HuFNSC02* *HuFNSC03* *HuFNSC04*; do
    DE=$(echo $file | cut -d'.' -f1)
    donor=$(echo $file | cut -d'.' -f2 | cut -d'-' -f2 | cut -d'_' -f1 | sed -e 's/HuFNSC0/Subject/g')
    echo $DE $donor $file
    less $file | awk '{if("'$DE'"=="UP"){de="UP.in.Cortex"}else{de="UP.in.GE"}; print $0"\t"de"\t""'$donor'"}' >> $dirOut/TableS5_NPC.txt
done
echo -e "ID\tRPKM1\tRPKM2\tp-value\tFDR\tDE\tSample" > $dirOut/TableS5_GW.txt
cd $dirIn/RNAseq/DEfine/gene/GW/
for file in *HuFNSC03*HuFNSC04*; do
    DE=$(echo $file | cut -d'.' -f1)
    cell=$(echo $file | cut -d'.' -f2 | cut -d'-' -f1)
    echo $DE $cell $file
    less $file | awk '{if("'$DE'"=="UP"){de="UP.in.GW15"}else{de="UP.in.GW13"}; print $0"\t"de"\t""'$cell'"}' >> $dirOut/TableS5_GW.txt
done
for file in *HuFNSC01*HuFNSC04* *HuFNSC02*HuFNSC04*; do
    DE=$(echo $file | cut -d'.' -f1)
    cell=$(echo $file | cut -d'.' -f2 | cut -d'-' -f1)
    echo $DE $cell $file
    less $file | awk '{if("'$DE'"=="UP"){de="UP.in.GW17"}else{de="UP.in.GW13"}; print $0"\t"de"\t""'$cell'"}' >> $dirOut/TableS5_GW.txt
done
for file in *HuFNSC01*HuFNSC03* *HuFNSC02*HuFNSC03*; do
    DE=$(echo $file | cut -d'.' -f1)
    cell=$(echo $file | cut -d'.' -f2 | cut -d'-' -f1)
    echo $DE $cell $file
    less $file | awk '{if("'$DE'"=="UP"){de="UP.in.GW17"}else{de="UP.in.GW15"}; print $0"\t"de"\t""'$cell'"}' >> $dirOut/TableS5_GW.txt
done

# Table S6. GWAS sites overlapped with NPC core enhancers.      
echo -e "enhancerChr\tenhancerStart\tenhancerEnd\tenhancerID\tgwasChr\tgwasStart\tgwasEnd\tgwasID\ttrait\tgenes" > $dirOut/TableS6.txt
less /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/core_enhancers.GWAS.txt >> $dirOut/TableS6.txt
