#!/bin/sh

# link RPKM files
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/RPKM/
mkdir -p $dirOut
for lib in 19 21 22 23 47; do
    echo $lib;
    ln -s $dirIn/CEMT_$lib/bams/RNA-Seq/qca/*/coverage/*.G.A.rpkm.pc $dirOut/CEMT_$lib.G.A.rpkm.pc
done

# IDH1/2 expression
dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/RPKM/
> $dirOut/IDH.RPKM
for lib in 19 21 22 23 47; do
    echo $lib;
    less $dirOut/CEMT_$lib.G.A.rpkm.pc | awk '$1 ~ /ENSG00000138413/ {print $1"\tIDH1\tCEMT_""'$lib'""\t"$3}' >> $dirOut/IDH.RPKM
    less $dirOut/CEMT_$lib.G.A.rpkm.pc | awk '$1 ~ /ENSG00000182054/ {print $1"\tIDH2\tCEMT_""'$lib'""\t"$3}' >> $dirOut/IDH.RPKM
done


