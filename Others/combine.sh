#!/bin/sh

less AML_A08912/coverage/AML_A08912.G.A.rpkm.pc | awk '{print $1}' > AML_master_RNAseq.rpkm.pc
header=("GeneID\t")
for folder in AML_A*; do
    echo $folder;
    less $folder/coverage/$folder.G.A.rpkm.pc | awk '{print $3}' > temp.rpkm;
    paste --delimiters="\t" AML_master_RNAseq.rpkm.pc temp.rpkm > temp.master;
    cp temp.master AML_master_RNAseq.rpkm.pc
    header+=$folder"\t"
done
rm temp*
echo -e ${header[@]} > header
cp AML_master_RNAseq.rpkm.pc AML_master_RNAseq.rpkm.pc1
cat header AML_master_RNAseq.rpkm.pc1 > AML_master_RNAseq.rpkm.pc
rm AML_master_RNAseq.rpkm.pc1 header
