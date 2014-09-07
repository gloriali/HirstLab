#!/bin/sh

grep GATA3 /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/Overlap_with_TFs/DMRs.-1_TFs > /home/lli/REMC/UMR/GATA3.lum.UMR
grep FOXA1 /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/Overlap_with_TFs/DMRs.-1_TFs > /home/lli/REMC/UMR/FOXA1.lum.UMR
grep ZNF217 /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/Overlap_with_TFs/DMRs.-1_TFs > /home/lli/REMC/UMR/ZNF217.lum.UMR
grep EGR1 /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/Overlap_with_TFs/DMRs.1_TFs > /home/lli/REMC/UMR/EGR1.myo.UMR

less /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/DMRs.p0.0005.s200.c3.-1_vs_genes/DMRs.-1_genes.pc > /home/lli/REMC/UMR/lum.UMR.genes
less /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/DMRs.p0.0005.s200.c3.1_vs_genes/DMRs.1_genes.pc > /home/lli/REMC/UMR/myo.UMR.genes



