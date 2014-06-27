%/gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab

addpath /home/lli/bin/matlab/ -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirOut='/home/lli/REMC/breast/vHMEC_fibr/exon/';
lib1='A18760'; cell1='fibr', donor1='RM070';
lib2='HS2263'; cell2='vHMEC', donor2='RM035';

[coord,geneid,n1,r1]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.exn.A.rpkm'),'%s %s %f %f');
[coord,geneid,n2,r2]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.exn.A.rpkm'),'%s %s %f %f');
[id]=[strcat(coord,'_',geneid)]

X=intersect(id,intersect(idgc,idl))
[C,ix,i]=intersect(id,X);
[C,ixl,i]=intersect(idl,X);
[C,ixgc,i]=intersect(idgc,X);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ixgc)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirOut='/home/lli/REMC/breast/vHMEC_fibr/exon/';
lib1='A18760'; cell1='fibr', donor1='RM070';
lib2='A18761'; cell2='vHMEC', donor2='RM071';

[coord,geneid,n1,r1]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.exn.A.rpkm'),'%s %s %f %f');
[coord,geneid,n2,r2]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.exn.A.rpkm'),'%s %s %f %f');
[id]=[strcat(coord,'_',geneid)]

X=intersect(id,intersect(idgc,idl))
[C,ix,i]=intersect(id,X);
[C,ixl,i]=intersect(idl,X);
[C,ixgc,i]=intersect(idgc,X);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ixgc)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirOut='/home/lli/REMC/breast/vHMEC_fibr/exon/';
lib1='A18472'; cell1='fibr', donor1='RM071';
lib2='HS2263'; cell2='vHMEC', donor2='RM035';

[coord,geneid,n1,r1]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.exn.A.rpkm'),'%s %s %f %f');
[coord,geneid,n2,r2]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.exn.A.rpkm'),'%s %s %f %f');
[id]=[strcat(coord,'_',geneid)]

X=intersect(id,intersect(idgc,idl))
[C,ix,i]=intersect(id,X);
[C,ixl,i]=intersect(idl,X);
[C,ixgc,i]=intersect(idgc,X);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ixgc)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

