%/gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab

addpath /home/mbilenky/matlab -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/';
dirOut='/home/lli/REMC/hESC/';
lib1='hESC_Derived_CD56plus_Ectoderm_Cultured_Cells'; cell1='ecto', donor1='hESC';
lib2='CD184'; cell2='endo', donor2='hESC';

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
dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/';
dirOut='/home/lli/REMC/hESC/';
lib1='hESC_Derived_CD56plus_Ectoderm_Cultured_Cells'; cell1='ecto', donor1='hESC';
lib2='hESC_Derived_CD56plus_Mesoderm_Cultured_Cells'; cell2='meso', donor2='hESC';

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
dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/';
dirOut='/home/lli/REMC/hESC/';
lib1='CD184'; cell1='endo', donor1='hESC';
lib2='hESC_Derived_CD56plus_Mesoderm_Cultured_Cells'; cell2='meso', donor2='hESC';

[coord,geneid,n1,r1]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.exn.A.rpkm'),'%s %s %f %f');
[coord,geneid,n2,r2]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.exn.A.rpkm'),'%s %s %f %f');
[id]=[strcat(coord,'_',geneid)]

X=intersect(id,intersect(idgc,idl))
[C,ix,i]=intersect(id,X);
[C,ixl,i]=intersect(idl,X);
[C,ixgc,i]=intersect(idgc,X);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ixgc)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);


