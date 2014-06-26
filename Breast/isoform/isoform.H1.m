%/gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab

addpath /home/lli/bin/matlab/ -end

% DE exons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn1='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
dirIn2='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirOut='/home/lli/REMC/breast/H1/exon/';
lib1='A17919'; cell1='myo', donor1='RM084';
lib2='H1_r1a'; cell2='H1', donor2='H1_r1a';

[coord,geneid,n1,r1]=textread(strcat(dirIn1,lib1,'/coverage/',lib1,'.G.exn.A.rpkm'),'%s %s %f %f');
[coord,geneid,n2,r2]=textread(strcat(dirIn2,lib2,'/coverage/',lib2,'.G.exn.A.rpkm'),'%s %s %f %f');
[id]=[strcat(coord,'_',geneid)]

X=intersect(id,intersect(idgc,idl))
[C,ix,i]=intersect(id,X);
[C,ixl,i]=intersect(idl,X);
[C,ixgc,i]=intersect(idgc,X);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ixgc)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn1='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
dirIn2='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.015
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');

dirOut='/home/lli/REMC/breast/H1/exon/';
lib1='A17919'; cell1='myo', donor1='RM084';
lib2='H1_r2a'; cell2='H1', donor2='H1_r2a';

[coord,geneid,n1,r1]=textread(strcat(dirIn1,lib1,'/coverage/',lib1,'.G.exn.A.rpkm'),'%s %s %f %f');
[coord,geneid,n2,r2]=textread(strcat(dirIn2,lib2,'/coverage/',lib2,'.G.exn.A.rpkm'),'%s %s %f %f');
[id]=[strcat(coord,'_',geneid)]

X=intersect(id,intersect(idgc,idl))
[C,ix,i]=intersect(id,X);
[C,ixl,i]=intersect(idl,X);
[C,ixgc,i]=intersect(idgc,X);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ixgc)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

% DE genes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn1='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
dirIn2='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.015

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirOut='/home/lli/REMC/breast/H1/gene/';
lib1='A17919'; cell1='myo', donor1='RM084';
lib2='H1_r1a'; cell2='H1', donor2='H1_r1a';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn1,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn2,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn1='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
dirIn2='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.015

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirOut='/home/lli/REMC/breast/H1/gene/';
lib1='A17919'; cell1='myo', donor1='RM084';
lib2='H1_r2a'; cell2='H1', donor2='H1_r2a';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn1,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn2,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);
