%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab
addpath /home/mbilenky/matlab -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_19 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_19 vs Cortex02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_19/';
sample1='CEMT_19'; name1='CEMT_19';
sample2='A03475.Cortex02'; name2='Cortex02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_19 vs GE02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_19/';
sample1='CEMT_19'; name1='CEMT_19';
sample2='A03476.GE02'; name2='GE02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_19 vs Cortex04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_19/';
sample1='CEMT_19'; name1='CEMT_19';
sample2='A15298.Cortex04'; name2='Cortex04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_19 vs GE04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_19/';
sample1='CEMT_19'; name1='CEMT_19';
sample2='A15299.GE04'; name2='GE04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_21 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_21 vs Cortex02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_21/';
sample1='CEMT_21'; name1='CEMT_21';
sample2='A03475.Cortex02'; name2='Cortex02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_21 vs GE02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_21/';
sample1='CEMT_21'; name1='CEMT_21';
sample2='A03476.GE02'; name2='GE02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_21 vs Cortex04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_21/';
sample1='CEMT_21'; name1='CEMT_21';
sample2='A15298.Cortex04'; name2='Cortex04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_21 vs GE04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_21/';
sample1='CEMT_21'; name1='CEMT_21';
sample2='A15299.GE04'; name2='GE04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_22 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_22 vs Cortex02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_22/';
sample1='CEMT_22'; name1='CEMT_22';
sample2='A03475.Cortex02'; name2='Cortex02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_22 vs GE02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_22/';
sample1='CEMT_22'; name1='CEMT_22';
sample2='A03476.GE02'; name2='GE02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_22 vs Cortex04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_22/';
sample1='CEMT_22'; name1='CEMT_22';
sample2='A15298.Cortex04'; name2='Cortex04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_22 vs GE04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_22/';
sample1='CEMT_22'; name1='CEMT_22';
sample2='A15299.GE04'; name2='GE04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_23 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_23 vs Cortex02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_23/';
sample1='CEMT_23'; name1='CEMT_23';
sample2='A03475.Cortex02'; name2='Cortex02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_23 vs GE02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_23/';
sample1='CEMT_23'; name1='CEMT_23';
sample2='A03476.GE02'; name2='GE02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_23 vs Cortex04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_23/';
sample1='CEMT_23'; name1='CEMT_23';
sample2='A15298.Cortex04'; name2='Cortex04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_23 vs GE04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_23/';
sample1='CEMT_23'; name1='CEMT_23';
sample2='A15299.GE04'; name2='GE04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_47 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_47 vs Cortex02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_47/';
sample1='CEMT_47'; name1='CEMT_47';
sample2='A03475.Cortex02'; name2='Cortex02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_47 vs GE02 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_47/';
sample1='CEMT_47'; name1='CEMT_47';
sample2='A03476.GE02'; name2='GE02';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_47 vs Cortex04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_47/';
sample1='CEMT_47'; name1='CEMT_47';
sample2='A15298.Cortex04'; name2='Cortex04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%% CEMT_47 vs GE04 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine//CEMT_47/';
sample1='CEMT_47'; name1='CEMT_47';
sample2='A15299.GE04'; name2='GE04';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);

