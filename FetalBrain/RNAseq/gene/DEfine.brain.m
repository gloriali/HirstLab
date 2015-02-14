%/gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab

addpath /home/mbilenky/matlab -end

%%%%%%%%%%%%%%%%%%%%% Cortex vs GE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A03474'; cell2='GE', donor2='HuFNSC01';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03475'; cell1='Cortex', donor1='HuFNSC02';
lib2='A03476'; cell2='GE', donor2='HuFNSC02';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A04599'; cell1='Cortex', donor1='HuFNSC03';
lib2='A15295'; cell2='GE', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A15298'; cell1='Cortex', donor1='HuFNSC04';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab
addpath /home/lli/HirstLab/Pipeline/matlab -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A04599'; cell2='Cortex', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A15298'; cell2='Cortex', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03475'; cell1='Cortex', donor1='HuFNSC02';
lib2='A04599'; cell2='Cortex', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03475'; cell1='Cortex', donor1='HuFNSC02';
lib2='A15298'; cell2='Cortex', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03474'; cell1='GE', donor1='HuFNSC01';
lib2='A15295'; cell2='GE', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03474'; cell1='GE', donor1='HuFNSC01';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03476'; cell1='GE', donor1='HuFNSC02';
lib2='A15295'; cell2='GE', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/FetalBrain/RNAseq/DEfine/gene/';
lib1='A03476'; cell1='GE', donor1='HuFNSC02';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);


%output format
%The UP/DN files format:
%Gene id, rpkm1, rpkm2, p-value, Multiple testing corrected p-value
%RPKM.corrected file format:


