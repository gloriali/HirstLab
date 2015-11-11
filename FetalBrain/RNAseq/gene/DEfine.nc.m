%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab

addpath /home/mbilenky/matlab -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03484'; cell1='Brain', donor1='HuFNSC01';
lib2='A07825'; cell2='Brain', donor2='HuFNSC02';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A03475'; cell2='Cortex', donor2='HuFNSC02';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03474'; cell1='GE', donor1='HuFNSC01';
lib2='A03476'; cell2='GE', donor2='HuFNSC02';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%% Cortex vs GE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab

addpath /home/mbilenky/matlab -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A03474'; cell2='GE', donor2='HuFNSC01';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03475'; cell1='Cortex', donor1='HuFNSC02';
lib2='A03476'; cell2='GE', donor2='HuFNSC02';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A04599'; cell1='Cortex', donor1='HuFNSC03';
lib2='A15295'; cell2='GE', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A15298'; cell1='Cortex', donor1='HuFNSC04';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab
addpath /home/mbilenky/matlab -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A04599'; cell2='Cortex', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A15298'; cell2='Cortex', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03475'; cell1='Cortex', donor1='HuFNSC02';
lib2='A04599'; cell2='Cortex', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03475'; cell1='Cortex', donor1='HuFNSC02';
lib2='A15298'; cell2='Cortex', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A04599'; cell1='Cortex', donor1='HuFNSC03';
lib2='A15298'; cell2='Cortex', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03474'; cell1='GE', donor1='HuFNSC01';
lib2='A15295'; cell2='GE', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03474'; cell1='GE', donor1='HuFNSC01';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03476'; cell1='GE', donor1='HuFNSC02';
lib2='A15295'; cell2='GE', donor2='HuFNSC03';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A03476'; cell1='GE', donor1='HuFNSC02';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=true; corr2=true; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdr=0.01

[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.nc.EnsID.GC','%s %f');

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/projects/epigenomics/users/lli/FetalBrain/RNAseq/DEfine/nc/';
lib1='A15295'; cell1='GE', donor1='HuFNSC03';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn,lib1,'/coverage/',lib1,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn,lib2,'/coverage/',lib2,'.G.A.rpkm.nc'),'%s %f %f %f %f %f');

[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ix)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin, eps, maxLim);


%output format
%The UP/DN files format:
%Gene id, rpkm1, rpkm2, p-value, Multiple testing corrected p-value
%RPKM.corrected file format:


