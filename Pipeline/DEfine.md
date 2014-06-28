## Matlab path on xhost
`/gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab`

## DEfine function
* Misha's directory: `addpath /home/mbilenky/matlab -end`
* My copy of DEfinev.0.9.2: `addpath /home/lli/bin/matlab -end`        

## DEfine on genes
* Sample code        
```
addpath /home/lli/bin/matlab -end
clear all; close all;
fdr=0.015;    
RPKMmin=0.005; Nmin=25; out=true; corr1=true; corr2=true; rpkm=true; figs=true; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID.GC','%s %f');
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/REMC/breast/gene/';
% for each library
lib1='A17918'; cell1='lum', donor1='RM084';
lib2='A17919'; cell2='myo', donor2='RM084';
[id,n1,r1,rmi,ra,rma]=textread(strcat(dirIn1,lib1,'/coverage/',lib1,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn2,lib2,'/coverage/',lib2,'.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl, r1(ix), r2(ix), n1(ix), n2(ix), [gl,gc], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);
```       

## DEfine on exons
* Sample code        
```
addpath /home/lli/bin/matlab -end
clear all; close all;
fdr=0.015;    
RPKMmin=0.005; Nmin=25; out=true; corr1=true; corr2=true; rpkm=true; figs=true; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.length','%s %f');
[idgc,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.GC','%s %f');
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/REMC/breast/exon/';
% for each library
lib1='A17918'; cell1='lum', donor1='RM084';
lib2='A17919'; cell2='myo', donor2='RM084';
[coord,geneid,n1,r1]=textread(strcat(dirIn1,lib1,'/coverage/',lib1,'.G.exn.A.rpkm'),'%s %s %f %f');
[coord,geneid,n2,r2]=textread(strcat(dirIn2,lib2,'/coverage/',lib2,'.G.exn.A.rpkm'),'%s %s %f %f');
[id]=[strcat(coord,'_',geneid)]
X=intersect(id,intersect(idgc,idl))
[C,ix,i]=intersect(id,X);
[C,ixl,i]=intersect(idl,X);
[C,ixgc,i]=intersect(idgc,X);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ix), r2(ix), n1(ix), n2(ix), [gl(ixl),gc(ixgc)], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);
```

## Output 
* The UP/DN files format: `Gene id    rpkm1   rpkm2   p-value Multiple testing corrected p-value`          
* RPKM.corrected file format: `Gene id    RPKM1corrected  N1 corrected    RPKM2 corrected N2 corrected    p-val1  p-valcorr1  p-val2    p-valcorr2`         
* Boxplot_GC: boxplot for GC bias correction        
* Boxplot_L: boxplot for gene/exon length bias correction          
* Plot2D: result summary 2D plot         
