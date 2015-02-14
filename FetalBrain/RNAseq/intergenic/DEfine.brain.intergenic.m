% /gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab

addpath /home/mbilenky/matlab -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
%rcut=2
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A03484'; cell1='Brain', donor1='HuFNSC01', ntotal1=;
lib2='A07825'; cell2='Brain', donor2='HuFNSC02', ntotal2=;

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%iduni=vertcat(id1unique,id2unique)
%r1uni=vertcat(r1unique,zeros(size(id2unique)))
%r2uni=vertcat(zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]
%iduni=vertcat(iduni,id1unique,id2unique)
%r1uni=vertcat(r1uni,r1unique,zeros(size(id2unique)))
%r2uni=vertcat(r2uni,zeros(size(id1unique)),r2unique)
%
%% apply RPKM cut off on unique bins
%index=find(r1uni+r2uni>rcut)
%unique=[iduni(index) num2cell(r1uni(index)) num2cell(r2uni(index))]
%fileOut=fopen(strcat(dirOut,'Unique.',cell1,'-',donor1,'_',cell2,'-',donor2,'_rmin_',num2str(rcut)),'wt')
%fprintf(fileOut,'%s\t%f\t%f\n',unique{,1},unique{,2},unique{,3})
%fclose(fileOut)
%fileOut=strcat(dirOut,'Unique.',cell1,'-',donor1,'_',cell2,'-',donor2,'_rmin_',num2str(rcut))
%dlmcell(fileOut,unique)

% apply DEfine on overlapping regions
id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A03475'; cell2='Cortex', donor2='HuFNSC02';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A03474'; cell1='GE', donor1='HuFNSC01';
lib2='A03476'; cell2='GE', donor2='HuFNSC02';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A04599'; cell1='Cortex', donor1='HuFNSC03';
lib2='A15298'; cell2='Cortex', donor2='HuFNSC04';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A15295'; cell1='GE', donor1='HuFNSC03';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A03473'; cell1='Cortex', donor1='HuFNSC01';
lib2='A03474'; cell2='GE', donor2='HuFNSC01';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A03475'; cell1='Cortex', donor1='HuFNSC02';
lib2='A03476'; cell2='GE', donor2='HuFNSC02';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A04599'; cell1='Cortex', donor1='HuFNSC03';
lib2='A15295'; cell2='GE', donor2='HuFNSC03';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=true; corr1=false; corr2=false; rpkm=true;figs=true; RPKMmin=0.005; Nmin=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdr=0.01
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A15298'; cell1='Cortex', donor1='HuFNSC04';
lib2='A15299'; cell2='GE', donor2='HuFNSC04';

strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
%id=vertcat(id,id1unique,id2unique)
%r1=vertcat(r1,r1unique,zeros(size(id2unique)))
%r2=vertcat(r2,zeros(size(id1unique)),r2unique)

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
%[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
%[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
%[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
%[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]

id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)

n1=zeros(size(r1))
n2=zeros(size(r2))
[cc,nfup,nfdn]=DEfine(id, r1, r2, n1, n2, [], dirOut, [cell1,'-',donor1], [cell2,'-',donor2], out,figs, fdr,corr1,corr2,rpkm,RPKMmin,Nmin);



