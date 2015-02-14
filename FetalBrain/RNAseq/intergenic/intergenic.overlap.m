% /gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/';
dirOut='/home/lli/RNAseq/DEfine/intergenic/';
lib1='A03484'; cell1='Brain', donor1='HuFNSC01';
lib2='A07825'; cell2='Brain', donor2='HuFNSC02';
eps=0.001

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

strand='neg'
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')

% apply fold change on overlapping regions
id=vertcat(id,idneg)
r1=vertcat(r1,r1neg)
r2=vertcat(r2,r2neg)
LFC=cellfun(@(a,b) a/b,max(r1,eps),max(r2,eps), 'UniformOutput', 0)
LFC=log2(LFC)
cdfplot(LFC)

