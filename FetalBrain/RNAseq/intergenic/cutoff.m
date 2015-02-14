strand='pos'
% BedTools-intersectBed: find overlapping regions of two bed files & unique bins
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
% get id & r1,r2 for overlapping bins
%[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
%[id]=[strcat(chr,':',chrstart,'-',chrend,'<1')]
%[chr1,start1,end1,r1,chr2,start2,end2,r2]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
% get id & r1,r2 for unique bins
[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<1')]
[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<1')]
iduni=vertcat(id1unique,id2unique)
r1uni=vertcat(r1unique,zeros(size(id2unique)))
r2uni=vertcat(zeros(size(id1unique)),r2unique)

strand='neg'
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED >' dirOut lib1 '_' lib2 '.' strand '.intersect.bed'])
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -wa -wb >' dirOut lib1 '_' lib2 '.' strand '.bed'])
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib1 '_not_' lib2 '.' strand '.bed'])
%[status]=system(['/home/lli/bedtools-2.17.0/bin/intersectBed -a' ' ' dirIn lib2 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib2 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -b' ' ' dirIn lib1 '/intergenic/coverage/genome_200/genomeProfile.bin200.step200.' lib1 '.q5.F1540.' strand '.RPKM.0.5.filtered.mappable.BED -v >' dirOut lib2 '_not_' lib1 '.' strand '.bed'])
%[chr,chrstart,chrend,value]=textread([dirOut lib1 '_' lib2 '.' strand '.intersect.bed'],'%s %s %s %f')
%[idneg]=[strcat(chr,':',chrstart,'-',chrend,'<-1')]
%[chr1,start1,end1,r1neg,chr2,start2,end2,r2neg]=textread([dirOut lib1 '_' lib2 '.' strand '.bed'],'%s %f %f %f %s %f %f %f')
[chr1unique,start1unique,end1unique,r1unique]=textread([dirOut lib1 '_not_' lib2 '.' strand '.bed'],'%s %s %s %f')
[id1unique]=[strcat(chr1unique,':',start1unique,'-',end1unique,'<-1')]
[chr2unique,start2unique,end2unique,r2unique]=textread([dirOut lib2 '_not_' lib1 '.' strand '.bed'],'%s %s %s %f')
[id2unique]=[strcat(chr2unique,':',start2unique,'-',end2unique,'<-1')]
iduni=vertcat(iduni,id1unique,id2unique)
r1uni=vertcat(r1uni,r1unique,zeros(size(id2unique)))
r2uni=vertcat(r2uni,zeros(size(id1unique)),r2unique)
cdfplot(r1uni)
cdfplot(r2uni)

