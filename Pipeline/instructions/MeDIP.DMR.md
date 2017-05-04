## Step1. Generate MeDIP only fractional methylation calls
### Calculate MeDIP signal around CpGs and background 
+ Run `RegionsCoverageFromWigCalculator.jar` on MeDIP libraries for CpG +/- 25bp regions and background regions.
+ Input files:
    * CpG +/- 25bp regions for each chr bed files: `/projects/epigenomics/users/mbilenky/CpG/hg19/CG_25_around_chr/`
    * Background regions for each chr bed files: `/projects/epigenomics/users/mbilenky/CpG/hg19/CG_empty_500_chr/`   
    * CpG +/- 25bp regions and Background regions for mm10: `/projects/epigenomics/resources/CpG/mm10/`       
    * MeDIP wig file        
+ Output files:   
    * `<chr>.<name>.coverage` file for each chr in each library: `chr   start   end     normalized coverage     max coverage`
    * `<chr>.<name>.covDist` file for each chr in each library: `coverage  #regions #accumulated regions`
    * `<chr>.<name>.cov` simplified file for Matlab input: `chr_start   normalized coverage`
+ Sample code:
```
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
LIB=/home/mbilenky/bin/Solexa_Java
csizes=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes
dirw=/projects/epigenomics2/users/lli/glioma/Kongkham/wig/
for region in "CG_25_around_chr" "CG_empty_500_chr"; do
	dirr=/projects/epigenomics/users/mbilenky/CpG/hg19/$region/
	for file in $dirw/*.wig.gz; do
		name=$(basename $file | sed -e 's/.wig.gz//g')
		out=/projects/epigenomics2/users/lli/glioma/Kongkham/$region/$name
		mkdir -p $out
		for chr in {1..22} "X" "Y"; do
			chr="chr"$chr
			mkdir -p $out/$chr
			echo $region $name $chr
			$JAVA -jar -Xmx15G $LIB/RegionsCoverageFromWigCalculator.jar -w $file -r $dirr/$chr.gz -o $out/$chr -c $csizes -n $name
			less $out/$chr/*.coverage | awk '{gsub("chr", "", $1); print $1"_"$2"\t"$4}' > $out/$chr/$chr"."$name.cov
		done
	done
done
```

### Generate MeDIP only fractional methylation calls
+ Run Matlab function by Misha: `medip_score2`
+ Input files:    
    * Signal file: `<chr>.<name>.cov` file generated for CpG +/- 25bp    
    * Background file: `<chr>.<name>.covDist` file generated for background regions    
+ Output files:
    * CDF plot for fractional methylation calls per chr
    * Fractional methylation calls per chr: `<chr>.<name>.dip` file: `chr_start   fractional calls`
+ Sample code: (Matlab: `/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab`)
```
addpath /home/mbilenky/matlab/dmr -end

close all; clear all;
set(0,'defaultaxesfontsize',18,'defaultlinelinewidth',2);
names={'IDHmut_hMC_merge.q5.F1028.SET_175','IDHmut_MC_merge.q5.F1028.SET_175','IDHwt_hMC_merge.q5.F1028.SET_175','IDHwt_MC_merge.q5.F1028.SET_175'};
chrs={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'};

for i = 1:4
    name=names{1,i};
    for j = 1:24
        chr=chrs{1,j}
        close all; 
        [l,cc] = textread(['/projects/epigenomics2/users/lli/glioma/Kongkham/CG_25_around_chr/',name,'/',chr,'/',chr,'.',name,'.cov'],'%s %f');
        [c,n,cn]=textread(['/projects/epigenomics2/users/lli/glioma/Kongkham/CpG_empty_500_coverage/',name,'/',chr,'/',chr,'.gz.',name,'.covDist'],'%f %f %f');

        x=c;
        y=cn/max(cn);
        z=cc;
        dip=medip_score2(x,y,z);

        figure('visible','off');box;
        cdfplot(dip);
        xlabel('Fractional methylation');
        ylabel('Fraction of CpGs');
        title(strcat(name,'.',chr));
        dirOut='/projects/epigenomics2/users/lli/glioma/Kongkham/CDF_5mC_plots/';
        nameOut=strcat(dirOut, 'CDF_5mC_',name,'.',chr);
        print(gcf, '-dpdf', strcat(nameOut, '.pdf'));

        t=size(dip); n=t(2);
        fileOut = fopen(strcat('/projects/epigenomics2/users/lli/glioma/Kongkham/CG_25_around_chr/',name,'/',chr,'/',chr,'.',name,'.dip'),'w');
        for i=1:n
			fprintf(fileOut,'%s\t', l{i});
			fprintf(fileOut,'%7.3f\t',dip(i));
			fprintf(fileOut,'\n');
        end
        fclose(fileOut);
    end
end
```

### Combine MeDIP methylation calls from all chrs
+ For each sample, join all `<chr>.<name>.dip` files into one `<name>.dip` file    
+ Sample code:
```
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/CG_25_around_chr/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/
for name in "IDHmut_hMC_merge.q5.F1028.SET_175" "IDHmut_MC_merge.q5.F1028.SET_175" "IDHwt_hMC_merge.q5.F1028.SET_175" "IDHwt_MC_merge.q5.F1028.SET_175"; do
    echo $name
    cat $dirIn/$name/*/*.dip > $dirOut/$name.dip
done
```

## Step2. DMR identification for pairwise comparisons
+ Parameters:
    * m: fractional methylation of one sample need to >= m 
    * delta: minimum difference in fractional calls to call DM CpG    
    * size: max distance between two consecutive CpGs     
    * cut: minimum number of CpGs per DMR
+ Input files: `<name>.dip` file from each library
+ Output files:
    * `DM.summary.stats`: `samples, m, delta, No.of DM CpGs, No.of hypermethylated DM CpGs, No.of hypomethylated DM CpGs`
    * `DMR.summary.stats`: `samples, size, cut, Average length of DMRs, Average No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs`
    * `DM.<cell1>-<donor1>_<cell2>-<donor2>.m<m>.d<delta>.bed `: DM CpGs: `chr   start   end     DM (hyper:1; hypo:-1)   mC1     mC2`
    * `DMR.<cell1>-<donor1>_<cell2>-<donor2>.m<m>.d<delta>.s<size>.c<cut>`: DMRs: `chr   start   end     ID  DM (hyper:1; hypo:-1)   No.of CpGs     length`
    * `DMR.<cell1>-<donor1>_<cell2>-<donor2>.m<m>.d<delta>.s<size>.c<cut>.hyper`: hypermethylated DMRs     
    * `DMR.<cell1>-<donor1>_<cell2>-<donor2>.m<m>.d<delta>.s<size>.c<cut>.hypo`: hypomethylated DMRs
+ Dependency:
    * `/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh`: Shell script to collapse DM CpGs into DMRs        
+ Sample code:
```
################### initial set up ###################
m=0.75    # fractional methylation of one sample need to > m 
delta=0.6 # minimum difference in fractional calls to call DM CpG
size=300  # max distance between two consecutive CpGs
cut=4     # minimum number of CpGs
cd /projects/epigenomics/users/lli/FetalBrain/MeDIP/
dirIn='/projects/epigenomics/users/lli/FetalBrain/MeDIP/'
dirOut=$dirIn/DMR/
mkdir -p $dirOut
> $dirOut/DMR.summary.stats  
> $dirOut/DM.summary.stats       
################### library information ###################
lib1="HS2788"; cell1="Brain"; donor1="HuFNSC01"; name1="HS2788.MeDIP.Brain01.q5.F1028.SET_174";
lib2="HS2790"; cell2="Brain"; donor2="HuFNSC02"; name2="HS2790.MeDIP.Brain02.q5.F1028.SET_174";
################### fixed code for each pairwise comparison from here on ###################
name=$cell1"-"$donor1"_"$cell2"-"$donor2
dm=DM.$name.m$m.d$delta.bed
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: m=$m, delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/$dm
less $dirOut/$dm | grep 'Bad line'
Ndm=($(wc -l $dirOut/$dm)); Nhyper=($(less $dirOut/$dm | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $dm -n $name.m$m.d$delta -s $size -c $cut
```


