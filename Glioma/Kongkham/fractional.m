% /gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab

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
        [c,n,cn]=textread(['/projects/epigenomics2/users/lli/glioma/Kongkham/CG_empty_500_chr/',name,'/',chr,'/',chr,'.gz.',name,'.covDist'],'%f %f %f');

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

