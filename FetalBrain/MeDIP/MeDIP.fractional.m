% /gsc/software/linux-x86_64-centos5/matlab-2012b/bin/matlab

addpath /home/mbilenky/matlab/dmr -end

close all; clear all;
set(0,'defaultaxesfontsize',18,'defaultlinelinewidth',2);
names={'HS2788.MeDIP.Brain01.q5.F1028.SET_174','HS2790.MeDIP.Brain02.q5.F1028.SET_174','HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174','HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174','HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157','HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166'};
chrs={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'};

for i = 1:6
    name=names{1,i};name1 = regexprep(name, '.q5.F1028.SET_[0-9]+', '')
    for j = 1:24
        chr=chrs{1,j}
        close all; 
        [l,cc] = textread(['/projects/epigenomics/users/lli/FetalBrain/MeDIP/CG_25_around_chr/',name,'/',chr,'/',chr,'.',name,'.cov'],'%s %f');
        [c,n,cn]=textread(['/projects/mbilenky/REMC/brain/MeDIP/analysis/CpG_empty_500_coverage/',chr,'/',chr,'.gz.',name1,'.covDist'],'%f %f %f');
         
        x=c;
        y=cn/max(cn);
        z=cc;
        dip=medip_score2(x,y,z);
        
        figure('visible','off');box;
        cdfplot(dip);
        xlabel('Fractional methylation');
        ylabel('Fraction of CpGs');
        title(strcat(name1,'.',chr));
        dirOut='/projects/epigenomics/users/lli/FetalBrain/MeDIP/CDF_5mC_plots/';
        nameOut=strcat(dirOut, 'CDF_5mC_',name1,'.',chr);
        print(gcf, '-dpdf', strcat(nameOut, '.pdf'));
         
        t=size(dip); n=t(2);
        fileOut = fopen(strcat('/projects/epigenomics/users/lli/FetalBrain/MeDIP/CG_25_around_chr/',name,'/',chr,'/',chr,'.',name,'.dip'),'w');
        for i=1:n
        fprintf(fileOut,'%s\t', l{i});
        fprintf(fileOut,'%7.3f\t',dip(i));
        fprintf(fileOut,'\n');
        end
        fclose(fileOut);
    end
end


% Resulting files can be joined together into the matrix per chromosome.
