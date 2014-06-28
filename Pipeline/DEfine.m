function [cc,nfup,nfdn]=DEfine(id, ir1, ir2, in1, in2, extpar, dirOut, name1, name2, out, figs, fdrth, corr1, corr2, rpkm, rmin, nmin)

suffix=strcat(name1,'_',name2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RE-NORMALIZATION and PAIRWISE COMPARISON
%
% DEfine v 0.9.1 March 29, 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%out=true;
%nmin=30; rmin=0.005;	% thresholds in number of reads and RPKM, no input paramters
may=6;
xstep=1;
nint=10;
h=15;

% parameters
fitp=1;	% defines a range (percentile) for selecting genes to fit normal
nbin=30; % number of bins in the original histogram
frac=0.75; % fraction of the max y-bin used to select bins that are used to improve center of the distribution
part=10; % percentile of data (left and right) of center that is used in the final histogram
nbinf=9; % number of bins in the final histogram


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',20,'defaultlinelinewidth',2);

ngen=length(id); 

ra=[ir1, ir2];
na=[in1, in2];
cc0(1)=corr(ra(:,1),ra(:,2),'type','Pearson');
cc0(2)=corr(ra(:,1),ra(:,2),'type','Spearman');


%%%%% all genes
ix=find(ra(:,1)>0 & ra(:,2)>0 & ra(:,1) + ra(:,2)>=rmin & na(:,1)+na(:,2)>=nmin);  %%%%% only genes when both expressions are not zero
ix1=find((ra(:,1)==0 & ra(:,2)>=rmin & na(:,2)>=nmin) | (ra(:,2)==0 & ra(:,1)>=rmin & na(:,1)>=nmin)); %%%%% one of two genes has zero read coverage
ix2=[ix;ix1];
ix0=find(ra(:,1) + ra(:,2)<rmin | na(:,1)+na(:,2)<nmin); %%%%% both genes have zero coverage

% global ratio
C12=median(ra(ix,1)./na(ix,1)./(ra(ix,2)./na(ix,2)));

ixup=double.empty; ixdn=double.empty; pvup=double.empty; pvdn=double.empty; fdrup=double.empty; fdrdn=double.empty;

x0=-15:0.0005:15;

ixup=double.empty; ixdn=double.empty; aixup=double.empty; aixdn=double.empty; apvup=double.empty; apvdn=double.empty;


gc_suffix='';
if corr1 
% correcting on GC
gc_suffix='_with_GC_correction';
var=extpar(:,2);

%{
figure;box;
cdfplot(var(ix));
xlabel('GC fraction in gene');
ylabel('Fraction of (coding) genes');
title('Ensembl protein coding genes');
xlim([0.3,0.75]);
grid on;
set(gcf,'position',[0 0 800 350]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'CDF_GC');
print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
%}

p=prctile(var(ix),0:100/nint:101);p(1)=min(var)-1;p(nint+1)=max(var);


% LOOP ***** 1st correction GC **********
for k=2:nint+1
ixs=ix(find(var(ix)>p(k-1) & var(ix)<=p(k)));
ixsa=find(var>p(k-1) & var<=p(k));
z=log2(ra(ixs,1)./ra(ixs,2));
b(ixs)=k-1;
zz(ixs)=z;

pf=prctile(z,[fitp,100-fitp]);
[y,x]=hist(z(find(z>=pf(1) & z<=pf(2))),nbin);del=x(2)-x(1);

maxi=min(find(y==max(y))); maxx=x(maxi);  maxy=y(maxi);

% get stable max
fy=find(y>frac*maxy);
zcenter=z(find(z>=min(x(fy))-del/2 & z<=max(x(fy))+del/2));
xcenter=median(zcenter);

xmin=prctile(z(find(z<=xcenter)),[part]); xmax=prctile(z(find(z>=xcenter)),[100-part]);
% in case of assymetric distribution
xleft=xcenter-xmin; xright=xmax-xcenter;
spread=min(xleft,xright);
xmax=xcenter+spread;xmin=xcenter-spread;

% alternative for spread
if xleft < xright
fs=find(z>=xcenter & z<=pf(2));
xmax=prctile(z(fs),[min(100,max(0,100*length(find(z>=xmin & z<=xcenter))/length(fs)))]);
else
fs=find(z>=pf(1) & z<=xcenter);
xmin=prctile(z(fs),[min(100,max(0,100*(1-length(find(z>=xcenter & z<=xmax))/length(fs))))]);
end

% re-bin
%[yf,xf]=hist(z(find(z>=xcenter-spread & z<=xcenter+spread)),nbinf);
[yf,xf]=hist(z(find(z>=xmin & z<=xmax)),nbinf); delf=xf(2)-xf(1);

a0(1)=xcenter;
a0(2)=max(xmax-xcenter,xmin-xleft);
a0(3)=max(yf)*delf/(normcdf(xcenter+delf/2,a0(1),a0(2))-normcdf(xcenter-delf/2,a0(1),a0(2)))*2;
[a1,chi2] = fminsearch(@(a) myfnorm(xf,yf/delf,a),a0,optimset('MaxIter',100000,'MaxFunEvals',100000,'TolX',1e-16));

n=abs(a1(3));m=a1(1);s=abs(a1(2));

%{
figure;box;hold;
stairs([x-del/2,x(nbin)+del/2],[y/del,y(nbin)/del]);
stairs([xf-delf/2,xf(nbinf)+delf/2],[yf/delf,yf(nbinf)/delf],'k');
plot(xf,yf/delf,'ok','MarkerFaceColor','r');
plot(x0,n*normpdf(x0,m,s),'r.','MarkerSize',2);
xlim([m-8*s,m+8*s]);
title(strcat('Correction 1 (',num2str(k-1),')'));
%}

chi2GC(k-1)=chi2; aveGC(k-1)=m; sigmaGC(k-1)=s;
%c=2^(m/2); ra(ixsa,1)=ra(ixsa,1)/c; ra(ixsa,2)=ra(ixsa,2)*c; % NORMALIZATION on GC 

end

if figs
figure;box;
boxplot(zz,b);
xlim([1.5,11.5]);
xlabel('Quantiles of GC content in gene');
ylabel('log2(RPKM1/RPKM2)');
title(strcat(name1,{' '},'vs',{' '},name2));
ylim([-3,3]);
l=line([0.5,11.5],[0,0]); set(l,'LineWidth',1,'Color','black');
set(gcf,'position',[0 0 800 350]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Boxplot_GC_',suffix);
print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
end

xc=prctile(var(ix),100/nint/2:100/nint:100);
%{
figure;hold;
plot(xc,chi2GC,'o-r'); 
xlabel('Quantile of GC');
ylabel('Chi2');
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Chi2_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));

figure;hold;
plot(xc,sigmaGC,'o-');
xlabel('Quantile of GC');
ylabel('Variation')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Sigma_Variation_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
%}

% fit with a polynomial

npol=3;
P = polyfit(xc,aveGC,npol);

if figs
figure;hold;
plot(xc,aveGC,'o-');
yc=P(npol+1);
for i=1:npol
yc=yc+xc.^(npol+1-i)*P(i);
end
end
%{
plot(xc,yc,'s-k');
xlabel('Quantile of GC');
ylabel('Average')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Average_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arg=P(npol+1);
for i=1:npol
arg=arg+var.^(npol+1-i)*P(i);
end
c=2.^(arg/2); ra(:,1)=ra(:,1)./c; ra(:,2)=ra(:,2).*c; % NORMALIZATION on GC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

L_suffix='';
if corr2
L_suffix='_with_L_correction';
% correcting on L
var=extpar(:,1);
nint=5;
p=prctile(var(ix),0:100/nint:101);p(1)=min(var)-1;p(nint+1)=max(var);

% LOOP ***************************************
for k=2:nint+1

ixs=ix(find(var(ix)>p(k-1) & var(ix)<=p(k)));
ixsa=find(var>p(k-1) & var<=p(k));

z=log2(ra(ixs,1)./ra(ixs,2));
b(ixs)=k-1;
zz(ixs)=z;

pf=prctile(z,[fitp,100-fitp]);
[y,x]=hist(z(find(z>=pf(1) & z<=pf(2))),nbin); del=x(2)-x(1);

maxi=min(find(y==max(y))); maxx=x(maxi);  maxy=y(maxi);

% get stable max
fy=find(y>frac*maxy);
zcenter=z(find(z>=min(x(fy))-del/2 & z<=max(x(fy))+del/2));
xcenter=median(zcenter);

xmin=prctile(z(find(z<=xcenter)),[part]); xmax=prctile(z(find(z>=xcenter)),[100-part]);
% in case of assymetric distribution
xleft=xcenter-xmin; xright=xmax-xcenter;
spread=min(xleft,xright);
xmax=xcenter+spread;xmin=xcenter-spread;

% alternative for spread
if xleft < xright
fs=find(z>=xcenter & z<=pf(2));
xmax=prctile(z(fs),[min(100,max(0,100*length(find(z>=xmin & z<=xcenter))/length(fs)))]);
else
fs=find(z>=pf(1) & z<=xcenter);
xmin=prctile(z(fs),[min(100,max(0,100*(1-length(find(z>=xcenter & z<=xmax))/length(fs))))]);
end

% re-bin
%[yf,xf]=hist(z(find(z>=xcenter-spread & z<=xcenter+spread)),nbinf);
[yf,xf]=hist(z(find(z>=xmin & z<=xmax)),nbinf); delf=xf(2)-xf(1);

a0(1)=xcenter;a0(2)=max(xmax-xcenter,xmin-xleft);a0(3)=max(yf)*delf/(normcdf(xcenter+delf/2,a0(1),a0(2))-normcdf(xcenter-delf/2,a0(1),a0(2)))*2;

[a1,chi2] = fminsearch(@(a) myfnorm(xf,yf/delf,a),a0,optimset('MaxIter',100000,'MaxFunEvals',100000,'TolX',1e-16));

n=abs(a1(3)); m=a1(1); s=abs(a1(2));

%{
figure;box;hold;
stairs([x-del/2,x(nbin)+del/2],[y/del,y(nbin)/del]);
stairs([xf-delf/2,xf(nbinf)+delf/2],[yf/delf,yf(nbinf)/delf],'k');
plot(xf,yf/delf,'ok','MarkerFaceColor','r');
plot(x0,n*normpdf(x0,m,s),'r.','MarkerSize',2);
xlim([m-8*s,m+8*s]);
title(strcat('Correction 2 (',num2str(k-1),')'));
%}

chi2L(k-1)=chi2; aveL(k-1)=m; sigmaL(k-1)=s;
c=2^(m/2); ra(ixsa,1)=ra(ixsa,1)/c; ra(ixsa,2)=ra(ixsa,2)*c;

end

if figs
figure;box;
boxplot(zz,b);
xlim([1.5,6.5]);
xlabel('Quantiles of coding gene length');
ylabel('log2(RPKM1/RPKM2)');
title(strcat(name1,{' '},'vs',{' '},name2));
ylim([-3,3]);
l=line([0.5,11.5],[0,0]); set(l,'LineWidth',1,'Color','black');
set(gcf,'position',[0 0 800 350]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Boxplot_L_',suffix,gc_suffix);
print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
end

xc=prctile(var(ix),100/nint/2:100/nint:100);

%{
figure;hold;
plot(xc,chi2L,'o-r'); 
xlabel('Quantile of L');
ylabel('Chi2')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Chi2_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));

figure;hold;
plot(xc,sigmaL,'o-');
xlabel('Quantile of L');
ylabel('Variation')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Sigma_Variation_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
%}

% fit with a line
%{
npol=3;
P = polyfit(xc,aveL,npol);
%}
if figs
figure;hold;
plot(xc,aveL,'o-');
end
%{
yc=P(npol+1);
for i=1:npol
yc=yc+xc.^(npol+1-i)*P(i);
end
plot(xc,yc,'s-k');
%}
%{
xlabel('Quantile of L');
ylabel('Average')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Average_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
arg=P(npol+1);
for i=1:npol
arg=arg+var.^(npol+1-i)*P(i);
end
c=2.^(arg/2); ra(:,1)=ra(:,1)./c; ra(:,2)=ra(:,2).*c; % NORMALIZATION on Length
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% Doing the job
nint=20;
%var=na(:,1).*na(:,2)./(na(:,1)+na(:,2)+1)./extpar(:,1);
%var=na(:,1).*na(:,2)./(na(:,1)+na(:,2)+1);
%var=(na(:,1)*C12+na(:,2))/2;
%var=sqrt(na(:,1).*na(:,2));
%var=max(na(:,1)*C12,na(:,2));
var=(ra(:,1)+ra(:,2))/2;
%var=sqrt(ra(:,1).*ra(:,2));
p=prctile(var(ix),0:100/nint:101);p(1)=min(var)-1;p(nint+1)=max(var);

% LOOP ***************************************
a1(1)=0;
for k=2:nint+1

ixs=ix(find(var(ix)>p(k-1) & var(ix)<=p(k)));
ixs1=ix1(find(var(ix1)>p(k-1) & var(ix1)<=p(k)));
ixsa=find(var>p(k-1) & var<=p(k));

z=log2(ra(ixs,1)./ra(ixs,2));

pf=prctile(z,[fitp,100-fitp]);
[y,x]=hist(z(find(z>=pf(1) & z<=pf(2))),nbin);
del=x(2)-x(1);
if figs
figure;box;hold;
stairs([x-del/2,x(nbin)+del/2],[y/del,y(nbin)/del]);
end

maxi=min(find(y==max(y))); maxx=x(maxi);  maxy=y(maxi);

% get stable max
fy=find(y>frac*maxy);
zcenter=z(find(z>=min(x(fy))-del/2 & z<=max(x(fy))+del/2));
xcenter=median(zcenter);

xmin=prctile(z(find(z<=xcenter)),[part]); xmax=prctile(z(find(z>=xcenter)),[100-part]);
% in case of assymetric distribution
xleft=xcenter-xmin; xright=xmax-xcenter;
spread=min(xleft,xright);
%xmax=xcenter+spread;xmin=xcenter-spread;

% alternative for spread
if xleft < xright
fs=find(z>=xcenter & z<=pf(2));
xmax=prctile(z(fs),[min(100,max(0,100*length(find(z>=xmin & z<=xcenter))/length(fs)))]);
else
fs=find(z>=pf(1) & z<=xcenter);
xmin=prctile(z(fs),[min(100,max(0,100*(1-length(find(z>=xcenter & z<=xmax))/length(fs))))]);
end

% re-bin
%[yf,xf]=hist(z(find(z>=xcenter-spread & z<=xcenter+spread)),nbinf);
[yf,xf]=hist(z(find(z>=xmin & z<=xmax)),nbinf);
delf=xf(2)-xf(1);
if figs
stairs([xf-delf/2,xf(nbinf)+delf/2],[yf/delf,yf(nbinf)/delf],'k');
plot(xf,yf/delf,'ok','MarkerFaceColor','r');
end

a0(1)=xcenter;a0(2)=max(xmax-xcenter,xmin-xleft);a0(3)=max(yf)*delf/(normcdf(xcenter+delf/2,a0(1),a0(2))-normcdf(xcenter-delf/2,a0(1),a0(2)))*2;
[a1,chi2] = fminsearch(@(a) myfnorm(xf,yf/delf,a),a0,optimset('MaxIter',100000,'MaxFunEvals',100000,'TolX',1e-16));

n=abs(a1(3)); m=a1(1); s=abs(a1(2));

if figs
plot(x0,n*normpdf(x0,m,s),'r.','MarkerSize',2);
xlim([m-8*s,m+8*s]);
title(strcat('Job (',num2str(k-1),')'));
end


%%%%%%%%%%%%%%%

chi2fin(k-1)=chi2; ave(k-1)=m; sigma(k-1)=s;
%{
c=2^(m/2); rca(ixsa,1)=ra(ixsa,1)/c; rca(ixsa,2)=ra(ixsa,2)*c;  % JOB normalization

%%%%%%%%%%%%%%%% adding those genes that have zero reads in one library
z1=log2(C12*max(na(ixs1,1),1)./max(na(ixs1,2),1));
z2=[z; z1]; 
ixs2=[ixs;ixs1];
%}
%{
% assign p-value for every z-value
dn=find(z2<=m); up=find(z2>m);
pvdn=  normcdf(z2(dn),m,s);  
pvup=(1-normcdf(z2(up),m,s));

apvup=[apvup; pvup];
apvdn=[apvdn; pvdn];

aixdn=[aixdn; ixs2(dn)];
aixup=[aixup; ixs2(up)];

% apply FDR-threshold
[dum,sixup]=sort(pvup);
fixup=up(find(pvup<fdrth*sixup/length(up)));
ixup=[ixup; ixs2(fixup)];

[dum,sixdn]=sort(pvdn); 
fixdn=dn(find(pvdn<fdrth*sixdn/length(dn)));
ixdn=[ixdn; ixs2(fixdn)];
%}

end
% LOOP ***************************************
%{
figure;hold;
plot(chi2fin,'o-r'); 
xlabel('Quantile of variable');
ylabel('Chi2')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Chi2_',suffix,'_',num2str(k-1));
print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
%}

% fit both with 3d 
% fit with a polynomial

xc=100/nint/2:100/nint:100;

% fit of 'sigma' with 3d degree polynomial
npolS=4;
PS = polyfit(xc,sigma,npolS);
yc=PS(npolS+1);
for i=1:npolS
yc=yc+xc.^(npolS+1-i)*PS(i);
end

if figs
figure;hold;
plot(xc,sigma,'o-');
plot(xc,yc,'xk-');
xlabel('Quantile of variable');
ylabel('Variation')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Sigma_Variation_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
end

% fit of 'average' with 3d degree polynomial
npolA=4;
PA = polyfit(xc,ave,npolA);
yc=PA(npolA+1);
for i=1:npolA
yc=yc+xc.^(npolA+1-i)*PA(i);
end

if figs
figure;hold;
plot(xc,ave,'o-');
plot(xc,yc,'xk-');
xlabel('Quantile of variable');
ylabel('Average')';
title(strcat(name1,'-vs- ',name2));
set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Average_',suffix,'_',num2str(k-1));
%print(gcf, '-dpdf',strcat(fileOut,'.pdf'));
end

% Introducing correction
aixup=double.empty; aixdn=double.empty; apvup=double.empty; apvdn=double.empty;
bound(1)=-1;
for i=1:1:100

A=PA(npolA+1);
for ii=1:npolA
A=A+i.^(npolA+1-ii)*PA(ii);
end

S=PS(npolS+1);
for ii=1:npolS
S=S+i.^(npolS+1-ii)*PS(ii);
end

bound(i+1)=prctile(var(ix),[i]);
ixs=ix(find(var(ix)>bound(i) & var(ix)<=bound(i+1)));
ixs1=ix1(find(var(ix1)>bound(i) & var(ix1)<=bound(i+1)));
ixsa=find(var>bound(i) & var<=bound(i+1));

z=log2(ra(ixs,1)./ra(ixs,2));

c=2^(A/2);
rca(ixsa,1)=ra(ixsa,1)/c; rca(ixsa,2)=ra(ixsa,2)*c;

% p-value

z1=log2(C12*max(na(ixs1,1),1)./max(na(ixs1,2),1));
z2=[z; z1]; 
ixs2=[ixs;ixs1];

dn=find(z2<=A); up=find(z2>A);
% note factor 2
pvdn=  normcdf(z2(dn),A,S)*2;  
pvup=(1-normcdf(z2(up),A,S))*2;

apvup=[apvup; pvup];
apvdn=[apvdn; pvdn];

aixdn=[aixdn; ixs2(dn)];
aixup=[aixup; ixs2(up)];

end

cc(1)=corr(rca(:,1),rca(:,2),'type','Pearson'); 
cc(2)=corr(rca(:,1),rca(:,2),'type','Spearman');

%%%%%%%%%%%%%%%%%%%%%%%%%
% new p-value calculation
%%%%%%%%%%%%%%%%%%%%%%%%%

%nfup=length(ixup);
%nfdn=length(ixdn);

clear ixup, ixdn;

% apply FDR-threshold globaly

LUP=length(aixup);
LDN=length(aixdn);

[dum,sixup]=sort(apvup);
for i=1:LUP
fdrup(sixup(i))=apvup(sixup(i))*LUP/i;
end
fixup=find(fdrup<fdrth);
ixup=aixup(fixup);


[dum,sixdn]=sort(apvdn); 
for i=1:LDN
fdrdn(sixdn(i))=apvdn(sixdn(i))*LDN/i;
end
fixdn=find(fdrdn<fdrth);
ixdn=aixdn(fixdn);

nfup=length(ixup);
nfdn=length(ixdn);

eps=0.0001;
if out
% plot p-value/FDR
%%%%%%%%%%%%%%%%%

maxv=max(max(ir1),max(ir2))*2;

col1=[0.85 0.85 0.85]; col2=[0.65 0.65 0.65];
colup1=[1 0.7 0.7]; coldn1=[0.35 0.7 1];
colup2=[0.8 0.4 0.4]; coldn2=[0.2 0.4 0.8];
colup3=[0.7 0.2 0.2]; coldn3=[0.1 0.2 0.7];
colup4=[1 0 0]; coldn4=[0 0 1];


if figs
figure;
set(0,'defaultaxesfontsize',16,'defaultlinelinewidth',1);
sh=scatterhist(log10(max(eps,rca(ix,2))),log10(max(eps,rca(ix,1))),'NBins',[50 50],'Direction','out');
so=get(sh(1),'children'); % handle for plot inside scaterplot axes
%set(so,'Marker','.','MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','w');
objs=findobj('type','patch');
set(objs(:),'FaceColor',[0.3 0.3 0.6],'EdgeColor',[0.5 0.5 0.5]);

hold;
plot(log10(max(eps,rca(:,2))),log10(max(eps,rca(:,1))),'s','MarkerSize',6,'MarkerEdgeColor',col1,'MarkerFaceColor',col1);
plot(log10(max(eps,rca(:,2))),log10(max(eps,rca(:,1))),'s','MarkerSize',2,'MarkerEdgeColor',col2,'MarkerFaceColor',col2);

plot(log10(max(eps,rca(ixup,2))),log10(max(eps,rca(ixup,1))),'s','MarkerSize',6,'MarkerEdgeColor',colup1,'MarkerFaceColor',colup1);
plot(log10(max(eps,rca(ixdn,2))),log10(max(eps,rca(ixdn,1))),'s','MarkerSize',6,'MarkerEdgeColor',coldn1,'MarkerFaceColor',coldn1);

plot(log10(max(eps,rca(ixup,2))),log10(max(eps,rca(ixup,1))),'s','MarkerSize',2,'MarkerEdgeColor',colup3,'MarkerFaceColor',colup3);
plot(log10(max(eps,rca(ixdn,2))),log10(max(eps,rca(ixdn,1))),'s','MarkerSize',2,'MarkerEdgeColor',coldn3,'MarkerFaceColor',coldn3);

plot(log10(max(eps,rca(:,2))),log10(max(eps,rca(:,1))),'.k','MarkerSize',2);


ll=log10(eps);ul=log10(5000);
xlim([ll,ul]);ylim([ll,ul]);diagl=line([ll,ul],[ll,ul]);set(diagl,'Color','green','LineWidth',1);

xlabel(strcat('log10( RPKM[',name2,'] )'));
ylabel(strcat('log10( RPKM[',name1,'] )'));

title( ['Total: ', num2str(length(id)), ' (R=', num2str(round(cc(1)*1000+0.5)/1000), ' RC=', num2str(round(cc(2)*1000+0.5)/1000),')',' FDR=',num2str(fdrth)]);
%,' (x=', num2str(find(ir2==0))) );
%,', y=',num2str(find(ir1==0)),', in=',num2str(find(ir1>0 & ir2>0)),')'));

% for the legend only

ax=axis();dx=(ax(2)-ax(1));dy=(ax(4)-ax(3));

hpd=plot(ax(1)-1,ax(2)+1,'s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',coldn1);
hpu=plot(ax(1)-1,ax(2)+1,'s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',colup1);
legend([hpd,hpu],strcat('DN:',num2str(nfdn)),strcat('UP:',num2str(nfup)),'Location','NorthWest');
legend('boxoff');

%{
yt=maxv/9;
tt=text(xt,yt,strcat('RPKM>',num2str(er),' and N-reads>',num2str(el)));
set(tt,'FontSize',20);

yt=maxv/27;
fc=(acut+1)/(1-acut);
%tt=text(xt,yt,strcat('A-cut: ',num2str(acutd),',',num2str(acutu)));
tt=text(xt,yt,strcat('Fold change: ',num2str(fc)));

set(gca,'XTick',[0.01 1 100 10000],'YTick',[0.01 1 100 10000]);
%}

set(gcf,'position',[0 0 800 750]);set(gcf, 'PaperPositionMode', 'auto');
fileOut=strcat(dirOut,'Plot2D_',suffix,gc_suffix,L_suffix);
print(gcf, '-dpdf',strcat(fileOut,'.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MA plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvar=log10(var);
yvar=log2(max(eps,rca(:,1))./max(eps,rca(:,2)));
figure;box;hold;
plot(xvar,yvar,'s','MarkerSize',6,'MarkerEdgeColor',col1,'MarkerFaceColor',col1);
plot(xvar,yvar,'s','MarkerSize',2,'MarkerEdgeColor',col2,'MarkerFaceColor',col2);

plot(xvar(ixup),log2(max(eps,rca(ixup,1))./max(eps,rca(ixup,2))),'s','MarkerSize',6,'MarkerEdgeColor',colup1,'MarkerFaceColor',colup1);
plot(xvar(ixup),log2(max(eps,rca(ixup,1))./max(eps,rca(ixup,2))),'s','MarkerSize',2,'MarkerEdgeColor',colup2,'MarkerFaceColor',colup2);
plot(xvar(ixdn),log2(max(eps,rca(ixdn,1))./max(eps,rca(ixdn,2))),'s','MarkerSize',6,'MarkerEdgeColor',coldn1,'MarkerFaceColor',coldn1);
plot(xvar(ixdn),log2(max(eps,rca(ixdn,1))./max(eps,rca(ixdn,2))),'s','MarkerSize',2,'MarkerEdgeColor',coldn2,'MarkerFaceColor',coldn2);

plot(xvar,log2(max(eps,rca(:,1))./max(eps,rca(:,2))),'.k','MarkerSize',2);

ymm=prctile(yvar,[0.1,99.9]); ymma=max(abs(ymm(1)),abs(ymm(2)));
xlim([-4,3.5]);
ylim([-ymma,ymma])
for k=2:nint
qul=line([log10(p(k)),log10(p(k))],[-ymma,ymma]);
set(qul,'Color','black','LineWidth',1);
end
bl=line([-4,3.5],[0,0]);
set(bl,'Color','green','LineWidth',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MA plot 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure;box;hold;

plot(log2(max(eps,rca(:,1)).*max(eps,rca(:,2))),log2(max(eps,rca(:,1))./max(eps,rca(:,2))),'s','MarkerSize',6,'MarkerEdgeColor',col1,'MarkerFaceColor',col1)
plot(log2(max(eps,rca(:,1)).*max(eps,rca(:,2))),log2(max(eps,rca(:,1))./max(eps,rca(:,2))),'s','MarkerSize',2,'MarkerEdgeColor',col2,'MarkerFaceColor',col2)

plot(log2(max(eps,rca(ix2(ixup),1)).*max(eps,rca(ix2(ixup),2))),log2(max(eps,rca(ix2(ixup),1))./max(eps,rca(ix2(ixup),2))),'s','MarkerSize',6,'MarkerEdgeColor',colup1,'MarkerFaceColor',colup1)
plot(log2(max(eps,rca(ix2(ixup),1)).*max(eps,rca(ix2(ixup),2))),log2(max(eps,rca(ix2(ixup),1))./max(eps,rca(ix2(ixup),2))),'s','MarkerSize',2,'MarkerEdgeColor',colup2,'MarkerFaceColor',colup2)
plot(log2(max(eps,rca(ix2(ixdn),1)).*max(eps,rca(ix2(ixdn),2))),log2(max(eps,rca(ix2(ixdn),1))./max(eps,rca(ix2(ixdn),2))),'s','MarkerSize',6,'MarkerEdgeColor',coldn1,'MarkerFaceColor',coldn1)
plot(log2(max(eps,rca(ix2(ixdn),1)).*max(eps,rca(ix2(ixdn),2))),log2(max(eps,rca(ix2(ixdn),1))./max(eps,rca(ix2(ixdn),2))),'s','MarkerSize',2,'MarkerEdgeColor',coldn2,'MarkerFaceColor',coldn2)

plot(log2(max(eps,rca(:,1)).*max(eps,rca(:,2))),log2(max(eps,rca(:,1))./max(eps,rca(:,2))),'.k','MarkerSize',2)

%xlim([p(1),p(11)]);
ylim([-may,may]);
bl=line([p(1),p(11)],[0,0]);
set(bl,'Color','black','LineWidth',1);
%}

% report genes
%%%%%%%%%%%%%%%%%
fileOut = fopen(strcat(dirOut,'UP.',suffix,'.FDR_',num2str(fdrth),'.rmin_',num2str(rmin),'.Nmin_',num2str(nmin)),'w');
for i=1:length(ixup)
fprintf(fileOut,'%s\t%f\t%f\t%f\t%f\n',id{ixup(i)},rca(ixup(i),:),apvup(fixup(i)),fdrup(fixup(i)));
end
fclose(fileOut);
fileOut = fopen(strcat(dirOut,'DN.',suffix,'.FDR_',num2str(fdrth),'.rmin_',num2str(rmin),'.Nmin_',num2str(nmin)),'w');
for i=1:length(ixdn)
fprintf(fileOut,'%s\t%f\t%f\t%f\t%f\n',id{ixdn(i)},rca(ixdn(i),:),apvdn(fixdn(i)),fdrdn(fixdn(i)));
end
fclose(fileOut);

end

% output re-normalized RPKM values
if rpkm
fileOut = fopen(strcat(dirOut,name1,'_',name2,'.RPKM.corrected'),'w');
for i=1:length(id)

[C,AU,BU]=intersect(i,aixup);
[C,AD,BD]=intersect(i,aixdn);
if AU==1
  t1u=apvup(BU)/2; t1d=1-t1u; % uncorrected 
  t2u=fdrup(BU)/2; t2d=1-t2u; % corrected
elseif AD==1
  t1d=apvdn(BD)/2; t1u=1-t1d;
  t2d=fdrdn(BD)/2; t2u=1-t2d;
else
  t1u=1;t1d=1;
  t2u=1;t2d=1;
end
cn=sqrt(rca(i,1)/rca(i,2)*na(i,2)/na(i,1));
fprintf(fileOut,'%s\t%i\t%f\t%i\t%f\t%f\t%f\t%f\t%f\n',id{i},na(i,1)*cn,rca(i,1),na(i,2)/cn,rca(i,2),t1u,t1d,t2u,t2d);
end
fclose(fileOut);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of DIFFERENTIAL EXPRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
