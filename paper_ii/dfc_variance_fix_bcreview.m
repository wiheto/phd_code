clear all

%% This scrip includes the all the code that was used in the final submission to brain connectivity 
%% THe one exception is the script dfc_variance_lambda.m, which estimates the lambda parameter for multiple windows. 


figdir='/data/william/Projects/dfc_contrasts/figures_bcresubmit/'

%% Simulation of normal distributions around normal distributions

rng(1)
ts=[]
clear rkeep cvkeep 
for n=1:1000
    cv=normrnd(.5,.1,1,1)
    R=mvnrnd([0.6 0.6],[1 cv; cv 1],100);
    cvkeep(n)=cv; 
    ts(end+1:end+100,1:2)=R;
    rraw=corr(R);
    rkeep(n)=rraw(1,2);
end

zrkeep=fisherz(rkeep)

lamba=boxcoxlm(rkeep',[ones(length(rkeep),1)],0,-5:0.01:5);
if lamba==0
    bcsim=log(rkeep);
else
    bcsim=bsxfun(@rdivide,(bsxfun(@power,rkeep,lamba)-1),lamba');
end

lamba=boxcoxlm(zrkeep,[ones(length(zrkeep),1)],0,-5:0.01:5);
if lamba==0
    zbcsim=log(zrkeep);
else
    zbcsim=bsxfun(@rdivide,(bsxfun(@power,zrkeep,lamba)-1),lamba');
end

figure
subplot(6,3,[1 4])
plot(ts(1:400,:))
axis square
subplot(6,3,[2 5]) 
hold on
hist(cvkeep-mean(cvkeep),-1:0.01:1,'k')
axis square
axis([-1 1 0 60])
subplot(6,3,[3 6]) 
hist(rkeep-mean(rkeep),-1:0.01:1)
axis square
axis([-1 1 0 60])
subplot(6,3,[7 10]) 
hist(zrkeep-mean(zrkeep),-1:0.01:1)
axis square
axis([-1 1 0 60])
subplot(6,3,[8 11]) 
hist(bcsim-mean(bcsim),-1:0.01:1)
axis square
axis([-1 1 0 60])
subplot(6,3,[9 12]) 
hist(zbcsim-mean(zbcsim),-1:0.01:1)
axis square
axis([-1 1 0 60])

skcv=skewness(cvkeep);
skr=skewness(rkeep);
skz=skewness(zrkeep);
skbc=skewness(bcsim);
skzbc=skewness(zbcsim);

vcv=var(cvkeep);
vr=var(rkeep);
vz=var(zrkeep);
vbc=var(bcsim);
vzbc=var(zbcsim);

[~, pcv, ~, swcv]=swtest_forced(cvkeep);
[~, psr, ~, swr]=swtest_forced(rkeep);
[~, pz, ~, swz]=swtest_forced(zrkeep);
[~, pbc, ~, swbc]=swtest_forced(bcsim);
[~, pzbc, ~, swzbc]=swtest_forced(zbcsim);


subplot(6,3,[13 16])
bar([skcv,skr,skz,skbc,skzbc])
axis square
axis([0.5 5.5 -.5 .5])         

subplot(6,3,[14 17])
bar([swcv,swr,swz,swbc,swzbc])
axis([0.5 5.5 -5 5])  
axis square


subplot(6,3,[15 18])
bar([vcv,vr,vz,vbc,vzbc])
axis([0.5 5.5 0 0.05])  
axis square

mkdir /data/william/Projects/dfc_contrasts/figures_bcresubmit
print(gcf,'/data/william/Projects/dfc_contrasts/figures_bcresubmit/simV2','-depsc','-r300')


%% EMPIRCAL PART


%Directory paramete
ldir = '/data/william/Projects/eoec_dfctime/' %load directory (data)
mdir = '/data/william/Projects/dfc_contrasts/' %workding directory 
figdir = [mdir 'figures_bcresubmit/'] %figures directory
mkdir([mdir 'figures_bcresubmit/'])
%Network information (node number and assigned networks)
cd /data/william/Projects/dfc_time/
powernetclassclusters
cd(mdir)

%Data parameters
% Allsubject approach
S = 48; %NumberOfSubjects
L = 240 %DataLength
C = 1   %Condition
sid=0
c=1
clear vox voxtmp
for s=1:S
    if s==8 || s== 22 %<- subjects that were excluded due to imcomplete data
    else
    s
    sid=sid+1;
    load([ldir 'timeseries_2016/eo_' num2str(s) '.mat'])
    dat = eo'; 
    for roi=1:264
         voxtmp(roi,:)=dat(roi,:)-mean(dat(roi,:));
    end
    tmp = voxtmp(:,:);
    tmp1 = reshape(tmp,size(tmp,1)*size(tmp,2),1); 
    tmp=(tmp-mean(tmp1))./std(tmp1);
    vox(:,:,c,sid)=tmp; 
    end
end


% Create static connectivity matrix
clear RStatic
for s=1:S-2
    for c=1:C
    RStatic(:,:,c,s)=corr(vox(:,:,c,s)');
    end
end
RStatic=reshape(RStatic,264*264,C,S-2);
RStatic=squeeze(RStatic(:,:,:,:)); 
ZRStatic=fisherz(RStatic); 


%Do calculate sliding window.
sw=31
s=1
rsw=[]; 
rsw=zeros(264,264,length(sw+1:L-sw),1,S-2);
for s=1:S-2
    s
tid=0
for n=sw+1:L-sw
    tid=tid+1; 
    for c=1:C
        rsw(1:264,1:264,tid,c,s)=corr(vox(:,n-sw:n+sw,c,s)');
    end
end
end

rsw=reshape(rsw,264*264,size(rsw,3),C,S-2);
rsw=squeeze(rsw(:,:,:)); 

zrsw=fisherz(rsw); 
zrsw=reshape(zrsw,size(rsw)); 


%% Load lambda, precalculate
slide=1
load(['/data/william/Projects/dfc_contrasts/boxcoxDMN_all_sw' num2str(sw) '_slide' num2str(slide)])

%% Show lambda per subject (reviewer request)

figure; 
subplot(2,2,1:2)
errorbar(1:46,mean(l),std(l),'kx')
subplot(2,2,3:4)
errorbar(1:46,mean(lz),std(lz),'kx')
print([figdir 'lambda_subjects'],'-depsc','-r300')

%% Plot distributions of lambda

figure; 
subplot(1,2,1) 
a=hist(reshape(squeeze(l),size(l,1)*46,1),-5:1:5)
bar(-5:5,a./length(reshape(squeeze(l),size(l,1)*46,1)))
axis square
axis([-5.5 5.5 0 .3])
subplot(1,2,2)
a=hist(reshape(squeeze(lz),size(lz,1)*46,1),-5:1:5)
bar(-5:5,a./length(reshape(squeeze(lz),size(lz,1)*46,1)))
axis square
axis([-5.5 5.5 0 .3])
print([figdir 'lambdadistribution'],'-depsc','-r300')

%% Transform the data

A=zeros(264);
ind=find(~tril(ones(264,264)));

rsw=rsw(ind,:,:);
zrsw=zrsw(ind,:,:);
RStatic=RStatic(ind,:,:);

%Preallocate
tsBC=zeros(size(rsw));
ztsBC=zeros(size(rsw));
clear tsBCtmp
for s=1:S-2
    %Preallocate
    tsBCtmp=zeros(size(rsw,1),size(rsw,2),C)*NaN; 
    %Scale data so that the minimum value is equal to 1
    dattmp=bsxfun(@plus,rsw(:,:,s),1-(min(rsw(:,:,s),[],2)));
    %perform boxcox
    lambda = l(:,s)'; 
    tsBCtmp(lambda~=0,:)=bsxfun(@rdivide,(bsxfun(@power,dattmp(lambda~=0,:),lambda(lambda~=0)')-1),lambda(lambda~=0)');
    tsBCtmp(lambda==0,:)=log(dattmp(lambda==0,:));
    %rescale so that the mean of the transformed data == mean of the
    %original data
    tsBC(:,:,s)=bsxfun(@minus,tsBCtmp,mean(tsBCtmp,2)-mean(rsw(:,:,s),2));

    
    %Preallocate
    tsBCtmp=zeros(size(zrsw,1),size(zrsw,2),C)*NaN; 
    %Scale data so that the minimum value is equal to 1
    dattmp=bsxfun(@plus,zrsw(:,:,s),1-(min(zrsw(:,:,s),[],2)));
    %perform boxcox
    lambda = lz(:,s)'; 
    tsBCtmp(lambda~=0,:)=bsxfun(@rdivide,(bsxfun(@power,dattmp(lambda~=0,:),lambda(lambda~=0)')-1),lambda(lambda~=0)');
    tsBCtmp(lambda==0,:)=log(dattmp(lambda==0,:));
    %rescale so that the mean of the transformed data == mean of the
    %original data
    ztsBC(:,:,s)=bsxfun(@minus,tsBCtmp,mean(tsBCtmp,2)-mean(zrsw(:,:,s),2));
    
end



%% Create supplementary figure correlation of mean of transforms vs sFC
d1=mean(abs(RStatic),2);
d2=mean(mean(abs(rsw),2),3);
figure; 

subplot(2,2,1)
scatter(d1,d2,'k.')
axis square
axis([0 1 0 1])
set(gca,'xTick',[0 0.5 1],'xTickLabel','')
set(gca,'yTick',[0 0.5 1],'yTickLabel','')
% print('/data/william/Projects/dfc_contrasts/meantransformation_vs_sFC1','-dpng','-r600')
% close all
% figure

subplot(2,2,2)
d2=mean(mean(abs(zrsw),2),3);
scatter(d1,d2,'k.')
axis square
axis([0 1 0 2])
set(gca,'xTick',[0 0.5 1],'xTickLabel','')
set(gca,'yTick',[0 0.5 1 1.5 2],'yTickLabel','')% close all

% figure

subplot(2,2,3)
d2=mean(mean(abs(tsBC),2),3);
scatter(d1,d2,'k.')
axis square
axis([0 1 0 1])
set(gca,'xTick',[0 0.5 1],'xTickLabel','')
set(gca,'yTick',[0 0.5 1],'yTickLabel','')% close all

% figure
subplot(2,2,4)
d2=mean(mean(abs(ztsBC),2),3);
scatter(d1,d2,'k.')
axis square
axis([0 1 0 2])
set(gca,'xTick',[0 0.5 1],'xTickLabel','')
set(gca,'yTick',[0 0.5 1 1.5 2],'yTickLabel','')
%Save as bitmap due to number of points 
print([figdir 'meantransformation_vs_sFC'],'-dpng','-r600')
%Create empty axis in eps format and can line up with png if nececesary
figure; 
for n=1:4
subplot(2,2,n)
scatter(1.5,1.5,'b.')
axis square
    if n==2 || n==4
        axis([0 1 0 2])
    else
        axis([0 1 0 1])
    end
end
print([figdir 'meantransformation_vs_sFC_axis2'],'-deps','-r600')


%% Show main example 

EOI=8083
s1=2
s2=7
figure
r=21
subplot(4,4,1)
plot(rsw(EOI,:,s1))
axis([1 200 0.5 1])
axis square
subplot(4,4,5)
plot(zrsw(EOI,:,s1))
axis([1 200 0 2])
axis square
subplot(4,4,9)
plot(tsBC(EOI,:,s1))
axis([1 200 0 1.5])
axis square
subplot(4,4,13)
plot(ztsBC(EOI,:,s1))
axis square
axis([1 200 -.5 3])


subplot(4,4,2)
v=rsw(EOI,:,s1);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])
axis square
subplot(4,4,6)
v=zrsw(EOI,:,s1);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])
axis square
subplot(4,4,10)
v=tsBC(EOI,:,s1);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])

axis square
subplot(4,4,14)
v=ztsBC(EOI,:,s1);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])



subplot(4,4,3)
plot(rsw(EOI,:,s2))
axis([1 200 0.5 1])
axis square
subplot(4,4,7)
plot(zrsw(EOI,:,s2))
axis([1 200 0.5 1.5])
axis square
subplot(4,4,11)
plot(tsBC(EOI,:,s2))
axis([1 200 0.5 .8])
axis square
subplot(4,4,15)
plot(ztsBC(EOI,:,s2))
axis([1 200 0.5 1])
axis square

subplot(4,4,4)
v=rsw(EOI,:,s2);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])
axis square
subplot(4,4,8)
v=zrsw(EOI,:,s2);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])
axis square
subplot(4,4,12)
v=tsBC(EOI,:,s2);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])

axis square
subplot(4,4,16)
v=ztsBC(EOI,:,s2);
histfit(v,r)
axis square
mu=mean(v); 
sigma=std(v); 
axis([mu-3*sigma mu+3*sigma 0 40])
set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
set(gca,'XTickLabel',[1 2 3 4 5])


print([figdir 'transform_examples'],'-depsc','-r300')



%% New example figure (due to revision) 


figure; 

h=hist(reshape(abs(RStatic),size(RStatic,1)*size(RStatic,2),1),0:0.025:1)
subplot(5,5,1:5)
bar(0:0.025:1,h/(size(RStatic,1)*size(RStatic,2)))
axis([0 1 0 .1])
s=2
ind=0
r=11
ylim=60
for n=0.1:0.2:0.9
    EOI=nearest(RStatic(:,10),n)
    subplot(5,5,6+ind)
    v=rsw(EOI,:,s);
    histfit(v,r)
    axis square
    mu=mean(v); 
    sigma=std(v); 
    axis([mu-3*sigma mu+3*sigma 0 ylim])
    set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
    set(gca,'XTickLabel',[1 2 3 4 5])

    subplot(5,5,11+ind)
    v=zrsw(EOI,:,s);
    histfit(v,r)
    axis square
    mu=mean(v); 
    sigma=std(v); 
    axis([mu-3*sigma mu+3*sigma 0 ylim])
    set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
    set(gca,'XTickLabel',[1 2 3 4 5])


    subplot(5,5,16+ind)
    v=tsBC(EOI,:,s);
    histfit(v,r)
    axis square
    mu=mean(v); 
    sigma=std(v); 
    axis([mu-3*sigma mu+3*sigma 0 ylim])
    set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
    set(gca,'XTickLabel',[1 2 3 4 5])


    subplot(5,5,21+ind)
    v=ztsBC(EOI,:,s);
    histfit(v,r)
    axis square
    mu=mean(v); 
    sigma=std(v); 
    axis([mu-3*sigma mu+3*sigma 0 ylim])
    set(gca,'XTick',[mu-2*sigma mu-sigma mu mu+1*sigma mu+2*sigma])
    set(gca,'XTickLabel',[1 2 3 4 5])

    ind=ind+1
end

print([figdir 'transform_overspectrum'],'-depsc','-r300')


%% Calculate SW statistic


clear p_BC p_raw p_z stat_BCz p_BCz stat_raw  stat_BCz stat_z bic_z bic_raw bic_bc stat_BC
for n=1:size(rsw,1)
    n
    for s=1:size(rsw,3) 
        [h p_BC(n,s), nonnormstat_BC(n,s), stat_BC(n,s)]=swtest_forced(tsBC(n,:,s));
        [h p_raw(n,s), nonnormstat_raw(n,s),  stat_raw(n,s)]=swtest_forced(rsw(n,:,s));
        [h p_z(n,s), nonnormstat_z(n,s),  stat_z(n,s)]=swtest_forced(zrsw(n,:,s));
        [h p_BCz(n,s), nonnormstat_BCz(n,s),  stat_BCz(n,s)]=swtest_forced(ztsBC(n,:,s));
    end
end



p_BC=reshape(p_BC,size(p_BC,1)*size(p_BC,2),1);
p_raw=reshape(p_raw,size(p_raw,1)*size(p_raw,2),1);
p_z=reshape(p_z,size(p_z,1)*size(p_z,2),1);
p_BCz=reshape(p_BCz,size(p_BCz,1)*size(p_BCz,2),1);

stat_BC=reshape(stat_BC,size(p_BC,1)*size(p_BC,2),1);
stat_raw=reshape(stat_raw,size(p_raw,1)*size(p_raw,2),1);
stat_z=reshape(stat_z,size(p_z,1)*size(p_z,2),1);
stat_BCz=reshape(stat_BCz,size(p_z,1)*size(p_z,2),1);

save('/data/william/Projects/dfc_contrasts/statsave','p_BC','p_raw','p_z','p_BCz','stat_raw','stat_BC','stat_BCz','stat_z')

load('/data/william/Projects/dfc_contrasts/statsave')


%% Non normal edges (fig5A)

ptest=[0.01 0.001 0.0001 0.00001]
clear NumSig
for p=1:length(ptest)
  
    NumSig(p,:)=[length(find(p_raw<ptest(p)))./length(p_raw) length(find(p_z<ptest(p)))./length(p_z) length(find(p_BC<ptest(p)))./length(p_BC)  length(find(p_BCz<ptest(p)))./length(p_BCz)]
  
end 
figure; 
hold on
col='krgb'
for n=1:4
   plot(NumSig(:,n),['x-' col(n)]) 
end
axis([0.5 4.5 0 1])
axis([0.5 4.5 0 1])
print([figdir 'nonnormaledges'],'-depsc','-r300')

%% percent of edges that had the lowest SW stat    

o=[]
for n=1:size(stat_raw,1)
        [~, o(end+1)]=min([stat_raw(n) stat_z(n)  stat_BC(n) stat_BCz(n) ]);
end
figure; 
a=hist(o,1:4)
bar(a./length(o))
print([figdir 'SWstat_lowest'],'-depsc','-r300')


%% Mean variance plot of ztsBC

clear tmp mtmp
for s=1:S-2
    s
tmp(:,s)=var(squeeze(ztsBC(:,:,s)),[],2);
mtmp(:,s)=mean(squeeze(ztsBC(:,:,s)),2);
end
t1=squeeze(mean(tmp,2));
t2=squeeze(mean(mtmp,2));
% t2=mean(reshape(RStatic,264*264,S-2),2);
figure; 
scatter(t2,t1,10,'ko','filled')
axis square; 
print([figdir 'mean_variance_zBC'],'-depsc','-r300')
print([figdir 'mean_variance_zBC'],'-dpng','-r300')

% [r p]=corr((t2(isnan(t1)==0)),(t1(isnan(t1)==0)),'type','Spearman')


%% Plots against static connectivity spectrum

vec=0:.025:1;
[a b]=histc(reshape(abs(RStatic),34716*46,1),vec);

clear stat statn stdstat
for n=1:length(a)
%     if sum(b==n)>0
        stat(n,1)=mean(stat_raw(b==n));
        stat(n,2)=mean(stat_z(b==n));
        stat(n,3)=mean(stat_BC(b==n));
        stat(n,4)=mean(stat_BCz(b==n));       
        
        stdstat(n,1)=std((stat_raw(b==n)));
        stdstat(n,2)=std((stat_z(b==n)));
        stdstat(n,3)=std((stat_BC(b==n)));
        stdstat(n,4)=std((stat_BCz(b==n)));

        statn(n)=sum(b==n)
%     end
end
figure; 
hold on
plot(vec,stat(:,1),'k'); 
plot(vec,stat(:,2),'r');
plot(vec,stat(:,3),'g');
plot(vec,stat(:,4),'b');
axis([0 1 0 6])
% shadedErrorBar(vec,stat(:,1),stdstat(:,1),'k')
% shadedErrorBar(vec,stat(:,2),stdstat(:,2),'r')
% shadedErrorBar(vec,stat(:,3),stdstat(:,3),'g')
% shadedErrorBar(vec,stat(:,4),stdstat(:,4),'b')
% legend('raw','z','Box Cox','Z-Box Cox')
print([figdir 'swstat_4conditions'],'-depsc','-r300')


%% Skewness

for s=1:S-2
    sBC(:,s)=skewness(tsBC(:,:,s)');
    sraw(:,s)=skewness(rsw(:,:,s)');
    sz(:,s)=skewness(zrsw(:,:,s)');
    szBC(:,s)=skewness(ztsBC(:,:,s)');
end


clear skew stdskew skewn
for n=1:length(a)
%     if sum(b==n)>0
        skew(n,1)=mean((sraw(b==n)));
        skew(n,2)=mean((sz(b==n)));
        skew(n,3)=mean((sBC(b==n)));
        skew(n,4)=mean((szBC(b==n)));
        
        stdskew(n,1)=std((sraw(b==n)));
        stdskew(n,2)=std((sz(b==n)));
        stdskew(n,3)=std((sBC(b==n)));
        stdskew(n,4)=std((szBC(b==n)));

        skewn(n)=sum(b==n)
%     end
end
figure; 
hold on
subplot(2,2,1)
shadedErrorBar(vec,skew(:,1),stdskew(:,1),'k')
axis([vec(1) vec(end) -2 1])
axis square
subplot(2,2,2)
shadedErrorBar(vec,skew(:,2),stdskew(:,2),'r')
axis([vec(1) vec(end) -2 1])
axis square
subplot(2,2,3)
shadedErrorBar(vec,skew(:,3),stdskew(:,3),'g')
axis([vec(1) vec(end) -2 1])
axis square
subplot(2,2,4)
shadedErrorBar(vec,skew(:,4),stdskew(:,4),'b')
axis([vec(1) vec(end) -2 1])
axis square
% legend('raw','z','Box Cox','Z-Box Cox')
print([figdir 'abskewness_4conditions_std'],'-depsc','-r300')




%% Median split on (could be bsxfuned for speed)

for s=1:S-2
    s
    for n=1:size(rsw,1)
        tmp=squeeze((rsw(n,:,s))); 
        tmp=(tmp-min(tmp))./(max(tmp)-min(tmp));
        rsw_medsplitvardif(n,s)=abs(var(tmp(tmp>median(tmp)))-var(tmp(tmp<median(tmp))));
  
        tmp=squeeze((zrsw(n,:,s))); 
        tmp=(tmp-min(tmp))./(max(tmp)-min(tmp));
        zrsw_medsplitvardif(n,s)=abs(var(tmp(tmp>median(tmp)))-var(tmp(tmp<median(tmp))));
        
        tmp=squeeze((tsBC(n,:,s))); 
        tmp=(tmp-min(tmp))./(max(tmp)-min(tmp));
        tsBC_medsplitvardif(n,s)=abs(var(tmp(tmp>median(tmp)))-var(tmp(tmp<median(tmp))));
        
        tmp=squeeze((ztsBC(n,:,s))); 
        tmp=(tmp-min(tmp))./(max(tmp)-min(tmp));
        ztsBC_medsplitvardif(n,s)=abs(var(tmp(tmp>median(tmp)))-var(tmp(tmp<median(tmp))));
    end
end


vec=0:.025:1;
[a b]=histc(reshape(abs(RStatic),34716*46,1),vec);

clear medsplitvardif medsplitvardifn stdmedsplitvardif
for n=1:length(a)
%     if sum(b==n)>0
        medsplitvardif(n,1)=mean(rsw_medsplitvardif(b==n));
        medsplitvardif(n,2)=mean(zrsw_medsplitvardif(b==n));
        medsplitvardif(n,3)=mean(tsBC_medsplitvardif(b==n));
        medsplitvardif(n,4)=mean(ztsBC_medsplitvardif(b==n));       
        
        stdmedsplitvardif(n,1)=std((rsw_medsplitvardif(b==n)));
        stdmedsplitvardif(n,2)=std((zrsw_medsplitvardif(b==n)));
        stdmedsplitvardif(n,3)=std((tsBC_medsplitvardif(b==n)));
        stdmedsplitvardif(n,4)=std((ztsBC_medsplitvardif(b==n)));

        medsplitvardifn(n)=sum(b==n)
%     end
end

figure; 
subplot(2,2,1)
shadedErrorBar(vec,medsplitvardif(:,1),stdmedsplitvardif(:,1),'k')
axis square
axis([0 1 -0.02 0.06])
subplot(2,2,2)
shadedErrorBar(vec,medsplitvardif(:,2),stdmedsplitvardif(:,2),'r')
axis square
axis([0 1 -0.02 0.06])
subplot(2,2,3)
shadedErrorBar(vec,medsplitvardif(:,3),stdmedsplitvardif(:,3),'g')
axis square
axis([0 1 -0.02 0.06])
subplot(2,2,4)
shadedErrorBar(vec,medsplitvardif(:,4),stdmedsplitvardif(:,4),'b')
axis square
axis([0 1 -0.02 0.06])
% legend('raw','z','Box Cox','Z-Box Cox')
print([figdir 'mediansplit_vardif_correlation'],'-depsc','-r300')






%% Variance difference 


vBC=squeeze(var(tsBC,[],2));
vraw=squeeze(var(rsw,[],2));
vz=squeeze(var(zrsw,[],2));
vzBC=squeeze(var(ztsBC,[],2));

EOI=1425
r=zeros(4,4); 
[r(1,2) p(1,2)]=corr(vraw(EOI,:)',vz(EOI,:)','type','Spearman');
[r(1,3) p(1,3)]=corr(vraw(EOI,:)',vBC(EOI,:)','type','Spearman');
[r(1,4) p(1,4)]=corr(vraw(EOI,:)',vzBC(EOI,:)','type','Spearman');
[r(2,3) p(2,3)]=corr(vz(EOI,:)',vBC(EOI,:)','type','Spearman');
[r(2,4) p(2,4)]=corr(vz(EOI,:)',vzBC(EOI,:)','type','Spearman');
[r(3,4) p(3,4)]=corr(vBC(EOI,:)',vzBC(EOI,:)','type','Spearman');



addpath('/data/william/functions/cbrewer/')
c=cbrewer('seq','YlOrBr',300); c(1,:)=[1 1 1]
figure; imagesc(1:4,1:4,triu(r(:,:))+triu(r(:,:))'); axis off; caxis([0 1]); axis square; colormap(c) 

print([figdir 'variance_correlation'],'-dpng','-r600')
figure; imagesc(1:4,1:4,r); colorbar; caxis([0 1]); axis square; colormap(c) 
print([figdir 'variance_correlation_colorbar'],'-dpng','-r600')







        
        
