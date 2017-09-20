  % Data for this project is taken from the EOEC data.
% It was preprocessed using CONN.
% Then 264 ROIs were created. They are saved in: [lDir]/'timeseries/eo_[subjectnumber].mat'

clear all
addpath(genpath('/data/william/Projects/dfc_time/'))

addpath('/data/william/Projects/dfc_timeseries/')
addpath(genpath('/data/william/toolbox/bc/'))


% These are founds in general functions folder
powernetclassclusters
roi_264_powerneuron

% Set paths
% Where data is loaded
ldir = '/data/william/Projects/eoec_dfctime/'
% Where figures are saved
sdir = '/data/william/Projects/Taylor/'


% Allsubject approach
S = 48; %NumberOfSubjects
L = 240 %DataLength
C = 1   %Condition
sid=0
c=1
clear vox mv mv_ec
for s=1:48

    if s==8 || s== 22
      % Skip subejcts with incomplete data
    else
      sid=sid+1;
      % Load eyes open data
      load([ldir 'timeseries/eo_' num2str(s) '.mat'])
      dat = eo';
      for roi=1:264
          voxtmp(roi,:,sid)=dat(roi,:)-mean(dat(roi,:));
      end
      tmp = voxtmp(:,:,sid);
      tmp1 = reshape(tmp,size(tmp,1)*size(tmp,2),1);
      tmp=(tmp-mean(tmp1))./std(tmp1);
      vox(:,:,c,sid)=tmp;

    end

end


% Static connectivity

clear StatConn
for s=1:S-2
StatConn(:,:,s)=corr(squeeze(vox(:,:,:,s)'),'type','Spearman');
end

cmap=cbrewer('div','RdBu',15)
cmap=flipud(cmap);

figure
imagesc(mean(StatConn(PlotOrder,PlotOrder,:),3))
colormap(cmap)
caxis([-.7 .7])
axis square
colorbar
print(gcf,[sdir '/data/william/Projects/Taylor/staticconn'],'-r300','-depsc')


% Sliding window

clear swdfc swdfc_pearson
window=22 %windowsize = 45 -> 22*2+1 (*2 = 90 seconds)
w=1
for s=1:S-2
    s
    parfor w=window+1:size(vox,2)-window
        swdfc(:,:,w-window,s)=corr(vox(:,w-window:w+window-1,s)','type','Spearman');
        swdfc_pearson(:,:,w-window,s)=corr(vox(:,w-window:w+window-1,s)');
    end
end
window=32 %windowsize = 65 -> *32+1 (*2 = 130 seconds)
w=1
for s=1:S-2
    s
    parfor w=window+1:size(vox,2)-window
        swdfc130(:,:,w-window,s)=corr(vox(:,w-window:w+window-1,s)','type','Spearman');
    end
end
window=12 %windowsize = 25 -> 12*2+1 (*2 = 50 seconds)
w=1
for s=1:S-2
    s
    parfor w=window+1:size(vox,2)-window
        swdfc50(:,:,w-window,s)=corr(vox(:,w-window:w+window-1,s)','type','Spearman');
    end
end
delete(gcp)

% Calculate mean and variance for each sliding window type
for s=1:S-2
    s
    c=0;
for n=1:264,
    for m=1:264,
        c=c+1;
        v(c,s)=var(swdfc(n,m,:,s));
        mu(c,s)=mean(swdfc(n,m,:,s));
        vp(c,s)=var(swdfc_pearson(n,m,:,s));
        mup(c,s)=mean(swdfc_pearson(n,m,:,s));
        v130(c,s)=var(swdfc130(n,m,:,s));
        mu130(c,s)=mean(swdfc130(n,m,:,s));
        v50(c,s)=var(swdfc50(n,m,:,s));
        mu50(c,s)=mean(swdfc50(n,m,:,s));
    end;
end
end
save([sdir '/data/william/Projects/Taylor/data_varmean'],'v','m','vp','mup','v130','mu130','v50','mu50')
load([sdir '/data/william/Projects/Taylor/data_varmean'])


%% Plot mean and variance

figure;
scatter(mean(mu(:,:),2),mean(v(:,:),2),2,'ko','filled')
axis([-1 1 0 .1])
set(gca,'FontName','Arial','FontSize',12)
xlabel('Mean')
ylabel('Variance')
print(gcf,[sdir '/data/william/Projects/Taylor/varmeanscatter'],'-r300','-depsc')
print(gcf,[sdir '/data/william/Projects/Taylor/fig2'],'-r300','-depsc')


%% Plot mean and variance  (supplementary figures)


figure('Position',[100 100 500 1000]);
subplot(3,1,2)
scatter(mean(mu130(:,:),2),mean(v130(:,:),2),2,'ko','filled')
axis([-1 1 0 .07])
set(gca,'FontName','Arial','FontSize',12)
xlabel('Mean')
ylabel('Variance')
title('Window Length = 130seconds')
axis square
subplot(3,1,1)
scatter(mean(mu50(:,:),2),mean(v50(:,:),2),2,'ko','filled')
axis([-1 1 0 .2])
set(gca,'FontName','Arial','FontSize',12)
xlabel('Mean')
ylabel('Variance')
title('Window Length = 50seconds')
axis square

subplot(3,1,3)
scatter(mean(mup(:,:),2),mean(vp(:,:),2),2,'ko','filled')
axis([-1 1 0 .1])
set(gca,'FontName','Arial','FontSize',12)
xlabel('Mean')
ylabel('Variance')
title('Window Length = 90seconds, Pearson Correlation')
axis square

print(gcf,[sdir '/data/william/Projects/Taylor/figs2'],'-r300','-depsc')
print(gcf,[sdir '/data/william/Projects/Taylor/figs2'],'-r300','-dpng')




%% Different ways to try and plot the main result.


figure
scatter(mean(mu(:,28),2),mean(v(:,28),2),'ko')
axis([-.8 1 0 .23])
print(gcf,[sdir '/data/william/Projects/Taylor/varmeanscatter_examplesubject'],'-r300','-depsc')


figure;
loglog(abs(mean(mu(:,:),2)),mean(v(:,:),2),'ko')
axis([0 1 0 .1])
print(gcf,[sdir '/data/william/Projects/Taylor/varmeanscatter_loglog'],'-r300','-depsc')



figure
loglog(abs(mean(mu(:,28),2)),mean(v(:,28),2),'ko')
axis([0 1 0 .23])
print(gcf,[sdir '/data/william/Projects/Taylor/varmeanscatter_examplesubject_loglog'],'-r300','-depsc')


%% Create figure 1

s1=.05.*sind(1:800)+.4;
s2=.2.*sind(180:979);
figure
hold on
plot(s1,'k')
plot(s2,'k')
axis([0 800 -.5 .5])
print(gcf,[sdir '/data/william/Projects/Taylor/var_mean_schema'],'-r300','-depsc')

figure
subplot(4,1,1)
th1=zeros(length(s1),1);
th2=zeros(length(s2),1);
th1(s1>0.4)=1;
th2(s2>0.4)=1;
imagesc([th1 th2]')
colormap('hot')

subplot(4,1,2)
th1=zeros(length(s1),1);
th2=zeros(length(s2),1);
th1(s1>0.3)=1;
th2(s2>0.3)=1;
imagesc([th1 th2]')
colormap('hot')

subplot(4,1,3)
th1=zeros(length(s1),1);
th2=zeros(length(s2),1);
th1(s1>0.15)=1;
th2(s2>0.15)=1;
imagesc([th1 th2]')
colormap('hot')

subplot(4,1,4)
th1=zeros(length(s1),1);
th2=zeros(length(s2),1);
th1(s1>mean(s1)+std(s1))=1;
th2(s2>mean(s2)+std(s2))=1;
imagesc([th1 th2]')
colormap('hot')

print(gcf,[sdir '/data/william/Projects/Taylor/var_mean_schema_thresholds'],'-r300','-depsc')

addpath(genpath('/data/william/functions/'))
cmap=cbrewer('div','RdBu',15)
cmap=flipud(cmap);

varmat=mean(var(swdfc,[],3),4);

figure

subplot(1,2,1)

imagesc(varmat(PlotOrder,PlotOrder))
colorbar
colormap(cmap)


%% Plot within and between connectivity (for each window type)

within=[]; between=[];
for n=1:264,
    for m=1:264
        if n>=m && n~=m
            if PowerNetClass(n)==PowerNetClass(m)
               within(end+1)=varmat(n,m);
            else
                between(end+1)=varmat(n,m);
            end
        end
    end
end

subplot(1,2,2)
hold on
errorbar([mean(within) mean(between)],[std(within) std(between)],'kx')
bar([mean(within) mean(between)],'k')
[h p]=ttest2(within,between)

print(gcf,[sdir '/data/william/Projects/Taylor/varnetwork'],'-r300','-depsc')


cmap=cbrewer('div','RdBu',15)
cmap=flipud(cmap);
varmat_pearson=mean(var(swdfc_pearson,[],3),4);


figure
subplot(1,2,1)
imagesc(varmat_pearson(PlotOrder,PlotOrder))
colorbar
colormap(cmap)

within=[]; between=[];
for n=1:264,
    for m=1:264
        if n>=m && n~=m
            if PowerNetClass(n)==PowerNetClass(m)
               within(end+1)=varmat_pearson(n,m);
            else
                between(end+1)=varmat_pearson(n,m);
            end
        end
    end
end

subplot(1,2,2)
hold on
errorbar([mean(within) mean(between)],[std(within) std(between)],'kx')
bar([mean(within) mean(between)],'k')
[h p]=ttest2(within,between)
axis square

print(gcf,[sdir '/data/william/Projects/Taylor/varnetwork_pearson'],'-r300','-depsc')



varmat_130=mean(var(swdfc130,[],3),4);

figure
subplot(1,2,1)
imagesc(varmat_130(PlotOrder,PlotOrder))
colorbar
colormap(cmap)

within=[]; between=[];
for n=1:264,
    for m=1:264
        if n>=m && n~=m
            if PowerNetClass(n)==PowerNetClass(m)
               within(end+1)=varmat_130(n,m);
            else
                between(end+1)=varmat_130(n,m);
            end
        end
    end
end

subplot(1,2,2)
hold on
errorbar([mean(within) mean(between)],[std(within) std(between)],'kx')
bar([mean(within) mean(between)],'k')
[h p]=ttest2(within,between)
axis square

print(gcf,[sdir '/data/william/Projects/Taylor/varnetwork_130'],'-r300','-depsc')



varmat_50=mean(var(swdfc50,[],3),4);

figure
subplot(1,2,1)
imagesc(varmat_50(PlotOrder,PlotOrder))
colorbar
colormap(cmap)

within=[]; between=[];
for n=1:264,
    for m=1:264
        if n>=m && n~=m
            if PowerNetClass(n)==PowerNetClass(m)
               within(end+1)=varmat_50(n,m);
            else
                between(end+1)=varmat_50(n,m);
            end
        end
    end
end

subplot(1,2,2)
hold on
errorbar([mean(within) mean(between)],[std(within) std(between)],'kx')
bar([mean(within) mean(between)],'k')
[h p]=ttest2(within,between)
axis square
print(gcf,[sdir '/data/william/Projects/Taylor/varnetwork_50'],'-r300','-depsc')
