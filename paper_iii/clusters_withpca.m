%% This script is very uncommented 

%% It clusters the fMRI data

%% Then plots the descriptives of the data / evaluation of the cluster

%% Then the sgraphlets 

%% Then distributions of the sgraphlets

%% Then correlatiing sgraphlet distance and s-graphlet transition probability 

%% Then performing clustering and sliding windowed data	



clear all
addpath(genpath('/data/william/Projects/dfc_time/'))

addpath('/data/william/Projects/dfc_timeseries/')
powernetclassclusters
roi_264_powerneuron
mdir = '/data/william/Projects/dfc_time/'
ldir = '/data/william/Projects/dfc_timeseries/'
addpath('/usr/local/fieldtrip-20141009/')
ft_defaults
addpath(genpath('/data/william/toolbox/NBS1.2/'))
%Create 264 roi from Power et al
%roi_264_powerneuron

sub = '1';
dir = ''


for r=1:size(coord,1)
    fid=find(coord(:,1)<coord(r,1)+15 & coord(:,1)>coord(r,1)-15);
    if isempty(fid)==0
        fid2=find(coord(fid,2)<coord(r,2)+15 & coord(fid,2)>coord(r,2)-15);
        if isempty(fid2)==0
            fid3=find(coord(fid(fid2),3)<coord(r,3)+15 & coord(fid(fid2),3)>coord(r,3)-15);
            if isempty(fid3)==0
                dellocal{r}=fid(fid2(fid3));
            end
        end
    end
end


sid=0
graph=[]
for s=1:100
    s
    sid=sid+1;
    load([ldir 'timeseries/' dir 'fwd_ts' num2str(s) '.mat'])

    clear tmp
    tmp=ft_preproc_bandpassfilter(timeseries',1/0.72,[0.01 0.1]);
    vox(:,:,sid) = reshape(tmp,size(tmp,1),size(tmp,2)*size(tmp,3));

end
for s=1:100
    for roi=1:264
        vox(roi,:,s)=vox(roi,:,s)-mean(vox(roi,:,s));
    end
    tmp = vox(:,:,s);
    tmp1 = reshape(tmp,size(tmp,1)*size(tmp,2),1); 
    tmp=(tmp-mean(tmp1))./std(tmp1);
%     tmp = (tmp-min(min(tmp)))/(max(max(tmp))-(min(min(tmp)))); %Scale between 0 and 1.     
    vox(:,:,s)=tmp; 
end


clear kmeangroups cost vol pow
voxall = reshape(vox,264,1200*100);

[coeff vox_pcascore latent a extent] = pca(voxall(:,:)');
for n=1:264
if sum(extent(1:n))<85
else
    break
end
end
numberofpca=n;

[a, pcaedat, c d extent]=pca(voxall(:,:)');


%% Cluster PCA with

clear tmp tmp1
%clear kmeangroups centroid sumD D
pool=parpool(20); stream=RandStream('mlfg6331_64'); options=statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
for k=1:20
   k
   [kmeangroups(:,k), centroid{k} sumD{k} D{k}]=kmeans(vox_pcascore(:,1:numberofpca),k,'options',options,'MaxIter',1000,'replicates',20);
end
delete(gcp)
save([mdir '/cost_of_kmean_run3' dir],'centroid','kmeangroups','sumD','D')



%% Plot evalation of clustering 


load([mdir '/cost_of_kmean_run3' dir])


clear s
for k=2:20
    k
    for i=1:120000
        nk = 1:k;
        nk(kmeangroups(i,k)) = []; 
        s(i) = (min(D{k}(i,nk))-D{k}(i,kmeangroups(i,k)))/max([D{k}(i,k) min(D{k}(i,nk))]);
    end
    sil(k) = mean(s)  
end
figure
plot(sil)
print(gcf,[mdir 'figures/silloutte'],'-depsc','-r300') 

eva = evalclusters(vox_pcascore(:,1:numberofpca),kmeangroups,'DaviesBouldin')


%% Clusters in time series

for n=1:8, 
    fid=find(kmeangroups(:,8)==n); 
    voxmean(n,:)=mean(voxall(:,fid),2); 
    voxmax(n,:)=max(voxall(:,fid),[],2); 
    voxmin(n,:)=min(voxall(:,fid),[],2); 
    voxvar(n,:)=std(voxall(:,fid),[],2); 
end


h=figure('Position',[200,200,600,200]); 
imagesc(voxmean(:,PlotOrder))
colormap(cmap)
set(h,'renderer','painters')
colorbar
addpath('/data/william/toolbox/exportfig/')
% caxis([-2.5 2.5])
export_fig([mdir '/figures/avg_ROI_per_state_ef.eps'], '-r300') 

h=figure('Position',[200,200,600,200]); 
imagesc(voxvar(:,PlotOrder))
colormap(cmap)
set(h,'renderer','painters')
colorbar
addpath('/data/william/toolbox/exportfig/')
% caxis([-2.5 2.5])
export_fig([mdir '/figures/var_ROI_per_state_ef.eps'], '-r300') 


h=figure('Position',[200,200,600,250]); 
subplot(2,1,1)
imagesc(voxmax(:,PlotOrder))
colormap(cmap)
caxis([-10 10])
set(h,'renderer','painters')
colorbar
subplot(2,1,2)
imagesc(voxmin(:,PlotOrder))
colormap(cmap)
set(h,'renderer','painters')
colorbar
caxis([-10 10])
addpath('/data/william/toolbox/exportfig/')
export_fig([mdir '/figures/maxmin_ROI_per_state_ef.eps'], '-r300') 



%% Examples of cluster assignment

load([mdir '/cost_of_kmean_run3'])
kchoice = 8

addpath(genpath('/data/william/toolbox/exportfig/'))
addpath('/data/william/akalla/william/Projects/dfc_timeseries/')
powernetclassclusters
clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a NetOrder]=sort(fid,'descend') 
NetOrder = [NetOrder -1];

addpath(genpath('/data/william/functions/'))
cmap=cbrewer('div','RdBu',200)
cmap=flipud(cmap);


h=figure
hist(kmeangroups(:,kchoice),1:kchoice)
axis([0.5 kchoice+.5 0 30000])
print(h,[mdir '/figures/hist_clusterassignments'],'-depsc','-r300') 


s=5
kt=kmeangroups(1+(s-1)*1200:s*1200-800,kchoice)
h=figure('Position',[200,200,900,1000]); 
c=hsv(kchoice)

subplot(3,1,1)
r=10
plot(vox(r,1:400,s),'k','LineWidth',2)
hold on
for n=1:8
    fid=find(kt==n)
    scatter(fid,vox(r,fid,s),30,'MarkerFaceColor',c(n,:),'MarkerEdgeColor','k');
end

subplot(3,1,2)
r=118
plot(vox(r,1:400,s),'k','LineWidth',2)
hold on
for n=1:8
    fid=find(kt==n)
    scatter(fid,vox(r,fid,s),30,'MarkerFaceColor',c(n,:),'MarkerEdgeColor','k');
end

subplot(3,1,3)
r=190
plot(vox(r,1:400,s),'k','LineWidth',2)
hold on
for n=1:8
    fid=find(kt==n)
    scatter(fid,vox(r,fid,s),30,'MarkerFaceColor',c(n,:),'MarkerEdgeColor','k');
end
export_fig([mdir '/figures/timeseries_vs_states.eps'], '-r300') 


figure
plot(kmeangroups(1:1200,kchoice))
figure
hist(kmeangroups(1:1200,kchoice),kchoice)
h=figure
s=5
subplot(3,1,1)
hist(kmeangroups(1200*(s-1)+1:1200*s,kchoice),1:kchoice)
axis([0.5 kchoice+.5 0 300])
set(gca,'YTick',[0 300]) 

subplot(3,1,2)
s=15
hist(kmeangroups(1200*(s-1)+1:1200*s,kchoice),1:kchoice)
axis([0.5 kchoice+.5 0 300])
set(gca,'YTick',[0 300]) 

s=46
subplot(3,1,3)
hist(kmeangroups(1200*(s-1)+1:1200*s,kchoice),1:kchoice)
set(gca,'YTick',[0 300]) 
axis([.5 kchoice+.5 0 300])
print(h,[mdir '/figures/example_subjectclusterassignment_hist'],'-depsc','-r300') 


figure
clear subject_ex
s=5
subject_ex(:,1) = kmeangroups(1200*(s-1)+1:1200*s,kchoice);
s=15
subject_ex(:,2) = kmeangroups(1200*(s-1)+1:1200*s,kchoice);
s=46
subject_ex(:,3) = kmeangroups(1200*(s-1)+1:1200*s,kchoice);
h=figure
subplot(3,1,1)
imagesc(subject_ex(:,1)')
colorbar
colormap(hsv(kchoice))
subplot(3,1,2)
imagesc(subject_ex(:,2)')
colorbar
colormap(hsv(kchoice))
subplot(3,1,3)
imagesc(subject_ex(:,3)')
colorbar
colormap(hsv(kchoice))
print(h,[mdir '/figures/example_subjectclusterassignment'],'-depsc','-r300') 


%% PLOT SGRAGPHLETS

addpath(genpath('/data/william/toolbox/bc/'))
addpath(genpath('/data/william/functions/'))
h=figure('Position',[100 0 150 900])
clear g graph
ii=0
clear keptnetid bingraph graph 
cmap=cbrewer('div','RdBu',15)
cmap=flipud(cmap);

for n=1:kchoice
    if sum(kmeangroups(:,kchoice)==n)>500
        ii=ii+1;
        keptnetid(ii)=n
        g(:,:,n)=corr(voxall(:,kmeangroups(:,kchoice)==n)');
        
        subplot(12,1,ii)

        for r=1:length(dellocal)
            g(r,dellocal{r},n)=0;
            g(r,r,n)=0;
        end
         graph(:,:,ii) = g(:,:,n); 
        imagesc(graph(PlotOrder,PlotOrder,ii))
      
    colormap(cmap)  
        caxis([-.5 .5])
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        tmp = threshold_proportional(graph(:,:,ii),0.05);
        bingraph(:,:,ii) = weight_conversion(tmp,'binarize');
        tmp = threshold_proportional(graph(:,:,ii),0.05);

    end
end

save([mdir 'sgraphlets.mat'],'graph')

print(h,[mdir 'figures/sgraphlets'],'-depsc','-r300') 
h=figure
colorbar
colormap(cmap) 
print(h,[mdir '/figures/sgraphlets_colorbar'],'-depsc','-r300') 

clear keptnetid pgraph

for perm=1:1000
    perm
    ii=0;
    [a r] = sort(rand(size(voxall,2),1));
    pvoxall = voxall(:,r); 
for n=1:kchoice
    if sum(kmeangroups(:,kchoice)==n)>500
        ii=ii+1;
        keptnetid(ii)=n;
        g(:,:,n)=corr(pvoxall(:,kmeangroups(:,kchoice)==n)');        

        pgraph(:,:,ii,perm) = g(:,:,n); 
        
    end
end
end

clear siggraph
siggraph=zeros(264,264,ii);
for n=1:ii
    n
    s=sort(squeeze(pgraph(:,:,n,:)),3);
    for i=1:264
        for j=1:264
            if graph(i,j,n)>s(i,j,975) && graph(i,j,n)>0
               siggraph(i,j,n)= graph(i,j,n);
            end
        end
    end
end
figure
for n=1:ii
    subplot(4,3,n)
    imagesc(siggraph(:,:,n))
    caxis([-.1 .1])
end
save([mdir '/figures/difference_to_randomizedtime'],'siggraph')    



%% Plot sgraphlet distributions


addpath('/data/william/toolbox/exportfig/')
clear kmeanfreq
for s=1:100
    ks=kmeangroups((s-1)*1200+1:s*1200,kchoice)
    ds=1:length(ks); 
    dsfid=find(diff(ks)==0); 
    ksds = ks(dsfid) 
    ii=0; 
    for n=length(dsfid):-1:1
        ds(dsfid(n)+1)=[]; 
    end
    kmeanfreq{s}(:,1)=ks(ds); %State in at switch
    kmeanfreq{s}(1:end-1,2)=diff(ds); %Length to next State
    kmeanfreq{s}(end,2)=diff([ds(end) 1200])+1; %Fill in length of last state
end
NextPair=[];
for s=1:100
    for n=1:length(kmeanfreq{s})-1
        NextPair(end+1,1)=kmeanfreq{s}(n,1); 
        NextPair(end,2)=kmeanfreq{s}(n+1,1); 
    end
end
clear prob
for n=1:max(NextPair(:,1),[],1)
    fid=find(NextPair(:,1)==n);
    counts=histc(NextPair(fid,2),1:max(NextPair(:,1),[],1));
    prob(n,:)=counts./length(fid); 
end
clear AvgNetLength  StdNetLength NNetLength
kdis=cell2mat(kmeanfreq');
for n=1:kchoice
    fid=find(kdis(:,1)==n);
    AvgNetLength(n)=mean(kdis(fid,2));
    StdNetLength(n)=std(kdis(fid,2));
    NNetLength(n)=length(kdis(fid,2));
end
h=figure
hold on
bar(AvgNetLength)
errorbar(AvgNetLength,StdNetLength,'xk')
box on
print(h,[mdir '/figures/lengthinstate'],'-depsc','-r300')

statesandlength=cell2mat(kmeanfreq')
h=figure
for n=1:kchoice
    fid=find(statesandlength(:,1)==n)
    subplot(2,4,n)
    hist(statesandlength(fid,2),1:50)
    title(num2str(n))
    axis([0 50 0 650])
end
export_fig([mdir '/figures/hist_length_in_each_state.eps'], '-r300') 




%% Distance vs transition figure 


for n=1:8
    for nn=1:8
        distgraph(n,nn)=sum(sum(abs(graph(:,:,n)-graph(:,:,nn)))) 
    end
end

clear d p
c=0
for n=1:8
    for nn=1:8
        if n==nn
        else
            c=c+1
            d(c)=distgraph(n,nn)./(264^2-264); 
            p(c)=prob(n,nn); 
        end
    end
end
figure('Position',[100 100 600 300])
scatter(d,p,50,'filled')
xlabel('Distance')
ylabel('Probability of transition')
lsline
[r pp]=corr(d',p','type','Spearman')
box on
addpath('/data/william/toolbox/exportfig/')
export_fig '/data/william/Projects/dfc_time/figures/distance_vs_transition' -r300 -eps



[ges s]=sort(ge)
figure
for n=1:kchoice
    subplot(3,3,n) 
    fidH=find(ges>ges(n))
    fidL=find(ges<ges(n))
    bar([sum(prob(s(n),s(fidH)),2) sum(prob(s(n),s(fidL)),2)])
end



%DOes one sgrpahlet correlate with movement? 

ldir = '/data/william/Projects/dfc_timeseries/'
addpath('/data/william/toolbox/Bramila/')
rmpath(genpath('/data/william/toolbox/spatialeconomics/'))
for s=1:100
    
    s   
    mv=load([ldir '/movement/' num2str(s) '.txt']);
    cfg.motionparam = mv(:,[1:6]);
    fwd=bramila_framewiseDisplacement(cfg);  
    load(([ldir '/timeseries/ts' num2str(s)]))    
    clear b
    badframes = find(fwd>0.5) 
    bf{s}=badframes; 
end



figure; 
kid=1; 
bfcol=[]; 
bf_all = []; 
for s=1:100
    tmp=kmeangroups(bf{s}+kid,kchoice);
    bfcol(end+1:end+length(tmp))=tmp;
    kid = kid+1200;
end
bf_all(end+1:end+length(bfcol))=bfcol; 
subplot(2,2,2) 
hist(bfcol,1:kchoice)

mv1=hist(bfcol,kchoice) 

title('Frames where fwd>0.5')
kid=1; 
bfcol=[]; 
for s=1:100
    tmp=kmeangroups(bf{s}+kid+1,kchoice);
    bfcol(end+1:end+length(tmp))=tmp;
    kid = kid+1200;
end
bf_all(end+1:end+length(bfcol))=bfcol; 
subplot(2,2,3) 
hist(bfcol,1:kchoice) 

mv2=hist(bfcol,kchoice) 
title('One volume after fwd>0.5')
kid=1; 
bfcol=[]; 
for s=1:100
    tmp=kmeangroups(bf{s}+kid+2,kchoice);
    bfcol(end+1:end+length(tmp))=tmp;
    kid = kid+1200;
end
%bf_all(end+1:end+length(bfcol))=bfcol; 
subplot(2,2,4) 
hist(bfcol,1:kchoice) 
mv3=hist(bfcol,kchoice) 

title('Two volumes after fwd>0.5')

kid=1; 
bfcol=[]; 
for s=1:100
    tmp=kmeangroups(bf{s}+kid-1,kchoice);
    bfcol(end+1:end+length(tmp))=tmp;
    kid = kid+1200;
end
bf_all(end+1:end+length(bfcol))=bfcol; 
subplot(2,2,1) 
hist(bfcol,1:kchoice) 
mv4=hist(bfcol,kchoice) 

addpath('/data/william/toolbox/exportfig/')
export_fig figures/sgraphlets_movement -r300 -eps

for n=1:kchoice
    percentwheremove(n)=sum(bf_all==n)./sum(kmeangroups(:,kchoice)==n)
    numwheremove(n)=sum(bf_all==n)
end



groupassignments=hist(kmeangroups(:,kchoice),kchoice)
figure
subplot(2,2,1)
scatter(groupassignments,mv1)
subplot(2,2,2)
scatter(groupassignments,mv2)
subplot(2,2,3)
scatter(groupassignments,mv3)
subplot(2,2,4)
scatter(groupassignments,mv4)










%Calculate different between subgraphs

clear pg1 pg2 pvals
permutations=10000
mycluster = parpool(8) 
for n=1:kchoice
    n
for nn=1:kchoice
    if n==nn
    else
        difgraph = zeros(264,264,permutations); 
        fidn=length(find(kmeangroups(:,kchoice)==n));
        fid=find(kmeangroups(:,kchoice)==n | kmeangroups(:,kchoice)==nn);
        parfor p=1:permutations
            [a r] = sort(rand(length(fid),1));
            pg1=corr(voxall(:,fid(r(1:fidn)))'); 
            pg2=corr(voxall(:,fid(r(fidn+1:end)))'); 
            difgraph(:,:,p)=pg1-pg2;
        end 
        d=graph(:,:,n)-graph(:,:,nn);
        for i=1:264
            for j=1:264
                L=length(find(d(i,j)<difgraph(i,j,:)))+1; 
                pvals(i,j,n,nn)=L/(permutations+1);
            end
        end
    end
    end
end
for n=1:kchoice
    pvals(:,:,n,n) = NaN; 
end
save([mdir 'pvalderived_betweensgraphs'],'pvals')

load([mdir 'pvalderived_betweensgraphs'])
p=reshape(pvals,264*264*kchoice*kchoice,1);
fid=find(isnan(p)==0);
[fdrp fdrco fdradj]= fdr(p(fid),0.025);
fdrcoall = zeros(length(p),1)*NaN; 
fdrcoall(fid) = fdradj; 
fdrcors = reshape(fdrcoall,264,264,kchoice,kchoice);
smat=zeros(264,264,kchoice,kchoice);
for n=1:kchoice
    n
    for nn=1:kchoice
        if n==nn
        else
        dd=graph(:,:,n);
        for i=1:264
            for j=1:264
                
                if fdrcors(i,j,n,nn)<0.0005 && dd(i,j)>0;
              
                        smat(i,j,n,nn)= dd(i,j);
                end
                
            end
        end
        end
    end
end


figure; 
c=0; 
for n=1:8, for nn=1:8, c=c+1;
subplot(8,8,c)
imagesc(smat(PlotOrder,PlotOrder,n,nn))
caxis([0 1])
colormap hot  
set(gca,'XTickLabel',[],'YTickLabel',[])
axis square
    end; end;
     
ni=0; clear GraphDif
PowerNetLabel{11,1}=-1
NetOrder(end)=11
for n=NetOrder
    f1 = find(PowerNetClass==PowerNetLabel{n,1});
    ni=ni+1; mi=0; 
    for m=NetOrder 
        f2 = find(PowerNetClass==PowerNetLabel{m,1});
        mi=mi+1;
        for k=1:kchoice
            for kk=1:kchoice 
                GraphDif(ni,mi,k,kk)=sum(sum(smat(f1,f2,k,kk)>0))/(length(f1)*length(f2));                   
            end
        end
    end
end
figure
kid=0; 
for k=1:kchoice,
    for kk=1:kchoice, 
        kid=kid+1; 
        subplot(kchoice,kchoice,kid) 
        imagesc(GraphDif(:,:,k,kk))
        caxis([0 1])
        colormap hot
        axis square 
        set(gca,'XTick',[])
        set(gca,'YTick',[])
    end
end

addpath('/data/william/toolbox/exportfig/')
export_fig figures/percent_differences_per_subgraph_state -r300 -eps









%% Sliding window part

%covariance sliding window (rudamentary)
clear swdfc
windowsize=120
w=1
swdfc=zeros(264,264,1080,100); 
for s=1:100
    s
    for w=1:size(vox,2)-windowsize
        swdfc(:,:,w,s)=cov(vox(:,w:w+windowsize-1,s)');
    end
end
swdfc=reshape(swdfc,264*264,size(swdfc,3)*100);


sw_pcanum=30;

s=5


[coeff sw_pcascore latent a extent] = pca(swdfc(:,1080*(s-1)+1:s*1080)');
sw_subject_ex(:,1) = kmeans(sw_pcascore(:,1:sw_pcanum),8,'replicates',20);
s=15
[coeff sw_pcascore latent a extent] = pca(swdfc(:,1080*(s-1)+1:s*1080)');
sw_subject_ex(:,2) = kmeans(sw_pcascore(:,1:sw_pcanum),8,'replicates',20);
s=46
[coeff sw_pcascore latent a extent] = pca(swdfc(:,1080*(s-1)+1:s*1080)');
sw_subject_ex(:,3) = kmeans(sw_pcascore(:,1:sw_pcanum),8,'replicates',20);

h=figure
subplot(3,1,1)
imagesc(sw_subject_ex(:,1)')
colorbar
colormap(hsv(8))
subplot(3,1,2)
imagesc(sw_subject_ex(:,2)')
colorbar
colormap(hsv(8))
subplot(3,1,3)
imagesc(sw_subject_ex(:,3)')
colorbar
colormap(hsv(8))
print(h,[mdir '/figures/example_slidingwindow_subjectclusterassignment'],'-depsc','-r300') 


