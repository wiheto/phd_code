clear all
addpath(genpath('/data/william/Projects/dfc_time/'))

addpath('/data/william/toolbox/hungarianalgorithem/')
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
dir = 'LR'


% Allsubject approach
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

clear tmp tmp1
%clear kmeangroups centroid sumD D
pool=parpool(20); stream=RandStream('mlfg6331_64'); options=statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);

for n=1:3
voxkeep=[]
if n==1
for s=1:100
voxkeep(end+1:end+400)=(1:400)+(s-1)*1200;
end
elseif n==2
for s=1:100
voxkeep(end+1:end+400)=([401:800])+(s-1)*1200;
end      
elseif n==3
for s=1:100
voxkeep(end+1:end+400)=(801:1200)+(s-1)*1200;
end   
end

for k=1:20
   k
   [k_kfold(:,k,n)]=kmeans(vox_pcascore(voxkeep,1:numberofpca),k,'options',options,'MaxIter',1000,'replicates',20);
end
end

% save([mdir 'cross_val_save'],'k_kfold');

for n=1:3
    voxkeep=[]
if n==1
for s=1:100
voxkeep(end+1:end+400)=(1:400)+(s-1)*1200;
end
elseif n==2
for s=1:100
voxkeep(end+1:end+400)=([401:800])+(s-1)*1200;
end      
elseif n==3
for s=1:100
voxkeep(end+1:end+400)=(801:1200)+(s-1)*1200;
end   
end
vk(:,n)=voxkeep;
    for k=1:20
        for kk=1:k
        g(:,:,kk,k,n)=corr(voxall(:,voxkeep(k_kfold(:,k,n)==kk))');
        end
    end
end

clear vm D
for k=1:20, 
    %Get mean voxel patterns per k. 
    for n=1:3, 
        for kk=1:k 
            vm(:,kk,n)=mean(voxall(:,vk(k_kfold(:,k,n)==kk,n)),2);
        end
    end  
    % Get the distance between each "set par" for each k pair 
    for kk=1:k
        for kkk=1:k
    for n=1:3
        for m=1:3
            D{n,m,k}(kk,kkk)=mean(abs(vm(:,kk,n)-vm(:,kkk,m)));
        end
    end
        end
    end
end

clear H
%Using Hungarian Algorithem, find optimal assignment of K. 
for k=1:20
for n=1:3
    for m=1:3
        if m<=n
        else
        H{n,m,k} = assignmentoptimal(D{n,m,k})
        end
    end
end
end

for k=1:20

for n=1:k
kdim{k}(n,1)=n 
kdim{k}(n,2)=H{1,2,k}(n)
kdim{k}(n,3)=H{2,3,k}(H{1,2,k}(n))
end
if sum(kdim{k}(:,3)==H{1,3,k})==k
    conv(k)=1;
else 
    conv(k)=0;
end

end

clear GD MI
for k=1:20
    if conv(k)==1
        i=0
        for n=1:3
            for m=1:3
                if m<=n
                else
                    i=i+1;
                    for kk=1:k
                        tmp = reshape(g(:,:,kk,k,n),264*264,1);
                        tmp2 = reshape(g(:,:,kdim{k}(kk,m),k,m),264*264,1);
%                         g1=threshold_proportional(g(:,:,kk,k,n),.05);
%                         g2=threshold_proportional(g(:,:,kdim{k}(kk,m),k,m),.05);
%                         [~, MI{k}(i,kk)]=partition_distance(g1,g2)
                        GD{k}(i,kk)= sum(abs(tmp-tmp2));
                    end
                end
            end
        end
        GGD(k)=mean(mean(GD{k}));
        GMI(k)=mean(mean(MI{k}));
    end
end


% gdiff(:,conv==0)=NaN; 
% gdiff(:,1)=NaN
GGD(:,1)=NaN
figure; plot(1:12,mean(GGD(:,:),1)./(264*264-264),'k','LineWidth',4); 
axis([0 20 0 0.1]); 
print(gcf,'/data/william/Projects/dfc_time/figures/Distance_kchoice_cv','-r300','-depsc')
d=mean(gdiff,1)

d(2:end)=(d(2:end)-min(d(2:end)))/(max(d(2:end))-min(d(2:end)))


clear vm_ex 
k=8
    for n=1:3, 
        for kk=1:k 
            vm_ex(:,kk,n)=mean(voxall(:,vk(k_kfold(:,k,n)==kdim{k}(kk,n),n)),2);
        end
    end
    
addpath(genpath('/data/william/functions/'))
cmap=cbrewer('div','RdBu',200)
cmap=flipud(cmap);
    
figure
subplot(3,1,1)
imagesc(vm_ex(PlotOrder,:,1)')
colormap(cmap)
caxis([-1.5 1.5])
subplot(3,1,2)
imagesc(vm_ex(PlotOrder,:,2)')
colormap(cmap)

caxis([-1.5 1.5])
subplot(3,1,3)
imagesc(vm_ex(PlotOrder,:,3)')
colormap(cmap)

caxis([-1.5 1.5])

addpath('/data/william/toolbox/exportfig/')
export_fig '/data/william/Projects/dfc_time/figures/meanvox_replicateK' '-r300' '-eps'

figure
imagesc(vm_ex(PlotOrder,:,3)')
colormap(cmap)

caxis([-1.5 1.5])
colorbar
export_fig '/data/william/Projects/dfc_time/figures/meanvox_replicateK_colorbar' '-r300' '-eps'



