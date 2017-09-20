
% The point of this code is to calculate the mutual information between two
% thresholded f-graphlets from two different datasets (RL (used in main
% text)) and LR - both from session 1. 

% NOTE - this code uses parittion distance which is NOT the optimal/correct way to do this.
% Since no graphlet has unique non-zero values, this is alright. 

% I've added taxicab distance of binary matrices in the code to show
% identical results (BinD) (except reversed). 

addpath('/data/william/toolbox/bc/')

th=0.09 %Chose proporitional threshold
ds='' % where data is kept (blank for RL, LR for LR)


load([mdir '/' ds '/graphdata2'])
clear thgraph binthgraph
for f=1:size(graphfreqsub,3)
    thgraph(:,:,f)=threshold_proportional(graphfreqsub(:,:,f),th);
    binthgraph(:,:,f) = weight_conversion(thgraph(:,:,f),'binarize');
end

ds='LR'
load([mdir '/' ds '/graphdata2'])
clear thgraphlr
for f=1:size(graphfreqsub,3)
    thgraphlr(:,:,f)=threshold_proportional(graphfreqsub(:,:,f),th);
    binthgraphlr(:,:,f) = weight_conversion(thgraphlr(:,:,f),'binarize');    
end

%Calculate MI between paritions (and taxicab distance)
for f=1:78
    f
    parfor ff=1:78
        [V9(f,ff), MI9(f,ff)]=partition_distance(thgraph(:,:,f),thgraphlr(:,:,ff));
        BinD9(f,ff) = sum(sum(abs(binthgraph(:,:,f)-binthgraphlr(:,:,ff))));
    end
end

%Repeat everything with th=0.05
th=0.05
ds=''
load([mdir '/' ds '/graphdata2'])
clear thgraph binthgraph
for f=1:size(graphfreqsub,3)
    thgraph(:,:,f)=threshold_proportional(graphfreqsub(:,:,f),th);
    binthgraph(:,:,f) = weight_conversion(thgraph(:,:,f),'binarize');
end

ds='LR'
load([mdir '/' ds '/graphdata2'])
clear thgraphlr
for f=1:size(graphfreqsub,3)
    thgraphlr(:,:,f)=threshold_proportional(graphfreqsub(:,:,f),th);
    binthgraphlr(:,:,f) = weight_conversion(thgraphlr(:,:,f),'binarize');    
end
for f=1:78
    f
    parfor ff=1:78
        [V5(f,ff), MI5(f,ff)]=partition_distance(thgraph(:,:,f),thgraphlr(:,:,ff));
        BinD5(f,ff) = sum(sum(abs(binthgraph(:,:,f)-binthgraphlr(:,:,ff))));
    end
end

%Get distance between delta freq
c=0; 
for f=1:78
    for ff=1:78
        c=c+1
        freqdist(c)=abs(f-ff);
        MIcol9(c)=MI9(f,ff);
    end
end
clear  avgMIc9
for c=0:77
    fid=find(freqdist==c)
    avgMIc9(c+1)=mean(MIcol9(fid)); 
    stdMIc9(c+1)=std(MIcol9(fid)); 
end

%Repeat with 5
c=0; 
for f=1:78
    for ff=1:78
        c=c+1
        freqdist(c)=abs(f-ff);
        MIcol5(c)=MI5(f,ff);
    end
end
for c=0:77
    fid=find(freqdist==c)
    avgMIc5(c+1)=mean(MIcol5(fid)); 
    stdMIc5(c+1)=std(MIcol5(fid)); 
end

%Plot distance delta freq
figure
subplot(1,2,1)
shadedErrorBar([0:freq(2)-freq(1):(freq(2)-freq(1))*77],avgMIc5,stdMIc5)
axis([0 .1 .5 .9])
axis square
subplot(1,2,2)
shadedErrorBar([0:freq(2)-freq(1):(freq(2)-freq(1))*77],avgMIc9,stdMIc9)
axis([0 .1 .5 .9])
axis square
print([mdir 'SupF2cd_replication'],'-depsc','-r300')

th=0.05
ds=''
load([mdir '/' ds '/graphdata2'])
clear thgraph
for f=1:size(graphfreqsub,3)
    thgraph(:,:,f)=threshold_proportional(graphfreqsub(:,:,f),th);
end

ds='LR'
load([mdir '/' ds '/graphdata2'])
clear thgraphlr
for f=1:size(graphfreqsub,3)
    thgraphlr(:,:,f)=threshold_proportional(graphfreqsub(:,:,f),th);
end

for f=1:78
    f
    parfor ff=1:78
        [V5(f,ff), MI5(f,ff)]=partition_distance(thgraph(:,:,f),thgraphlr(:,:,ff));
    end
end

% plot frequency-frequency MI (as noted in introduction, this is almost identical for taxicab of binary)


addpath(genpath('/data/william/functions/cbrewer/'))
cmap=cbrewer('div','RdBu',101)
cmap=flipud(cmap);

figure;
subplot(1,2,1)
imagesc(freq,freq,MI5)
axis square
caxis([0.7 .9])
colormap(cmap((length(cmap)+1)/2:end,:))

subplot(1,2,2)
imagesc(freq,freq,MI9)
axis square
colormap(cmap((length(cmap)+1)/2:end,:))

caxis([0.7 .9])
print([mdir 'SupF2_replication'],'-depsc','-r300')


figure
imagesc([])
caxis([0.7 .9])
colorbar
colormap(cmap((length(cmap)+1)/2:end,:))
print([mdir 'SupF2_replication_colorbar'],'-depsc','-r300')
