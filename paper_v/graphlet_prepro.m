%% Graphlet Preprocessing

mdir = '/data/william/Projects/dfc_timeseries/'
ds = '' %LR='LR' RL=''

%% Procedure: 
%Load, standardize, delete local connections, average over subjects
%Repeat for empirical, permutation (amplitude), permutation (imag/phase)


%% Prelim
addpath(mdir)
roi_264_powerneuron %Give "dellocal" which is the nearby connection to zero out

%% Empirical 

c=0;
clear graphfreqsub
for s=1:100
    c=c+1
    load([mdir '/graph/' ds '/s' num2str(s)])
    graphfreqsub(:,:,:,s)=graph(:,:,1:78); 
end
clear graph 

%standerdize per f-graphlet
clear z_graphfreqsub 
for s=1:size(graphfreqsub,4)
    s
    tmp1 = graphfreqsub(:,:,:,s);
    for f=1:78
        tmp = reshape(tmp1(:,:,f),264*264,1);
        ztmp=(tmp1(:,:,f)-mean(tmp))./std(tmp);
        ztmp= reshape(ztmp,264,264); 
        z_graphfreqsub(:,:,f,s)=ztmp;
    end
end
clear graphfreqsub tmp tmp1 tmpz

%Zero connections
for r=1:length(dellocal)
    z_graphfreqsub(r,dellocal{r},:,:)=0;
    z_graphfreqsub(r,r,:,:)=0;
end

% Take average over groups
graphfreq=mean(z_graphfreqsub,4); 
clear z_graphfreqsub


permgraphfreq=zeros([size(graphfreq) 100]);
%Permutation load
clear graphfreqsub
for perm=1:100
    %wavelet
    perm
    c=0;
    for s=1:100
        c=c+1;
        load(['/data/william/Projects/dfc_timeseries/' ds '/permgraph/s' num2str(s) '_' num2str(perm)])
        graphfreqsub(:,:,1:78,s)=graph(:,:,1:78); 
    end

    %Stadnardise using mean and std
    clear z_graphfreqsub norm_graphfreqsub
    for s=1:size(graphfreqsub,4)
        tmp1 = graphfreqsub(:,:,:,s);
        for f=1:78
            tmp = reshape(tmp1(:,:,f),264*264,1);
            ztmp=(tmp1(:,:,f)-mean(tmp))./std(tmp);
            ztmp= reshape(ztmp,264,264); 
            z_graphfreqsub(:,:,f,s)=ztmp;
        end
    end

    for r=1:length(dellocal)
        z_graphfreqsub(r,dellocal{r},:,:)=0;
        z_graphfreqsub(r,r,:,:)=0;
    end

    permgraphfreq(:,:,:,perm)=mean(z_graphfreqsub,4); 
end


%% ImagPhase Perm

phpermgraphfreq=zeros([size(graphfreq) 100]);
%Permutation load
clear graphfreqsub
for perm=1:100
    %wavelet
    perm
    graphfreqsub=zeros(264,264,78,100);
    for s=1:100
        load(['/data/william/Projects/dfc_timeseries/' ds '/permgraph_phase/s' num2str(s) '_' num2str(perm)])
        graphfreqsub(:,:,1:78,s)=graph(:,:,1:78); 
    end

    %Stadnardise using mean and std
    clear z_graphfreqsub norm_graphfreqsub
    parfor s=1:size(graphfreqsub,4)
        tmp1 = graphfreqsub(:,:,:,s);
        for f=1:78
            tmp = reshape(tmp1(:,:,f),264*264,1);
            ztmp=(tmp1(:,:,f)-mean(tmp))./std(tmp);
            ztmp= reshape(ztmp,264,264); 
            z_graphfreqsub(:,:,f,s)=ztmp;
        end
    end

    for r=1:length(dellocal)
        z_graphfreqsub(r,dellocal{r},:,:)=0;
        z_graphfreqsub(r,r,:,:)=0;
    end

    phpermgraphfreq(:,:,:,perm)=mean(z_graphfreqsub,4); 
end
    
%% Scale each group matrix between 0 and 1
clear graphfreqsub permgraphfreqsub
for f=1:78
    for p=1:100
    permgraphfreqsub(:,:,f,p) = permgraphfreq(:,:,f,p)./(max(max(abs(permgraphfreq(:,:,f,p)))));
    phpermgraphfreqsub(:,:,f,p) = phpermgraphfreq(:,:,f,p)./(max(max(abs(phpermgraphfreq(:,:,f,p)))));
    end
    graphfreqsub(:,:,f) = graphfreq(:,:,f)./(max(max(abs(graphfreq(:,:,f))))); %Scale between -1 and 1. 
end



save([mdir '/' ds '/graphdata1'],'permgraphfreqsub','-v7.3')
save([mdir '/' ds '/graphdata2'],'graphfreqsub','-v7.3')
save([mdir '/' ds '/graphdata3'],'phpermgraphfreqsub','-v7.3')