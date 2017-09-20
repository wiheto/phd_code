
%This code calculates the efficiency and betweenness centrality per frequency bin


%% Parameters
%Requires mdir
th=0.05
ds = 'LR' 
addpath('/data/william/toolbox/bc')


%Load data from other files
load([mdir '/graph/foi'])
freq=freq(1:78)

load([mdir '/' ds '/graphdata1'])
load([mdir '/' ds '/graphdata2'])
load([mdir '/' ds '/graphdata3'])

%Threshold at 0.05 or 0.09
clear thgraph ththgraph phthgraph

%Threshold models 
%Null models
for p=1:100
    p
    for f=1:size(graphfreqsub,3)
        ththgraph(:,:,f,p)=threshold_proportional(permgraphfreqsub(:,:,f,p),th);
        phththgraph(:,:,f,p)=threshold_proportional(phpermgraphfreqsub(:,:,f,p),th);
    end
end

%Empirical data
clear thgraph
for f=1:size(graphfreqsub,3)
    thgraph(:,:,f)=threshold_proportional(graphfreqsub(:,:,f),th);
end

%% Efficiency
clear nmge ge
% Null model
for p=1:100
    p
parfor f=1:78
nmge(f,p) = efficiency_wei(ththgraph(:,:,f,p));
phnmge(f,p) = efficiency_wei(phththgraph(:,:,f,p));
end
end
phnmge=sort(phnmge,2);
nmge=sort(nmge,2);

%Empirical data
parfor f=1:78
    f
    ge(f)=efficiency_wei(thgraph(:,:,f));
end


%% Betweeness 
clear thbet phthbet
%Null Models
for p=1:200
    p
    parfor f=1:78
        thbet(:,f,p)=betweenness_wei(ththgraph(:,:,f,p));       
        phthbet(:,f,p)=betweenness_wei(phththgraph(:,:,f,p));       
    end
end
phthbet=sort(phnmge,3);
thbet=sort(nmge,3);

%Empirical 
clear bet
for f=1:78
    f
    bet(:,f)=betweenness_wei(thgraph(:,:,f));       
end

if th==0.05
    save(['/data/william/Projects/dfc_timeseries/' ds 'powgraphbetl_individualscaling' ],'bet','thbet','ge','nmge','thgraph','ththgraph','graphfreqsub','permgraphfreqsub','phththgraph','-v7.3')    
else
    save(['/data/william/Projects/dfc_timeseries/' ds 'powgraphbetl_individualscalingph_9' ],'bet','phthbet','ge','thbet','ththgraph','nmge','phnmge','thgraph','phththgraph','-v7.3')
end







%Whole spectra done by request of reviewers
for s=1:100
    
    
    load([mdir 'timeseries/fwd_ts' num2str(s)])
    for roi=1:264
        dat.label{roi,1}=num2str(roi);
    end

    TR=0.72;
    Fs=1/TR;

    dat.time{1} = 1/Fs:1/Fs:(size(timeseries,1)/Fs); 
    for r=1:264
    dat.trial{1}(r,:)=timeseries(:,r)'-mean(timeseries(:,r));
    end
    dat.fsample = Fs;
    
    cfg=[]; 
    cfg.method = 'wavelet'
    cfg.output = 'pow'; 
    cfg.toi=dat.time{1}(1):dat.time{1}(1):dat.time{1}(end);
    cfg.width =  6
    cfg.gwidth = 3
    cfg.polyremovel = 0
    pow = ft_freqanalysis(cfg,dat); 
    sdir=[mdir '/power_tfr/s' num2str(s)]
    mkdir(sdir)
%     save(sdir,'pow')

    %Wavelet
%     load([mdir '/power_tfr/s' num2str(s)])
    for f=1:length(pow.freq)
        fid=find(isnan(pow.powspctrm(2,f,:))==0);
        if isempty(fid)==1
            graph(1:264,1:264,f,s)=NaN;
        else
        graph(:,:,f,s)=corr(squeeze(pow.powspctrm(:,f,fid))',squeeze(pow.powspctrm(:,f,fid))','type','Spearman');
        end
    end
 

    
end

