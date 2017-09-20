
addpath('/data/william/Projects/dfc_timeseries/')

%% Parameters

foilim = [0.01 0.1] %Frequencies of interest
TR = 0.72

mdir=['/data/william/Projects/dfc_timeseries/']

ds = 'LR' %Blank for RL, LR for LR.

%% F-graphlet derivation.


powernetclassclusters

for s=1:100


    load([mdir 'timeseries/' ds 'fwd_ts' num2str(s)])
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
    cfg.foilim=foilim;
    cfg.width =  6
    cfg.gwidth = 3
    cfg.polyremovel = 0
    pow = ft_freqanalysis(cfg,dat);
    sdir=[mdir ds '/power_tfr/s' num2str(s)]
    mkdir(sdir)
%     save(sdir,'pow')

    %Wavelet
    clear graph
%     load([mdir '/power_tfr/s' num2str(s)])
    for f=1:length(pow.freq)
        fid=find(isnan(pow.powspctrm(2,f,:))==0);
        graph(:,:,f)=corr(squeeze(pow.powspctrm(:,f,fid))',squeeze(pow.powspctrm(:,f,fid))','type','Spearman');
    end

    save([mdir '/graph/' ds '/s' num2str(s)],'graph')

%     Below code is used to print a figure used in figure 1
    if s==10
        addpath(genpath('/data/william/functions/cbrewer/'))
        cmap=cbrewer('div','RdBu',15)
        cmap=flipud(cmap);
        h=figure('Position',[100 700 800 300])
        imagesc(pow.time,pow.freq,squeeze(pow.powspctrm(120,:,:)))
        colormap([repmat([1 1 1],300,1); hot(300)]);
        axis xy
        caxis([-15000 15000])
        colorbar
        print(h,[mdir '/figures/tfr_ex1'],'-depsc','-r300')

        h=figure('Position',[100 700 800 300])
        imagesc(pow.time,pow.freq,squeeze(pow.powspctrm(180,:,:)))
        colormap([repmat([1 1 1],300,1); hot(300)]);
        axis xy
        caxis([-15000 15000])
        colorbar
        print(h,[mdir '/figures/tfr_ex2'],'-depsc','-r300')

         f=70
                    h=figure

        imagesc(graph(PlotOrder,PlotOrder,f))
        caxis([-.8 .8])
        colormap(cmap)

        print(h,[mdir '/figures/conmat_ex'],'-depsc','-r300')

    end


end


%% permute imaginary part of fourier


clear dat pow
for s=1:100
    s

    load([mdir 'timeseries/' ds 'fwd_ts' num2str(s)])
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


        s

        dat.trial{1} = dat.trial{1}(:,:);

        cfg=[];
        cfg.method = 'wavelet'
        cfg.toi=dat.time{1}(1):dat.time{1}(1):dat.time{1}(end);
        cfg.foilim=[0.01 .1];
        cfg.width =  6;
        cfg.gwidth = 3;
        cfg.polyremovel = 0;
        cfg.output = 'fourier'
        freq = ft_freqanalysis(cfg,dat);


        %Wavelet
    for perm=1:100
        clear graph
        fq = freq;
        for f=1:78

            fid=find(isnan(fq.fourierspctrm(1,2,f,:))==0);
            [~, r] = sort(rand(length(fid),1));
            irand=imag(fq.fourierspctrm(:,:,f,fid(r)))*i;
            realnonrand=real(fq.fourierspctrm(:,:,f,fid));
            fq.powspctrm(:,f,fid)=squeeze(abs(irand+realnonrand).^2);
%             irande=imag(fq.fourierspctrm(:,:,f,fid))*i;
%             fq.emppowspctrm(:,f,fid)=squeeze(abs(irande+realnonrand).^2);
%      R=abs(Z);
%             P=angle(Z);
%             Z2=R.*exp(i*P(r));
%             R2=abs(Z2);
        end
        for f=1:78
            fid=find(isnan(fq.fourierspctrm(1,2,f,:))==0);
            graph(:,:,f)=corr(squeeze(fq.powspctrm(:,f,fid))',squeeze(fq.powspctrm(:,f,fid))','type','Spearman');
        end

        save(['/data/william/Projects/dfc_timeseries/' ds '/permgraph_phase/s' num2str(s) '_' num2str(perm)],'graph')

    %     Below code is used to print corresponding to figure 1, but not used
%         in text
% load(['/data/william/Projects/dfc_timeseries/' ds '/permgraph_phase/s' num2str(s) '_' num2str(perm)])
        if s==10
            addpath(genpath('/data/william/functions/cbrewer/'))
            cmap=cbrewer('div','RdBu',15)
            cmap=flipud(cmap);
            h=figure('Position',[100 700 800 300])
            imagesc(fq.time,fq.freq,squeeze(fq.powspctrm(120,:,:)))
            colormap([repmat([1 1 1],300,1); hot(300)]);
            axis xy
            caxis([-15000 15000])
            colorbar
            print(h,[mdir '/figures/tfr_ex1_imagsshuffle'],'-depsc','-r300')

            h=figure('Position',[100 700 800 300])
            imagesc(fq.time,fq.freq,squeeze(fq.powspctrm(180,:,:)))
            colormap([repmat([1 1 1],300,1); hot(300)]);
            axis xy
            caxis([-15000 15000])
            colorbar
            print(h,[mdir '/figures/tfr_ex2_imagsshuffle'],'-depsc','-r300')

             f=70
                        h=figure
            imagesc(graph(PlotOrder,PlotOrder,f))
            caxis([-.8 .8])
            colormap(cmap)
            axis square
            colorbar
            print(h,[mdir '/figures/conmat_ex_imagsshuffle'],'-depsc','-r300')

        end
    end
end



%%Ampltiude based putmutations



clear dat pow
for s=1:100
    s

    load([mdir 'timeseries/' ds 'fwd_ts' num2str(s)])
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


    for perm=1:100
        s
        [o RndOdr]=sort(rand(length(dat.trial{1}(1,:)),1));

        dat.trial{1} = dat.trial{1}(:,RndOdr);

        cfg=[];
        cfg.method = 'wavelet'
        cfg.output = 'pow';
        cfg.toi=dat.time{1}(1):dat.time{1}(1):dat.time{1}(end);
        cfg.foilim=[0.01 .1];
        cfg.width =  6;
        cfg.gwidth = 3;
        cfg.polyremovel = 0;
        pow = ft_freqanalysis(cfg,dat);

        %Wavelet
        for f=1:78
            fid=find(isnan(pow.powspctrm(2,f,:))==0);
            graph(:,:,f)=corr(squeeze(pow.powspctrm(:,f,fid))',squeeze(pow.powspctrm(:,f,fid))','type','Spearman');
        end

        save(['/data/william/Projects/dfc_timeseries/' ds '/permgraph/s' num2str(s) '_' num2str(perm)],'graph')

    %     Below code is used to print corresponding to figure 1, but not used
%         in text
load(['/data/william/Projects/dfc_timeseries/' ds '/permgraph/s' num2str(s) '_' num2str(perm)])
        if s==10
            addpath(genpath('/data/william/functions/cbrewer/'))
            cmap=cbrewer('div','RdBu',15)
            cmap=flipud(cmap);
            h=figure('Position',[100 700 800 300])
            imagesc(pow.time,pow.freq,squeeze(pow.powspctrm(120,:,:)))
            colormap([repmat([1 1 1],300,1); hot(300)]);
            axis xy
            caxis([-15000 15000])
            colorbar
            print(h,[mdir '/figures/tfr_ex1_timesshuffle'],'-depsc','-r300')

            h=figure('Position',[100 700 800 300])
            imagesc(pow.time,pow.freq,squeeze(pow.powspctrm(180,:,:)))
            colormap([repmat([1 1 1],300,1); hot(300)]);
            axis xy
            caxis([-15000 15000])
            colorbar
            print(h,[mdir '/figures/tfr_ex2_timesshuffle'],'-depsc','-r300')

             f=70
                        h=figure
            imagesc(graph(PlotOrder,PlotOrder,f))
            caxis([-.8 .8])
            colormap(cmap)
            axis square
            colorbar
            print(h,[mdir '/figures/conmat_ex_timesshuffle'],'-depsc','-r300')
        end

    end
end
