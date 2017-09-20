%% This file takes human connectome files, loads them - detrends and then performs
%% extracts the ROI information (in roi_264_powerneuron.m)
%% scrubbing using the framewise displacement
%% require marsbar and bramila_framewiseDisplacement() which can be found https://users.aalto.fi/~eglerean/bramila.html

%% Parameters

clear all
mdir = '/data/william/Projects/dfc_timeseries/' %Working directory
cd(mdir)
roi_264_powerneuron
sdir = [mdir 'timeseries/LR'];
ldir = ['/data/william/datasets/LR'];
sphere_radius=10;
FDThreshold = .5; %Threshold for scrubbing
S = 100 % NumberOfSubjects




%% Set up

if exist(sdir)==0, mkdir(sdir); end
cd(mdir)
roi_264_powerneuron
addpath('/data/william/Projects/dfc_timeseries/')
addpath(genpath('/data/william//toolbox/marsbar-0.43'))
addpath('/data/william//toolbox/Bramila/')
rmpath('/data/william//toolbox/marsbar-0.43/spm99/')
addpath('/usr/local/spm8/')

clear sphere
for n=1:size(coord,1)
    cfg.centre = coord(n,:)
    cfg.radius = sphere_radius;
    sphere{n} = maroi_sphere(cfg);
end

%% Go through each subject, detrend and then extract rois
%45 63
for s=1:100

    fname = [ldir '/' num2str(s) '.nii']
    dat=read_nifti(fname);

    cfg=[];
    cfg.infile = fname;
    cfg.vol = dat;
    cfg.write = 0;
    vol=bramila_detrend(cfg);

    ainfo=spm_vol([ldir '/' num2str(s) '.nii']);

    nn=0;
    timeseries=[];
    %Done in steps of 400 due to hardcoded SPM problems
    for c=1:400:801
        clear volinfo

        for n=1:400
            nn=nn+1;
            volinfo=ainfo(n);
            volinfo=rmfield(volinfo,'pinfo');
            volinfo.fname=[pwd '/tmp.nii'];
            volinfo.private.dat.fname=[pwd '/tmp.nii'];
%             volinfo.dim(4)=1;
            spm_write_vol(volinfo,vol(:,:,:,nn));
        end
        fdat = spm_select('ExtList', pwd,['tmp.nii'],1:400) ;
        %Extract the spheres from marsbar
        tmp = get_marsy(sphere{:},fdat,'mean','v');
        %Build the timeseries variable in steps of 400
        timeseries(end+1:end+400,:) = summary_data(tmp,'mean');
        !rm tmp.nii
    end


    save([sdir '/ts' num2str(s)],'timeseries')

end




%% load the timeseries from the previous cell and perform scrubing at FDThreshold
%% INterpolate with cubic spline

for s=1:S

    s

    mv=load([ldir '/movement/' num2str(s) '.txt']); %Orginally named Motion_Regressors.txt
    cfg=[];
    cfg.motionparam = mv(:,[1:6]); %Use first 6 motion regressors.
    fwd=bramila_framewiseDisplacement(cfg); 

    load(([sdir '/ts' num2str(s)]))

    clear b
    badframes = find(fwd>FDThreshold)
    timeseres(badframes,:) = NaN;
    numbadframes(s)=length(badframes);
    for r=1:264
        t=timeseries(:,r);
        timeseries(:,r)=spline(1:length(t),t,1:length(t)); %Interpolate bad frames with a 3order spline (similar to Allen et al 2014)
    end
    save([sdir 'fwd_ts' num2str(s)],'timeseries')

end
