%% Little outdate file that I once used on HCP data

clear all
mdir = '/data/william/Projects/dfc_timeseries/'
cd(mdir)

addpath(genpath('/data/william//toolbox/bc/'))

addpath(genpath('/data/william//toolbox/marsbar-0.43'))
mkdir('roi/')
%Create 264 roi from Power et al
roi_264_powerneuron
addpath('/data/william//toolbox/Bramila/')
rmpath('/data/william//toolbox/marsbar-0.43/spm99/')    
clear sphere
sphere_radius=5; 
for n=1:size(coord,1)   
    cfg.centre = coord(n,:)
    cfg.radius = sphere_radius; 
    sphere{n} = maroi_sphere(cfg);
end
for b=2:10

    addpath(['1RL/b' num2str(b)])
   
    d=dir(['1RL/b' num2str(b)])
    
for s=3:length(d) 
    
    ss=10*(b-1)+s-2;
      
    fname = [d(s).name]
    dat=read_nifti(fname); 
    
    cfg=[];
    cfg.infile = fname; 
    cfg.vol = dat; 
    cfg.write = 0; 
    vol=bramila_detrend(cfg);
    
    ainfo=spm_vol([mdir '1RL/b' num2str(b) '/' fname]);

    nn=0; 
    timeseries=[];
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
        tmp = get_marsy(sphere{:},fdat,'mean','v');
        timeseries(end+1:end+400,:) = summary_data(tmp,'mean');     
        !rm tmp.nii
    end
    
    cfg=[];
    TR=0.72
    cfg.ts = timeseries; 
    cfg.plot = 1; 
    dvars = bramila_dvars(cfg);
    
    save([mdir 'timeseries/ts' num2str(ss)],'timeseries','dvars')
    
    
    setenv('RMFILE',[mdir '1RL/b' num2str(b) '/' fname]);
    !rm $RMFILE
    
end

end