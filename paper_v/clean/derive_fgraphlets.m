%Start Frieldtrip

ft_defaults

%%Parameters 

foilim = [0.01 0.1] %Frequencies of interest
TR = 0.72
Fs=1/TR; #Sampling rate
roiNumber = 264 

%% Load demeaned/standerdized data in roi x time matrix
    
load([YOURDATA]) 
data=YOURDATA; 

%% Create fieldtrip structure
for roi=1:roiNumber
    dat.label{roi,1}=num2str(roi);
end
dat.time{1} = 1/Fs:1/Fs:(size(data,1)/Fs); 
dat.trial{1}=data; 
dat.fsample = Fs;

%% Use Fieldtip function 
cfg=[]; 
cfg.method = 'wavelet' %I prefer mtmconvol but then the parameters are a little different. But this is what I used in the paper 
cfg.output = 'pow'; 
cfg.toi=dat.time{1}(1):dat.time{1}(1):dat.time{1}(end);
cfg.foilim=foilim;
cfg.width =  6 % Number of cycles
cfg.gwidth = 3 % SIze of wavelet 
cfg.polyremovel = 0 
pow = ft_freqanalysis(cfg,dat); 
