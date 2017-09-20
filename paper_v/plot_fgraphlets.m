%With the pre and postprocessing of f-graphlets in mainscript, this script. Also incldues the str contribution calculations. 
% plots the figures used in the paper. 
% This code is quite long and poorly commented (Sorry). 
% Was first phd project and additions asked in reviewer got added a non-linearly into the code. 

% Requires: schemaball, sigstar, brain connectivity toolbox, brainnet viewer 

clear all
mdir = '/data/william/Projects/dfc_timeseries/'
cd(mdir)

addpath('/data/william/functions/sigstar/')
addpath(genpath('/data/william/toolbox/bc/'))

addpath(genpath('/data/william/toolbox/marsbar-0.43'))
addpath(genpath(['/data/william/akalla/william/functions/']))
addpath(genpath('/data/william/toolbox/brainviewer'))
addpath('/data/william/toolbox/Bramila/')
addpath('/data/william/toolbox/exportfig/')

%Create 264 roi from Power et al

roi_264_powerneuron
powernetclassclusters



load(['/data/william/Projects/dfc_timeseries/powgraphbetl_individualscaling' ])    
load(['/data/william/Projects/dfc_timeseries/powgraphbetl_individualscalingph' ])    

load([mdir '/graph/foi'])
freq=freq(1:78)


%% Figure 1


h=figure
plot(freq,squeeze(sum(sum(thgraph))),'k')

print(h,[mdir '/figures/strength'],'-depsc','-r300')

nmge = sort(nmge,2) 
h=figure
plot(freq(1:size(ge,2)),ge,'k')
    hold on
plot(freq(1:size(ge,2)),nmge(:,99),'k--')
plot(freq(1:size(ge,2)),phnmge(:,99),'k:')
print(h,[mdir '/figures/global_efficiency'],'-depsc','-r300')


h=figure
plot(freq(1:size(ge,2)),ge./mean(nmge,2)','k')
hold on
plot(freq,nmge(:,99)./mean(nmge,2),'k--')
mc(1:78)=max(nmge(:,99)./mean(nmge,2)); 
plot(freq,mc,'k-.')


print(h,[mdir '/figures/normalized_global_efficiency'],'-depsc','-r300')

h=figure
plot(freq(1:size(ge,2)),ge./mean(phnmge,2)','k')
hold on
plot(freq,phnmge(:,99)./mean(phnmge,2),'k--')
mc(1:78)=max(phnmge(:,99)./mean(phnmge,2)); 
plot(freq,mc,'k-.')

%[6 21 57]
sigstar([freq(5) freq(7)])
sigstar([freq(16) freq(24)])
sigstar([freq(53) freq(62)])

print(h,[mdir '/figures/normalized_global_phefficiency'],'-depsc','-r300')



corder = [4 6 9 10 7 2 1 5 3 8];
u=(1/255)
cm=[u*128 0 0;
    0 0 u*128; 
    0 u*255 u*255; 
    0 u*128 u;
    u u*255 0;
    0 0 0; 
    u*255 u*128 0; 
    u*255 u*128 u*128; 
    u*128 0 u*128; 
    u*255 0 u*255]



sththgraph = mean(ththgraph,4); 
    
c=0;
foi=1:78;
th = 195
h=figure('Position',[100 600 1000 400])
for n=corder(1:end-1)
    c=c+1
    subplot(2,9,c) 
    hold on
    fid=find(PowerNetClass==PowerNetLabel{n,1})
    NET{n} = fid
    plot(freq(foi),squeeze(sum(sum(thgraph(NET{n},NET{n},foi)))),'Color',cm(c,:))
    %plot(freq(foi),squeeze(sum(sum(null_outpers(NET{n},NET{n},foi,th)))),'--','Color',cm(c,:))
    subplot(2,9,c+9)
    NOTNET{n} = 1:264; NOTNET{n}(fid)=[]; 
    hold on
    plot(freq(foi),squeeze(sum(sum(thgraph(NET{n},NOTNET{n},foi)))),'Color',cm(c,:))
%     plot(freq(foi),squeeze(sum(sum(sththgraph(NET{n},NOTNET{n},foi,th)))),'--','Color',cm(c,:))

end

print(h,[mdir 'figures/degree_connections_5'],'-r300','-depsc')

%Example group connectivity matrices
for f=[1 2 78]
    h=figure
    imagesc(graphfreqsub(:,:,f)) 
    caxis([-1 1])
    print(h,[mdir '/figures/example_groupconnectivity_' num2str(f)],'-r300','-depsc')
end



%degree differences

h=figure('Position',[100 100 1000 200])

c=0;
for n=corder(1:end-1)
    c=c+1
    subplot(1,9,c) 
    NOTNET{n} = 1:264; NOTNET{n}(fid)=[]; 
    plot(freq(foi),squeeze(sum(sum(thgraph(NET{6},NET{n},foi)))),'Color',cm(c,:))
end

print(h,[mdir 'figures/fp_degree_connections_5'],'-r300','-depsc')

h=figure('Position',[100 100 1000 200])

c=0;
for n=corder(1:end-1)
    c=c+1
    subplot(1,9,c) 
    NOTNET{n} = 1:264; NOTNET{n}(fid)=[]; 
    plot(freq(foi),squeeze(sum(sum(thgraph(NET{5},NET{n},foi)))),'Color',cm(c,:))
end


print(h,[mdir 'figures/vis_degree_connections_5'],'-r300','-depsc')

h=figure('Position',[100 100 1000 200])

c=0;
for n=corder(1:end-1)
    c=c+1
    subplot(1,9,c) 
    NOTNET{n} = 1:264; NOTNET{n}(fid)=[]; 
    plot(freq(foi),squeeze(sum(sum(thgraph(NET{4},NET{n},foi)))),'Color',cm(c,:))
end


print(h,[mdir 'figures/dmn_degree_connections_5'],'-r300','-depsc')


%Betweenness 

nbet=bet./((size(bet,1)-1)*(size(bet,1)-2)/2);
nthbet=(thbet)./(((size((thbet),1)-1)*(size((thbet),1)-2))/2);

for i=1:264
    for f=1:78
s=sort(squeeze(((nthbet(i,f,:))))); 
fid=find(nbet(i,f)>s(99));
if length(fid)==1
sbet(i,f)=nbet(i,f);
else 
    sbet(i,f)=0; 
    end
    end
end
figure
imagesc(sbet)

clear meanbet betpercluster nbetpercluster
for n=1:10
    fid=find( PowerNetClass==PowerNetLabel{n,1})
    betpercluster{n}=sbet(fid,:);
end

for n=1:10
    figure    
    imagesc(betpercluster{n})%-nbetpercluster{n})
    caxis([0 0.2])
    colormap('hot')
    print([mdir '/figures/bet_net' num2str(n)],'-r300','-depsc')
end
close all

foi=[6 21 57];
id = 4; 

bet_foi(:,1)=mean(sbet(:,foi(1)-id:foi(1)+id),2);
bet_foi(:,2)=mean(sbet(:,foi(2)-id:foi(2)+id),2);
bet_foi(:,3)=mean(sbet(:,foi(3)-id:foi(3)+id),2);
nbet_foi(:,1)=mean(nbet(:,foi(1)-id:foi(1)+id),2);
nbet_foi(:,2)=mean(nbet(:,foi(2)-id:foi(2)+id),2);
nbet_foi(:,3)=mean(nbet(:,foi(3)-id:foi(3)+id),2);

save([mdir '/figures/bet4blobs'],'bet_foi')
save([mdir '/figures/nbetsave'],'nbet_foi')


clear p
for f=1:78
   P(:,f)=participation_coef(thgraph(:,:,f),PowerNetClass);
end
for f=1:3
    for n=1:10
        fid=find(PowerNetClass==PowerNetLabel{n,1});
        Pm(n,f)=mean(mean(P(fid,foi(f)-id:foi(f)+id)));
        Pn{n,f} = mean(P(fid,foi(f)-id:foi(f)+id),2) 
    end
end
clear pv1 pv2 pPm v1 v2 emp
for f=1:3
    f
    for ff=1:3
        emp(:,f,ff) = Pm(:,f)-Pm(:,ff); 
    for p=1:1000
        [a s]=sort(rand(264,2),2);    
        for n=1:264
            if s(n,1)==1
                v1(n,1)=mean(P(n,foi(f)-id:foi(f)+id));
                v2(n,1)=mean(P(n,foi(ff)-id:foi(ff)+id));
            elseif s(n,1)==2
                v1(n,1)=mean(P(n,foi(ff)-id:foi(ff)+id));
                v2(n,1)=mean(P(n,foi(f)-id:foi(f)+id));
            end
        end   
        for n=1:10
            fid=find(PowerNetClass==PowerNetLabel{n,1});
            pv1(n)=(mean(v1(fid)));
            pv2(n)=(mean(v2(fid)));
        end
        pPm(:,f,ff,p)=pv1-pv2; 
    end
    end
end
sig = zeros(3,3,10)
pPm = sort(pPm,4); 
for f=1:3
    for ff=1:3
        for n=1:10
            if emp(n,f,ff)>pPm(n,f,ff,992) || emp(n,f,ff)<pPm(n,f,ff,8)
                sig(f,ff,n) = 1
            end
        end
    end
end
clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a o]=sort(fid,'descend') 
h=figure('Position',[100,700,800,400])
Pm(:,4) = NaN; 
e=reshape(Pm(o,:)', size(Pm,2)*size(Pm,1),1);
bar(e)
x=1; 


for n=1:10
    fid=find(sig(:,:,o(n))==1);
    stars = {}; 
    if sum(ismember(fid,4))==1;
        sigstar([x x+1]);
    end
    if sum(ismember(fid,8))==1;
        sigstar([x+1, x+2]);
    end
    if sum(ismember(fid,7))==1;
        sigstar([x, x+2]);
    end    
    x=x+4;
end
       



%PhBetweenness 

nbet=bet./((size(bet,1)-1)*(size(bet,1)-2)/2);
nthbet=(phthbet)./(((size((phthbet),1)-1)*(size((phthbet),1)-2))/2);

for i=1:264
    for f=1:78
s=sort(squeeze(((nthbet(i,f,:))))); 
fid=find(nbet(i,f)>s(99));
if length(fid)==1
sbet(i,f)=nbet(i,f);
else 
    sbet(i,f)=0; 
    end
    end
end
figure
imagesc(sbet)

clear meanbet betpercluster nbetpercluster
for n=1:10
    fid=find( PowerNetClass==PowerNetLabel{n,1})
    betpercluster{n}=sbet(fid,:);
end

for n=1:10
    figure    
    imagesc(betpercluster{n})%-nbetpercluster{n})
    caxis([0 0.2])
    colormap('hot')
    print([mdir '/figures/bet_net' num2str(n)],'-r300','-depsc')
end
close all

bi{1}=6;
bi{2}=16:24;
bi{3}=53:61;

bet_foi(:,1)=mean(sbet(:,bi{1}),2);
bet_foi(:,2)=mean(sbet(:,bi{2}),2);
bet_foi(:,3)=mean(sbet(:,bi{3}),2);
nbet_foi(:,1)=mean(nbet(:,bi{1}),2);
nbet_foi(:,2)=mean(nbet(:,bi{2}),2);
nbet_foi(:,3)=mean(nbet(:,bi{3}),2);

save([mdir '/figures/bet4blobs'],'bet_foi')
save([mdir '/figures/nbetsave'],'nbet_foi')


cmap=cbrewer('div','RdBu',15)
cmap=flipud(cmap);

figure
c=0;
for n=[1 2 3 78]
    c=c+1
subplot(4,1,c)
    imagesc(graphfreqsub(PlotOrder,PlotOrder,n))
caxis([-.5 .5])
colormap(cmap)
end
print([mdir '/figures/fgraphlet_example_group'],'-r300','-depsc')


%% Network-Network Interactions 





th = 0.05;
thtext = '5'
clear dat

bi{1}=6;
bi{2}=16:24;
bi{3}=53:61;

for f=1:78
    dat(:,:,f)=threshold_proportional(mean(graphfreqsub(:,:,f),3),th); 
%     dat(:,:,f)=thgraph(:,:,f);
    dat(:,:,f) = 2.*(dat(:,:,f)./sum(sum(dat(:,:,f))));
end


clear permdif pv1 pv2 
f1id=0;
for f1=1:3
    f1id=f1id+1
    f2id=0;
for f2=1:3
    f2id=f2id+1; 
    if f1>=f2
    else
    for p=1:1000
        clear perm1 perm2
        for n=1:10
            fid=find(PowerNetClass==PowerNetLabel{n,1});
            
            v1=mean(dat(fid,fid,bi{f1}),3); 
            v2=mean(dat(fid,fid,bi{f2}),3);
                    
            [a r] = sort(rand(length(fid),length(fid),2),3);
            clear pv1 pv2
            for a=1:length(fid) 
            for b=1:length(fid)
                if a<=b
                if r(a,b,1)==1
                    pv1(a,b) = v1(a,b);
                    pv2(a,b) = v2(a,b);
                    pv1(b,a) = v1(a,b);
                    pv2(b,a) = v2(a,b);
                else
                    pv1(a,b) = v2(a,b);
                    pv2(a,b) = v1(a,b);
                    pv1(b,a) = v2(a,b);
                    pv2(b,a) = v1(a,b);
                end
                end
            end
            end  
            perm1(n) = sum(sum(pv1))./2; 
            perm2(n) = sum(sum(pv2))./2;
        end
        permdif(f1id,f2id,:,p) = (perm1-perm2);  
    end
    end
end
end

clear emp emps emp1 emp2
f1id=0; 
for f1=1:3
    f1id = f1id+1; 
    f2id=0;     
for f2=1:3
    f2id = f2id+1;     
    c=c+1
    for n=1:10
        fid=find(PowerNetClass==PowerNetLabel{n,1});
        emp1(n) = sum(sum(mean(dat(fid,fid,bi{f1}),3)))./2;
        emp2(n) = sum(sum(mean(dat(fid,fid,bi{f2}),3)))./2;        
    end
    emp(f1id,f2id,:) = emp1-emp2; 
    emps(:,f1id) = emp1;    
end
end

sig=zeros(3,3,10);
for f1id=1:3
for f2id=1:3
    if f1id>=f2id
    else
        for n=1:10
            s=sort(squeeze(permdif(f1id,f2id,n,:)));
            if s(992) < emp(f1id,f2id,n) || s(8) > emp(f1id,f2id,n)      
               sig(f1id,f2id,n) = 1;  
            end
        end
    end
end
end

clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a o]=sort(fid,'descend') 

h=figure('Position',[100,700,800,400])
emps(:,4) = NaN; 
e=reshape(emps(o,:)', size(emps,2)*size(emps,1),1);
bar(e)
x=1; 
for n=1:10
    fid=find(sig(:,:,o(n)));
    stars = {}; 
    if sum(ismember(fid,4))==1;
        sigstar([x x+1]);
    end
    if sum(ismember(fid,8))==1;
        sigstar([x+1, x+2]);
    end
    if sum(ismember(fid,7))==1;
        sigstar([x, x+2]);
    end    
    x=x+4;
end
       
        
sigstar([1.1,1.2])
set(gca,'XTick',2:4:40) 
set(gca,'XTickLabel',{PowerNetLabel{o,2}}') 


print(h,[mdir '/figures/withinwork_' thtext '_.eps'],'-depsc','-r300')

clear permdif perm1 perm2


clear permdif pv1 pv2 
f1id=0;
for f1=1:3
    f1id=f1id+1
    f2id=0;
for f2=1:3
    f2id=f2id+1; 
    if f1>=f2
    else
    for p=1:1000
        clear perm1 perm2
        for n=1:10
            fid=find(PowerNetClass==PowerNetLabel{n,1});
            
            nfid = 1:264;
            nfid(fid) = []; 
            v1=mean(dat(fid,nfid,bi{f1}),3); 
            v2=mean(dat(fid,nfid,bi{f2}),3);
            
            [a r] = sort(rand(length(fid),length(nfid),2),3);
            clear pv1 pv2
            for a=1:length(fid) 
            for b=1:length(nfid) 
                if r(a,b,1)==1
                    pv1(a,b) = v1(a,b);
                    pv2(a,b) = v2(a,b);
                else
                    pv1(a,b) = v2(a,b);
                    pv2(a,b) = v1(a,b);
                end
            end
            end  
            perm1(n) = sum(sum(pv1)); 
            perm2(n) = sum(sum(pv2));
        end
        permdif(f1id,f2id,:,p) = perm1-perm2;  
    end
    end
end
end

clear emp emps emp1 emp2
f1id=0; 
for f1=1:3
    f1id = f1id+1; 
    f2id=0;     
for f2=1:3
    f2id = f2id+1;     
    for n=1:10
        fid=find(PowerNetClass==PowerNetLabel{n,1});
        nfid = 1:264;
        nfid(fid) = []; 
        emp1(n) = sum(sum(mean(dat(fid,nfid,bi{f1}),3))); 
        emp2(n) = sum(sum(mean(dat(fid,nfid,bi{f2}),3)));         
    end
    emp(f1id,f2id,:) = emp1-emp2; 
    emps(:,f1id) = emp1;    
end
end

sig=zeros(3,3,10);
for f1id=1:3
for f2id=1:3
    if f1id>=f2id
    else
        for n=1:10   
            s=sort(squeeze(permdif(f1id,f2id,n,:)));
            if s(992) < emp(f1id,f2id,n) || s(8) > emp(f1id,f2id,n)
               sig(f1id,f2id,n) = 1;  
            end  
        end
    end
end
end

clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a o]=sort(fid,'descend') 



h=figure('Position',[100,700,800,400])
emps(:,4) = NaN; 
e=reshape(emps(NetOrder(1:10),:)', size(emps,2)*size(emps,1),1);
bar(e)
x=1; 

for n=1:10
    fid=find(sig(:,:,o(n)));
    stars = {}; 
    if sum(ismember(fid,4))==1;
        sigstar([x x+1]);
    end
    if sum(ismember(fid,8))==1;
        sigstar([x+1, x+2]);
    end
    if sum(ismember(fid,7))==1;
        sigstar([x, x+2]);
    end    
    x=x+4
end
       
        
sigstar([1.1,1.2])
set(gca,'XTick',2:4:40) 
set(gca,'XTickLabel',{PowerNetLabel{o,2}}') 

print(h,[mdir '/figures/betweenwork_' thtext '_.eps'],'-depsc','-r300')




clear permdif perm1 perm2 emps emp1 emp2 e


clear emp emps emp1 emp2
f1id=0; 
for f1=1:3
    f1id = f1id+1; 
    f2id=0;     
for f2=1:3
    f2id = f2id+1;     
    for n=1:10
    for nn=1:10
        fid=find(PowerNetClass==PowerNetLabel{n,1});
        fid2=find(PowerNetClass==PowerNetLabel{nn,1});
        emp1(n,nn) = sum(sum(mean(dat(fid,fid2,bi{f1}),3))); 
        emp2(n,nn) = sum(sum(mean(dat(fid,fid2,bi{f2}),3)));
        num(n,nn,f1id,f2id)=length(find((mean(dat(fid,fid2,bi{f1}),3))));        
    end
    end
    emp(f1id,f2id,:,:) = emp1-emp2; 
    emps(:,:,f1id) = emp1;    
end
end


clear permdif pv1 pv2 
f1id=0;
for f1=1:3
    f1id=f1id+1
    f2id=0;
for f2=1:3
    f2id=f2id+1; 
    for p=1:1000
        clear perm1 perm2 
        for n=1:10
            for nn=1:10
                
            fid=find(PowerNetClass==PowerNetLabel{n,1});
            fid2=find(PowerNetClass==PowerNetLabel{nn,1});
            
            v1=mean(dat(fid,fid2,bi{f1}),3); 
            v2=mean(dat(fid,fid2,bi{f2}),3);
            
            [a r] = sort(rand(length(fid),length(fid2),2),3);
            clear pv1 pv2
            for a=1:length(fid) 
            for b=1:length(fid2) 
                if r(a,b,1)==1
                    pv1(a,b) = v1(a,b);
                    pv2(a,b) = v2(a,b);
                else
                    pv1(a,b) = v2(a,b);
                    pv2(a,b) = v1(a,b);
                end
            end
            end  
            perm1(n,nn) = sum(sum(pv1)); 
            perm2(n,nn) = sum(sum(pv2));
            end
        end
        permdif(f1id,f2id,:,:,p) = perm1-perm2;  
    end
end
end



sig=zeros(3,3,10,10);
for f1id=1:3
for f2id=1:3
    for n=1:10  
        for nn=1:10
            s=sort(squeeze(permdif(f1id,f2id,n,nn,:)));
            if s(992) < emp(f1id,f2id,n,nn) 
               sig(f1id,f2id,n,nn) = 1; 
            end

        end
    end
end
end



clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a o]=sort(fid,'descend') 

conmat=zeros(10,10,3)
for n=1:10 
    for nn=1:10
        if sig(1,2,n,nn)==1 && sig(1,3,n,nn)==1 && max(max(num(n,nn,:,:)))>=1
            conmat(n,nn,1) = 4;

        end
        if sig(2,1,n,nn)==1 && sig(2,3,n,nn)==1 && max(max(num(n,nn,:,:)))>=1
            conmat(n,nn,2) = 4; 

        end
        if sig(3,1,n,nn)==1 && sig(3,2,n,nn)==1 && max(max(num(n,nn,:,:)))>=1
            conmat(n,nn,3) = 4;     
         end
    end
end
      
schemaball({PowerNetLabel{o,2}},conmat(o,o,1)) 
print(gcf,[mdir '/figures/network_to_network_' thtext '_1n.eps'],'-depsc','-r300')

schemaball({PowerNetLabel{o,2}},conmat(o,o,2)) 
print(gcf,[mdir '/figures/network_to_network_' thtext '_2n.eps'],'-depsc','-r300')

schemaball({PowerNetLabel{o,2}},conmat(o,o,3)) 
print(gcf,[mdir '/figures/network_to_network_' thtext '_3n.eps'],'-depsc','-r300')















%% NEGATIVE NETWORKS










thtext = 'abs0'
clear dat

bi{1}=6;
bi{2}=16:24;
bi{3}=53:61;

for f=1:3
    dat(:,:,f)=threshold_absolute(-1.*mean(graphfreqsub(:,:,bi{f}),3),0); 
    dat(:,:,f) = 2.*(dat(:,:,f)./sum(sum(dat(:,:,f))));
end




cmap=cbrewer('div','RdBu',15)
cmap=flipud(cmap);

figure;
subplot(1,3,1)
imagesc(mean(graphfreqsub(PlotOrder,PlotOrder,bi{1}),3))
caxis([-.5 .5])
colormap(cmap)
subplot(1,3,2)
imagesc(mean(graphfreqsub(PlotOrder,PlotOrder,bi{1}),3).*-1)
caxis([-.5 .5])
colormap(cmap)
subplot(1,3,3)
imagesc(threshold_absolute(-1.*graphfreqsub(PlotOrder,PlotOrder,bi{1}),0))
caxis([-.5 .5])
colormap(cmap)
print(gcf,[mdir '/figures/negcorr_example'],'-r300','-depsc')

clear permdif pv1 pv2 
f1id=0;
for f1=1:3
    f1id=f1id+1
    f2id=0;
for f2=1:3
    f2id=f2id+1; 
    if f1>=f2
    else
    for p=1:1000
        clear perm1 perm2
        for n=1:10
            fid=find(PowerNetClass==PowerNetLabel{n,1});
            
            v1=mean(dat(fid,fid,bi{f1}),3); 
            v2=mean(dat(fid,fid,bi{f2}),3);
                    
            [a r] = sort(rand(length(fid),length(fid),2),3);
            clear pv1 pv2
            for a=1:length(fid) 
            for b=1:length(fid)
                if a<=b
                if r(a,b,1)==1
                    pv1(a,b) = v1(a,b);
                    pv2(a,b) = v2(a,b);
                    pv1(b,a) = v1(a,b);
                    pv2(b,a) = v2(a,b);
                else
                    pv1(a,b) = v2(a,b);
                    pv2(a,b) = v1(a,b);
                    pv1(b,a) = v2(a,b);
                    pv2(b,a) = v1(a,b);
                end
                end
            end
            end  
            perm1(n) = sum(sum(pv1))./2; 
            perm2(n) = sum(sum(pv2))./2;
        end
        permdif(f1id,f2id,:,p) = (perm1-perm2);  
    end
    end
end
end

clear emp emps emp1 emp2
f1id=0; 
for f1=1:3
    f1id = f1id+1; 
    f2id=0;     
for f2=1:3
    f2id = f2id+1;     
    c=c+1
    for n=1:10
        fid=find(PowerNetClass==PowerNetLabel{n,1});
        emp1(n) = sum(sum(mean(dat(fid,fid,bi{f1}),3)))./2;
        emp2(n) = sum(sum(mean(dat(fid,fid,bi{f2}),3)))./2;        
    end
    emp(f1id,f2id,:) = emp1-emp2; 
    emps(:,f1id) = emp1;    
end
end

sig=zeros(3,3,10);
for f1id=1:3
for f2id=1:3
    if f1id>=f2id
    else
        for n=1:10
            s=sort(squeeze(permdif(f1id,f2id,n,:)));
            if s(992) < emp(f1id,f2id,n) || s(8) > emp(f1id,f2id,n)      
               sig(f1id,f2id,n) = 1;  
            end
        end
    end
end
end

clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a o]=sort(fid,'descend') 

h=figure('Position',[100,700,800,400])
emps(:,4) = NaN; 
e=reshape(emps(o,:)', size(emps,2)*size(emps,1),1);
bar(e)
x=1; 
for n=1:10
    fid=find(sig(:,:,o(n)));
    stars = {}; 
    if sum(ismember(fid,4))==1;
        sigstar([x x+1]);
    end
    if sum(ismember(fid,8))==1;
        sigstar([x+1, x+2]);
    end
    if sum(ismember(fid,7))==1;
        sigstar([x, x+2]);
    end    
    x=x+4;
end
       
  
set(gca,'XTick',2:4:40)
set(gca,'XTickLabel',{PowerNetLabel{o,2}}') 
legend('0.016', '0.035', '0.075')

print(h,[mdir '/figures/NEGwithinwork_' thtext '_.eps'],'-depsc','-r300')

clear permdif perm1 perm2


clear permdif pv1 pv2 
f1id=0;
for f1=1:3
    f1id=f1id+1
    f2id=0;
for f2=1:3
    f2id=f2id+1; 
    if f1>=f2
    else
    for p=1:1000
        clear perm1 perm2
        for n=1:10
            fid=find(PowerNetClass==PowerNetLabel{n,1});
            
            nfid = 1:264;
            nfid(fid) = []; 
            v1=mean(dat(fid,nfid,bi{f1}),3); 
            v2=mean(dat(fid,nfid,bi{f2}),3);
            
            [a r] = sort(rand(length(fid),length(nfid),2),3);
            clear pv1 pv2
            for a=1:length(fid) 
            for b=1:length(nfid) 
                if r(a,b,1)==1
                    pv1(a,b) = v1(a,b);
                    pv2(a,b) = v2(a,b);
                else
                    pv1(a,b) = v2(a,b);
                    pv2(a,b) = v1(a,b);
                end
            end
            end  
            perm1(n) = sum(sum(pv1)); 
            perm2(n) = sum(sum(pv2));
        end
        permdif(f1id,f2id,:,p) = perm1-perm2;  
    end
    end
end
end

clear emp emps emp1 emp2
f1id=0; 
for f1=1:3
    f1id = f1id+1; 
    f2id=0;     
for f2=1:3
    f2id = f2id+1;     
    for n=1:10
        fid=find(PowerNetClass==PowerNetLabel{n,1});
        nfid = 1:264;
        nfid(fid) = []; 
        emp1(n) = sum(sum(mean(dat(fid,nfid,bi{f1}),3))); 
        emp2(n) = sum(sum(mean(dat(fid,nfid,bi{f2}),3)));         
    end
    emp(f1id,f2id,:) = emp1-emp2; 
    emps(:,f1id) = emp1;    
end
end

sig=zeros(3,3,10);
for f1id=1:3
for f2id=1:3
    if f1id>=f2id
    else
        for n=1:10   
            s=sort(squeeze(permdif(f1id,f2id,n,:)));
            if s(992) < emp(f1id,f2id,n) || s(8) > emp(f1id,f2id,n)
               sig(f1id,f2id,n) = 1;  
            end  
        end
    end
end
end

clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a o]=sort(fid,'descend') 



h=figure('Position',[100,700,800,400])
emps(:,4) = NaN; 
e=reshape(emps(o,:)', size(emps,2)*size(emps,1),1);
bar(e)
x=1; 

for n=1:10
    fid=find(sig(:,:,o(n)));
    stars = {}; 
    if sum(ismember(fid,4))==1;
        sigstar([x x+1]);
    end
    if sum(ismember(fid,8))==1;
        sigstar([x+1, x+2]);
    end
    if sum(ismember(fid,7))==1;
        sigstar([x, x+2]);
    end    
    x=x+4
end
       
        
sigstar([1.1,1.2])
set(gca,'XTick',2:4:40)
set(gca,'XTickLabel',{PowerNetLabel{o,2}}') 
legend('0.016', '0.035', '0.075')

print(h,[mdir '/figures/NEGbetweenwork_' thtext '_.eps'],'-depsc','-r300')







clear permdif perm1 perm2 emps emp1 emp2 e


clear emp emps emp1 emp2
f1id=0; 
for f1=1:3
    f1id = f1id+1; 
    f2id=0;     
for f2=1:3
    f2id = f2id+1;     
    for n=1:10
    for nn=1:10
        fid=find(PowerNetClass==PowerNetLabel{n,1});
        fid2=find(PowerNetClass==PowerNetLabel{nn,1});
        emp1(n,nn) = sum(sum(mean(dat(fid,fid2,bi{f1}),3))); 
        emp2(n,nn) = sum(sum(mean(dat(fid,fid2,bi{f2}),3)));
        num(n,nn,f1id,f2id)=length(find((mean(dat(fid,fid2,bi{f1}),3))));        
    end
    end
    emp(f1id,f2id,:,:) = emp1-emp2; 
    emps(:,:,f1id) = emp1;    
end
end


clear permdif pv1 pv2 
f1id=0;
for f1=1:3
    f1id=f1id+1
    f2id=0;
for f2=1:3
    f2id=f2id+1; 

    for p=1:1000
        clear perm1 perm2 
        for n=1:10
            for nn=1:10
                
            fid=find(PowerNetClass==PowerNetLabel{n,1});
            fid2=find(PowerNetClass==PowerNetLabel{nn,1});
            
            v1=mean(dat(fid,fid2,bi{f1}),3); 
            v2=mean(dat(fid,fid2,bi{f2}),3);
            
            [a r] = sort(rand(length(fid),length(fid2),2),3);
            clear pv1 pv2
            for a=1:length(fid) 
            for b=1:length(fid2) 
                if r(a,b,1)==1
                    pv1(a,b) = v1(a,b);
                    pv2(a,b) = v2(a,b);
                else
                    pv1(a,b) = v2(a,b);
                    pv2(a,b) = v1(a,b);
                end
            end
            end  
            perm1(n,nn) = sum(sum(pv1)); 
            perm2(n,nn) = sum(sum(pv2));
            end
        end
        permdif(f1id,f2id,:,:,p) = perm1-perm2;  

    end
end
end



sig=zeros(3,3,10,10);
for f1id=1:3
for f2id=1:3
        for n=1:10  
        for nn=1:10
            s=sort(squeeze(permdif(f1id,f2id,n,nn,:)));
            if s(992) < emp(f1id,f2id,n,nn) 
               sig(f1id,f2id,n,nn) = 1; 
            end

        end
        end
end
end



clear fid 
for c=1:10
     fid(c)=length(find(PowerNetClass==PowerNetLabel{c,1}));
end
[a o]=sort(fid,'descend') 

conmat=zeros(10,10,3)
for n=1:10 
    for nn=1:10
        if sig(1,2,n,nn)==1 && sig(1,3,n,nn)==1 && max(max(num(n,nn,:,:)))>=1
            conmat(n,nn,1) = 4;

        end
        if sig(2,1,n,nn)==1 && sig(2,3,n,nn)==1 && max(max(num(n,nn,:,:)))>=1
            conmat(n,nn,2) = 4; 

        end
        if sig(3,1,n,nn)==1 && sig(3,2,n,nn)==1 && max(max(num(n,nn,:,:)))>=1
            conmat(n,nn,3) = 4; 
     
         end
    end
end
      
schemaball({PowerNetLabel{o,2}},conmat(o,o,1)) 
print(gcf,[mdir '/figures/NEGnetwork_to_network_' thtext '_1n.eps'],'-depsc','-r300')

schemaball({PowerNetLabel{o,2}},conmat(o,o,2)) 
print(gcf,[mdir '/figures/NEGnetwork_to_network_' thtext '_2n.eps'],'-depsc','-r300')

schemaball({PowerNetLabel{o,2}},conmat(o,o,3)) 
print(gcf,[mdir '/figures/NEGnetwork_to_network_' thtext '_3n.eps'],'-depsc','-r300')



%% BETWEENESS RESIDUALS AND HUBS

bet_foi_mult = bet_foi*100
[a o] = sort(bet_foi_mult,1,'descend')
topbet=o; 

pa=40
betkeep = zeros(264,3) 
betkeep(o(1:pa,1),1) = bet_foi_mult(o(1:pa,1),1);
betkeep(o(1:pa,2),2) = bet_foi_mult(o(1:pa,2),2);
betkeep(o(1:pa,3),3) = bet_foi_mult(o(1:pa,3),3); 



for f=1:3   
    node=betkeep(:,f)
    export_power264_graph(node,1,PowerNetClass,['resub_hubs' num2str(f)])
end
node(1:264)=1
export_power264_graph(node,1,PowerNetClass,['resub_allnondes'])



close all
BrainNet_MapCfg(['resub_hubs1.node'],'settings.mat',['hubs_f1.jpg'])
close all
BrainNet_MapCfg(['resub_hubs2.node'],'settings.mat',['hubs_f2.jpg'])
close all
BrainNet_MapCfg(['resub_hubs3.node'],'settings.mat',['hubs_f3.jpg'])


c=0;
clear rho p
h=figure('Position',[200 200 500 300])
for f=1:3
    for ff=1:3
        if f<ff
c=c+1
U=union(find(nbet_foi(:,f)>0),find(nbet_foi(:,ff)>0))
[rho(f,ff) p(f,ff)] = corr(nbet_foi(:,f),nbet_foi(:,ff),'type','Spearman')
r=regstats(nbet_foi(:,f),nbet_foi(:,ff), 'linear')
hettest(f,ff) = TestHet(r.r,r.yhat,'-BPK')

subplot(2,3,c)
scatter(nbet_foi(U,f),nbet_foi(U,ff),15,'k','filled');
lsline
xlabel(num2str(f))
ylabel(num2str(ff))
box on
axis equal
axis([0 .2 0 .2])
set(gcf, 'PaperUnits', 'inches','PaperPosition',[0 0 10 10]);
subplot(2,3,c+3)
 scatter(r.yhat, r.r,15,'k','filled')
box on
 axis([0 .15 -0.1 .1])
axis square

        end
%topbet
    end
end

export_fig 'hetroskadcicity' '-eps' '-q100'


clear num
for foi=1:3
for n=1:11, 
    if n~= 11, 
        fid=find(PowerNetClass(topbet(1:40,foi))==PowerNetLabel{n,1}), 
        num(n,foi)=length(fid); 
    else, 
        fid=find(PowerNetClass(topbet(1:40,foi))==-1), 
        num(n,foi)=length(fid); 
    end; 
end
end

NO=NetOrder
NO(end)=11
figure
bar(num(NO,:))
print(gcf,'/data/william/Projects/dfc_timeseries/betweenness_subgraphpropertion','-r300','-depsc')


tb=topbet(1:40,:)
a=length(intersect(tb(:,1),tb(:,2)))./40
b=length(intersect(tb(:,1),tb(:,3)))./40
c=length(intersect(tb(:,2),tb(:,3)))./40
d=length(intersect(intersect(tb(:,1),tb(:,2)),tb(:,3)))./40
figure
venn([1 1 1],[a b c d],'EdgeColor','black')
axis off
set(gcf,'render','painters')
print(gcf,'/data/william/Projects/dfc_timeseries/betweenness_venn','-r300','-depsc')

figure
[H,S] = venn([1 1 1],[a b c d],'FaceAlpha', 0.6);

  %Now label each zone 
  for i = 1:3
      text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), ['Zone ' num2str(i)])
  end

