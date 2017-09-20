%% T - graphlet analysis with bursts. 


%% Load data
%Requires voxall to be loaded (see clusters_withpca.m)


load('/data/william/Projects/dfc_time/cost_of_kmean_run3')
cd(mdir)

parpool(20)
for k=[5 10 12 8]

kchoice = k


cd /data/william/Projects/dfc_timeseries/
roi_264_powerneuron
cd /data/william/Projects/dfc_time

clear g graph

addpath(genpath('/data/william/toolbox/bc/'))

% Make s graphlets

for n=1:kchoice
    if sum(kmeangroups(:,kchoice)==n)>1

        g(:,:,n)=corr(voxall(:,kmeangroups(:,kchoice)==n)');

        for r=1:length(dellocal)
            g(r,dellocal{r},n)=0;
            g(r,r,n)=0;
        end
         graph(:,:,n) = g(:,:,n); 
    end
end

clear   bingraph thgraph 
for th=[5 10]
    ii=0;
    for n=1:kchoice
            ii=ii+1
            thgraph{th}(:,:,ii) = threshold_proportional(graph(:,:,ii),th/100);
            bingraph{th}(:,:,ii) = weight_conversion(thgraph{th}(:,:,ii),'binarize');
    end
end

%% Calculate burstiness

clear B ts
for th=[5 10]
B=zeros(264,264)*NaN; 
dcol=[]
gr=bingraph{th};
parfor n=1:264
    for m=1:264
        dcol=[];
        for s=1:100
            ts=squeeze(gr(n,m,kmeangroups((s-1)*1200+1:s*1200,kchoice)));
            fid=find(ts==1);
            if ~isempty(fid)
                d=diff(fid);
                dcol(end+1:end+length(d))=d; 
            end
        end
        B(n,m)=(std(dcol)-mean(dcol))./(std(dcol)+mean(dcol)); 
    end
end
save(['/data/william/Projects/dfc_time/B' num2str(k) '_' num2str(th)],'B')



%% Permution test on burstiness 

clear PB ts
PB=zeros(264,264)*NaN; 
dcol=[]
parfor p=1:200    
[~, r]=sort(rand(1200,1));
for n=1:264
    for m=1:264
        dcol=[];
        for s=1:100
            ts=squeeze(gr(n,m,kmeangroups((s-1)*1200+1:s*1200,kchoice)));
            fid=find(ts(r)==1);
            if ~isempty(fid)
                d=diff(fid);
                dcol(end+1:end+length(d))=d; 
            end
        end
        PB(n,m,p)=(std(dcol)-mean(dcol))./(std(dcol)+mean(dcol)); 
    end
end
end
save(['/data/william/Projects/dfc_time/PB_k' num2str(k) '_' num2str(th)],'PB')
    


end
end


%% For loop over clusters and thresholds
clear grcount; 
for k=[5 8 12]
for th=[5 10]



	load(['/data/william/Projects/dfc_time/PB_k' num2str(k) '_' num2str(th)])
	load(['/data/william/Projects/dfc_time/B' num2str(k) '_' num2str(th)])

	% Caclulate average burstiness per network-network interaction 

	for n=1:10
	    for m=1:10
		fid1=find(PowerNetClass==PowerNetLabel{n,1});
		fid2=find(PowerNetClass==PowerNetLabel{m,1});
		Bnet(n,m)=nanmean(nanmean(B(fid1,fid2)));
	    end 
	end


	Bbet=[]; 
	Bbet_e{k}{th} = [];  
	Bdiag = []; 
	Bdiag_e{k}{th} = [];  

	% Within or between burstiness 

	for n=1:10
	    for m=1:10
		fid1=find(PowerNetClass==PowerNetLabel{n,1});
		fid2=find(PowerNetClass==PowerNetLabel{m,1});
		if n<m
		    Bbet(end+1)=Bnet(n,m); 
		    Btmp=reshape(B(fid1,fid2),length(fid1)*length(fid2),1);
		    Btmp(isnan(Btmp)==1)=[];
		    Bbet_e{k}{th}(end+1:end+length(Btmp))=Btmp;
		    
		elseif n==m
		    Bdiag(end+1)=Bnet(n,n);
		    q = find(triu(true(size(B(fid1,fid2))))==1);
		    Btmp=B(q);
		    Btmp(isnan(Btmp)==1)=[];
		    Bdiag_e{k}{th}(end+1:end+length(Btmp))=Btmp;
		end
	    end
	end


	% Plot differences in between and within network burstiness 

	figure; bar([nanmean(Bbet) nanmean(Bdiag)])

	v1=Bbet
	v2=Bdiag
	for n=1:1000
	    [~, ro] = sort(rand(length(v1)+length(v2),1));
	    v=[v1 v2];
	    pv1=v(ro(1:length(v1)));
	    pv2=v(ro(length(v1)+1:end)); 
	    p(n)=nanmean(pv1)-nanmean(pv2); 
	end
	p=sort(p); 
	    
	hold on
	errorbar([nanmean(Bbet),nanmean(Bdiag)],[nanstd(Bbet),nanstd(Bdiag)],'kx')
	if p(25) > nanmean(Bbet)-nanmean(Bdiag) || p(975) < nanmean(Bbet)-nanmean(Bdiag)
	    sigstar([1,2]);
	end 
	axis([0.5 2.5 -.8 .6 ])
	print(gcf,[mdir '/figures/burstyresults_diag_vs_upper_' num2str(th) '_k' num2str(k)],'-r300','-depsc')


	% Find if single edge is significant

	spermBmax = sort(squeeze(max(max(PB)))); 
	spermBmin = sort(squeeze(min(min(PB)))); 
	spermB=sort(PB,3);
	SB=zeros(264,264) 
	for n=1:264
	    for nn=1:264
		if B(n,nn)<0 
		    if B(n,nn) <= spermBmin(5) %spermB(n,nn,5)
		        SB(n,nn) = B(n,nn); 
		    end
		elseif B(n,nn)>0 
		    if B(n,nn) >= spermBmax(195) %spermB(n,nn,195)
		        SB(n,nn) = B(n,nn); 
		    end
		end
	    end
	end

	
	% Additional plot when k=8. The idea of this plot is to show the distribution of thresholded connections are not biased over states 
	if k==8 
	    grcount{th} = zeros(kchoice,1)
	    for s=1:8
	        Bperstate_within{s}=[]; 
	        Bperstate_bet{s} = []; 
	    end
	    for i=1:264
		for j=1:264
		    if i<j
		        if SB(i,j)>0
		            grcount{th}=grcount{th}+squeeze(bingraph{th}(i,j,:));
		        end
		        for s=1:8 
		            if bingraph{th}(i,j,s)==1
		                if PowerNetClass(i)==PowerNetClass(j) && PowerNetClass(i)>0
		                   Bperstate_within{s}(end+1) = B(i,j);
		                else
		                   Bperstate_bet{s}(end+1) = B(i,j);
		                end
		            end
		        end
		    end
		end
            end
		
	    if th==10
		figure
		hold on
		plot(grcount{5},'x-')
		plot(grcount{10},'x-')
		axis([0 9 0 1000])
		box on
		print(gcf,'/data/william/Projects/dfc_time/figures/bursty_states','-r300','-depsc')
	    end
	end


	% Plots significant periodic and bursty edges

	for n=1:10
	    for nn=1:10 
		f1 = find(PowerNetClass==PowerNetLabel{n,1});
		f2 = find(PowerNetClass==PowerNetLabel{nn,1});
		Perodic(n,nn)=length(find(SB(f1,f2)<0))./(length(f1)*length(f2));
		Bursty(n,nn)=length(find(SB(f1,f2)>0))./(length(f1)*length(f2));
	    end
	end
	figure
	subplot(2,1,1)
	imagesc(Perodic(NetOrder(1:10),NetOrder(1:10)))
	cm = hot(2000)
	cm(1,:) = [0 0 0]
	colormap(cm)
	colorbar
	if th==5
	caxis([0 .4])
	else
	    caxis([0 .6])
	end
	subplot(2,1,2)
	imagesc(Bursty(NetOrder(1:10),NetOrder(1:10)))
	if th==5
	caxis([0 .4])
	else
	    caxis([0 .6])
	end
	colormap(cm)
	colorbar
	title([num2str(k) ' - ' num2str(th)])
	print(gcf,[mdir '/figures/g_burstyresults_mccorr_' num2str(th) '_k' num2str(k)],'-r300','-depsc')

end
end


%% plot distributions of B valutes

figure 
subplot(2,2,1); 
hist(Bbet_e{8}{5},-1:0.1:1); 
axis([-1.2 1.2 0 250])
axis square
subplot(2,2,2); 
hist(Bdiag_e{8}{5},-1:0.1:1)
axis([-1.2 1.2 0 50])
axis square
subplot(2,2,3); 
hist(Bbet_e{8}{10},-1:0.1:1); 
axis([-1.2 1.2 0 650])
axis square
subplot(2,2,4); 
hist(Bdiag_e{8}{10},-1:0.1:1)
axis([-1.2 1.2 0 70])
axis square
print(gcf,[mdir '/figures/distribution_of_burstycoeefs_k8_th55_1010'],'-depsc','-r300')



%% Intercontact times distribution example
d=[]
for s=1:100
n=192, m=134,
bingraphstime(:,:,1:1200)=bingraph{th}(:,:,kmeangroups((s-1)*1200+1:s*1200,kchoice));
fid=find(squeeze(bingraphstime(n,m,:))==1);
d(end+1:end+length(fid)-1)=diff(fid); 
end

f=poissfit(d)
g=poisspdf(d,f)
f=gpfit(d)
h1=gppdf(1:max(d),f(1),f(2),0)
h2=gppdf(1:max(d),f(1),f(2),-1)
[i, o] = sort(d);

t=i.^-1.75
t2=i.^-exp(1)
figure
semilogy(1:max(d),hist(d,1:max(d))./sum(hist(d,1:max(d))),'ko')
hold on
semilogy(d(o),g(o),'k--'); 
hold on
%semilogy(d(o),t2,'r--'); 
semilogy(d(o),t,'k'); 
% semilogy(1:max(d),h1,'c--'); 
% semilogy(1:max(d),h2,'r--'); 
hold on; 
set(gca,'XTick',[1 10:10:max(d)+10],'XTickLabel',[1 10:10:max(d)+10]) 
axis([0 300 0.000001 1])
xlabel('Inter Edge Interval (scans)')
ylabel('Probability')
legend('Emperical Data','Poissen distribution','t^-a')
print(gcf,'/data/william/Projects/dfc_time/figbursty_example_distribution','-r300','-depsc')
% end



