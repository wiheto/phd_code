
%% Parameters

S = 46; %NumberOfSubjects
L = 240 %DataLength
C = 2   %Condition
N = 264 % Nuber of ROIs 
sdir = '/data/directory/to/save/'

%% Load data
data=get_eoecNdata(); 
%data should be roi x time x conditions x subjects
data=reshape(data,N,L,C,S);

%% Standerdize data
%not really necessary unless PCA is going to be performed
for s=1:S
    for c=1:C
        data(:,:,c,s)=bsxfun(@rdivide,bsxfun(@minus,data(:,:,c,s),mean(data(:,:,c,s),2)),std(data(:,:,c,s),[],2));
    end
end


%% Weighted Pearson correlation

rg=zeros(N,N,L,C,S); 
for c=1:C
for s=1:S
    s
    % Global weights, (i.e. scaled inverse euclidean distance)
    Ggt=squareform(pdist(data(:,:,c,s)'));
    timew=1./Ggt;
    timew(timew==inf)=NaN; 
    timew=bsxfun(@rdivide,bsxfun(@minus,timew,min(min(timew))),max(max(timew))-min(min(timew))); 
    timew(isnan(timew)==1)=1;
        
    
    parfor r1=1:N
        rvec{r1}=zeros(N,L); 
        for r2=1:N 
            if r1<r2
                data2corr=data([r1 r2],:,c,s)';
                for t=1:L
                    tmp=weightedcorrs(data2corr,timew(t,:)+0.000001); 
                    rvec{r1}(r2,t)=tmp(1,2);
                end
            end
        end
    end
    rg(:,:,:,c,s)=permute(reshape(cell2mat(rvec),N,L,N),[3 1 2]); 
    rg(:,:,:,c,s)=rg(:,:,:,c,s)+permute(rg(:,:,:,c,s),[2 1 3]);
end
end

%% Fisher and boxcox

clear tsBCtmp
for c=1:C
    for s=1:S
        l={};
        s
        tsBCtmp{s}{c}=zeros(N*N,L)*NaN; 
        dattmp=fisherz(rg(:,:,:,c,s)); 
        dattmp=reshape(dattmp,size(tsBCtmp{s}{c},1),L);
        % Scale so min(dattmp)=1
        dattmp=bsxfun(@plus,dattmp,1-min(dattmp,[],2));
        % Find lambda
        parfor n=1:size(dattmp,1)
            l{n}=boxcoxlm(dattmp(n,:)',[ones(L,1)],0,-5:0.1:5);
        end
        l=cell2mat(l); 
        lamba{s}=l; 
        % perform boxcox tansform
        tsBCtmp{s}{c}(l~=0,:)=bsxfun(@rdivide,(bsxfun(@power,dattmp(l~=0,:),l(l~=0)')-1),l(l~=0)');
        tsBCtmp{s}{c}(l==0,:)=log(dattmp(l==0,:));
    end
end
clear gtsBC
%Rescale values back to original mean
for c=1:C
for s=1:S
    s
    gtsBC(:,:,c,s)=bsxfun(@minus,tsBCtmp{s}{c},mean(tsBCtmp{s}{c},2)-mean(reshape(rg(:,:,:,c,s),N*N,L),2));
end
end
% Convert to z values (not necessary either, can be useful depending on thresholding techniques)
for s=1:S
    for c=1:C
        gtsBC(:,:,c,s)=bsxfun(@rdivide,bsxfun(@minus,gtsBC(:,:,c,s),mean(gtsBC(:,:,c,s),2)),std(gtsBC(:,:,c,s),[],2));
    end
end

%% rescale to roi x roi x L x C x S

gtsBC=reshape(gtsBC,N,N,L,C,S); 

%% Save 

save([sdir 'wcorrTSBCg'],'gtsBC','-v7.3')





























