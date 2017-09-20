clear all

ldir = '/data/william/Projects/eoec_dfctime/'
mdir = '/data/william/Projects/dfc_jk/'
cd /data/william/Projects/dfc_time/
powernetclassclusters
dellocal=[];
for r=1:size(coord,1)
    fid=find(coord(:,1)<coord(r,1)+18 & coord(:,1)>coord(r,1)-18);
    if isempty(fid)==0
        fid2=find(coord(fid,2)<coord(r,2)+18 & coord(fid,2)>coord(r,2)-18);
        if isempty(fid2)==0
            fid3=find(coord(fid(fid2),3)<coord(r,3)+18 & coord(fid(fid2),3)>coord(r,3)-18);
            if isempty(fid3)==0
                dellocal{r}=fid(fid2(fid3));
            end
        end
    end
end
cd(mdir)
% Allsubject approach
S = 48; %NumberOfSubjects
L = 240 %DataLength
C = 1   %Condition
sid=0
c=1
clear vox voxtmp
for s=1:S
    if s==8 || s== 22
    else
    s
    sid=sid+1;
    load([ldir 'timeseries_2016/eo_' num2str(s) '.mat'])
    dat = eo'; 
    for roi=1:264
         voxtmp(roi,:)=dat(roi,:)-mean(dat(roi,:));
    end
    tmp = voxtmp(:,:);
    tmp1 = reshape(tmp,size(tmp,1)*size(tmp,2),1); 
    tmp=(tmp-mean(tmp1))./std(tmp1);
    vox(:,:,c,sid)=tmp; 
    end
end



clear RStatic
for s=1:S-2
    for c=1:C
    RStatic(:,:,c,s)=corr(vox(:,:,c,s)');
    end
end




% rperm=zeros(264,264,size(dat,1),32);
% Marked=zer


for sw=[31]

for slide = [1 round(sw/2)] 

s=1
rsw=[]; 
rsw=zeros(264,264,length(sw+1:slide:L-sw),1,S-2);
for s=1:S-2
    s
tid=0
for n=sw+1:slide:L-sw
    tid=tid+1; 
    for c=1:C
        rsw(1:264,1:264,tid,c,s)=corr(vox(:,n-sw:n+sw,c,s)');
    end
end
end

A=zeros(264,264);
ind=find(~tril(ones(size(A)))); 

rsw=reshape(rsw,264*264,size(rsw,3),S-2);
rsw=squeeze(rsw(ind,:,:)); 

zrsw=fisherz(rsw); 
zrsw=reshape(zrsw,size(rsw)); 

l=zeros(size(dattmp,1),S-2);
for s=1:S-2
    s
    dattmp=bsxfun(@plus,rsw(:,:,s),1-(min(rsw(:,:,s),[],2)));
    for e=1:size(dattmp,1)
        l(e,s)=boxcoxlm(dattmp(e,:)',[ones(size(rsw,2),1)],0,-5:0.01:5);
    end
end


lz=zeros(size(dattmp,1),S-2);
for s=1:S-2
    s
    dattmp=bsxfun(@plus,zrsw(:,:,s),1-(min(zrsw(:,:,s),[],2)));
    %dattmp=real(dattmp); 
    for e=1:size(dattmp,1)
        lz(e,s)=boxcoxlm(dattmp(e,:)',[ones(size(rsw,2),1)],0,-5:0.01:5);
    end
end

save(['/data/william/Projects/dfc_contrasts/boxcoxDMN_all_sw' num2str(sw) '_slide' num2str(slide)],'l','lz','-v7.3')

end
end




