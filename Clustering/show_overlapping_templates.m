function [nclusters,newname] = show_overlapping_templates(clusterdir,fileprefix,forcemerge)

%allparameters; 
global onewaveformlength; global numberchanns;

if nargin < 3
    forcemerge = [];
end

d = dir(clusterdir);
filenames = {d.name};
inds = strmatch(fileprefix,filenames);
nclusters = length(inds);

T = zeros(nclusters,onewaveformlength*numberchanns);
for i=1:nclusters
    load([clusterdir,fileprefix,num2str(i),'.mat']);
    T(i,:) = median(C);
end

TT = reshape(T',[onewaveformlength,numberchanns,size(T,1)]);
TT = cat(1,zeros(5,numberchanns,size(T,1)),TT,zeros(5,numberchanns,size(T,1)));
TT = reshape(TT,[(onewaveformlength+10)*numberchanns,size(T,1)]);

TT_shift = zeros((onewaveformlength+10)*numberchanns,size(T,1)*10);
for i=1:size(T,1)
    for j=1:10
        TT_shift(:,(i-1)*10+j) = circshift(TT(:,i),[j-5,0]);
    end
end
CC = corrcoef(TT_shift);
C = zeros(size(T,1));
shifts = zeros(size(T,1));
for i=1:size(T,1)
    for j=1:size(T,1)
        crosscorr = CC((i-1)*10+1:i*10,(j-1)*10+1:j*10);
        [m,maxind] = max(crosscorr(:));
        [t0,t1] = ind2sub(size(crosscorr),maxind);
        C(i,j) = max(crosscorr(:));
        shifts(i,j) = t1-t0;
    end
end

flags = sparse(false(size(C)));
for i=1:size(forcemerge,1)
    flags(forcemerge(i,1),forcemerge(i,2)) = true;
    flags(forcemerge(i,2),forcemerge(i,1)) = true;
end

for i=1:size(T,1)
    for j=(i+1):size(T,1)
        if C(i,j) > 0.5
            t1 = reshape(TT(:,i),onewaveformlength+10,numberchanns);
            t2 = circshift(reshape(TT(:,j),onewaveformlength+10,numberchanns),[shifts(i,j),0]);
            for ch=1:numberchanns
                subplot(6,5,ch); plot([t1(:,ch),t2(:,ch)]); ylim([min(min(t1(:)),min(t2(:))),max(max(t1(:)),max(t2(:)))]);
            end
            query = input('Merge (y/n)? ','s');
            if ~flags(i,j)
                flags(i,j) = strcmp(query,'y');
                flags(j,i) = flags(i,j);
            end
        end
    end
end

unprocessed = 1:size(T,1);
dups = {};
while ~isempty(unprocessed)
    queue = [unprocessed(1)];
    conncomp = [];
    while ~isempty(queue)
        conncomp = [conncomp, queue(1)];
        unprocessed = unprocessed(unprocessed ~= queue(1));
        neighbors = find(flags(queue(1),:));
        neighbors = intersect(neighbors,unprocessed);
        queue = queue(2:end);
        queue = [queue,setdiff(neighbors,queue)];
    end
    dups = [dups {conncomp}];
end

newname = ['m',fileprefix];

for i=1:length(dups)
    mC = [];
    mIndexOrigin = [];    
    for j=dups{i}
        load([clusterdir,fileprefix,num2str(j),'.mat']);
        Cpad = reshape(C',[onewaveformlength,numberchanns,size(C,1)]);
        Cpad = cat(1,zeros(5,numberchanns,size(C,1)),Cpad,zeros(5,numberchanns,size(C,1)));        
        Cpad = circshift(Cpad,[shifts(dups{i}(1),j),0,0]);
        C = reshape(Cpad(6:(onewaveformlength+5),:,:),[onewaveformlength*numberchanns,size(C,1)])';
        mC = [mC; C];
        mIndexOrigin = [mIndexOrigin; IndexOrigin];
    end
    C = mC; IndexOrigin = mIndexOrigin;
    
    save([clusterdir,newname,num2str(i),'.mat'],'C','IndexOrigin');
end

nclusters = length(dups);
