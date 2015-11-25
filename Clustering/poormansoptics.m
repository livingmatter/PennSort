function [ind,reach] = poormansoptics(r,type)

% type = 'dis' if r is a distance, 'sim' if r is a similarity metric
if nargin < 2
    type = 'dis';
end
reach = Inf(1,size(r,1));
unprocessed = true(1,size(r,1));
ind = zeros(1,size(r,1));

%i = ceil(rand*size(r,1));
i=1;
ind(1) = i;
%reach(1) = Inf;
unprocessed(i) = false;
reach(unprocessed) = min(reach(unprocessed),r(ind(i),unprocessed));
for i=2:size(r,1)
    %[m,ix] = min(r(ind((i-min(tail,i-1)):(i-1)),unprocessed),[],2);    
    %[reach(i),ix2] = min(m);
    if strcmp(type,'dis')
        [m,ix] = min(reach(unprocessed));
    elseif strcmp(type,'sim')
        [m,ix] = max(reach(unprocessed));
    end
    %[m,ix] = min(r(ind(i-1),unprocessed));
    tmp = find(unprocessed);
    ind(i) = tmp(ix);
    unprocessed(ind(i))=false;
    if strcmp(type,'dis')
        reach(unprocessed) = min(reach(unprocessed),r(ind(i),unprocessed));
    elseif strcmp(type,'sim')
        reach(unprocessed) = max(reach(unprocessed),r(ind(i),unprocessed));
    end
end

