infileprefix = 'trial7g2';
inputdir = '/Users/jason/MEA/MEAData/09_25_09/trial7/SortCompare/SpikeTimes_trial7_s11/';
clusterdir = '/Users/jason/MEA/MEAData/09_25_09/trial7/SortCompare/trial7_group2clusters_s7/';
clusterhead = 'mmCluster';

global onewaveformlength;
global numberchanns;

onewaveformlength = 32;
numberchanns = 30;

fitfiles = 1:241;

maxnum = 100000;

nincomp = 0;
for i=fitfiles
    load([inputdir,'incomplete_group2_',num2str(i),'.mat']);
    nincomp = nincomp + length(incompFits);
end

S = zeros(min(maxnum,nincomp),numberchanns*onewaveformlength);
peakch = zeros(min(maxnum,nincomp),1);
chmask = zeros(min(maxnum,nincomp),numberchanns);

ix = 0;
for i=fitfiles(randperm(length(fitfiles)))
    display(['Processing ', num2str(i), '...']);
    load([inputdir,'incomplete_group2_',num2str(i),'.mat']);
    
    if ix+length(incompFits) <= maxnum
        S(ix+1:ix+length(incompFits),:) = reshape([incompFits(:).S],onewaveformlength*numberchanns,length(incompFits))';       
        peakch(ix+1:ix+length(incompFits)) = [incompFits(:).peakchan]';        
        chmask(ix+1:ix+length(incompFits),:) = reshape([incompFits(:).chmask],numberchanns,length(incompFits))';        

        ix = ix + length(incompFits);
    else
        break;
    end

end

S = S(1:ix,:);
paekch = peakch(1:ix);
chmask = chmask(1:ix,:);

save([inputdir,'incompleteall.mat'],'S','peakch','chmask');
%%
ix = find(min(S,[],2) < -45);

length(ix)
chmask = logical(chmask);
ind=order_points(S(ix,:),chmask(ix,:),peakch(ix),clusterdir);
save([clusterdir, 'DataToCluster_',infileprefix,'_2ndpass.mat'],'S','ind');
display('Data sorted, ready for clustering.');

clu = plots;

nclusters = length(clu);

d = dir(clusterdir);
[startind,endind] = regexp({d.name},[clusterhead,'[0-9]+']);
clustind = 0;
for i=1:length(d)
    if ~isempty(startind{i})
        newclustind = str2double(d(i).name(length(clusterhead)+startind{i}:endind{i}));
        if newclustind > clustind
            clustind = newclustind;
        end
    end
end

for i=1:nclusters
    IndexOrigin = ind(clu{i});
    C = S(IndexOrigin,:);
    save([clusterdir,clusterhead,num2str(clustind + i),'.mat'],'IndexOrigin','C');
end
cd(clusterdir)
