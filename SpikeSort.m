%% Spike Sorting code written by JSP and PCN

%% Parameters.


global onewaveformlength;
global numberchanns;
% Length of each spike clip, in samples.
onewaveformlength = 32;                                                    
numberchanns = 30;

% Upsampling factor 
ndecimate = 5;                                                             

% Electrode group to analyze
group = 2;                                                                 

% Data directory.
datadir = '/Users/jason/MEA/MEAData/02_01_12/trial2/';                     
outputdir = '/Users/jason/MEA/MEAData/02_01_12/SpikeTimes_trial2/';
% Directory in which to store clustered data
clusterdir = ['/Users/jason/MEA/MEAData/02_01_12/group',num2str(group),'clusters'];
% Prefix of files to process
infileprefix = 'trial2_';
% Prefix of files to write spike times to
outfileprefix = 'ST_'; 
% List of files to cluster
clustfiles = [30];                                                   

% File from which to extract noise
noisefile = [datadir,infileprefix,num2str(clustfiles(1))]; 
% List of file numbers to fit
fitfiles = 1:100;                                                           

% Length of time per data file (samples)
filelen = 30e4;                                                            

% Data loading function: its first argument should be a file name and it should return
% [y,time], where each row of y is the time series from one channel and time contains the
% time stamps for each sample.
loadfun = @(name) new_LoadAndFilterPenn(name,filelen);                     
                                                                            

% Specify array type
array = MEA(@geom_2F30EL);
%array = MEA(@geom_tetrode);                                               

% Threshold for spike extraction (uV)
threshold = -40;
% Threshold for noise extraction 
%(uV, will be compared to absolute value of voltage so must be positive)
noisethresh = 30;                                                          

whiten = true;



%% Clustering and template building.
% Pre-Clustering
[S,T,chmask,peakch]=extractSpikes([datadir,infileprefix],...
                    clustfiles,'clustering',group,threshold,loadfun,array);
                
display([num2str(size(S,1)),' spikes extracted. Getting noise...']);

[K,noise] = getNoiseAndWhiteningFilter([datadir,infileprefix],...
                clustfiles(1),group,noisethresh,loadfun,array);

if whiten, S = S*K; end   % Whiten

allspikesN = size(S,1);

% THIS IS THE MAIN CALL FOR OPTICS CLUSTERING
ind = order_points(S,chmask,peakch,datadir);


if ~exist(clusterdir,'dir')
    mkdir(clusterdir);
end
save([clusterdir, '/DataToCluster_',infileprefix,'.mat'],'S','ind');
display('Data sorted, ready for clustering.');

% Manual clustering
clu = plots;
nclusters = length(clu);

for i=1:nclusters
    IndexOrigin = ind(clu{i});
    C = S(IndexOrigin,:);
    save([clusterdir,'/Cluster',num2str(i),'.mat'],'IndexOrigin','C');
end
cd(clusterdir)


%% Manual decisions on cluster merging & splitting.

[nclusters,newname] = show_overlapping_templates([clusterdir,'/'],'Cluster');

nclusters = clustersplit;
%%
[nclusters,newname] = show_overlapping_templates([clusterdir,'/'],newname);


%%
% Build templates.
[TemplateMatrix,clusterind] = makeTemplates(newname,nclusters,array,group,ndecimate);
nclusters = size(TemplateMatrix,2);
%%
% Set templates to ignore.
%ignorelist = [4,5,11];     % G1
%ignorelist = [1,8,9];       % G2
ignorelist = [];
keeplist = setdiff(1:nclusters,ignorelist);
save('TemplateMatrix.mat','TemplateMatrix','keeplist');

% Collect statistics.

TemplateMatrix_reshaped = reshape(TemplateMatrix,...
                size(TemplateMatrix,1)/numberchanns,numberchanns,size(TemplateMatrix,2));           
[tempamp,minch] = min(min(TemplateMatrix_reshaped));
amplist = cell(1,nclusters);
avgrate = 0;
allspikesN = 0;

for i=1:length(keeplist)
    load([newname,num2str(keeplist(i)),'.mat']);

    amplist{keeplist(i)} = zeros(1,size(C,1));
    C_reshaped = reshape(C',size(C,2)/numberchanns,numberchanns,size(C,1));
    for j=1:size(C,1)
        amplist{keeplist(i)}(j) = min(C_reshaped(:,minch(keeplist(i)),j));
    end
    amplist{keeplist(i)} = amplist{keeplist(i)} / tempamp(keeplist(i)); 
    avgrate = avgrate+size(C,1);
    allspikesN = allspikesN + size(C,1);
end

totalElapsedTime = length(clustfiles)*filelen;
avgrate = avgrate/totalElapsedTime / ndecimate;

% Doesn't do anything but remember what the different columns of stats mean...
statnames = {'sumamp','sumampsquare','amplmean','amplvar','nspikes','rate'};
stats = zeros(nclusters,6);
stats = updateClusterStats(stats,amplist,totalElapsedTime*ndecimate);

% Parameters of exponential fit to noise autocorrelation function.
[fitCorMag,fitCorLen] = getCorrFit(noise);  

xout = -50:0.1:0;
h = histc(min(noise,[],2),xout) / size(noise,1);
Cnoise = [xout',cumsum(h)];
save('Stats.mat','stats','fitCorMag','fitCorLen','Cnoise','keeplist','allspikesN','totalElapsedTime','K');

%%
% Restrict to kept templates.
TemplateMatrix = TemplateMatrix(:,keeplist);
stats = stats(keeplist,:);
nclusters = length(keeplist);

cd(datadir);
%% Fitting.
% Big loop:
%  a. Extract spikes from one 30 s. file.
%  b. Fit using cmultifit.
%  c. (optionally) update priors.
nclusters = size(TemplateMatrix,2);
allamps = cell(1,nclusters);
[shiftadj,shifttemp,candidatelist,shifts] = preFit(fitCorLen,fitCorMag,TemplateMatrix,ndecimate);
for i=1:length(fitfiles)
    % Get spikes and times
    [S,T] = extractSpikes([datadir,infileprefix],fitfiles(i),'fitting',group,threshold,loadfun,array);
    if whiten, S = S*K; end    % Whiten
    
    % Fit to templates
    [SpikeTimes,amplist,incompFits] = cmultifit(S,T,TemplateMatrix,shiftadj,shifttemp,candidatelist,shifts,stats,Cnoise,array,group);
    
    % Update cluster statistics
    for j=1:nclusters
         allspikesN = allspikesN + length(SpikeTimes{j});
         allamps{j} = [allamps{j},amplist{j}];
    end
    totalElapsedTime = totalElapsedTime + filelen;
    
%   %Update priors       
%    stats = updateClusterStats(stats,amplist,totalElapsedTime*ndecimate);

    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
    save([outputdir,outfileprefix,'_group',num2str(group),'_',num2str(fitfiles(i)),'.mat'],'SpikeTimes','stats','amplist');
    save([outputdir,'incomplete_group',num2str(group),'_',num2str(fitfiles(i)),'.mat'],'incompFits');
end


