%function [SpikeTimes, amplist] = dev_cmultifit(S,T,TemplateMatrix,stats,fitCorLen,fitCorMag,Cnoise,array,group)%,realID,realTau,realAmp)
function [SpikeTimes, amplist,incompFits] = cmultifit(S,T,TemplateMatrix,shiftadj,shifttemp,candidatelist,shifts,stats,Cnoise,array,group,realID,realTau)

% pcn 7/09 
% jsp 3/10

% Subfunctions:
% -----------------------------------------------
% function [spike1,t1,n] = concatenate(S,T,n)
%   Concatenate overlapping triggers
%

% exitcode= 'AI' "all templates ineligible"; 'LA' "low amplitude"; 'LR' "low probability ratio"; 'TM' "too many spikes found"; 'AM' "ambiguous"
                                
nclock=clock; disp([date 'CMULTIFIT.M begin ' num2str(nclock(4:5))])

% behavior
makeplots=false;
adjmethod = 'fast';  % 'fast' = fast but memory intensive, 'slow' = slow but memory efficient.
%watchlist=[];   % always display these events, identified by nsequence
%plotpause=true;
       
% params
global onewaveformlength; global numberchanns; 
midwaveform = floor(onewaveformlength/2);   % Location of peak in downsampled templates

timewindow = min(10,midwaveform-1);
interestingVThresh = -25;
interestingRThresh = 0; 
Nmaxfitspike=16;
spikenbd=21;

% vars
allrightscorelist =[];% list of likelihood ratios for best fits
madeAplot=false;
totalspikes=0;  % how many we fit
spikeTimesAllEvents=[]; % collect results
spikeIDsAllEvents=[];
spikeAmplAllEvents=[];
nclusters = size(TemplateMatrix,2);


% Find optimal thresholds to separate spikes from noise.
phi = zeros(size(Cnoise,1),nclusters);
for i=1:nclusters
    phi(:,i) = (1-Cnoise(:,2)).*normcdf(Cnoise(:,1),stats(i,3)*min(TemplateMatrix(:,i)),sqrt(stats(i,4))*abs(min(TemplateMatrix(:,i))));    
end
[separability,noisethresh] = max(phi);
noisethresh = Cnoise(noisethresh,1);


nshifts = size(shifttemp,4);
nwidths = size(shifttemp,5);

% Store neighborhoods of each channel; precompile template norms

nbds = zeros(spikenbd,numberchanns);
templnormbank = zeros(numberchanns,nclusters*nshifts*nwidths);
euclidtemplnormbank = zeros(numberchanns,nclusters*nshifts*nwidths);

if strcmp(adjmethod,'fast')
    shiftadjbank = zeros(numberchanns, (2*timewindow+1)*spikenbd,nclusters,nshifts*nwidths);
    shiftbank = zeros(numberchanns, (2*timewindow+1)*spikenbd,nclusters,nshifts*nwidths);
end
for nc=1:numberchanns,
    [nbrgrp,nbrch] = array.Neighborhood(group,nc,'channels',spikenbd);
    nbds(:,nc) = nbrch;
    resshift = reshape(shifttemp(midwaveform-timewindow:midwaveform+timewindow,nbrch,:,:,:), (2*timewindow+1)*spikenbd,nclusters*nshifts*nwidths);
    restemp = reshape(shiftadj(midwaveform-timewindow:midwaveform+timewindow,nbrch,:,:,:), (2*timewindow+1)*spikenbd,nclusters*nshifts*nwidths);
    templnormbank(nc,:) = sum(resshift .* restemp);
    euclidtemplnormbank(nc,:) = sum(resshift.^2);
    if strcmp(adjmethod,'fast')
        shiftadjbank(nc,:,:,:) = reshape(restemp, size(restemp,1), nclusters, nshifts*nwidths);
        shiftbank(nc,:,:,:) = reshape(resshift, size(resshift,1), nclusters, nshifts*nwidths);
    end
end

neighbors = cell(1,numberchanns);
for nc=1:numberchanns
     [nbrgrp,nbrch] = array.Neighborhood(1,nc,'distance',30*sqrt(2));
     neighbors{nc} = nbrch;
end

    
templnormbank = reshape(templnormbank,numberchanns,nclusters,nshifts*nwidths);
euclidtemplnormbank = reshape(euclidtemplnormbank,numberchanns,nclusters,nshifts*nwidths);


gamma = repmat(stats(:,3), [1,nshifts*nwidths]);
sigmasq = repmat(stats(:,4), [1,nshifts*nwidths]);
lnKmua = repmat(log(stats(:,6)), [1,nshifts*nwidths]);

nevents = size(S,1);
nsequence=1;


npresent = zeros(1,nclusters);
nfit = zeros(1,nclusters);
incompFits = [];

% MAIN LOOP:
tic;disp('go!');
while nsequence <= nevents   % at the end of this loop nsequence hops ahead over any concatenated triggers
    [spike1,t1,nsequence] = concatenate(S,T,nsequence);        % t1 is absolute time of first sample of spike1

    originalspike = spike1;                                    % spike1 will get templates subtracted; save original for eventual display

    fittedspikes = zeros(size(spike1,1),numberchanns);         % as we find spikes, place them here
    spikeIDsThisEvent=[];
    spikeTimesThisEvent=[];
    spikeAmplThisEvent=[];          % initialize lists of IDs and absolute times
    exitcode=''; 
    refractReason=''; % keep track of why the fit terminated
    
    questionablespike = zeros(size(spike1,1),numberchanns); % if we fit a spike but with poor R, we'll save it here.
    
    % SECONDARY LOOP        
    for secondloop=1:Nmaxfitspike,   % keep fitting spikes until done        
        [tmp,tmptimes] = min(spike1);                      % for each channel get the time when it hits peak
        [tmpampl,peakchan]=min(tmp);                       % get absolute peak
        
        if tmpampl>interestingVThresh, 
            exitcode='LA'; break; end                     % nothing left to fit: exit secondary loop        

        peaktime = tmptimes(peakchan);                     % time of abs peak
        peaknbd=nbds(:,peakchan);
                
        
        % Truncate spike to exactly 2*timewindow+1 samples near peak, and
        % to channels in neighborhood of peak:
        % The peak will be at position timewindow+1 in spiketrunc.
        
        spiketrunc = zeros(2*timewindow+1,length(peaknbd));
        timeR = min(peaktime+timewindow, size(spike1,1)) - peaktime;
        timeL = peaktime - max(peaktime-timewindow, 1);
        spiketrunc(timewindow+1-timeL:timewindow+1+timeR,:) = spike1(peaktime-timeL:peaktime+timeR,peaknbd);
        
                        
        % examine candidate templates:        
        candidates = candidatelist{peakchan};
        ncandidates = length(candidates);
                     
        if ~isempty(candidates)            
            % Compute log-posterior:
            
            if strcmp(adjmethod,'fast')
                adj = reshape(shiftadjbank(peakchan,:,candidates,:), size(shiftadjbank,2), ncandidates*nshifts*nwidths);            
                euclidadj = reshape(shiftbank(peakchan,:,candidates,:), size(shiftbank,2), ncandidates*nshifts*nwidths);            
            else
                [nbrgrp,nbrch] = array.Neighborhood(group,peakchan,'channels',spikenbd);
                adj = reshape(shiftadj(midwaveform-timewindow:midwaveform+timewindow,nbrch,candidates,:,:), (2*timewindow+1)*spikenbd,ncandidates*nshifts*nwidths);
                euclidadj = reshape(shifttemp(midwaveform-timewindow:midwaveform+timewindow,nbrch,candidates,:,:), (2*timewindow+1)*spikenbd,ncandidates*nshifts*nwidths);                
            end
            crossprod = spiketrunc(:)' * adj;
            euclidcrossprod = spiketrunc(:)' * euclidadj;
            crossprod = reshape(crossprod, ncandidates, nshifts*nwidths);
            euclidcrossprod = reshape(euclidcrossprod, ncandidates, nshifts*nwidths);
            
            templnorms = reshape(templnormbank(peakchan,candidates,:), length(candidates), nshifts*nwidths);
            euclidtemplnorms = reshape(euclidtemplnormbank(peakchan,candidates,:), length(candidates), nshifts*nwidths);           
            
            Astarnum = gamma(candidates,:) + sigmasq(candidates,:).*crossprod;          % Numerator
            Astarden = 1 + sigmasq(candidates,:).*templnorms;                           % Denominator                            

            lnPmut = lnKmua(candidates,:) + ...
                     (1./(2*sigmasq(candidates,:))) .* (-gamma(candidates,:).^2 + (Astarnum.^2 ./ Astarden)) - ...
                     0.5*log(Astarden);
                        
            % Marginalize over widths:
            ref = max(lnPmut(:));
            lnPmut = log(mean(exp(reshape(lnPmut,length(candidates),nshifts,nwidths)-ref),3))+ref;

            % Rank templates & find best shifts
            [Pbests,bestshifts]=max(lnPmut,[],2);  % best shift for each template, relative to peak
            [rankorder,ranklist]=sort(Pbests,'descend');

            Pbest = rankorder(1);               
            bestmu = candidates(ranklist(1));
            
            % "Second fit" to identify best amplitude, shift, and width which minimizes residual:            
             
             Astar = euclidcrossprod(ranklist(1), :) ./ euclidtemplnorms(ranklist(1), :);
             resid = -euclidcrossprod(ranklist(1), :).^2 ./ euclidtemplnorms(ranklist(1),:);
             [r,bestind] = min(resid);
             bestA = Astar(bestind);
             [bestshift,bestwidth] = ind2sub([nshifts,nwidths],bestind);
            
            
                
            bestTemp = shifttemp(:,:,bestmu,bestshift,bestwidth);
            Pmarg=log(sum(exp(lnPmut(ranklist(1),:)-Pbest)))+Pbest;   % marginalize over timeshifts
            
            if Pmarg < interestingRThresh
                exitcode='LR';         % this spike has "Low R" but save it anyway for scrutiny
                questionablespike = reconstructSpike(bestA,peaktime,bestTemp,1:size(questionablespike,1),onewaveformlength,numberchanns);
                    
                break; % terminate secondary loop; don't enter the questionable spike in fittedspikes
            end
            
 

            % Only count spike if amplitude exceeds threshold. Subtract
            % it either way, though.
            
            if bestA*min(TemplateMatrix(:,bestmu)) < noisethresh(bestmu)            
                totalspikes = totalspikes+1;  % this counts all the spikes found in the entire dataset

                spikeTimesThisEvent = [spikeTimesThisEvent, shifts(bestshift)+peaktime+t1];
                spikeIDsThisEvent=[spikeIDsThisEvent bestmu];
                spikeAmplThisEvent=[spikeAmplThisEvent bestA];
                allrightscorelist=[allrightscorelist Pmarg];
                           
            end
            
            % add the fitted spike to fittedspikes and subtract it from waveform:
            fitspike = reconstructSpike(bestA,peaktime,bestTemp,1:size(spike1,1),onewaveformlength,numberchanns);                                    
            
            
                 
            if makeplots 
                tmax = size(originalspike,1);
                nscalep = max(originalspike(:)); nscalem = min(originalspike(:));
                figure(1);
                for ni=1:numberchanns
                    subplot(6,5,ni);
                    plot(1:tmax,spike1(:,ni),'k',1:tmax,fitspike(:,ni),'r');%,1:tmax,fittedspikes(:,ni),'b');
                    xlim([1,tmax]); ylim([nscalem,nscalep]);
                end  
                pause;       

            end

            fittedspikes = fittedspikes + fitspike;
            spike1 = spike1 - fitspike;
            
            if strcmp(exitcode,'AM'), break;end;        % don't fit any more spikes if this one was ambiguous
        else

            exitcode='AI'; break;                       % terminate secondary loop "All Ineligible"
        end
    end % end secondary loop; try fitting another spike to residual

    % Catch and log incomplete fits
    if strcmp(exitcode,'AI') || strcmp(exitcode, 'LR')
        newIncompFit.peakchan = peakchan;
        
        mint = max(1,peaktime-floor(onewaveformlength/2));
        maxt = min(size(spike1,1), peaktime + ceil(onewaveformlength/2)-1);
        newIncompFit.S = zeros(onewaveformlength,numberchanns);
        newIncompFit.S(floor(onewaveformlength/2)+1-(peaktime-mint) : floor(onewaveformlength/2) + 1 + (maxt - peaktime), :) = spike1(mint:maxt,:);
        newIncompFit.S = reshape(newIncompFit.S,1,onewaveformlength*numberchanns);
                
        newIncompFit.chmask = false(1,numberchanns);
        newIncompFit.chmask(neighbors{peakchan}) = true;
        incompFits = [incompFits, newIncompFit];
    end
    
    if secondloop==Nmaxfitspike, exitcode='TM'; end                         % upper bound on number of spikes was reached
    Nfoundspikes = length(spikeIDsThisEvent);  % this counts the spikes found in this one event
    
    fitTau = spikeTimesThisEvent - t1 - 16;
    fitIDs = spikeIDsThisEvent;
    
    [spikeTimesThisEvent,tmpreorder]=sort(spikeTimesThisEvent);
    spikeIDsThisEvent = spikeIDsThisEvent(tmpreorder);
    spikeAmplThisEvent = spikeAmplThisEvent(tmpreorder);

    spikeTimesAllEvents=[spikeTimesAllEvents spikeTimesThisEvent];
    spikeIDsAllEvents=[spikeIDsAllEvents spikeIDsThisEvent];
    spikeAmplAllEvents=[spikeAmplAllEvents spikeAmplThisEvent];


   
end % end main loop; get another trigger

tmp=toc;
disp(['total fit spikes = ' num2str(totalspikes) '; time per fit spike = ' num2str(tmp/totalspikes)]);
amplist = cell(1,nclusters);
SpikeTimes = cell(1,nclusters);

for i=1:nclusters    
    SpikeTimes{i} = spikeTimesAllEvents(spikeIDsAllEvents == i);
    amplist{i} = spikeAmplAllEvents(spikeIDsAllEvents == i);
end

end




function [spike1,t1,n] = concatenate(S,T,n)
% concatenate overlapping triggers
    global onewaveformlength; global numberchanns;
    spike1 = reshape(S(n,:), onewaveformlength, numberchanns);    
    t1 = T(n) - floor(onewaveformlength/2);
    iter=1;
    while (n+1)<=length(T) && T(n+1)-T(n)<onewaveformlength % does next trigger overlap?                
        overlappingspike=reshape(S(n+1,:), onewaveformlength, numberchanns);  
        spike1 = [spike1(1:T(n+1)-T(n),:); overlappingspike];
        n = n+1;
        %display(['iter = ',num2str(iter)]);
        iter=iter+1;
    end
    n = n+1;
end







