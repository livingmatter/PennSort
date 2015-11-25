function [TemplateMatrix,clusterind] = makeTemplates(fileprefix,nclusters,array,group,ndecimate)


%allparameters;
global onewaveformlength; global numberchanns;
deltat = 1/ndecimate;
mid=floor((onewaveformlength-1)/2)*ndecimate + 1;
spikenbd=21;
smoothedtime=length([1:deltat:onewaveformlength]); 
ccwidth=13; % one half the zone of allowed shifts
trimwidth=1; trimup=trimwidth/deltat;
TemplateMatrix=zeros(smoothedtime*numberchanns,nclusters); %place all templates here
% main loop
clusterind = [];

if length(nclusters)==1, goodclusters = 1:nclusters; else goodclusters=nclusters; end
for i=1:length(goodclusters)
%for examinethis = goodclusters
    examinethis = goodclusters(i);
%for examinethis=1:nclusters%1:min(clusterstoshow,nclusters),     
    
    disp(['examining ' num2str(examinethis)]);
%    if (whichset==105)&&ismember(examinethis,masklist), disp('masking');end 
    besttime=[]; %accumulate time shifts here
    
    load([fileprefix,num2str(examinethis),'.mat']);
    if isempty(C), disp(['skipping #' num2str(examinethis)]); continue;end  % skip this cluster if empty
   % Assign a peak channel to each cluster. Find some neighborhood of it, then
   % take median of waveforms within this neighborhood. This is the draft
   % template.
    numevents = size(C,1);
    C_reshaped = reshape(C',onewaveformlength,numberchanns,numevents);
    [m,nleader] = min(min(median(C_reshaped,3)));
    [nbrgrp,includedchannels] = array.Neighborhood(group,nleader,'channels',spikenbd);
    alltemplates = zeros(onewaveformlength,numberchanns,numevents);
    temp = reshape(C',onewaveformlength,numberchanns,numevents);
    alltemplates(:,includedchannels,:) = temp(:,includedchannels,:);
    drafttemplate = median(alltemplates,3);
    
 %   load([fprefix 'clusters/Templates/alltemplate' num2str(examinethis) '.mat']);% obtain 'drafttemplate', 'includedchannels', 'alltemplates', 'nleader', 'nmean','allpeakchs'
    
%     nsizeall=size(alltemplates);
%     if length(nsizeall)==2, numevents=1; disp('only one exemplar'); else, numevents=nsizeall(3) 
%     end
%     if (whichset==105)&&ismember(examinethis,masklist), % is this a cluster that we have chosen to mask?
%         for np=1:numberchanns,
%             tmpmask=1;
%             for nm=1:length(masklist),
%                 if masklist(nm)==examinethis,
%                     nleader=maskleaders(nm);       % after masking generally we need to choose the leader channel manually
%                     tmpmask=masks(nm,np);break; end% get a 1 or 0
%             end;
%             if tmpmask~=1,
%                 drafttemplate(:,np)=drafttemplate(:,np)*tmpmask;
%                 for nq=1:length(includedchannels),
%                     if includedchannels(nq)==np, alltemplates(:,nq,:)=alltemplates(:,nq,:)*tmpmask; end
%                 end
%             end
%         end;
%     end % end if whichset...
%    save([fprefix 'clusters/Templates/allmaskedtemplate' num2str(examinethis) '.mat'],'alltemplates', 'drafttemplate')
    
   % draft=drafttemplate(:,includedchannels); % restrict attn to those channels with content
    
    smoothdraft=interp1([1:onewaveformlength],drafttemplate,[1:deltat:onewaveformlength],'spline');
%    smoothdraft=interp1([1:onewaveformlength],draft,[1:deltat:onewaveformlength],'linear');
    smoothdraft([1:trimup,end-trimup+1:end],:)=0;  % remove first and last 0.1us
    %smoothdraft(end-trimup+1:end,:)=zeros(trimup,nnumchann);
% center the smoothed draft template to its absolute peak 
    [tmp,tmptimes]=min(smoothdraft);  % for each channel get the time when it hits peak
    [ampl,peakchan]=min(tmp);         % get absolute peak
    peaktime=tmptimes(peakchan);
    midpoint=floor((onewaveformlength-1)/2)/deltat + 1; % where we want abs peak
    shift=peaktime-midpoint;
    tempsmooth=zeros(smoothedtime,numberchanns);
    if shift>0,
        tempsmooth((1:(smoothedtime-shift)),:)=smoothdraft(((1+shift):smoothedtime),:);
        smoothdraft(:,:)=tempsmooth(:,:);end;
    if shift<0,    %of course do nothing if shift==0
        tempsmooth(((1-shift):smoothedtime),:)=smoothdraft((1:(smoothedtime+shift)),:);
        smoothdraft(:,:)=tempsmooth(:,:);end;
% %! show selected draft template:
%     nscalep=max(smoothdraft(:));nscalem=min(smoothdraft(:));
%     figure(3);
%     for panel=[1:nnumchann],
%         subplot(6,5,includedchannels(panel));plot([1:deltat:onewaveformlength],smoothdraft(:,panel));ylim([nscalem,nscalep]);xlim([1,onewaveformlength]);
%     end;
%     subplot(6,5,1);text(4,0,['draft of cluster' num2str(examinethis)],'FontSize',17);
% increase time resolution (upsample)
    smoothall=interp1([1:onewaveformlength],alltemplates,[1:deltat:onewaveformlength],'spline');
    % smoothall=interp1([1:onewaveformlength],alltemplates,[1:deltat:onewaveformlength],'linear');
    smoothall(1:trimup,:,:)=zeros(trimup,numberchanns,numevents);
    smoothall(end-trimup+1:end,:,:)=zeros(trimup,numberchanns,numevents);% remove first and last 0.1us
% cross-correlate each waveform to the draft template and shift
    nmaxmax=-10000;nminmax=0;keeperlist=[];
% SECONDARY LOOP
    for nevent=1:numevents,
          nscalep=max(max(smoothall(:,:,nevent)));nscalem=min(min(smoothall(:,:,nevent)));
          nmaxmax=max(nmaxmax,nscalem);nminmax=min(nminmax,nscalem);
%!       show all waveforms in selected cluster (prior to aligning)
%         if (examinethis==40),
%             figure;
%             for panel=[1:nnumchann],
%                 subplot(6,5,includedchannels(panel));plot([1:deltat:onewaveformlength],smoothall(:,panel,nevent));ylim([nscalem,nscalep]);xlim([1,onewaveformlength]);
%             end;
%             text(4,0,[num2str(examinethis) ':' num2str(nevent)],'FontSize',17);
%         end;%
% find crosscorrelation in leader channel only
        crosscorr=zeros(1,1+(2*ccwidth/deltat));        
        for tau=(-ccwidth/deltat):(ccwidth/deltat),
            sumindex=[max(1,1-tau):min(smoothedtime,smoothedtime-tau)];sumindextau=sumindex+tau;
            crosscorr(tau+(ccwidth/deltat)+1)=sum(smoothdraft(sumindex,nleader).*smoothall(sumindextau,nleader,nevent));
        end; %end for tau loop
%         if (examinethis==3)&&(nevent==70),figure;plot([1:1+16/deltat],crosscorr);title([num2str(examinethis) ':' num2str(nevent)]);end
        [maxcor,best]=max(crosscorr);
        best=best-((ccwidth/deltat)+1);  %relative to zero
        if (best<(-(ccwidth-1)/deltat))||(best>(ccwidth-1)/deltat),
        % (We can't keep any event that can't be aligned. If we get this error need to discard some events.)
            disp(['event ' num2str(nevent) ' twiddle out of range=' num2str(best)]);
        else;
            keeperlist=[keeperlist nevent];
        end
        besttime(nevent)=best; %save for later
% shift the entire event to align:
        tempsmooth=zeros(smoothedtime,numberchanns);
        if best>0,
            tempsmooth([1:(smoothedtime-best)],:)=smoothall([(1+best):smoothedtime],:,nevent);
            smoothall(:,:,nevent)=tempsmooth(:,:);end;
        if best<0,    %of course do nothing if best==0
            tempsmooth([(1-best):smoothedtime],:)=smoothall([1:(smoothedtime+best)],:,nevent);
            smoothall(:,:,nevent)=tempsmooth(:,:);end;
    end; %end for nevent (SECONDARY LOOP)
    if length(keeperlist)==0
        display(['Skipping ', num2str(i),': no keepers']);
    else
        varbesttime=var(besttime); %besttime now has the shifts for each event in this cluster

    % average the aligned traces 
        %TemplateMatrix=zeros(smoothedtime,numberchanns);
    %    consensus=sum(smoothall,3)/numevents;
        consensus=median(smoothall(:,:,keeperlist),3);
        [tmp,tmptimes]=min(consensus(:,:));  % for each channel get the time when it hits peak
        [ampl,peakchan]=min(tmp);         % get absolute peak
        peaktime=tmptimes(peakchan);      % time of abs peak

        % Center template on abs. min
        consensus = [zeros(ndecimate,numberchanns); consensus; zeros(ndecimate,numberchanns)];
        consensus = circshift(consensus, [mid-peaktime,0]);
        consensus = consensus((ndecimate+1):(smoothedtime+ndecimate),:);

        TemplateMatrix(:,i)=consensus(:);
        clusterind = [clusterind, i];
    end
    
%     disp([' template ' num2str(examinethis) ' template peak time is ' num2str(peaktime-midpoint)]);
%     nscalep=max(consensus(:));nscalem=min(consensus(:));
%     if (whichset==105)&&ismember(examinethis,masklist), disp('masking');end %MOVED
%     for panel=1:length(includedchannels),
% %         tmpmask=1;
% %         if (whichset==105)&&ismember(examinethis,masklist), %MOVED
% %             for nm=1:length(masklist),
% %                 if masklist(nm)==examinethis,
% %                 tmpmask=masks(nm,includedchannels(panel));break; % get a 1 or 0
% %                 end;end;
% %         end;
%         TemplateMatrix(:,includedchannels(panel))=consensus(:,panel);%*tmpmask;
% %         if makeplots||(ismember(examinethis,watchlist)),figure(2);
% %             subplot(6,5,includedchannels(panel));plot([1:deltat:onewaveformlength],consensus(:,panel) ... %*tmpmask
% %                 );ylim([nscalem,nscalep]);xlim([1,onewaveformlength]);end;
%     end;
%     if makeplots||(ismember(examinethis,watchlist)),subplot(6,5,1);text(4,0,['consensus' num2str(examinethis)],'FontSize',9);
%         text(4,-30,[num2str(numevents) ' exemplar'],'FontSize',9);
%         pause;end
%     for panel=1:numberchanns,subplot(6,5,includedchannels(panel),'replace');    end
%     if (examinethis==1000)||makemultiplots, figure(1);
%         for panel=1:length(includedchannels), subplot(6,5,includedchannels(panel),'replace'); end
%     end

% repackage as all numberchanns traces (not just the nnumchann interesting ones); write the average to a file
    %save([fprefix 'clusters/Templates/template' num2str(examinethis) '.mat'],'TemplateMatrix', 'includedchannels', 'besttime');
    %save([fprefix 'clusters/Templates/exemplars' num2str(examinethis) '.mat'],'smoothall');
end %end examinethis (MAIN LOOP)

TemplateMatrix = TemplateMatrix(:,clusterind);