function [bigS,bigT,bigchmask,bigpeakch] = extractSpikes(fileprefix,...
                                filelist,mode,groups,threshold,loadfun,array)
    
    global onewaveformlength; global numberchanns;
    
    filelen = 30*1e4;    
    %array = MEA(2,@geom_2F30EL);
    neighbors = cell(1,numberchanns);
    for i=1:numberchanns
         [nbrgrp,nbrch] = array.Neighborhood(1,i,'distance',30*sqrt(2));
         neighbors{i} = nbrch;
    end

%     ChannelGroups = cell(1,2);
%     ChannelGroups{1} = [45,38,16,17,25, ...
%                         55,37,46,18,26, ...
%                         53,54,24,19,27, ...
%                         52,59,48,08,00, ...
%                         51,58,49,09,01, ...
%                         50,57,56,47,39] + 1;
% 
%     ChannelGroups{2} = [33,41,42,43,36, ...
%                         32,40,35,31,23, ...
%                         20,28,34,44,30, ...
%                         02,10,05,29,22, ...
%                         03,11,06,15,21, ...
%                         04,12,13,14,07] + 1;
    ChannelGroups = array.ChannelGroups;
    
    bigS = [];
    bigT = [];
    bigchmask = [];
    bigpeakch = [];
    for f=filelist        
        % Load data file
        display(['Processing ',fileprefix,num2str(f),'...']);
        %y = new_LoadAndFilterPenn([fileprefix,num2str(f)]);
        [y,time] = loadfun([fileprefix,num2str(f)]);
        
        for group=groups
            data = y(ChannelGroups{group},:);
            clear y;
            
            % Extract spikes
            eventlengths = [];
            switch mode
                case 'fitting'
                    crossings = diff(min(data)<threshold);
                    RisingCrossings = find(crossings == 1);                    
                    FallingCrossings = find(crossings == -1);
                    
                    if FallingCrossings(1) < RisingCrossings(1)
                        FallingCrossings = FallingCrossings(2:end);
                    end
                    T = [];
                    peakch = [];
                    lastT = 0;
                    %peakch = zeros(size(RisingCrossings));
                    for i=1:min(length(RisingCrossings),length(FallingCrossings))
                        [minV,dT] = min(min(data(:,RisingCrossings(i):FallingCrossings(i))));
                        nextT = RisingCrossings(i) + (dT-1);
                        %if nextT > lastT + onewaveformlength
                            T = [T,nextT];
                            [minV,minch] = min(data(:,nextT)');
                            peakch = [peakch,minch];
                            lastT = nextT;
                        %end
                        eventlengths = [eventlengths, FallingCrossings(i)-RisingCrossings(i)];
                    end
                    
                case 'clustering'
                    % Breadth-first search to find connected components:
                    % While there are unprocessed elements:
                    %  Add first unprocessed threshold crossing to queue.
                    %  While queue is nonempty do:
                    %   Pop top element of queue; remove from unprocessed; move to queue_out.
                    %   Enqueue all unprocessed, threshold-crossing, spatiotemporal neighbors.
                    %  return
                    %  queue_out now stores a full connected component. Find absolute peak,
                    %  empty queue_out
                    % return.
                    %fig = figure;
                    T = [];
                    peakch = [];
                    unprocessed = find(data<threshold);
                    lastT = zeros(1,numberchanns);
                    
                    tprev = -Inf;
                    chprev = [];
                    count=0;
                    while ~isempty(unprocessed)
                        queue = [];
                        queue(1) = unprocessed(1);
                        unprocessed = unprocessed(2:end);
                        queue_out = [];
                        while ~isempty(queue)
                            [ch0,t0] = ind2sub(size(data),queue(1));
                            queue_out = [queue_out; queue(1)];
                            % Find spatiotemporal neighbors
                            inds = sub2ind(size(data),neighbors{ch0},t0*ones(size(neighbors{ch0})));
                            if t0==1
                                inds = [inds; sub2ind(size(data),ch0,t0+1)];
                            elseif t0==size(data,2)
                                inds = [inds; sub2ind(size(data),ch0,t0-1)];
                            else
                                inds = [inds; sub2ind(size(data),[ch0;ch0],[t0-1;t0+1])];
                            end
                            % Restrict to unprocessed threshold crossings
                            [ism,loc] = ismember(inds,unprocessed);
                            unprocessed(loc(loc~=0)) = [];
                            % Pop first element and enqueue inds
                            queue = [queue(2:end); inds(ism)];
                        end

                        [m,minind] = min(data(queue_out));
                        [clustch, clustt] = ind2sub(size(data), queue_out);
                        
                        [chmin,tmin] = ind2sub(size(data),queue_out(minind));
                        midwaveform = floor(onewaveformlength/2);
%                         Stoshow = data(:,tmin-midwaveform:tmin+((onewaveformlength-midwaveform)-1));
%                         Stoshow = reshape(Stoshow', onewaveformlength*size(data,1),1);
%                         if tprev < tmin && tmin < tprev + 10
%                             count = count+1;
%                             if count==27
%                                 show_waveforms(array,Stoshow,fig,unique(clustch),chprev); pause;
%                             end
%                             
%                         end
                        tprev = tmin; chprev = unique(clustch);
                        if tmin > lastT(chmin)+8
                            T = [T,tmin];
                            peakch = [peakch, chmin];                        
                            lastT(chmin) = tmin;
                        end
                        
                    end               
            end

            midwaveform = floor(onewaveformlength/2);
            peakch = peakch(T>midwaveform & T<=size(data,2)-(midwaveform-1));
            T = T(T>midwaveform & T<=size(data,2)-(midwaveform-1));
            [T,sortind] = sort(T,'ascend');
            peakch = peakch(sortind);

            S = zeros(length(T),onewaveformlength*numberchanns);
            for i=1:length(T)
                 S(i,:) = reshape(data(:,T(i)-midwaveform:T(i)+...
                     ((onewaveformlength-midwaveform)-1))',1,onewaveformlength*numberchanns);
            end


            T = T + (f-1)*filelen;
            switch mode
                case 'clustering'
                    chmask = false(length(T),numberchanns);
                    for i=1:length(T)
                        nbrs = neighbors{peakch(i)}; 
                        chmask(i,nbrs) = true;
                    end
            end
        end
        bigT = [bigT, T];
        bigpeakch = [bigpeakch, peakch];
        bigS = [bigS; S];
        if strcmp(mode,'clustering')
            bigchmask = logical([bigchmask; chmask]);
        end
    end
    display('*** Finished! ***');


  