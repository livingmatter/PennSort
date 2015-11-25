function [K,noise,C,data] = getNoiseAndWhiteningFilter(noisefileprefix,filenums,group,thresh,loadfun,array)

    global onewaveformlength; global numberchanns;
    
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
                    
    %data = new_LoadAndFilterPenn([noisefileprefix,num2str(filenums(1))]);
    data = loadfun([noisefileprefix,num2str(filenums(1))]);
    data = data(ChannelGroups{group},:);  
    spatialnoise = data(:,min(data)>-thresh & max(data)<thresh);
    C = cov(spatialnoise');
    K = inv(sqrtm(C));
    K = K / min(diag(K));
    K = kron(K,eye(onewaveformlength));
    
    noise = [];
    for f = filenums        
        %data = new_LoadAndFilterPenn([noisefileprefix,num2str(f)]);
        data = loadfun([noisefileprefix,num2str(filenums(1))]);
        data = data(ChannelGroups{group},:);  
    
        % Get noise matrix
        subthresh = max(abs(data)) < thresh;
        markers = find(abs(diff(subthresh))==1)+1;
        t = markers(diff(markers)>onewaveformlength+10);
        t = t(subthresh(t));

        noisetmp = zeros(length(t),onewaveformlength*numberchanns);
        for i=1:length(t)
            temp = data(:,t(i):t(i)+onewaveformlength+9);
            temp = temp(:,11:end);
            noisetmp(i,:) = reshape(temp',1,onewaveformlength*numberchanns);
            %noise(i,:) = reshape(data(:,t(i):t(i)+onewaveformlength+9)',1,onewaveformlength*numberchanns);
        end
        noise = [noise; noisetmp*K];
    end
end
