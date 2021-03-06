function [G, GroupNum, ChannelGroups] = geom_2F30EL(n)
    GroupNum = 2;
    
    chnum=(1:30)';
    
%     G = 30.*([rem((chnum-1),5),floor((chnum-1)./5)]);
%     if n==2
%         G(:,1) = G(:,1) + 320;
%     end
    
    G = cell(1,2);
    for i=1:2
        G{i} = 30.*([rem((chnum-1),5),floor((chnum-1)./5)]);       
    end
    G{2}(:,1) = G{2}(:,1) + 320;
    
    ChannelGroups = cell(1,2);
    ChannelGroups{1} = [45,38,16,17,25, ...
                        55,37,46,18,26, ...
                        53,54,24,19,27, ...
                        52,59,48,08,00, ...
                        51,58,49,09,01, ...
                        50,57,56,47,39] + 1;

    ChannelGroups{2} = [33,41,42,43,36, ...
                        32,40,35,31,23, ...
                        20,28,34,44,30, ...
                        02,10,05,29,22, ...
                        03,11,06,15,21, ...
                        04,12,13,14,07] + 1;