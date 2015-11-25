function [G,GroupNum,ChannelGroups] = geom_tetrode(n)

GroupNum = 1;
G{1} = [-1, 0;
         0,-1;
         0, 1;
         1, 0];

ChannelGroups{1} = 1:4;