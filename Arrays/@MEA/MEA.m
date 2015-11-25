classdef MEA
    properties (SetAccess = private)
        GroupNum;       % Number of electrode groups.
        ChNum;          % Array with GroupNum entries, giving number of channels in each group.
        ChannelGroups;  % Cell array with GroupNum entries, each containing list of channel numbers belonging to that group.
        Geometry;       % Cell array with GroupNum entries. Each entry has two columns, the x and y coordinates of each electrode.
    end
    properties (SetAccess = private, GetAccess = private)
        ChIndex;        
        
    end
    
    methods
        function obj = MEA(geomfun)
            % GroupNum = number of channel groups
            % geomfun = function handle to a function returning a matrix
            % giving electrode x-y coordinates.
            
            [obj.Geometry, obj.GroupNum, obj.ChannelGroups] = geomfun();
            %obj.Geometry = cell(1,GroupNum);
            obj.ChNum = zeros(1,obj.GroupNum);
            for n=1:obj.GroupNum
            %    obj.Geometry{n} = geomfun(n);
                obj.ChNum(n) = length(obj.ChannelGroups{n});
            end
            
            obj.ChIndex = zeros(obj.GroupNum, max(obj.ChNum));
            index = 1;
            for i=1:obj.GroupNum
                for j=1:obj.ChNum(i)
                    obj.ChIndex(i,j) = index;
                    index = index+1;
                end
            end           
            
        end
        
        function [r] = XY(obj,group,ch)
            r = obj.Geometry{group}(ch,:);
        end
        
        function [d] = Distance(obj,group1,ch1,group2,ch2)
            r1 = obj.XY(group1,ch1);
            r2 = obj.XY(group2,ch2);            
            d = norm(r1-r2);
        end
        
        function [out] = isNeighbor(obj,group1,ch1,group2,ch2,thresh)
            if nargin < 6
                thresh = 1.5*30;
            end
            out = obj.Distance(group1,ch1,group2,ch2) <= thresh;
        end
        
        function [nbrgrp,nbrch] = Neighborhood(obj,group,ch,type,thresh)
            
            d = Inf(obj.GroupNum, max(obj.ChNum));
            for i=1:obj.GroupNum
                for j=1:obj.ChNum(i)
                    d(i,j) = obj.Distance(i,j,group,ch);
                    
                end
            end
            switch type
                case 'distance'                    
                    if nargin<5
                        thresh = 1.5*30;
                    end
                    %[nbrgrp,nbrch] = find(d <= thresh);
                    ix = find(d<=thresh);
                case 'channels'
                    if nargin < 5
                        thresh = 9;
                    end
                    [d_sort,order] = sort(d(:));
                    ix = order(1:thresh,:);
                    %[nbrgrp,nbrch] = ind2sub(size(d),order(1:thresh,:));
            end
            [nbrgrp,nbrch] = ind2sub(size(d),ix);
            nbrgrp = nbrgrp(:); nbrch = nbrch(:);
        end
                
        function Show(obj,ch)
            C = [];
            R = [];
            for i=1:obj.GroupNum
                for j=1:obj.ChNum(i)
                    R = [R; obj.XY(i,j)];
                    C = [C; 0,0,1];
                end
            end           
            
            if nargin==2
                for i=1:size(ch,1)
                    C(obj.ChIndex(ch(i,1),ch(i,2)),:) = [1,0,0];
                end
            end
            
            scatter(R(:,1),R(:,2),[],C,'MarkerEdgeColor','k','MarkerFaceColor','k');
            xlim([min(R(:,1))-30, max(R(:,1))+30]);
            ylim([min(R(:,2))-30, max(R(:,2))+30]);
        end
        
        function [handles] = makeAxes(obj,parent)
             group = 1;

            x = obj.Geometry{group}(:,1);
            y = obj.Geometry{group}(:,2);
            x = (x - min(x)) ./ range(x);
            y = (y - min(y)) ./ range(y);
            dx_min = zeros(size(x));
            dy_min = zeros(size(y));
            for i=1:length(x)
                if ~isempty(x(x>x(i)))
                    dx_min(i) = min(x(x > x(i))) - x(i);
                else
                    dx_min(i) = Inf;
                end
                if ~isempty(y(y>y(i)))
                    dy_min(i) = min(y(y > y(i))) - y(i);
                else
                    dy_min(i) = Inf;
                end
            end
            w = min(dx_min);
            h = min(dy_min);
            x = x ./ (1+w);
            y = y ./ (1+h);
            w = w / (1+w);
            h = h / (1+h);
            handles = zeros(1,obj.ChNum(group));
            for ch=1:obj.ChNum(group)
                handles(ch) = axes('Parent',parent,'Position',[x(ch), max(y) - y(ch), w, h]);
            end
        end
    end
end