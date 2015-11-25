function [fig] = show_waveforms(array, S, fig, color1, color2)

if nargin<3
    fig = figure;
end
figure(fig);

if nargin < 4    
    color1 = []; color2 = [];
end

handles = array.makeAxes(fig);
S = reshape(S,size(S,1)/length(handles),length(handles),size(S,2));
for ch=1:length(handles)
    
    if ismember(ch,color1)
        plot(handles(ch),squeeze(S(:,ch,:)),'color','r');%,'LineWidth',2);        
    elseif ismember(ch,color2)
        plot(handles(ch),squeeze(S(:,ch,:)),'color','b');%,'LineWidth',2);        
    else
        plot(handles(ch),squeeze(S(:,ch,:)),'color','k');%,'LineWidth',2);
    end
%     if ch==1
%         xlabel('Time(ms)'); ylabel('Voltage (\mu V)');
%         set(handles(ch),'XTick',0:0.5:3.2,'YTick',-500:50:50,'XLim',[1,size(S,1)],'YLim',[-500,50]);
%     else
    set(handles(ch),'XLim',[1,size(S,1)],'YLim',[min(S(:))-50,max(S(:))+50]);
        %set(handles(ch),'TickDir','in','TickLength',[0.01,0.01],'LineWidth',1,'XTick',10:10:30,'YTick',-300:100:100,'XLim',[1,size(S,1)],'YLim',[-400,100],'XTickLabel',[],'YTickLabel',[]);
    end
end
