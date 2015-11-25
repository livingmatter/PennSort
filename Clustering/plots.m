function varargout = plots(varargin)
% PLOTS M-file for plots.fig
%      PLOTS, by itself, creates a new PLOTS or raises the existing
%      singleton*.
%
%      H = PLOTS returns the handle to a new PLOTS or the handle to
%      the existing singleton*.
%
%      PLOTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTS.M with the given input arguments.
%
%      PLOTS('Property','Value',...) creates a new PLOTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plots_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plots_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plots

% Last Modified by GUIDE v2.5 19-Aug-2009 14:45:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plots_OpeningFcn, ...
                   'gui_OutputFcn',  @plots_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%guidata(hObject,handles)

% --- Executes just before plots is made visible.
function plots_OpeningFcn(hObject, eventdata, handles, varargin)
global modify_rectangles;
global draw_average;
global plot_position;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plots (see VARARGIN)

%Create the data to plot

%ishandle(hObject)

%{
str = [];
for i = 1:1050 / 5
str = [str, '%*f %*f %*f %*f %f ' ];
end
str = [str, '%f %f \"%d\"'];
fid = fopen('big-10000-15.cri','r');
disp 'getting data'
handles.spikedata = fscanf(fid,str,[1050 /5 + 3,100]);   %1050
%}


% Choose default command line output for plots
handles.output = hObject;
%handles.spikedata = varargin{1};
%handles.datalength = size(handles.spikedata,2);
%handles.startpoint = 1;
%handles.stepsize = 10;
%handles.blocksize = handles.datalength;
handles.cluster = [];
handles.i = 0;
handles.active = [];
handles.nr_areas = 0;
handles.close = 0;
handles.h = [];             %contains the handles to the rectangles
handles.figure2 = 0;        %figure for the animation
handles.current_spike = 0;
handles.current_frame = 1;
handles.spikedata = [];
position = get(0,'ScreenSize');
plot_position = [(position(3) - position(3)/2-500) (position(4)-position(4)/3+500) 500 300];

modify_rectangles = 0;
draw_average = 0;
handles.t = timer;
set(handles.t,'TimerFcn',{@update_plot, handles.figure1})
set(handles.t,'ExecutionMode','fixedDelay','BusyMode','drop','Period',0.2)





%startpoint = handles.startpoint;
%blocksize = handles.blocksize;
%endpoint = startpoint + blocksize - 1;
%colormap(jet(200));
%axes(handles.axes1)
%image([1,size(handles.spikedata,2)],[1,1053],20 - handles.spikedata(:,1:size(handles.spikedata,2)));
%zoom xon
%pan xon
%uipushtool3_ClickedCallback(hObject, eventdata, handles)
guidata(hObject,handles);
% 
% if length(varargin)==1
%     load_data(hObject,varargin{1},handles);
% end
    
% UIWAIT makes plots wait for user response (see UIRESUME)
 %uiwait(handles.figure1);
uiwait;


% function load_data(hObject,filename,handles)
%     if (filename ~= 0)
%         %filename = sprintf('%s%s',pname,fname);
%         load(filename);
% 
%         for i = 1:length(handles.cluster)
%             delete(handles.h(i))
%         end
%         %clear handles.cluster
% 
%         handles.cluster = [];
%         handles.i = 0;
%         handles.active = [];
%         handles.nr_areas = 0;
%         handles.close = 0;
%         handles.h = [];             %contains the handles to the rectangles
%         handles.figure2 = 0;        %figure for the animation
%         handles.current_spike = 0;
%         handles.current_frame = 1;
%         modify_rectangles = 0;
%         draw_average = 0;
% 
%         handles.spikedata = [S(ind,:),ind]';
% 
%         handles.frames = size(S,2)/30;
% 
%         set(handles.figure1, 'Name',['SPIKE Tool', '    ',filename])
% 
%         handles.active = 0;
%         handles.figure2 = 0;
% 
% 
%         colormap(jet(250));
%         axes(handles.axes1)
%         image([1,size(handles.spikedata,2)],[1,size(S,2) + 1],20 - handles.spikedata(:,1:size(handles.spikedata,2)));
%         zoom xon
% 
%         guidata(hObject,handles)
%     end


% --- Outputs from this function are returned to the command line.
function varargout = plots_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
handles = guidata(hObject);

if ~isempty(handles.cluster) && ~isempty(handles.spikedata)
    for i = 1:length(handles.cluster)
        data{i} = round(handles.cluster(i).pos(1)):round(handles.cluster(i).pos(2));
%        data{i} = handles.spikedata(end,round(handles.cluster(i).pos(1)):round(handles.cluster(i).pos(2)));
    end
    varargout{1} = data;
else varargout{1} = [];
end
%h = gcf;
%saveas(gcf,'plotsGUI','pdf');
%savefig('draftGUIcmyk.tiff','tiff', '-cmyk','-r600')

delete(gcf)







% --- Executes during object creation, after setting all properties.
function popupmenue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% -gets executed when the rectangles should bee drawn
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
set(handles.uitoggletool7,'Enable','off')
set(handles.uitoggletool8,'Enable','off')
set(handles.uitoggletool12,'Enable','off')
but = 1;
num = 0;
i = 0;
while but == 1
    i = i + 1;

[x(i) y(i) but] = ginput(1);

num = num + 1;
if mod(num,2) == 0  
    
handles.i = handles.i + 1;
handles.cluster(handles.i).pos = x(num-1:num);

guidata(hObject,handles)
ii = handles.i;
handles.nr_areas = handles.nr_areas + 1;
handles.h(handles.nr_areas) = rectangle('Position',[handles.cluster(ii).pos(1),0,handles.cluster(ii).pos(2)-handles.cluster(ii).pos(1),1053], ...
    'EdgeColor','none','EraseMode','xor','FaceColor',[0.8 0.8 0.8]);

set(handles.h(handles.nr_areas),'ButtonDownFcn',{@button_down, handles.figure1})
%handles.figure1
guidata(hObject,handles)

%spikedata_pushbutton_Callback(hObject, eventdata, handles)
end

end
set(handles.uitoggletool7,'Enable','on')
set(handles.uitoggletool8,'Enable','on')
set(handles.uitoggletool12,'Enable','on')
set(handles.uitoggletool7,'State','off')
set(handles.uitoggletool8,'State','off')
set(handles.uitoggletool7,'State','off')
set(handles.uitoggletool12,'State','off')
pan xon









%-gets executed when a rectangle is marked
function button_down(src,eventdata,handle_figure1)
% src - the object that is the source of the event
% evnt - empty for this property
global modify_rectangles;
global plot_position;
   %mmm = get(src,parent)
   %get(mmm,parent)
   sel_typ = get(gcbf,'SelectionType');
   handles = guidata(handle_figure1);
   
   if modify_rectangles == 2
       %if handles.active ~= 0
       %  set(handles.active,'Selected','off','EraseMode','xor','FaceColor',[0.8,0.8,0.8])
       %end
       
       if handles.active == src
           handles.active = 0;
       end
       delete(src)
       ind = find(handles.h == src);
       handles.h(ind) = [];
       handles.cluster(ind) = [];
       handles.i = handles.i - 1;
       handles.nr_areas = handles.nr_areas -  1;
       
       guidata(handle_figure1,handles)
       return    
   end
       
       
   %gets executed when modification is not blocked
   if modify_rectangles == 1
   %switch sel_typ 
      %case 'normal'
         if handles.active ~= 0
         set(handles.active,'Selected','off','EraseMode','xor','FaceColor',[0.8,0.8,0.8])
         end
         %disp('User clicked left-mouse button')
         set(src,'Selected','on')
         set(src,'EraseMode','xor','FaceColor',[0.6 0.6 0.6])
         handles.active = src;
         
        handles.current_spike = 0;
        handles.current_frame = 1;

        
        stop(handles.t)

        handles.figure2 = figure;
        set(handles.figure2,'OuterPosition',plot_position)
        set(handles.figure2,'DeleteFcn',{@animation_deleted,handle_figure1} )

       
        colormap(jet(250));

        
        set(handles.uitoggletool1,'Enable','off')
        set(handles.uitoggletool2,'Enable','off')
        set(handles.uitoggletool3,'Enable','off')
        set(handles.uitoggletool4,'Enable','off')
        set(handles.uitoggletool7,'Enable','off')
        set(handles.uitoggletool8,'Enable','off')
        set(handles.uitoggletool12,'Enable','off')
        set(handles.uipushtool3,'Enable','off')
        %set(handles.figure1,'HitTest','off')
        modify_rectangles = 0;
        guidata(handle_figure1,handles)
        start(handles.t)
        %guidata(handle_figure1,handles)
         % guidata(handle_figure1,handles)
        %guidata(handle_figure1,handles)
         
         %set(src,'FaceColor','b')
         %guidata(handle_figure1,handles)
         %return
      %{
case 'extend'
         disp('User did a shift-click')
         set(src,'Selected','on')
         set(src,'FaceColor','r')
      case 'alt'
         disp('User did a control-click')
         delete(src)
         set(src,'Selected','on')
         %set(src,'FaceColor','g')
         %set(src,'SelectionHighlight','on')
   end
   %}
   end

   
   
   
   
% --------------------------------------------------------------------
function uipushtool4_ClickedCallback(hObject,eventdata,handles)
global modify_rectangles;
global draw_average

set(handles.uitoggletool7,'State','on')


pan off
zoom off
datacursormode off

set(handles.uitoggletool1,'State','off')
        set(handles.uitoggletool2,'State','off')
        set(handles.uitoggletool3,'State','off')
        set(handles.uitoggletool4,'State','off')
        set(handles.uitoggletool8,'State','off')
        set(handles.uitoggletool12,'State','off')
        %handles.button_h = hObject;
guidata(hObject,handles)

draw_average = 0;
modify_rectangles = 1;







function update_plot(hObject,eventdata,handle_figure1,handle_axes2)
global draw_average

handles = guidata(handle_figure1);
figure(handles.figure2)
%-part which is esecuted, when the average should be drawn
if draw_average == 1;
    
    if handles.current_frame == handles.frames + 1
        
    rang = handles.cluster(handles.h == handles.active).pos;
    %figure(handles.figure2)
    s1 = subplot(1,2,1);
    image(zeros(6*30,5*30));
    axis image
    
    %title(['Spike: ',int2str(handles.current_spike - round(rang(1))) , ' / ' ,int2str(round(rang(2)) - round(rang(1)) + 1)]);

%{    
a = colorbar('peer',s1,[0.48 0.2 0.02584 0.6]);
ff = get(a,'YTick');
ff = ff - 20;
set(a,'YTickLabel',ff);
%}
    
    
    subplot(1,2,2);
    image(zeros(6,5))
    axis image
    drawnow
    %handles.current_spike = handles.current_spike + 1;
    handles.current_frame = 1;
    guidata(handle_figure1,handles)
    return
    end
    
    rang = handles.cluster(handles.h == handles.active).pos;
    p = handles.spikedata(1:handles.frames*30,round(rang(1)):round(rang(2)));
    %p1 = mean(p,2);
    p1 = median(p,2);
%    p1(end-2:end) = [];
    pp = reshape(p1,[handles.frames,5,6]);
    
%figure(handles.figure2)


[X,Y] = meshgrid(0:30:4*30,0:30:5*30);
[XI,YI] = meshgrid(-15:(4*30+15),-15:(5*30+15));
frame = handles.current_frame;
pp1 = squeeze(pp(frame,:,:))';

ZI = interp2(X,Y,pp1,XI,YI,'cubic');
s1 =  subplot(1,2,1);
image(20 - ZI);
%title(['Spike: ',int2str(handles.current_spike - round(rang(1))) , ' / ' ,int2str(round(rang(2)) - round(rang(1)) + 1)]);


%{
a = colorbar('peer',s1,[0.48 0.2 0.02584 0.6]);
ff = get(a,'YTick');
ff = ff - 20;
set(a,'YTickLabel',ff);
%}

axis image

subplot(1,2,2);
image(20 - pp1)
axis image
%axis off

drawnow
  handles.current_frame = handles.current_frame + 1;  
  guidata(handle_figure1,handles)  
    return
end


%-part which should be executeed when every single spike should be drawn
if handles.current_frame == handles.frames + 1
    rang = handles.cluster(handles.h == handles.active).pos;
    %figure(handles.figure2)
    
    
    spike = handles.current_spike;
    frame = handles.current_frame - 1;


    p = handles.spikedata(:,spike);
    p(end-2:end) = [];
    pp = reshape(p,[handles.frames,5,6]);

    [X,Y] = meshgrid(0:30:4*30,0:30:5*30);
    [XI,YI] = meshgrid(-15:(4*30+15),-15:(5*30+15));
    pp1 = squeeze(pp(frame,:,:))';
    ZI = interp2(X,Y,pp1,XI,YI,'cubic');
    
    
    
    
    s1 = subplot(1,2,1);
    image(20 - ZI);
    axis image
    
title(['Spike: ',int2str(handles.current_spike - round(rang(1)) + 1) , ' / ' ,int2str(round(rang(2)) - round(rang(1)) + 1)]);

%{
a = colorbar('peer',s1,[0.48 0.2 0.02584 0.6]);
ff = get(a,'YTick');
ff = ff - 20;
set(a,'YTickLabel',ff);
%}

    
    subplot(1,2,2);
    image(zeros(6,5))
    axis image
    %axis off
    drawnow
    
    handles.current_spike = handles.current_spike + 1;
    handles.current_frame = 1;
    guidata(handle_figure1,handles)
    return
end


rang = handles.cluster(handles.h == handles.active).pos;
if handles.current_spike == 0
   handles.current_spike = round(rang(1));
elseif handles.current_spike > round(rang(2))
   stop(handles.t)
   close(handles.figure2)
   %handles.figure2 = 0;
   %guidata(hObject,handles)
   return;
end

spike = handles.current_spike;
frame = handles.current_frame;


p = handles.spikedata(:,spike);
p(end-2:end) = [];
pp = reshape(p,[handles.frames,5,6]);

[X,Y] = meshgrid(0:30:4*30,0:30:5*30);
[XI,YI] = meshgrid(-15:(4*30+15),-15:(5*30+15));
pp1 = squeeze(pp(frame,:,:))';
ZI = interp2(X,Y,pp1,XI,YI,'cubic');


%figure(handles.figure2)



s1 = subplot(1,2,1);
image(20 - ZI);
axis image
title(['Spike: ',int2str(handles.current_spike - round(rang(1)) + 1) , ' / ' ,int2str(round(rang(2)) - round(rang(1)) + 1)]);




s2 = subplot(1,2,2);
image(20 - pp1)
axis image
%axis off

%{
a = colorbar('peer',s2,[0.48 0.2 0.02584 0.6]);
ff = get(a,'YTick');
ff = ff - 20;
set(a,'YTickLabel',ff);
%}

%align([s1,a,s2],'distribute','bottom')
drawnow

handles.current_frame = handles.current_frame + 1;

guidata(handle_figure1,handles)














%-Executes when user closes animation window
function animation_deleted(hObject,eventdata,handle_figure1)
global modify_rectangles;
global plot_position

handles = guidata(handle_figure1);
stop(handles.t)



delete(handles.figure2)
figure(handles.figure1)

%set(handles.figure1,'HitTest','on')
handles.figure2 = [];

 set(handles.uitoggletool1,'Enable','on')
        set(handles.uitoggletool2,'Enable','on')
        set(handles.uitoggletool3,'Enable','on')
        set(handles.uitoggletool4,'Enable','on')
        set(handles.uitoggletool7,'Enable','on')
        set(handles.uitoggletool8,'Enable','on')
        set(handles.uitoggletool12,'Enable','on')
        %set(handles.uitoggletool12,'State','off')
        %set(handles.figure1,'HitTest','off')
        
        %set(handles.uitoggletool7,'State','on')
        set(handles.uipushtool3,'Enable','on')
        
plot_position = get(hObject,'OuterPosition');

modify_rectangles = 1;
figure(handles.figure1)
%guidata(handles.figure1,handles)






% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selection = questdlg('Do you want to close the GUI?',...
                     'Close Request Function',...
                     'Yes','No','Yes');

switch selection,
   case 'Yes',
    %stop(handles.t)
    delete(handles.t)
    %clear handles.t
    %handles = guidata(hObject);
    %handles.output = handles.cluster;
   %varargout{1} =  handles.output;
   %guidata(hObject,handles)
    %varargout{2} = 'jan';
    %handles.close = 1;
    %guidata(hObject,handles)
    uiresume;
    %plots_OutputFcn(hObject, eventdata, handles) 
    %close(handles.figure1)
    %cc = handles.cluster;
    %delete(gcf)
    %close(gcf)
     
   case 'No'
     return
     
end

%guidata(hObject,handles)
%uiresume;










% --------------------------------------------------------------------
function uitoggletool3_ClickedCallback(hObject, eventdata, handles)
set(handles.uitoggletool7,'State','off')
set(handles.uitoggletool8,'State','off')
set(handles.uitoggletool12,'State','off')
%modify_rectangles = 0;
guidata(hObject,handles)
pan xon
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function uitoggletool2_ClickedCallback(hObject, eventdata, handles)
set(handles.uitoggletool7,'State','off')
set(handles.uitoggletool8,'State','off')
set(handles.uitoggletool2,'State','off')
%modify_rectangles = 0;
%guidata(hObject,handles)
zoom(0.5)

pan xon
% hObject    handle to uitoggletool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function uitoggletool7_OffCallback(hObject, eventdata, handles)
set(handles.uitoggletool7,'State','off')
global modify_rectangles
modify_rectangles = 0;
guidata(hObject,handles)
% hObject    handle to uitoggletool7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function uitoggletool4_ClickedCallback(hObject, eventdata, handles)
global modify_rectangles

pan off
zoom off
set(handles.uitoggletool7,'State','off')
set(handles.uitoggletool8,'State','off')
set(handles.uitoggletool12,'State','off')
datacursormode(handles.figure1,'on')

modify_rectangles = 0;

% hObject    handle to uitoggletool4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function uitoggletool8_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global draw_average
global modify_rectangles

modify_rectangles = 0;
draw_average = 0;



% --------------------------------------------------------------------
function uitoggletool8_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global modify_rectangles;
global draw_average


set(handles.uitoggletool1,'State','off')
        set(handles.uitoggletool2,'State','off')
        set(handles.uitoggletool3,'State','off')
        set(handles.uitoggletool4,'State','off')
        set(handles.uitoggletool7,'State','off')
        set(handles.uitoggletool12,'State','off')
        %handles.button_h = hObject;        
        
        draw_average = 1;
modify_rectangles  = 1;

pan off
zoom off
datacursormode off
modify_rectangles = 1;

guidata(hObject,handles)




% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles)
set(handles.uitoggletool7,'State','off')
set(handles.uitoggletool8,'State','off')
set(handles.uitoggletool12,'State','off')
zoom xon
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool1_OffCallback(hObject, eventdata, handles)
set(handles.uitoggletool7,'State','off')
set(handles.uitoggletool12,'State','off')
zoom off
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function open_optics_Callback(hObject, eventdata, handles)
global modify_rectangles;
global draw_average;
global plot_position;

    [fname,pname] = uigetfile('*.opt','Enter data file');
if (fname ~= 0)
filename = sprintf('%s%s',pname,fname);



for i = 1:length(handles.cluster)
    delete(handles.h(i))
end
%clear handles.cluster

handles.cluster = [];
handles.i = 0;
handles.active = [];
handles.nr_areas = 0;
handles.close = 0;
handles.h = [];             %contains the handles to the rectangles
handles.figure2 = 0;        %figure for the animation
handles.current_spike = 0;
handles.current_frame = 1;
modify_rectangles = 0;
draw_average = 0;

fid = fopen(filename,'r');

line = fgetl(fid);
mark = find(line == '"');
line = line(1:mark(1)-1);
s_line = str2num(line);
dim = length(s_line) - 2;           %minus 2, because OPTICS adds to lines

fseek(fid,0,'bof');


str = [];
for i = 1:dim / 5                   %the nummber of dimensions should be divisible by 5!!
str = [str, '%*f %*f %*f %*f %f ' ];
end
str = [str, '%f %f \"%d\"'];


handles.frames = dim / 30 / 5;      %gets the number of frames drawn

axes(handles.axes1)
cla
axis([0 1 0 1])
axis xy
%axis off
text(0.05,0.92,['loading    ',filename, '  ...'],'FontSize',14,'Interpreter','none');
drawnow
handles.spikedata = fscanf(fid,str,[dim /5 + 3,inf]);  
fclose(fid);


set(handles.figure1, 'Name',['SPIKE Tool', '    ',filename])

handles.active = 0;
handles.figure2 = 0;


colormap(jet(250));
axes(handles.axes1)
image([1,size(handles.spikedata,2)],[1,dim + 3],20 - handles.spikedata(:,1:size(handles.spikedata,2)));
zoom xon

guidata(hObject,handles)

end



guidata(hObject,handles)
% hObject    handle to open_optics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function open_mat_Callback(hObject, eventdata, handles)
% hObject    handle to open_mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [fname,pname] = uigetfile('*.mat','Enter data file');    
    %load_data(hObject,filename,handles);    
    if (fname ~= 0)
        filename = sprintf('%s%s',pname,fname);
        %filename = sprintf('%s%s',pname,fname);
        load(filename);

        for i = 1:length(handles.cluster)
            delete(handles.h(i))
        end
        %clear handles.cluster

        handles.cluster = [];
        handles.i = 0;
        handles.active = [];
        handles.nr_areas = 0;
        handles.close = 0;
        handles.h = [];             %contains the handles to the rectangles
        handles.figure2 = 0;        %figure for the animation
        handles.current_spike = 0;
        handles.current_frame = 1;
        modify_rectangles = 0;
        draw_average = 0;

        handles.spikedata = [S(ind,:),ind]';

        handles.frames = size(S,2)/30;

        set(handles.figure1, 'Name',['SPIKE Tool', '    ',filename])

        handles.active = 0;
        handles.figure2 = 0;


        colormap(jet(250));
        axes(handles.axes1)
        image([1,size(handles.spikedata,2)],[1,size(S,2) + 1],20 - handles.spikedata(:,1:size(handles.spikedata,2)));
        zoom xon
    end

    guidata(hObject,handles)

% --------------------------------------------------------------------
function open_clusters_Callback(hObject, eventdata, handles)
[fname,pname] = uigetfile('*.clu','Enter data file');
if (fname ~= 0)
    filename = sprintf('%s%s',pname,fname);
    load(filename,'-mat')

    for i = 1:length(handles.cluster)
        delete(handles.h(i))
    end

    handles.cluster = [];
    handles.cluster = clust; 
    guidata(hObject,handles)

    for i = 1:length(clust)
        handles.i = length(clust);
        handles.nr_areas = length(clust);

        handles.h(i) = rectangle('Position',[handles.cluster(i).pos(1),0,handles.cluster(i).pos(2)-handles.cluster(i).pos(1),1053], ...
            'EdgeColor','none','EraseMode','xor','FaceColor',[0.8 0.8 0.8]);

        set(handles.h(i),'ButtonDownFcn',{@button_down, handles.figure1})
    end
    handles.active = 0;
    handles.figure2 = 0;
    guidata(hObject,handles)
end
% hObject    handle to open_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_clusters_Callback(hObject, eventdata, handles)
[fname,pname] = uiputfile('*.clu','Save cluster file');
if (fname ~= 0)
filename = sprintf('%s%s',pname,fname);
clust = handles.cluster; 
save(filename,'clust')
end
% hObject    handle to save_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool12_OnCallback(hObject, eventdata, handles)
global modify_rectangles;

set(gcf,'Pointer','cross');

set(handles.uitoggletool1,'State','off')
        set(handles.uitoggletool2,'State','off')
        set(handles.uitoggletool3,'State','off')
        set(handles.uitoggletool4,'State','off')
        set(handles.uitoggletool7,'State','off')
        set(handles.uitoggletool8,'State','off')
        %handles.button_h = hObject;        
        

pan off
zoom off
datacursormode off
modify_rectangles = 2;
% hObject    handle to uitoggletool12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool12_OffCallback(hObject, eventdata, handles)
global modify_rectangles;
modify_rectangles = 0;
set(handles.uitoggletool12,'State','off')
set(gcf,'Pointer','arrow');
% hObject    handle to uitoggletool12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





