function varargout = clustersplit(varargin)
% CLUSTERSPLIT M-file for clustersplit.fig
%      CLUSTERSPLIT, by itself, creates a new CLUSTERSPLIT or raises the existing
%      singleton*.
%
%      H = CLUSTERSPLIT returns the handle to a new CLUSTERSPLIT or the handle to
%      the existing singleton*.
%
%      CLUSTERSPLIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTERSPLIT.M with the given input arguments.
%
%      CLUSTERSPLIT('Property','Value',...) creates a new CLUSTERSPLIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before clustersplit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to clustersplit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help clustersplit

% Last Modified by GUIDE v2.5 14-Dec-2009 16:20:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @clustersplit_OpeningFcn, ...
                   'gui_OutputFcn',  @clustersplit_OutputFcn, ...
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


% --- Executes just before clustersplit is made visible.
function clustersplit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to clustersplit (see VARARGIN)

% Choose default command line output for clustersplit
handles.output = hObject;

[FileName,PathName] = uigetfile('mCluster*.mat');
%load([PathName,FileName]);
files = dir([PathName,'mCluster*.mat']);
clustnums = zeros(1,length(files));
for i=1:length(files)
    file = files(i).name;
    [startind,endind] = regexp(file,'[0123456789]+');
    clustnums(i) = str2num(file(startind:endind));
end
clustnums = sort(clustnums,'ascend');
clustnumstr = [];
for i=1:length(clustnums)
    clustnumstr = [clustnumstr, sprintf('%03d',clustnums(i)), '|'];
end    

set(handles.dirtext,'String',PathName);
set(handles.clusterlist,'String',clustnumstr);
princompstring = [];
for i=1:32
    princompstring = [princompstring, sprintf('%02d',i), '|'];
end
set(handles.princomplist,'String',princompstring);

handles.clusterpath = [PathName,'mCluster'];

handles.numclusters = max(clustnums);
handles.array = MEA(@geom_2F30EL);
handles.channelaxes = handles.array.makeAxes(handles.displaypanel);
for ch=1:length(handles.channelaxes)
    set(handles.channelaxes(ch),'ButtonDownFcn',{@channelaxes_Callback,handles});
    set(handles.channelaxes(ch),'XTick',[]);
    set(handles.channelaxes(ch),'HandleVisibility','callback');
end

handles.ClusterData = [];
handles.IndexOrigin = [];
handles.IndicesToDisplay = [];
handles.PCACoeff = [];
handles.PCAScore = [];
handles.CutLine = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes clustersplit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = clustersplit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.numclusters;


% --- Executes on selection change in clusterlist.
function clusterlist_Callback(hObject, eventdata, handles)
% hObject    handle to clusterlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns clusterlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clusterlist


% --- Executes during object creation, after setting all properties.
function clusterlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_ClusterData(handles)

if get(handles.allwaveformsradio,'Value')
    data = handles.ClusterData(:,:,handles.IndicesToDisplay);
else
    data = median(handles.ClusterData(:,:,handles.IndicesToDisplay),3);
end

for ch=1:length(handles.channelaxes)    
    plot(handles.channelaxes(ch),squeeze(data(:,ch,:)));
    set(handles.channelaxes(ch),'XTick',[],'XLim',[1,size(handles.ClusterData,1)],'YLim',[min(handles.ClusterData(:)),max(handles.ClusterData(:))]);
    set(handles.channelaxes(ch),'ButtonDownFcn',{@channelaxes_Callback,handles});
    set(handles.channelaxes(ch),'HandleVisibility','callback','DrawMode','fast');
end


% --- Executes on button press in updatebutton.
function updatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to updatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currentclust = get(handles.clusterlist,'Value');
load([handles.clusterpath, num2str(handles.currentclust)]);

handles.IndexOrigin = IndexOrigin;
handles.ClusterData = reshape(C',size(C,2)/length(handles.channelaxes),length(handles.channelaxes),size(C,1));
handles.IndicesToDisplay = 1:size(C,1);

template = median(handles.ClusterData,3);
[m,handles.pcchan] = min(min(template));
set(handles.pcchanneltext,'String',num2str(handles.pcchan));

set(handles.princompbutton,'Enable','on');
set(handles.pcaforward,'Enable','on');
set(handles.pcaback,'Enable','on');

guidata(hObject,handles);

plot_ClusterData(handles);

% --- Executes on button press in forwardbutton.
function forwardbutton_Callback(hObject, eventdata, handles)
% hObject    handle to forwardbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentclust = get(handles.clusterlist,'Value');
if currentclust < handles.numclusters
    set(handles.clusterlist,'Value',currentclust+1);
    drawnow;
    guidata(hObject,handles);
    updatebutton_Callback(hObject,eventdata,handles);
end

% --- Executes on button press in backbutton.
function backbutton_Callback(hObject, eventdata, handles)
% hObject    handle to backbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentclust = get(handles.clusterlist,'Value');
if currentclust > 1
    set(handles.clusterlist,'Value',currentclust-1);
    drawnow;
    guidata(hObject,handles);
    updatebutton_Callback(hObject,eventdata,handles);
end    

% --- Executes on selection change in princomplist.
function princomplist_Callback(hObject, eventdata, handles)
% hObject    handle to princomplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns princomplist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from princomplist


% --- Executes during object creation, after setting all properties.
function princomplist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to princomplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in princompbutton.
function princompbutton_Callback(hObject, eventdata, handles)
% hObject    handle to princompbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.PCACoeff,handles.PCAScore] = princomp(squeeze(handles.ClusterData(:,handles.pcchan,:))');
pctohist = get(handles.princomplist,'Value');
hist(handles.pcaxes,handles.PCAScore(:,pctohist),0.3*length(handles.PCAScore));
handles.CutLine = line('XData',[0,0],'YData',[0,0],'Color','r','Parent',handles.pcaxes);
guidata(hObject,handles);


% --- Executes on button press in pcaforward.
function pcaforward_Callback(hObject, eventdata, handles)
% hObject    handle to pcaforward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentpc = get(handles.princomplist,'Value');
if currentpc < 32
    set(handles.princomplist,'Value',currentpc+1);
    guidata(hObject,handles);
    princompbutton_Callback(hObject,eventdata,handles);
end

% --- Executes on button press in pcaback.
function pcaback_Callback(hObject, eventdata, handles)
% hObject    handle to pcaback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentpc = get(handles.princomplist,'Value');
if currentpc > 1
    set(handles.princomplist,'Value',currentpc-1);
    guidata(hObject,handles);
    princompbutton_Callback(hObject,eventdata,handles);
end


% --- Executes on button press in cutbutton.
function cutbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%axes(handles.pcaxes);
currentpc = get(handles.princomplist,'Value');
[x,y] = ginput(1);

set(handles.CutLine,'XData',[x,x],'YData',get(handles.pcaxes,'YLim'));
drawnow;
%line([x,x],get(handles.pcaxes,'YLim'),'Color','r'); drawnow;
handles.ixL = find(handles.PCAScore(:,currentpc) <= x);
handles.ixR = find(handles.PCAScore(:,currentpc) > x);

guidata(hObject,handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pcchanneltext_Callback(hObject, eventdata, handles)
% hObject    handle to pcchanneltext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcchanneltext as text
%        str2double(get(hObject,'String')) returns contents of pcchanneltext as a double


% --- Executes during object creation, after setting all properties.
function pcchanneltext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcchanneltext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function channelaxes_Callback(hObject, eventdata, handles)

pcchan = find(handles.channelaxes == gca);
set(handles.pcchanneltext,'String',num2str(pcchan));
handles.pcchan = pcchan;
guidata(hObject,handles);



% --- Executes on button press in showlbutton.
function showlbutton_Callback(hObject, eventdata, handles)
% hObject    handle to showlbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.IndicesToDisplay = handles.ixL;
guidata(hObject,handles);
plot_ClusterData(handles);


% --- Executes on button press in showrbutton.
function showrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to showrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.IndicesToDisplay = handles.ixR;
guidata(hObject,handles);
plot_ClusterData(handles);


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ix = handles.IndicesToDisplay;
C = reshape(handles.ClusterData(:,:,ix),size(handles.ClusterData,1)*size(handles.ClusterData,2),length(ix))';
IndexOrigin = handles.IndexOrigin(ix);
handles.numclusters = handles.numclusters + 1;
save([handles.clusterpath,num2str(handles.numclusters),'.mat'],'C','IndexOrigin');

allix = 1:size(handles.ClusterData,3);
ix = allix(~ismember(allix,ix));
C = reshape(handles.ClusterData(:,:,ix),size(handles.ClusterData,1)*size(handles.ClusterData,2),length(ix))';
IndexOrigin = handles.IndexOrigin(ix);
save([handles.clusterpath,num2str(handles.currentclust),'.mat'],'C','IndexOrigin');

clustnumstr = [];
for i=1:handles.numclusters
    clustnumstr = [clustnumstr, sprintf('%03d',i), '|'];
end    
set(handles.clusterlist,'String',clustnumstr);

guidata(hObject,handles);


% --- Executes on button press in allwaveformsradio.
function allwaveformsradio_Callback(hObject, eventdata, handles)
% hObject    handle to allwaveformsradio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allwaveformsradio

plot_ClusterData(handles);


% --- Executes on button press in medianradio.
function medianradio_Callback(hObject, eventdata, handles)
% hObject    handle to medianradio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of medianradio

plot_ClusterData(handles);
