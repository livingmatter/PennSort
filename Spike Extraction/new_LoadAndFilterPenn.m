%function [y,time,tmp,Channels, clock,myheader] = LoadAndFilterPenn(FileName, type)
% y=LoadAndFilter(FileName,Channels);
function [y,time] = new_LoadAndFilterPenn(FileName,filelen,type)

%if nargin < 2, FileLen = Inf; end;
if nargin < 3
    type = 'Penn';
end
%FileLen = Inf;

switch type
    case 'Penn'
        n = regexp(FileName,'_[0-9]*');
        n = n(end);
        seq_index = str2num(FileName(n+1:end));
        
        fid=fopen(FileName,'r','ieee-be');
        myheader=new_Penn_read_lv_header(fid);
        display(['Header read: Samples/iter = ', num2str(myheader.samplesperiteration), ' #Channels = ', num2str(myheader.numdatach)]);
        tmp=fread(fid,[myheader.samplesperiteration, myheader.numdatach],'int16');
        while ~feof(fid)
          tmp=[tmp; fread(fid,[myheader.samplesperiteration, myheader.numdatach],'int16')];
        end

        fclose(fid);

        [samples, numch] = size(tmp);

        Channels = eval(myheader.channel_list);
        Channels = Channels+1;
        Channels = intersect(Channels,1:60);
        tmp = tmp';
        
        time = (1:size(tmp,2)) + (n-1)*filelen;
    case 'Ronen'
        load(FileName);
        tmp = double(data(1:63,:));
        Channels = 1:63;
        numch = 64;
    otherwise
        display('Not a valid type.');
        tmp = 0;
end

%trigger = tmp(:,numch);
% tmp = tmp(:,Channels);

%timeline=0:samples-1;
%timeline=timeline./myheader.scanrate;


%  timeline is the time axis
%  y contains the recorded data


disp('***  Datafile read ***'); 

%trigger = tmp(numch,51:end);
%clock = tmp(numch-1,:);

tmp=tmp(Channels,:);

% for i=1:60
%     y(i,:)=HighPassFilter(tmp(i,:));
% end

% The filter is peaked at 51, so y ends up having an offset of 50 relative
% to the original data.
y=DataLowPassFilter(tmp);
y = y(:,51:end);
time = time(1:end-50);
disp('***  channels highpass-filtered   ***'); 

