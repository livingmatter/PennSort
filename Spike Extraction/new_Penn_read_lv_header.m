function header=Penn_read_lv_header(fid)

header.type=char(fread(fid, 1, 'char'));
header.size=fread(fid, 1, 'int32');
header.channel_list_length=fread(fid, 1, 'int32');
header.channel_list=[char(fread(fid, header.channel_list_length, 'char'))]';
header.channel_info_length=fread(fid, 1, 'int32');

header.numch = fread(fid,1,'int32');
for i=1:header.numch
    chconfig.channel = ReadLVString(fid);
    chconfig.uilim = fread(fid,1,'float32');
    chconfig.lilim = fread(fid,1,'float32');
    chconfig.range = fread(fid,1,'float32');
    chconfig.polarity = fread(fid,1,'uint16');
    chconfig.gain = fread(fid,1,'float32');
    chconfig.coupling = fread(fid,1,'uint16');
    chconfig.inputmode = fread(fid,1,'uint16');
    chconfig.scalemult = fread(fid,1,'float32');
    chconfig.scaleoff = fread(fid,1,'float32');
    header.channel(i)=chconfig;
end

header.scanrate = fread(fid,1,'float32');
header.channeldelay = fread(fid,1,'float32');

header.time='';
new_char=char(fread(fid, 1, 'char'));
while new_char~='	'
    header.time=[header.time new_char];
    new_char=char(fread(fid, 1, 'char'));
end

header.date='';
new_char=char(fread(fid, 1, 'char'));
while new_char~='	'
    header.date=[header.date new_char];
    new_char=char(fread(fid, 1, 'char'));
end

header.user_header='';
while ftell(fid)<(header.size+1+4)
% +1 for the file type character, +4 for the header size itself!!
header.user_header=[header.user_header char(fread(fid, 1, 'char'))];
% header.user_header=[header.user_header char(fread(fid, 300, 'char'))];
end

header.user_header;

% Read out samples per iteration and number of data channels (JP)
fseek(fid,header.size+1+4,'bof');
header.numdatach = fread(fid,1,'int32');
header.samplesperiteration = fread(fid,1,'int32');
position=ftell(fid);
%disp(['AT END OF HEADER, FILE POSITION IS: ' num2str(position) ]);


