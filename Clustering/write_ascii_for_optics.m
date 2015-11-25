function write_ascii_for_optics(files,s_interp)
%
% write_ascii_for_optics('folder + filename',interpolated spikes)

str = [];
for i = 1:(size(s_interp,1) * size(s_interp,2))
str = [str, '%f '];
end
str = [str, '\n'];

fid = fopen([files],'wt');          %change that line here!!

fprintf(fid,str,s_interp);                    %that line takes a while
