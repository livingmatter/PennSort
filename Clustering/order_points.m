function [ind,clu] = order_points(S,chmask,peakch,datadir)

%allparameters;
global onewaveformlength; global numberchanns;

% epsilon = 40;
% nmin = 10;

optics_minpts = 10;
optics_eps = 1e6;

optics_folder = '/Users/jason/MEA/SpikeSort/Clustering/';
if ~exist([datadir 'OPTICSfiles'],'dir')
    mkdir([datadir 'OPTICSfiles']);
end
optics_file_prefix = [datadir 'OPTICSfiles/temp'];



% [chmask_sorted,initind] = sortrows(chmask);
% markers = find(any(diff(chmask_sorted,[],1)~=0,2));
% if isempty(markers)
%     markers = [length(initind)];
% end
% markers = [1;markers+1];

interpfactor = 5;
ti = 1:(1/interpfactor):onewaveformlength;
S_masked = zeros(size(S,1),numberchanns*(length(ti)+2*interpfactor));

for i=1:size(S,1)
    spike = reshape(S(i,:),onewaveformlength,numberchanns);


    % Interpolate & pad
    spike = interp1(1:onewaveformlength,spike,ti,'spline');
    spike = [zeros(interpfactor,size(spike,2));...
             spike;...
             zeros(interpfactor,size(spike,2))];

    % Align to new minimum
    wavemid = floor(onewaveformlength/2);
    [m,offset] = min(spike([(wavemid*interpfactor)+1:(wavemid+2)*interpfactor],peakch(i)));
    spike = circshift(spike,[interpfactor-offset+1,0]);

    % Zero out distant channels
    spike2 = zeros(size(spike));
    spike2(:,chmask(i,:)') = spike(:,chmask(i,:)');    
    S_masked(i,:) = reshape(spike2,1,size(spike2,1)*numberchanns);

%     spike_masked = zeros(size(spike));
%     spike_masked(:,chmask(i,:)') = spike(:,chmask(i,:)');
%     S_masked(i,:) = reshape(spike_masked,1,onewaveformlength*numberchanns);
end


%r = cell(1,length(markers)-1);
%clu = {};
%for m=1:length(markers)-1
ind = [];
for m=unique(peakch)
    display(['Sorting neighborhood of channel ',num2str(m),'...']);

    %spikes = initind(markers(m):markers(m+1)-1);
    spikes = find(peakch == m)';
%     r_euclid = zeros(length(spikes));
%     for i=1:length(spikes)
%         for j=1:(i-1)
%             r_euclid(i,j) = norm(S_masked(spikes(i),:)-S_masked(spikes(j),:));
% %             r_euclid(i,j) = norm(S_masked(spikes(i),:)-S_masked(spikes(j),:),inf);
%             r_euclid(j,i) = r_euclid(i,j);
%         end
%     end
    
    % *** Run OPTICS ***
%    r_euclid = r_euclid + tril(r_euclid,-1)';
    dim = size(S_masked,2);
    num_points = length(spikes);
    ascii_file_name = [optics_file_prefix num2str(m) '.ascii'];
    dump_file_name = [optics_file_prefix num2str(m) '.dump'];
    optics_file_name = [optics_file_prefix num2str(m) '.opt'];
    
    % Save spikes out in ascii format.
    write_ascii_for_optics(ascii_file_name,reshape(S_masked(spikes,:)',dim/numberchanns,numberchanns,num_points));      
    %convert .ascii to .dump file
    command = [[optics_folder,'ascii2dump'],' ', ascii_file_name,' ', num2str(dim),' ', dump_file_name];
    system(command);    
    %run optics on the .dump file -> .opt file as result
    command = [optics_folder,'optics',' ', '-s',' ', num2str(optics_eps,'%2.0e'),' ', num2str(optics_minpts), ' '...
        , num2str(dim),' ', num2str(num_points),' ',dump_file_name, ' ', optics_file_name];
    system(command);
    
    % Load .opt output to get ordering.
    str = [];
    for i = 1:dim
        str = [str, '%f ' ];
    end
    str = [str, '%f %f \"%d\"'];

    fid = fopen(optics_file_name,'r');
    OPTICS_output = fscanf(fid,str,[dim + 3,inf]);  
    ind_m = OPTICS_output(dim+3,:)+1;
    clear OPTICS_output;
    fclose(fid);

  % ind_m = poormansoptics(r_euclid);
%    clu_local = dbscan(r_euclid,epsilon,nmin);
    
%     ind_m = [];
%     for i=1:length(clu_local)
%         clu = [clu,{spikes(clu_local{i})}];
%         ind_m = [ind_m,clu_local{i}];
%     end
    %subplot(1,2,1);
    %imagesc(-S(spikes(ind_m),:)',[-50,500]); drawnow;
    %subplot(1,2,2); 
    %imagesc(r_euclid(ind_m,ind_m));drawnow; pause;
%    hist(r_euclid(:),500); drawnow; pause;
%    plot(S_masked(spikes,:)'); drawnow; pause;
    
    ind = [ind; spikes(ind_m)];
%    r{m} = r_euclid;
end
