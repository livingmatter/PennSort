function [shiftadj,shifttemp,candidatelist,shifts] = preFit(fitCorLen,fitCorMag,TemplateMatrix,ndecimate)
    
    % Build up bank of stretched and shifted templates, both decorrelated
    % and raw.

    global onewaveformlength; global numberchanns;
    midwaveform = floor(onewaveformlength/2);   % Location of peak in downsampled templates
    relevantthresh=0.7;                % cutoff to determine which channels are relevant to a template
    
    correlationrange = 3;              % explore this range around the expected peak location
    deltat = 1/ndecimate;
    shifts = -correlationrange:deltat:correlationrange;
    widths = [0.9,1,1.1];
    %widths = [1];
    
    nclusters = size(TemplateMatrix,2);
    n = size(TemplateMatrix,1)/numberchanns;
    fulltemplatetable = reshape(TemplateMatrix,n,numberchanns, nclusters);
   
    % Generate all shifted, scaled, and downsampled templates.    
    shifttemp = zeros(onewaveformlength, numberchanns, nclusters, length(shifts),length(widths));
    for i=1:length(shifts)
        for j=1:length(widths)
            t = 1:deltat:onewaveformlength;            
            tprime = ((1:onewaveformlength)-midwaveform-shifts(i)) / widths(j) + midwaveform;
                 
            shifttemp(:,:,:,i,j) = interp1(t,fulltemplatetable,tprime,'spline',0);
        end
    end
        
    % Generate inverse correlation matrix:
    covX = exp(-1/fitCorLen); 
    %covX=0;
    covA = (1/fitCorMag)*(1+covX^2)/(1-covX^2); % diagonal entries of inverse covariance
    covB = -(1/fitCorMag)*covX/(1-covX^2);      % offdiagonal entries    
    Chati = spdiags([covB*ones(onewaveformlength,1),covA*ones(onewaveformlength,1),covB*ones(onewaveformlength,1)],[-1,0,1],onewaveformlength,onewaveformlength);
    
    % Compute the adjoints
    shifttemp = reshape(shifttemp, onewaveformlength, numberchanns*nclusters*length(shifts)*length(widths));
    shiftadj = Chati * shifttemp;
    shifttemp = reshape(shifttemp, onewaveformlength, numberchanns, nclusters, length(shifts), length(widths));
    shiftadj = reshape(shiftadj, onewaveformlength, numberchanns, nclusters, length(shifts), length(widths));
    
    % Make candidatelist: for each channel, list of templates which have
    % significant amplitude on that channel.
    candidatelist = cell(1,numberchanns);
    TemplateAmps = min(TemplateMatrix);
    TemplateMins = squeeze(min(reshape(TemplateMatrix,n,numberchanns,nclusters)));
    for i=1:numberchanns
        candidatelist{i} = find(TemplateMins(i,:) < relevantthresh*TemplateAmps);
    end
end
