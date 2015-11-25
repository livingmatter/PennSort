function [stats] = updateClusterStats(stats, amplist, N)

% Hyperparameters:
% Conjugate to Gaussian amplitude distribution:
alpha = 5;
beta = alpha * 0.001;
nu = 5;
lambda = 1;

% Conjugate to Binomial:
a = 0.1; b = 0.1;

for i=1:size(stats,1)
    stats(i,5) = stats(i,5) + length(amplist{i});           % Number of spikes
    
    stats(i,1) = stats(i,1) + sum(amplist{i});              % Sum of amplitudes
    stats(i,2) = stats(i,2) + sum(amplist{i}.^2);           % Sum of squared amplitudes

    % Mean amplitude            
    stats(i,3) = (stats(i,1) + lambda*nu) / (stats(i,5)+nu);

    % Amplitude variance
    stats(i,4) = (2*beta + stats(i,2) - 2*stats(i,3)*stats(i,1) + stats(i,5)*stats(i,3)^2 + nu*(stats(i,3)-lambda)^2)/(stats(i,5)+2*alpha);
end

% p
stats(:,6) = (a + stats(:,5)) / (a + b + N);
