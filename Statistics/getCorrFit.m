function [fitCorMag, fitCorLen,acorr] = getCorrFit(noise,epsilon)

% Find autocorrelation function of noise, fit to exponential, then return
% scale and correlation length.

% noise = n x 960 matrix of noise traces.

global onewaveformlength; global numberchanns;
%allparameters;
if nargin==1
    epsilon = 0.01;
end

fitCorMag = var(noise(:));

noise = reshape(noise',onewaveformlength,numberchanns*size(noise,1));
noise = [zeros(floor(onewaveformlength/2),size(noise,2)); noise; zeros(floor(onewaveformlength/2),size(noise,2))];
ftnoise = fft(noise);
pwr = mean(ftnoise.*conj(ftnoise),2);
acorr = ifft(pwr);

% Fit an exponential to initial part of autocorr. function greater than
% epsilon*acorr(1).
firstneg = find(acorr<epsilon*acorr(1),1,'first');
x = 1:(firstneg-1);
y = log(acorr(x));
fitcoeffs = [ones(size(x))', x'] \ y;

fitCorLen = -(1/fitcoeffs(2));
