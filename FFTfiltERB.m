function [output,cfFilt,ERBcenter,nERB] = FFTfiltERB(input,fs,range,chanWidthERB)
%%
% TAKE THE FOURIER TRANSFORM OF THE SIGNAL
% SET ALL FREQUENCIES OUTSIDE THE PASSBAND TO 0
% INVERSE FOURIER TRANSFORM THE RESULT
% CUT-OFF AND CENTER FREQUENCIES ARE ESTIMATED ON an ERB SCALE
% (GLASBERG & MOORE 1990) 
% 
%
%       [output,cfFilt,ERBcenter,nERB] = FFTfiltERB(input,fs,range,chanWidthERB)
%
%       INPUT
%   input:             signal to be filtered => m-by-1 array
%   fs:                sampling frequency => integer
%   range:             frequency range of the filterbank => 2-by-1 array
%   chanWidthERB:      width of the bands in ERB => integer
% 
%       OUTPUT
%   output:            output signal => m-by-1 array
%   cfFilt:            cutoffs of the bands => (n+1)-by-1 array
%   ERBcenter:         center frequencies of the bands => n-by-1 array
%   nERB:              # of channels => integer (upper frequency of the highest channel <= max(range))
%% ============================== VARIABLES ===============================
input = input-mean(input);
limLow = range(1);
limHigh = range(2);
N = size(input,1);
dF = fs/N;
f = (-fs/2:dF:fs/2-dF)';
% ====== Estimate # of ERBs within range ======
maxERB = 21.4*log10(0.00437*limHigh+1)-21.4*log10(0.00437*limLow+1);
nERB = floor(maxERB/chanWidthERB);
% ====== Set lower cut-off frequency ======
cfFilt = limLow;
% ====== Estimate cut-offs & center frequencies ======
for cf = 1:nERB
    ERBlow = 21.4*log10(0.00437*cfFilt(cf)+1);
    ERBcenter(cf) = (exp(log(10)*((ERBlow+0.5*chanWidthERB)/21.4))-1)/0.00437;
    cfFilt(cf+1) = (10^((ERBlow+chanWidthERB)/21.4)-1)/0.00437;
end
output = [];
%% ========================= APPLY FFT FILTERING ==========================
parfor band = 1:nERB
    % ====== Set cutoffs
	filtArray = ((cfFilt(band)  < abs(f)) & (abs(f) < cfFilt(band+1)));
	% ====== FFT/IFFT
    spektrum = fftshift(fft(input))/N;
    spektrum = filtArray.*spektrum;
    output(:,band) = real(ifft(ifftshift(spektrum)));
end