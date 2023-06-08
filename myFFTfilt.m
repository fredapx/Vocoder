function output = myFFTfilt(input,fs,cutoff,fType)
%%
% TAKE THE FOURIER TRANSFORM OF THE SIGNAL
% SET ALL FREQUENCIES OUTSIDE THE PASSBAND TO 0
% INVERSE FOURIER TRANSFORM THE RESULT
%
%
%       output = myFFTfilt(input,fs,cutoff,fType)
%
%       INPUT
%   input:             signal to be filtered => n-by-1 array
%   fs:                sampling frequency => integer
%   cutoff:            frequency range of the filterbank => integer or 2-by-1 array
%   fType:             filter type, specified as one of the following:
%                      'bp' specifies a bandpass filter
%                      'lo' specifies a lowpass filter
%                      'hi' specifies a highpass filter
% 
%       OUTPUT
%   output:            filtered output signal => n-by-1 array
%% ============================== VARIABLES ===============================
N = size(input,1);
dF = fs/N;
f = (-fs/2:dF:fs/2-dF)';
%% ========================= APPLY FFT FILTERING ==========================
input = input-mean(input);
% ====== Set cutoff
switch fType
    case "bp"
        filtArray = ((cutoff(1)  < abs(f)) & (abs(f) < cutoff(2)));
    case "lo"
        filtArray = abs(f) < cutoff;        
    case "hi"
        filtArray = abs(f) > cutoff;
        cutoff
end
% ====== FFT/IFFT
spektrum = fftshift(fft(input))/N;
spektrum = filtArray.*spektrum;
output = real(ifft(ifftshift(spektrum)));
% figure;
% subplot(2,1,1);
% plot (f,abs(spektrum));
% subplot(2,1,2);
% plot (f,abs(spektrum));
end