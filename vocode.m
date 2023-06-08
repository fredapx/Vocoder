function vocodedSig = vocode(input,fs,range,chanWidth,envCf,carrierType)
%%
% ENVELOPE VOCODER Extract the temporal envelope in a series of frequency bands
% and use these envelopes to modulate a series of carriers in
% the corresponding frequency bands.
% Bands are regularly spaced on an ERB scale
%
%
%       vocodedSig = vocode(input,fs,range,chanWidth,envCf,carrierType)
%
%       INPUT
%   input:              the sound file to be vocoded => m-by-1 array
%   fs:                 sample frequency of the input => integer
%   range:              frequency range of the vocoder => 2-by-1 array
%   chanWidth:          the width of the frequency bands (in ERB) => integer
%   envCf:              envelope cut-off frequency => integer (1=1/2 ERB)
%   carrierType:        carrier type => 1=noise, 2=tone
%
%       OUTPUT
%   vocodedSig:         envelope vocoded signal => m-by-1 array

%%
    Nyq = fs/2;
    duree = length(input);
    [filtSig,chanCF,chanCenter,nChan] = FFTfiltERB(input,fs,range,chanWidth);
    
    
    if envCf == 1 % Estimate envelope cutoff as half the bandwidth of the auditory filter tuned to the center freq of the analysis band
        envCf = zeros(nChan,1);
        parfor ii = 1:nChan
            ERB = 24.7*(0.00437*chanCenter(ii)+1);
            envCf(ii) = round(ERB/2);
        end
    else
        [b,a] = butter(3, envCf/Nyq);
    end

    voc = zeros(duree,nChan);
    %%
    switch carrierType
        case 1 % % Gaussian noise
            noise = randn(duree,1);
            [carrier,TRASH,IGNORE] = FFTfiltERB(noise,fs,range,chanWidth);            
        case 2 % Tone
            carrier = voc;
            t = (1:duree)/fs;
            for ii = 1:nChan
                theta = rand(1,1)*2*pi;
                carrier(:,ii) = sin(2*pi*chanCenter(ii)*t+theta);
            end        
    end

    parfor chan = 1:nChan
        band = filtSig(:,chan);
        bandLevel = rms(band);
        if envCf == 1 % Estimate envelope cutoff as half the bandwidth of the auditory filter tuned to the center freq of the analysis band
            ERB = 24.7*(4.37*(chanCenter(chan)/1000)+1);
            cfe = round(ERB/2);
        else
            cfe = envCf;
        end
        [b,a] = butter(3, cfe/Nyq);
        env = abs(hilbert(band));
        env = filtfilt(b,a,env);
    %     TFS = cos(angle(hilbert(band)));
        bandSig = env.*carrier(:,chan);
        bandSig = myFFTfilt(bandSig,fs,[chanCF(chan) chanCF(chan+1)],"bp")
        voc(:,chan) = bandSig*bandLevel;
    end

    vocodedSig = sum(voc,2);

end