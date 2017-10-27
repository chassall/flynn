function [fftResults,frequencies] = doFFT(trimmedEEG)
%doFFT Does an FFT
%   C. Hassall
%   October, 2017

% run the fft
fftOutput = fft(trimmedEEG.data,[],2);

% determine frequency resolution
frequencyResolution = trimmedEEG.srate / trimmedEEG.pnts;

%We want to scale the output of the fft function to be divided by the length of the data submitted because The FFT computation involves a lot of summing over time points, and so dividing by N puts the data back to the scale of the original input.
fftOutput = fftOutput/trimmedEEG.pnts;

% we want to remove the negative half of the fft output because it is redundant
fftOutput = fftOutput(:,1:trimmedEEG.pnts/2,:);

% you take abs of the output of the fft because you want the resultant of the real and imaginary output
fftOutput = abs(fftOutput);

% you want to you multiply by 2 to ensure that the positive aspect of the fft output has to be corrected for the loss of the negative aspect
fftOutput = 2 * fftOutput;

% finally, remove the very first value as it is the power at 0 Hz (the DC power)
fftOutput(:,1,:) = [];

% Compute the average FFT output for this condition
fftResults = squeeze(mean(fftOutput,3));

% create a frequency bin
frequencies = frequencyResolution:frequencyResolution:length(fftOutput)*frequencyResolution;

end