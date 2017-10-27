close all; clear all;

% FFT Test Signal
load('cognitive_assessment_01.mat');
[nChannels,nTimepoints,nTrials] = size(EEG.data);

% Make a signal with 1 Hz, 6 Hz, and 20 Hz
timepoints = 0.004:0.004:nTimepoints*0.004;
testSignal = sin(2*pi*1*timepoints) + sin(2*pi*6*timepoints) + sin(2*pi*20*timepoints);
subplot(1,3,1);
plot(timepoints,testSignal);
ylim([-5 5]);

% Replace all trials with test signal
for t = 1:nTrials
    for c = 1:nChannels
        EEG.data(c,:,t) = testSignal;
    end
end

save('ffttestsignal_01.mat','EEG');

% Wavelet Test Signal (400 ms theta burst centered on 300 ms: 100 ms to 500 ms)
thetaTimes = 0.004:0.004:100*0.004;
thetaBurst = zeros(1,nTimepoints);
thetaBurst(251:350) = 2*sin(2*pi*6*thetaTimes);
subplot(1,3,2);
plot(timepoints,thetaBurst);
ylim([-5 5]);
testSignal2 = testSignal + thetaBurst;
subplot(1,3,3);
plot(timepoints,testSignal2);
ylim([-5 5]);

% Replace all trials with test signal
for t = 1:nTrials
    for c = 1:nChannels
        EEG.data(c,:,t) = testSignal2;
    end
end

save('wvlttestsignal_01.mat','EEG');