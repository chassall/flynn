function DISC = FLYNN( pathToConfigFile, pathToLocsFile )
%FLYNN 3.2.0 Takes a config file pathname and a locations file pathname, then loads, organizes, and
%analyzes EEG data.
%
% C. Hassall and O. Krigolson
% December, 2017
%
% FLYNN 3.0 .mat input (EEGLAB format), multiple .mat output
% FLYNN 2.0 trial eeg text file input, multiple .mat output
% FLYNN 1.0 average eeg text file input, single .mat output
% Requires: disc.wav, flynn.jpg, stats toolbox

% FLYNN version number (major, minor, revision)
version = '3.2.0';

% Load config file
if nargin == 0
    configFileId = fopen('FLYNNConfiguration.txt');
    userLocsFile = [];
elseif nargin == 1
    configFileId = fopen(pathToConfigFile);
    userLocsFile = [];
elseif nargin == 2
    configFileId = fopen(pathToConfigFile);
    userLocsFile = readlocs(pathToLocsFile);
end
C = textscan(configFileId, '%q','CommentStyle','%');
fclose(configFileId);
answer = C{1};

% Parse ANSWER
basefilename = answer{1};
subjectnumbers = strsplit(answer{2},',');
numberofsubjects = length(subjectnumbers);
baselinesettings = str2num(answer{3});
artifactsettings = str2num(answer{4});
outfile = answer{5};

% Determine which ERP/FFT/WVLT to do
numAnalyses = length(answer)-5;
if numAnalyses == 0
    disp('Error: No analysis specified');
    return;
end

% ERP variables
ERP.markers = {};
ERP.startTime = [];
ERP.endTime = [];
ERP.conditions = {};
numErpConditions = 0;
numErpMarkersByCondition = [];

% ALL variables
ALL.markers = {};
ALL.startTime = [];
ALL.endTime = [];
ALL.conditions = {};
numAllConditions = 0;
numAllMarkersByCondition = [];
ALL.whichMarker = {};
ALL.isArtifact = {};

% FFT variables
FFT.markers = {};
FFT.startTime = [];
FFT.endTime = [];
FFT.conditions = {};
numFftConditions = 0;
numFftMarkersByCondition = [];

% WAV variables
WAV.markers = {};
WAV.startTime = [];
WAV.endTime = [];
WAV.baselineStart = [];
WAV.baselineEnd = [];
WAV.frequencyStart = [];
WAV.frequencyEnd = [];
WAV.frequencySteps = [];
WAV.rangeCycles = [];
WAV.conditions = {};
numWavConditions = 0;
numWavMarkersByCondition = [];

for i = 1:length(answer)-5
    thisAnalysis = answer{5+i};
    temp = strsplit(thisAnalysis,',');
    if strcmp(temp{1},'ERP')
        numErpConditions = numErpConditions + 1;
        numMarkers = length(temp) - 4;
        numErpMarkersByCondition(numErpConditions) = numMarkers;
        for k = 1:numMarkers
            ERP.markers{k,numErpConditions} = temp{1+k};
        end
        ERP.startTime{numErpConditions} = temp{2 + numMarkers};
        ERP.endTime{numErpConditions} = temp{3 + numMarkers};
        ERP.conditions{numErpConditions} = temp{4+numMarkers};
    elseif strcmp(temp{1},'ALL')
        numAllConditions = numAllConditions + 1;
        numMarkers = length(temp) - 4;
        numAllMarkersByCondition(numAllConditions) = numMarkers;
        for k = 1:numMarkers
            ALL.markers{k,numAllConditions} = temp{1+k};
        end
        ALL.startTime{numAllConditions} = temp{2 + numMarkers};
        ALL.endTime{numAllConditions} = temp{3 + numMarkers};
        ALL.conditions{numAllConditions} = temp{4+numMarkers};
    elseif strcmp(temp{1},'FFT')
        
        numFftConditions = numFftConditions + 1;
        numFftMarkers = length(temp) - 4;
        numFftMarkersByCondition(numFftConditions) = numFftMarkers;
        for k = 1:numFftMarkers
            FFT.markers{k,numFftConditions} = temp{1+k};
        end
        FFT.startTime{numFftConditions} = temp{2 + numFftMarkers};
        FFT.endTime{numFftConditions} = temp{3 + numFftMarkers};
        FFT.conditions{numFftConditions} = temp{4+numFftMarkers};
        
    elseif strcmp(temp{1},'WAV')
        
        numWavConditions = numWavConditions + 1;
        numWavMarkers = length(temp) - 10;
        numWavMarkersByCondition(numWavConditions) = numWavMarkers;
        for k = 1:numWavMarkers
            WAV.markers{k,numWavConditions} = temp{1+k};
        end
        WAV.startTime{numWavConditions} = temp{2 + numWavMarkers};
        WAV.endTime{numWavConditions} = temp{3 + numWavMarkers};
        WAV.baselineStart{numWavConditions} = temp{4+numWavMarkers};
        WAV.baselineEnd{numWavConditions} = temp{5+numWavMarkers};
        WAV.frequencyStart{numWavConditions} = temp{6+numWavMarkers};
        WAV.frequencyEnd{numWavConditions} = temp{7+numWavMarkers};
        WAV.frequencySteps{numWavConditions} = temp{8+numWavMarkers};
        WAV.rangeCycles{numWavConditions} = temp{9+numWavMarkers};
        WAV.conditions{numWavConditions} = temp{10+numWavMarkers};
    else
        disp('Error: Unknown analysis');
        return;
    end
end

% DISC will hold participant summaries
DISC.EEGSum = []; % EEG Summary (participant, channels, datapoints)
DISC.ALLSum = []; % ALL Summary (participant, channels, datapoints)
DISC.ERPSum = []; % ERP Summary (participant, kept epochs, removed epochs)
DISC.FFTSum = []; % FFT Summary (participant, kept epochs, removed epochs)
DISC.WAVSum = []; % WAV Summary (participant, kept epochs, removed epochs)

firstLocsFile = [];

% Do analysis for each participant (ERP, FFT, WVLT)
for p = 1:numberofsubjects
    if isempty(subjectnumbers{p})
        disp('Error: No participants present');
        return;
    end
    % Data Import
    disp(['Current Subject Being Loaded: ' subjectnumbers{p}]);
    filename = [basefilename subjectnumbers{p} '.mat'];
    load(filename);
    
    %% Attempt to sort data if there is a user-defined locs file
    if ~isempty(userLocsFile)
        newOrder = nan(1,length(EEG.chanlocs)); % New channel order
        % Compare user channels to actual channels - if there is a match,
        % record in which position it was found
        for i = 1:length(userLocsFile)
            for k = 1:length(EEG.chanlocs)
                if strcmp(userLocsFile(i).labels,EEG.chanlocs(k).labels)
                   newOrder(i) = k; 
                end
            end
        end
        % Error checking
        if length(userLocsFile) ~= length(EEG.chanlocs) || any(isnan(newOrder))
            disp('Error: Locs file mismatch');
            return;
        else
           EEG.chanlocs = EEG.chanlocs(newOrder); % Reorder locs TODO: just use user-specified locs?
           EEG.data = EEG.data(newOrder,:,:); % Reorder data
        end
    end
    
    chanlocs = EEG.chanlocs;
    srate = EEG.srate;
    times = EEG.xmin*1000:1000/EEG.srate:EEG.xmax*1000;
    thisParticipantNumber = str2num(cell2mat(regexp(subjectnumbers{p},'\d','match'))); % Remove non-digits first
    if p == 1
       firstLocsFile = chanlocs; 
       DISC.EEGSum = [DISC.EEGSum; thisParticipantNumber EEG.nbchan EEG.pnts 1]; % First one is OK
    elseif isequal(firstLocsFile,chanlocs)
       DISC.EEGSum = [DISC.EEGSum; thisParticipantNumber EEG.nbchan EEG.pnts 1]; % Channel locs match
    else
       DISC.EEGSum = [DISC.EEGSum; thisParticipantNumber EEG.nbchan EEG.pnts 0]; % Channel locs don't match
    end
    
    %% Baseline Correction (if specified)
    if ~isempty(baselinesettings)
        baselinePoints = dsearchn(times',baselinesettings(:)); % Find the baseline indices
        baseline = mean(EEG.data(:,baselinePoints(1):baselinePoints(2) ,:),2);
        EEGb = EEG.data - repmat(baseline,[1,EEG.pnts,1]); % EEG data, with baseline correction applied
    else
        EEGb = EEG.data;
    end
    
    %% Epoching
    allMarkers = {EEG.epoch.eventtype}; % Markers within each epoch
    
    % Problem: epochs contain multiple markers - to know which one is at 0 ms, we need to check latencies
    
    latencies = {EEG.epoch.eventlatency}; % Latencies of all events within each epoch
    actualMarkers = {}; % Marker of interest for each epoch
    for m = 1:length(allMarkers)
        thisSetOfMarkers = allMarkers{m};
        theseLatencies = cell2mat(latencies{m});
        [~, whichOne] = min(abs(theseLatencies - abs(EEG.xmin)*1000000)); % Find the latency (in nanoseconds?) closest to 0 ms
        if isempty(whichOne)
            disp('Error: Timing error in EEGLAB file');
            return;
        end
        actualMarkers{m} = thisSetOfMarkers{whichOne};
    end
    
    %% ERP Analysis
    for c = 1:length(ERP.conditions)
        
        isThisCondition = false(1,length(actualMarkers));
        % Make a logical vector so that all relevant markers are inccluded
        for m = 1:numErpMarkersByCondition(c)
            isThisCondition = isThisCondition | strcmp(actualMarkers,ERP.markers{m,c});
        end
        
        if sum(isThisCondition) == 0
            disp('No epochs found');
            return;
        end
        
        ERP.timepoints{c} = str2num(ERP.startTime{c}):1000/EEG.srate:str2num(ERP.endTime{c});
        ERP.data{c} = nan(EEG.nbchan,length(ERP.timepoints{c}));
        erpPoints = dsearchn(times', [str2num(ERP.startTime{c}) str2num(ERP.endTime{c})]');
        erpEEG = EEGb(:,erpPoints(1):erpPoints(2),:);
        
        % ERP Artifact Rejection TODO: Make this a function
        % Artifact Rejection - Gradient
        maxAllowedStep = artifactsettings(1)*(1000/EEG.srate); % E.g. 10 uV/ms ~= 40 uV/4 ms... Equivalent to Analyzer?
        gradient = erpEEG(:,2:end,:) - erpEEG(:,1:end-1,:);
        gradientViolation = squeeze(any(gradient > maxAllowedStep,2));
        
        % Artifact Rejection - Difference
        maxAllowedDifference = artifactsettings(2);
        diffEEG = max(erpEEG,[],2) - min(erpEEG,[],2);
        differenceViolations = squeeze(diffEEG > maxAllowedDifference);
        
        allViolations = sum(gradientViolation) + sum(differenceViolations);
        isArtifact = allViolations ~= 0;
        
        ERP.nAccepted{c} = sum(~isArtifact & isThisCondition);
        ERP.nRejected{c} = sum(isArtifact & isThisCondition);
        
        thisAverage = mean(erpEEG(:,:,~isArtifact & isThisCondition),3);
        %         plot(thisAverage(34,:));
        %         hold on;
        ERP.data{c} = thisAverage;
        
        DISC.ERPSum = [DISC.ERPSum; thisParticipantNumber c ERP.nAccepted{c} ERP.nRejected{c}];
    end

    %% ALL Analysis (will store all trials of a certain type)
    for c = 1:length(ALL.conditions)
        
        isThisCondition = false(numAllMarkersByCondition(c),length(actualMarkers));
        % Make a logical vector so that all relevant markers are inccluded
        for m = 1:numAllMarkersByCondition(c)
            isThisCondition(m,:) = strcmp(actualMarkers,ALL.markers{m,c});
        end
        isAnyCondition = sum(isThisCondition) ~= 0;
        
        if sum(isAnyCondition) == 0
            disp('Error: No epochs found');
            return;
        end
        
        ALL.timepoints{c} = str2num(ALL.startTime{c}):1000/EEG.srate:str2num(ALL.endTime{c});
        ALL.data{c} = nan(EEG.nbchan,length(ALL.timepoints{c}));
        allPoints = dsearchn(times', [str2num(ALL.startTime{c}) str2num(ALL.endTime{c})]');
        allEEG = EEGb(:,allPoints(1):allPoints(2),:);
        
        % ERP Artifact Rejection TODO: Make this a function
        % Artifact Rejection - Gradient
        maxAllowedStep = artifactsettings(1)*(1000/EEG.srate); % E.g. 10 uV/ms ~= 40 uV/4 ms... Equivalent to Analyzer?
        gradient = abs(allEEG(:,2:end,:) - allEEG(:,1:end-1,:));
        gradientViolation = squeeze(any(gradient > maxAllowedStep,2));
        
        % Artifact Rejection - Difference
        maxAllowedDifference = artifactsettings(2);
        diffEEG = max(allEEG,[],2) - min(allEEG,[],2);
        differenceViolations = squeeze(diffEEG > maxAllowedDifference);
        
        allViolations = sum(gradientViolation) + sum(differenceViolations);
        isArtifact = allViolations ~= 0;
        
        ALL.nAccepted{c} = sum(~isArtifact & isAnyCondition);
        ALL.nRejected{c} = sum(isArtifact & isAnyCondition);
        % ALL.data{c} = erpEEG(:,:,~isArtifact & isThisCondition);
        ALL.data{c} = allEEG(:,:,isAnyCondition);
        ALL.whichMarker{c} = isThisCondition(:,isAnyCondition); % Marker for each trial
        ALL.isArtifact{c} = isArtifact;
        
        DISC.ALLSum = [DISC.ALLSum; thisParticipantNumber c ALL.nAccepted{c} ALL.nRejected{c}];
    end
  
    %% FFT Analysis
    for c = 1:length(FFT.conditions)
        FFT.timepoints{c} = str2num(FFT.startTime{c}):1000/EEG.srate:str2num(FFT.endTime{c});
        FFT.frequencyResolution{c} = EEG.srate / length(FFT.timepoints{c});
        
        fftPoints = dsearchn(times', [str2num(FFT.startTime{c}) str2num(FFT.endTime{c})]');
        fftEEG = EEGb(:,fftPoints(1):fftPoints(2),:);
        
        % ERP Artifact Rejection
        % Artifact Rejection - Gradient
        maxAllowedStep = artifactsettings(1)*(1000/EEG.srate); % E.g. 10 uV/ms ~= 40 uV/4 ms... Equivalent to Analyzer?
        gradient = abs(fftEEG(:,2:end,:) - fftEEG(:,1:end-1,:));
        gradientViolation = squeeze(any(gradient > maxAllowedStep,2));
        
        % Artifact Rejection - Difference
        maxAllowedDifference = artifactsettings(2);
        diffEEG = max(fftEEG,[],2) - min(fftEEG,[],2);
        differenceViolations = squeeze(diffEEG > maxAllowedDifference);
        
        allViolations = sum(gradientViolation) + sum(differenceViolations);
        isArtifact = allViolations ~= 0;
        
        % Contruct a boolean indicating if an epoch should be included
        isThisCondition = false(1,length(actualMarkers));
        % Make a logical vector so that all relevant markers are inccluded
        for m = 1:numFftMarkersByCondition(c)
            isThisCondition = isThisCondition | strcmp(actualMarkers,FFT.markers{m,c});
        end
        
        % Return if no epochs found
        if sum(isThisCondition) == 0
            disp('Error: No epochs found');
            return;
        end
        
        % Store the number of good epochs for this condition and the
        % proportion rejected
        FFT.nAccepted{c} = sum(~isArtifact & isThisCondition);
        FFT.nRejected{c} = sum(isArtifact & isThisCondition);
        
        DISC.FFTSum = [DISC.FFTSum; thisParticipantNumber c FFT.nAccepted{c} FFT.nRejected{c}];
        
        % Prepare the EEG on which the FFT will be run
        trimmedEEG.data = EEG.data(:,fftPoints(1):fftPoints(2),~isArtifact & isThisCondition);
        trimmedEEG.pnts = length(fftPoints(1):fftPoints(2));
        trimmedEEG.srate = EEG.srate;
        
        % Call doFFT
        [FFT.data{c},FFT.frequencies{c}] = doFFT(trimmedEEG);
    end
    
    %% Wavelet Analysis (TODO)
    for c = 1:length(WAV.conditions)
        
        WAV.timepoints{c} = str2num(WAV.startTime{c}):1000/EEG.srate:str2num(WAV.endTime{c});
        WAV.frequencyResolution{c} = EEG.srate / length(WAV.timepoints{c});
        
        wavPoints = dsearchn(times', [str2num(WAV.startTime{c}) str2num(WAV.endTime{c})]');
        wavEEG = EEGb(:,wavPoints(1):wavPoints(2),:);
        
        % ERP Artifact Rejection
        % Artifact Rejection - Gradient
        maxAllowedStep = artifactsettings(1)*(1000/EEG.srate); % E.g. 10 uV/ms ~= 40 uV/4 ms... Equivalent to Analyzer?
        gradient = abs(wavEEG(:,2:end,:) - wavEEG(:,1:end-1,:));
        gradientViolation = squeeze(any(gradient > maxAllowedStep,2));
        
        % Artifact Rejection - Difference
        maxAllowedDifference = artifactsettings(2);
        diffEEG = max(wavEEG,[],2) - min(wavEEG,[],2);
        differenceViolations = squeeze(diffEEG > maxAllowedDifference);
        
        allViolations = sum(gradientViolation) + sum(differenceViolations);
        isArtifact = allViolations ~= 0;
        
        % Contruct a boolean indicating if an epoch should be included
        isThisCondition = false(1,length(actualMarkers));
        % Make a logical vector so that all relevant markers are inccluded
        for m = 1:numWavMarkersByCondition(c)
            isThisCondition = isThisCondition | strcmp(actualMarkers,WAV.markers{m,c});
        end
        
        % Return if no epochs found
        if sum(isThisCondition) == 0
            disp('Error: No epochs found');
            return;
        end
        
        % Store the number of good epochs for this condition and the
        % proportion rejected
        WAV.nAccepted{c} = sum(~isArtifact & isThisCondition);
        WAV.nRejected{c} = sum(isArtifact & isThisCondition);
        
        DISC.WAVSum = [DISC.WAVSum; thisParticipantNumber c WAV.nAccepted{c} WAV.nRejected{c}];
        
        trimmedEEG.data = EEG.data(:,wavPoints(1):wavPoints(2),~isArtifact & isThisCondition);
        [~,~,trimmedEEG.trials] = size(trimmedEEG.data);
        trimmedEEG.times = WAV.timepoints{c};
        trimmedEEG.srate = EEG.srate;
        trimmedEEG.pnts =  length(wavPoints(1):wavPoints(2));
        
        baseline_windows = [str2num(WAV.baselineStart{c}) str2num(WAV.baselineEnd{c})];
        min_freq = str2num(WAV.frequencyStart{c});
        max_freq = str2num(WAV.frequencyEnd{c});
        num_frex = str2num(WAV.frequencySteps{c});
        range_cycles = str2num(WAV.rangeCycles{c});
        [WAV.data{c},WAV.dataPercent{c},WAV.frequencies{c}] = doWavelet(trimmedEEG,baseline_windows,min_freq,max_freq,num_frex,range_cycles);
        
    end
    %% Data Export
    outfilename = [outfile subjectnumbers{p} '.mat'];
    save(outfilename,'version','srate','chanlocs','ERP','ALL','FFT','WAV');
end

%% Visualization
scrsz = get(groot,'ScreenSize');
figure_loc = [1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/2]; % This is where figures will be drawn
%
% Flynn quotations
fquotes = {'On the other side of the screen, it all looks so easy.','How are you going to run the universe if you can''t answer a few unsolvable problems, huh?','Come on, you scuzzy data, be in there. Come on.','It''s time I leveled with you. I''m what you guys call a User.','I shouldn''t have written all of those tank programs.','It''s all in the wrists.','Greetings, programs!','Now for some real User power.','Did we make it? Hooray for our side.','Like the man says, there''s no problems, only solutions.','No sweat. I play video games better than anybody.','Come on, I''m - I''m scared of the dark. All this technology scares me.','Damn recognizer. Just go straight! I gotta get to that I/O tower.'};
this_fq = fquotes{randi(length(fquotes))};

fig = figure('Name',['FLYNN ' version],'NumberTitle','off','OuterPosition',figure_loc,'Toolbar','none','MenuBar','none');
whitebg(fig,'k');
set(gcf,'color','black');
flynn = imread('flynn.jpg');
subplot(2,4,[1 2 5 6]);
imshow(flynn,'InitialMagnification',33);
title(['\fontname{Courier}FLYNN ' version]);
if all(DISC.EEGSum(:,4))
    xlabel(['\fontname{Courier}' sprintf(['"' this_fq '" -Flynn'])]);
    disp('CONSISTENT CHANNELS');
else
    whichOnes = find(DISC.EEGSum(:,4) == 0);
    xlabel(['\fontname{Courier}USER ERROR!!! ' strjoin(subjectnumbers(whichOnes),', ') ' inconsistent with ' subjectnumbers{1}]);
    disp('USER ERROR!!! INCONSISTENT CHANNELS');
    disp([strjoin(subjectnumbers(whichOnes),', ') ' inconsistent with ' subjectnumbers{1}]);
end

% EEG Summary
subplot(2,4,3);
bar(DISC.EEGSum(:,2));
% lgd = legend({},'Location', 'northoutside','Orientation','horizontal');
title('Number of Channels');
xticklabels(num2str(DISC.EEGSum(:,1)));
xlabel('Participant');
ylabel('Number of Channels');

% ERP Summary
if ~isempty(DISC.ERPSum) || ~isempty(DISC.ALLSum)
    subplot(2,4,4);
    if ~isempty(DISC.ALLSum)
        bar(DISC.ALLSum(:,3:4),'stacked');
        lgd1 = legend('Accepted','Rejected','Location', 'northoutside','Orientation','horizontal');
        xticklabels(num2str(DISC.ALLSum(:,1:2)));
        title(lgd1,'ALL Artifacts');
    else
        bar(DISC.ERPSum(:,3:4),'stacked');   
        lgd1 = legend('Accepted','Rejected','Location', 'northoutside','Orientation','horizontal');
        xticklabels(num2str(DISC.ERPSum(:,1:2)));
        title(lgd1,'ERP Artifacts');
    end
    xlabel('Participant, Condition');
    ylabel('Number of Epochs');
    %lgd = legend(strrep(conditionnames,'_',''),'Location', 'north','Orientation','horizontal');
    %ylim([0 max(max(DISC(:,:,3))) + 20]);
end

% FFT Summary
if ~isempty(DISC.FFTSum)
    subplot(2,4,7);
    bar(DISC.FFTSum(:,3:4),'stacked');
    lgd2 = legend('Accepted','Rejected','Location', 'northoutside','Orientation','horizontal');
    xticklabels(num2str(DISC.FFTSum(:,1:2)));
    xlabel('Participant, Condition');
    ylabel('Number of Epochs');
    title(lgd2,'FFT Artifacts');
end

% WAV Summary
if ~isempty(DISC.WAVSum)
    subplot(2,4,8);
    bar(DISC.WAVSum(:,3:4),'stacked');
    lgd2 = legend('Accepted','Rejected','Location', 'northoutside','Orientation','horizontal');
    xticklabels(num2str(DISC.WAVSum(:,1:2)));
    xlabel('Participant, Condition');
    ylabel('Number of Epochs');
    title(lgd2,'WAV Artifacts');
end

end

% Use the commented-out code below to display results (condition 1, channel 1)
% plot(FFT.frequencies{1},FFT.data{1}(1,:));
% contourf(WAV.timepoints{1}, WAV.frequencies{1}, squeeze(WAV.data{1}(1,:,:)),'linecolor','none');