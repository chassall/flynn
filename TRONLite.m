% Sample TRON script
% C. Hassall
% March, 2018

% Load participant data
filePrefix = 'cognitive_assessment_flynn_';
ps = {'01','02','03','04','05','06','07','08','09','10'};
allWin = nan(length(ps),62,200);
allLose = nan(length(ps),62,200);
for p = 1:length(ps)
    load([filePrefix ps{p}]);
    allWin(p,:,:) = ERP.data{1};
    allLose(p,:,:) = ERP.data{2};
end

% Get timepoints
timePoints = ERP.timepoints{1};

% Compute the difference waves
allDiff = allWin - allLose;

% Pick a channel
whichChannel = 34; % FCz - check chanlocs

% Compute reward positivity
timeWindow = [232 336]; % 232 ms to 336 ms post feedback
pntsWindow = dsearchn(timePoints',timeWindow'); % Find corresponding points
rewPos = squeeze(mean(allDiff(:,whichChannel,pntsWindow(1):pntsWindow(2)),3));

%% Plot grand averages
figure;
aveWin = squeeze(mean(allWin(:,whichChannel,:),1));
aveLose = squeeze(mean(allLose(:,whichChannel,:),1));
aveDiff = aveWin - aveLose;
plot(timePoints,aveWin,'LineWidth',2);
hold on;
plot(timePoints,aveLose,'LineWidth',2);
plot(timePoints,aveDiff,'k--','LineWidth',2);
hold off;
legend('Win','Lose','Difference');
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.LineWidth = 1.5;
ax.YLabel.String = 'Voltage (\muV)';
ax.XLabel.String = 'Time (ms)';
ax.FontWeight = 'bold';
ax.Box = 'off';
ax.YDir = 'reverse';
ax.Legend.Location = 'southwest';
ax.Legend.Box = 'off';