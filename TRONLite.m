% Sample TRON script
% C. Hassall
% March, 2018
% Added a topo plot and image saving August, 2018

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

% Make a vector for the topo plot
topoVector = squeeze(mean(mean(allDiff(:,:,pntsWindow(1):pntsWindow(2)),3),1));

%% Plot grand averages

% Figure setup
fig = figure;
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
width = 14; % 3.0, 9.0, 14.0, or 19.0 aka minimum size, single column, 1.5-column, or double column
height = 7; % 24 cm max
padding = 0.1;
fig.PaperPosition = [0, 0, width, height];
fig.PaperSize = [width, height];
fig.Position = [padding, padding, width-padding, height-padding];
fig.Resize = 'off';
fig.InvertHardcopy = 'off';
fig.Color = 'white';

aveWin = squeeze(mean(allWin(:,whichChannel,:),1));
aveLose = squeeze(mean(allLose(:,whichChannel,:),1));
aveDiff = aveWin - aveLose;
plot(timePoints,aveWin,'k','LineWidth',2);
hold on;
plot(timePoints,aveLose,'k--','LineWidth',2);
plot(timePoints,aveDiff,'k:','LineWidth',2);
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

print('waveforms','-dtiff','-r600'); % Save to TIFF format

%% Topo

% Figure setup
fig = figure;
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
width = 9; % 3.0, 9.0, 14.0, or 19.0 aka minimum size, single column, 1.5-column, or double column
height = 10; % 24 cm max
padding = 0.1;
fig.PaperPosition = [0, 0, width, height];
fig.PaperSize = [width, height];
fig.Position = [padding, padding, width-padding, height-padding];
fig.Resize = 'off';
fig.InvertHardcopy = 'off';
fig.Color = 'white';

numContour = 6;
minTopo = min(topoVector);
maxTopo = max(topoVector);
t = topoplot(topoVector,'Standard-10-20-NEL-62.locs','numcontour',numContour,'maplimits',[minTopo, maxTopo],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','both');
t.Parent.XLim = [-0.6 0.6]; % For some reason this is necessary to make the nose/ears show up
t.Parent.YLim = [-0.6 0.6];
c = colorbar();
c.Label.String = 'Voltage (\muV)';
c.FontName = 'Arial';
c.FontWeight = 'bold';
c.FontSize = 12;
c.Location = 'southoutside';
colormap('parula'); % Don't use jet

fig.Color = 'white'; % Set figure background back to white in case topoplot changed it
print('topo','-dtiff','-r600'); % Save to TIFF format