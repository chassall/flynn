function fig = plotdisc(DISC)
%PLOTDISC Visualize the output of FLYNN
%   Plots artifact rejection summaries. Returns a figure handle because why
%   not.

scrsz = get(groot,'ScreenSize');
figure_loc = [1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/2]; % This is where figures will be drawn

% Flynn quotations
fquotes = {'On the other side of the screen, it all looks so easy.','How are you going to run the universe if you can''t answer a few unsolvable problems, huh?','Come on, you scuzzy data, be in there. Come on.','It''s time I leveled with you. I''m what you guys call a User.','I shouldn''t have written all of those tank programs.','It''s all in the wrists.','Greetings, programs!','Now for some real User power.','Did we make it? Hooray for our side.','Like the man says, there''s no problems, only solutions.','No sweat. I play video games better than anybody.','Come on, I''m - I''m scared of the dark. All this technology scares me.','Damn recognizer. Just go straight! I gotta get to that I/O tower.'};
this_fq = fquotes{randi(length(fquotes))};

fig = figure('Name',['FLYNN ' DISC.version],'NumberTitle','off','OuterPosition',figure_loc,'Toolbar','none','MenuBar','none');
whitebg(fig,'k');
set(gcf,'color','black');
flynn = imread('flynn.jpg');
subplot(2,4,[1 2 5 6]);
imshow(flynn,'InitialMagnification',33);
title(['\fontname{Courier}FLYNN ' DISC.version]);
xlabel(['\fontname{Courier}' sprintf(['"' this_fq '" -Flynn'])]); 

% ERP Summary
if ~isempty(DISC.ERPSum)
    subplot(2,4,3);
    thisProportion = DISC.ERPSum(:,4)./(DISC.ERPSum(:,3) + DISC.ERPSum(:,4));
    thisProportion = reshape(thisProportion,max(DISC.ERPSum(:,2)),DISC.N)';
    imagesc(thisProportion,[0 1]);
    c = colorbar;
    c.Label.String = 'Proportion Rejected';
    xticks(1:max(DISC.ERPSum(:,2)));
    xticklabels(DISC.ERPConditions);
    yticks(1:length(DISC.participants));
    yticklabels(DISC.participants);
    title('ERP Artifacts');
    xlabel('Condition');
    ylabel('Participant');
end

% FFT Summary
if ~isempty(DISC.FFTSum)
    subplot(2,4,4);
    thisProportion = DISC.FFTSum(:,4)./(DISC.FFTSum(:,3) + DISC.FFTSum(:,4));
    thisProportion = reshape(thisProportion,max(DISC.FFTSum(:,2)),DISC.N)';
    imagesc(thisProportion,[0 1]);
    c = colorbar;
    c.Label.String = 'Proportion Rejected';
    xticks(1:max(DISC.FFTSum(:,2)));
    xticklabels(DISC.FFTConditions);
    yticks(1:length(DISC.participants));
    yticklabels(DISC.participants);
    title('FFT Artifacts');
    xlabel('Condition');
    ylabel('Participant');
end

% WAV Summary
if ~isempty(DISC.WAVSum)
    subplot(2,4,7);
    thisProportion = DISC.WAVSum(:,4)./(DISC.WAVSum(:,3) + DISC.WAVSum(:,4));
    thisProportion = reshape(thisProportion,max(DISC.WAVSum(:,2)),DISC.N)';
    imagesc(thisProportion,[0 1]);
    c = colorbar;
    c.Label.String = 'Proportion Rejected';
    xticks(1:max(DISC.WAVSum(:,2)));
    xticklabels(DISC.WAVConditions);
    yticks(1:length(DISC.participants));
    yticklabels(DISC.participants);
    title('WAV Artifacts');
    xlabel('Condition');
    ylabel('Participant');
end

% ALL Summary
if ~isempty(DISC.ALLSum)
    subplot(2,4,8);
    thisProportion = DISC.ALLSum(:,4)./(DISC.ALLSum(:,3) + DISC.ALLSum(:,4));
    thisProportion = reshape(thisProportion,max(DISC.ALLSum(:,2)),DISC.N)';
    imagesc(thisProportion,[0 1]);
    c = colorbar;
    c.Label.String = 'Proportion Rejected';
    xticks(1:max(DISC.ALLSum(:,2)));
    xticklabels(DISC.ALLConditions);
    yticks(1:length(DISC.participants));
    yticklabels(DISC.participants);
    title('ALL Artifacts');
    xlabel('Condition');
    ylabel('Participant');
end

colormap('parula');

end

