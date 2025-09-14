%% Create Figure 3

% SET PROJECT DIRECTORY
project_root = 'D:\Backup\Itiel_BM\BM_article_project';
cd(project_root)

% Add all subfolders to path so MATLAB can find your functions
addpath(genpath(project_root))

% Choose e.g. session
egDateInx = 4;
dates = {'2009_09_16' '2009_09_23' '2009_10_07' '2009_10_14'};
egDate = dates{egDateInx};
vsd_folder = ['data\BM\Frodo\VSDI\', egDate];
fileDates = {'1609_' '2309_' '0610_' '1410_'};
load("data\BM\Frodo\VSDI\datesConsCell.mat");
load("data\BM\Frodo\VSDI\maskDatesCons.mat");
S = 1;

% Load matrices
load([vsd_folder, '\shani_ROIs'])
chamber = importdata([vsd_folder, '\ID_chamber.mat']);
V1_pix = rioPix_allV1;
ANALYSIS_ROI = V1_pix;


% Adjust as desire
full_time_course = ((1:80)-27) * 10;

% circle_contour_gui() % save a variable 'chamber'

%% Draw mimg of noMask e.g. session

condsn2 = importdata(fullfile(vsd_folder, datesConsCell{2, egDateInx}, 'condsn2_bl.mat'));
% condsn4 = importdata(fullfile(vsd_folder, datesConsCell{2, egDateInx}, 'condsn4_bl.mat'));
verCleanTrials = [1 2 4 6 7];

frames = 33:4:45;
pltMaps = ones(10000, numel(frames)) * 1000;
smoothMaps = mfilt2(nanmean(condsn2(:,frames,verCleanTrials), 3), 100, 100)-1;
pltMaps(chamber, :) = smoothMaps(chamber, :);

timeLabel = ([frames] - 27)*10;

figure;
mimg2(pltMaps, 100, 100, 0.00005, 0.002, timeLabel, 'k'); colormap mapgeog

% plot time-course
pltFrames = 25:62;
pltTime = (pltFrames - 27)*10;

noMask_signal = condsn2(ANALYSIS_ROI, pltFrames, verCleanTrials)-1;
noMask_mean = squeeze(mean(noMask_signal, 1));
noMask_sem = std(noMask_mean,[], 2) / sqrt(size(noMask_mean, 2));

figure; hold on
shadedErrorBar(pltTime, mean(noMask_mean, 2), noMask_sem);
% plot(pltTime, noMask_mean);


%% Draw mimg of mask150 e.g. session

condsn2 = importdata(strcat(vsd_folder, '\e\condsn2_bl.mat'));

frames = 47:5:62;
pltMaps = ones(10000, numel(frames)) * 1000;
smoothMaps = mfilt2(nanmean(condsn2(:,frames,:), 3), 100, 100)-1;
pltMaps(chamber, :) = smoothMaps(chamber, :);

timeLabel = ([frames] - 27)*10;

figure;
mimg2(pltMaps, 100, 100, 0.00005, 0.0025, timeLabel, 'k'); colormap mapgeog

% plot time-course
pltFrames = 25:67;
pltTime = (pltFrames - 27)*10;

mask150_signal = condsn2(ANALYSIS_ROI, pltFrames, :)-1;
mask150_mean = squeeze(mean(mask150_signal, 1));
mask150_sem = std(mask150_mean,[], 2) / sqrt(size(mask150_mean, 2));

figure; hold on
shadedErrorBar(pltTime, mean(mask150_mean, 2), mask150_sem);


%% Draw mimg of SOA150 e.g. session

condsn2 = importdata(fullfile(vsd_folder, datesConsCell{3, egDateInx}, 'condsn2_bl.mat'));
% condsn4 = importdata(fullfile(vsd_folder, datesConsCell{3, egDateInx}, 'condsn4_bl.mat'));

frames = 33:4:61;
pltMaps = ones(10000, numel(frames)) * 1000;
smoothMaps = mfilt2(nanmean(condsn2(:,frames,:), 3), 100, 100)-1;
pltMaps(chamber, :) = smoothMaps(chamber, :);

timeLabel = ([frames] - 27)*10;

figure;
mimg2(pltMaps, 100, 100, 0.00005, 0.0025, timeLabel, 'k'); colormap mapgeog

%% Create Time-Course of noMask vs. maskAlone- E.g session
% plot time-course
pltFrames = 25:62;
pltTime = (pltFrames - 27)*10;

% Colors palette
colors = [
    [80 0 80]/255;
    [120 20 120]/255;
    [160 60 160]/255;
    [200 100 200]/255;
];

% Load target+bg data
verCleanTrials = [1 2 4 6 7];

vertical = importdata(fullfile(vsd_folder, datesConsCell{2, egDateInx},'condsn2_bl.mat'));
vertical = nanmean(vertical(ANALYSIS_ROI,pltFrames, verCleanTrials), 1);
vertical = squeeze(vertical) -1;

norm_value = max(nanmean(vertical, 2));
norm_act = vertical / norm_value;
stimTC = nanmean(norm_act, 2);
targetSEM = nanstd(norm_act, [],  2) / sqrt(size(norm_act, 2));

figure(41);
shadedErrorBar(pltTime, mean(norm_act,2), targetSEM);


enclosedArea = cell(3,1);

for condInx = 1:3
    if maskDatesCons{condInx, egDateInx} == 'N' % omit 14/10 soa100?  | maskDatesCons{condInx, date_inx} == 'c'
        continue
    end
    if condInx < 3
        soaFolder = datesConsCell{condInx+2, egDateInx};
    else
        soaFolder = datesConsCell{condInx+3, egDateInx};
    end
    maskFolder = maskDatesCons{condInx, egDateInx};

    % Load BM data
    vertical = importdata(fullfile(vsd_folder, soaFolder, 'condsn2_bl.mat'));
    horizontal = importdata(fullfile(vsd_folder, soaFolder, 'condsn4_bl.mat'));
    common_act = cat(3, nanmean(vertical(ANALYSIS_ROI, : ,: ), 1), nanmean(horizontal(ANALYSIS_ROI, : ,: ), 1));
    common_act = squeeze(common_act(:,pltFrames,:)) -1;
    norm_value = max(nanmean(common_act, 2));
    norm_act = common_act / norm_value;
    bmTC = nanmean(norm_act, 2);
    bmSEM = nanstd(norm_act, [],  2) / sqrt(size(norm_act, 2));

    % Load mask alone data
    vertical = importdata(fullfile(vsd_folder, maskFolder, 'condsn2_bl.mat'));
    horizontal = importdata(fullfile(vsd_folder, maskFolder, 'condsn4_bl.mat'));
    common_act = cat(3, nanmean(vertical(ANALYSIS_ROI, : ,: ), 1), nanmean(horizontal(ANALYSIS_ROI, : ,: ), 1));
    common_act = squeeze(common_act(:,pltFrames,:)) -1;
    norm_value = max(nanmean(common_act, 2));
    norm_act = common_act / norm_value;
    maskTC = nanmean(norm_act, 2);
    maskSEM = nanstd(norm_act, [],  2) / sqrt(size(norm_act, 2));

    figure(41); hold on
    shadedErrorBar(pltTime, maskTC, maskSEM,...
        {'color', colors(condInx, :), 'linewidth', 1}, 0);
    hold off

    figure; hold on
    shadedErrorBar(pltTime, bmTC, bmSEM,...
        {'color', [0 0 128]/255, 'linewidth', 1}, 0);
    shadedErrorBar(pltTime, maskTC, maskSEM,...
        {'color', colors(condInx, :), 'linewidth', 1}, 0);
    shadedErrorBar(pltTime, stimTC, targetSEM,...
        {'color', [128 128 128]/255, 'linewidth', 1}, 0);
    xline(0, 'k')
    ylim([-0.1 1])
end

%% Create Time-Course of V1 activity in BM
% plot time-course
pltFrames = 25:65;
pltTime = (pltFrames - 27)*10;

% prepare soa150 activity
vertical = importdata('data\BM\Frodo\VSDI\2009_10_14\b\condsn2_bl.mat');
horizontal = importdata('data\BM\Frodo\VSDI\2009_10_14\b\condsn4_bl.mat');
soa150 = cat(3, vertical(ANALYSIS_ROI,pltFrames,:), horizontal(ANALYSIS_ROI,pltFrames,:));
soa150 = squeeze(mean(soa150, 1))-1;
sem150 = std(soa150,[], 2) / sqrt(size(soa150, 2));

% prepare soa100 activity
vertical = importdata('data\BM\Frodo\VSDI\2009_10_14\c\condsn2_bl.mat');
horizontal = importdata('data\BM\Frodo\VSDI\2009_10_14\c\condsn4_bl.mat');
soa100 = cat(3, vertical(ANALYSIS_ROI,pltFrames,:), horizontal(ANALYSIS_ROI,pltFrames,:));
soa100 = squeeze(mean(soa100, 1))-1;
sem100 = std(soa100,[], 2) / sqrt(size(soa100, 2));

% prepare soa50 activity
vertical = importdata('data\BM\Frodo\VSDI\2009_10_14\d\condsn2_bl.mat');
horizontal = importdata('data\BM\Frodo\VSDI\2009_10_14\d\condsn4_bl.mat');
soa50 = cat(3, vertical(ANALYSIS_ROI,pltFrames,:), horizontal(ANALYSIS_ROI,pltFrames,:));
soa50 = squeeze(mean(soa50, 1))-1;
sem50 = std(soa50,[], 2) / sqrt(size(soa50, 2));

% plot shaded errorbars
figure; hold on
shadedErrorBar(pltTime, mean(soa150, 2), sem150);
shadedErrorBar(pltTime, mean(soa100, 2), sem100);
shadedErrorBar(pltTime, mean(soa50, 2), sem50);
yline(0, 'k')
xline(0, 'k')
ylim([-0.0003 0.0025])
hold off

% Calculate stats
sem = @(x) nanstd(x)/sqrt(numel(x));
frame110ms = find(pltTime == 110);
frame130ms = find(pltTime == 130);
firstPeak150 = mean(soa150(frame110ms:frame130ms, :), 1);
firstPeak50 = mean(soa50(frame110ms:frame130ms, :), 1);
fprintf('1st peak respnse of soa150, 110-130ms: %.3d, SEM: %3.d, signrank(0): %.3d\n',...
    mean(firstPeak150), sem(firstPeak150), signrank(firstPeak150, 0))
fprintf('1st peak respnse of soa50, 110-130ms: %.3d, SEM: %3.d, signrank(0): %.3d\n',...
    mean(firstPeak50), sem(firstPeak50), signrank(firstPeak50, 0))
fprintf('1st peak ranksum of soa150&50 , 110-130ms: %.3d',...
    ranksum(firstPeak150, firstPeak50))

frame280ms = find(pltTime == 280);
frame300ms = find(pltTime == 300);
scndPeak150 = mean(soa150(frame280ms:frame300ms, :), 1);
scndPeak50 = mean(soa50(frame280ms:frame300ms, :), 1);
fprintf('2nd peak respnse of soa150, 280-300ms: %.3d, SEM: %3.d, signrank(0): %.3d\n',...
    mean(scndPeak150), sem(scndPeak150), signrank(scndPeak150, 0))
fprintf('2nd peak respnse of soa50, 280-300ms: %.3d, SEM: %3.d, signrank(0): %.3d\n',...
    mean(scndPeak50), sem(scndPeak50), signrank(scndPeak50, 0))
fprintf('2nd peak ranksum of soa150&50 , 280-300ms: %.3d',...
    ranksum(scndPeak150, scndPeak50))


%% Collect all data for Grand Analysis

allBmData = cell(3,1); % 150, 100, 50
allMaskData = cell(3,1); % 150, 100, 50, tg+bg
allStimData = [];
pltFrames = 25:65;
pltTime = (pltFrames - 27)*10;

for date_inx = 1:4
    date = dates{date_inx};
    vsd_folder = ['data\BM\Frodo\VSDI\', date, '\'];
    load([vsd_folder, '\shani_ROIs'])
    analysisROI = rioPix_allV1;

    % Load target+bg data
    stimFolder = maskDatesCons{4, date_inx};
    if stimFolder ~= 'N'
        vertical = importdata(fullfile(vsd_folder, stimFolder, 'condsn2_bl.mat'));
        horizontal = importdata(fullfile(vsd_folder, stimFolder, 'condsn4_bl.mat'));
        common_act = cat(3, nanmean(vertical(analysisROI, : ,: ), 1), nanmean(horizontal(analysisROI, : ,: ), 1));
        common_act = squeeze(common_act(:,pltFrames,:)) -1;
        norm_value = max(nanmean(common_act, 2));
        allStimData = cat(2, allStimData, common_act / norm_value);
    end


    % Load BM data
    for condInx = 1:3
        if condInx < 3
            bmFolder = datesConsCell{condInx+2, date_inx};
        else
            bmFolder = datesConsCell{condInx+3, date_inx};
        end
        if bmFolder == 'N'
            continue
        end
        vertical = importdata(fullfile(vsd_folder, bmFolder, 'condsn2_bl.mat'));
        horizontal = importdata(fullfile(vsd_folder, bmFolder, 'condsn4_bl.mat'));
        common_act = cat(3, nanmean(vertical(analysisROI, : ,: ), 1), nanmean(horizontal(analysisROI, : ,: ), 1));
        common_act = squeeze(common_act(:,pltFrames,:)) -1;
        norm_value = max(nanmean(common_act, 2));
        allBmData{condInx} = cat(2, allBmData{condInx}, common_act / norm_value);
    end

    % Load mask alone data
    for condInx = 1:3
        maskFolder = maskDatesCons{condInx, date_inx};
        if maskFolder == 'N'
            continue
        end
        vertical = importdata(fullfile(vsd_folder, maskFolder, 'condsn2_bl.mat'));
        horizontal = importdata(fullfile(vsd_folder, maskFolder, 'condsn4_bl.mat'));
        common_act = cat(3, nanmean(vertical(analysisROI, : ,: ), 1), nanmean(horizontal(analysisROI, : ,: ), 1));
        common_act = squeeze(common_act(:,pltFrames,:)) -1;
        norm_value = max(nanmean(common_act, 2));
        allMaskData{condInx} = cat(2, allMaskData{condInx}, common_act / norm_value);
    end
end


%% GA of mask modulation
enclosedArea = cell(3,1);

for condInx = 1:3
    maskTC = nanmean(allMaskData{condInx}, 2);
    bmValues = allBmData{condInx};

    % Specify the starting index and range
    switch condInx
        case 1
            startIndex = 18; % Starting index for x
            endIndex = 38;
        case 2
            startIndex = 13; % Starting index for x
            endIndex = 33;
        case 3
            startIndex = 8; % Starting index for x
            endIndex = 28;
    end
    startIndex = 15; % Starting index for x

    for j=1:size(bmValues, 2)   %pay attention to start index!!
        aboveIndices = maskTC > bmValues(:,j);
        area = maskTC.*aboveIndices - bmValues(:,j).*aboveIndices;
        enclosedArea{condInx}(j) = sum(area(startIndex:endIndex));
    end
end


meanArea = cellfun(@mean, enclosedArea);
sem = @(x) nanstd(x)/sqrt(numel(x));
semArea = cellfun(sem, enclosedArea);
sgrnkArea = cellfun(@(x) signrank(x, 0), enclosedArea);

% Plot bars of mean mask-modulation
figure; hold on
bar([1 2 3], meanArea)
errorbar([], meanArea, semArea, 'r', 'linestyle', 'none', 'linewidth', 2);
title('Mask Modulation')
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'soa150' 'soa100' 'soa50'});

% Get pValues
p150_100 = ranksum(enclosedArea{1}, enclosedArea{2});
p150_50 = ranksum(enclosedArea{1}, enclosedArea{3});
p100_50 = ranksum(enclosedArea{2}, enclosedArea{3});

% Plot the stars
groups = {[1 2] [1 3] [2 3]};
sigstar(groups,[(p150_100), (p150_50), (p100_50)], 1);

% sigstar([1 2], P)
hold off

%% GA of Interaction Index - according to Rev #2
% ((mask-alone + stimulus-alone) - BM) / (mask-alone + stimulus-alone).

interactionIndex = cell(3,1);
stimTC = nanmean(allStimData, 2);

figure; hold on
for condInx = 1:3
    maskTC = nanmean(allMaskData{condInx}, 2);
    bmValues = allBmData{condInx};

    % Specify the starting index and range
    switch condInx
        case 1
            startIndex = find(pltTime == 200); % Starting index for x
            endIndex = find(pltTime == 300);
        case 2
            startIndex = find(pltTime == 150); % Starting index for x
            endIndex = find(pltTime == 250);
        case 3
            startIndex = find(pltTime == 100); % Starting index for x
            endIndex = find(pltTime == 200);
    end
    nominatorTC = maskTC + stimTC - mean(bmValues, 2);
    denominatorTC = maskTC + stimTC;
    plot(pltTime, nominatorTC ./ denominatorTC)

    nominator = sum(maskTC(startIndex:endIndex)) + sum(stimTC(startIndex:endIndex)) - sum(bmValues(startIndex:endIndex,:), 1);
    denominator = sum(maskTC(startIndex:endIndex)) + sum(stimTC(startIndex:endIndex));
    interactionIndex{condInx} = nominator ./ denominator;
end
title('Interaction Index')
legend({'soa150' 'soa100' 'soa50'})
ylim([-0.5,1])
yline(0, 'k')

% Calculate Index and SEM
meanIndex = cellfun(@mean, interactionIndex);
sem = @(x) nanstd(x)/sqrt(numel(x));
semIndex = cellfun(sem, interactionIndex);
sgrnkIndex = cellfun(@(x) signrank(x, 0), interactionIndex);

% Plot bars of mean mask-modulation
figure; hold on
bar([1 2 3], meanIndex)
errorbar([], meanIndex, semIndex, 'r', 'linestyle', 'none', 'linewidth', 2);
title('Interaction Index')
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'soa150' 'soa100' 'soa50'});

% Get pValues
p150_100 = ranksum(interactionIndex{1}, interactionIndex{2});
p150_50 = ranksum(interactionIndex{1}, interactionIndex{3});
p100_50 = ranksum(interactionIndex{2}, interactionIndex{3});

% Plot the stars
groups = {[1 2] [1 3] [2 3]};
sigstar(groups,[(p150_100), (p150_50), (p100_50)], 1);

% sigstar([1 2], P)
hold off

%% GA of mask-stim intefearance
overlapArea = cell(3,1);

stimTC = nanmean(allStimData, 2);

startIndex = find(pltTime == 50);
endIndex = find(pltTime == 230);
stimArea = sum(stimTC(startIndex:endIndex));

for condInx = 1:3
    maskValues = allMaskData{condInx};
    for j=1:size(maskValues, 2)
        underMaskIndices = stimTC >= maskValues(:,j);
        underStimIndices = stimTC < maskValues(:,j);
                        
        underMaskArea = maskValues(:,j).*underMaskIndices;
        underStimArea = stimTC.*underStimIndices;

        area = underMaskArea' + underStimArea';
        overlapArea{condInx}(j) = sum(area(startIndex:endIndex));
    end
end


meanArea = cellfun(@mean, overlapArea)/stimArea;
sem = @(x) nanstd(x)/sqrt(numel(x));
semArea = cellfun(sem, overlapArea)/stimArea;
sgrnkArea = cellfun(@(x) signrank(x, 0), overlapArea);

% Plot bars of mean mask-modulation
figure; hold on
bar([1 2 3], meanArea)
errorbar([], meanArea, semArea, 'r', 'linestyle', 'none', 'linewidth', 2);
title('Mask-stim interference')
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'soa150' 'soa100' 'soa50'});

% Get pValues
p150_100 = ranksum(overlapArea{1}, overlapArea{2});
p150_50 = ranksum(overlapArea{1}, overlapArea{3});
p100_50 = ranksum(overlapArea{2}, overlapArea{3});

% Plot the stars
groups = {[1 2] [1 3] [2 3]};
sigstar(groups,[(p150_100), (p150_50), (p100_50)], 1);

% sigstar([1 2], P)
hold off