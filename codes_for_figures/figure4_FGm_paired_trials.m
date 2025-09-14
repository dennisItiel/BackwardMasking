%% Visualise FGm in noMask VSD map
% the ROIs through all the way is verOnly vs. horOnly (i.e, non-overlapping
% ROIs.In this script, I pair between cond2 and cond4 trials, one of each.
% The FGm is calculted, per each pair, as:
%   mean(cond2_verOnly - cond2_horOnly, cond4_horOnly - cond4_verOnly).
% Last modified: 07/10/2024

% SET PROJECT DIRECTORY
project_root = 'D:\Backup\Itiel_BM\BM_article_project';
cd(project_root)

% Add all subfolders to path so MATLAB can find your functions
addpath(genpath(project_root))

% Load data
eg_data = importdata('data\BM\Frodo\VSDI\2009_12_09\b\condsn4_bl.mat');
shaniRoi = importdata('data\BM\Frodo\VSDI\2009_12_09\shani_ROIs.mat');
chamber = importdata('data\BM\Frodo\VSDI\2009_12_09\ID_chamber.mat');

ver_only_roi = setdiff(shaniRoi.roiPix_fit_ver, shaniRoi.roiPix_fit_hor);
hor_only_roi = setdiff(shaniRoi.roiPix_fit_hor, shaniRoi.roiPix_fit_ver);
fig_roi = hor_only_roi;
bg_roi = ver_only_roi;

% Process maps
mapsTimes = 50:50:400;
mapsFrames = (mapsTimes / 10) + 27;
egMaps = mean(eg_data(:, mapsFrames, :), 3)-1;
smoothedMaps = mfilt2(egMaps, 100, 100);
chamberMaps = ones(10000, numel(mapsFrames)) * 1000;
chamberMaps(chamber,:) = smoothedMaps(chamber,:);

% display mimg of activity
figure; hold on
mimg2(chamberMaps, 100, 100, 0.00, 0.0025, mapsTimes); colormap mapgeog
cont = zeros(10000, 1); cont(shaniRoi.roiPix_fit_ver) = 1;
contour(reshape(cont, 100, 100)')

% plot shaded-errorbars of fig vs. bg
frames = 22:67;
pltTime = ((frames) - 27) * 10;
fig = squeeze(mean(eg_data(fig_roi, frames,:), 1))-1;
figMean = mean(fig, 2);
figSEM = std(fig, [], 2)/sqrt(size(fig,2));
bg = squeeze(mean(eg_data(bg_roi, frames,:), 1))-1;
bgMean = mean(bg, 2);
bgSEM = std(bg, [], 2)/sqrt(size(bg,2));

figure; hold on
shadedErrorBar(pltTime, figMean, figSEM)
shadedErrorBar(pltTime, bgMean, bgSEM)
xlabel('Time (ms)')
ylabel('DF/F')
legend({'' 'figure' '' 'background'})
hold off

% find statistics of fig and bg
tp1 = find(pltTime==80):find(pltTime==100);
tp2 = find(pltTime==180):find(pltTime==200);


stats(1,1) = mean(figMean(tp1)); % 1st peak mean
stats(1,2) = std(mean(fig(tp1,:), 1), [], 2)/sqrt(size(fig,2)); % 1st peak SEM
stats(1,3) = mean(figMean(tp2)); % 2nd peak mean
stats(1,4) = std(mean(fig(tp2,:), 1), [], 2)/sqrt(size(fig,2)); % 2nd peak SEM
stats(2,1) = mean(bgMean(tp1)); % 1st peak mean
stats(2,2) = std(mean(bg(tp1,:), 1), [], 2)/sqrt(size(bg,2)); % 1st peak SEM
stats(2,3) = mean(bgMean(tp2)); % 2nd peak mean
stats(2,4) = std(mean(bg(tp2,:), 1), [], 2)/sqrt(size(bg,2)); % 2nd peak SEM
stats(3,1) = mean(figMean(tp1) - bgMean(tp1)); % 1st peak mean
stats(3,2) = std(mean(fig(tp1,:) - bg(tp1,:), 1), [], 2)/sqrt(size(fig,2)); % 1st peak SEM
stats(3,3) = mean(figMean(tp2) - bgMean(tp2)); % 2nd peak mean
stats(3,4) = std(mean(fig(tp2,:) - bg(tp2,:), 1), [], 2)/sqrt(size(fig,2)); % 2nd peak SEM

rnksm(1) = ranksum(mean(fig(tp1,:), 1), mean(bg(tp1,:), 1));
rnksm(2) = ranksum(mean(fig(tp2,:), 1), mean(bg(tp2,:), 1));


%% Calculating the FGm of all trials

monkey = 'Frodo';
switch monkey
    case 'Frodo'
        dates = {'2009_09_16' '2009_09_23' '2009_10_07' '2009_10_14'};
        fileDates = {'1609_' '2309_' '0610_' '1410_'};
        load("data\BM\Frodo\VSDI\datesConsCell.mat");
    case 'Tolkin'
        dates = {'2012_05_24' '2012_05_30' '2012_06_27' '2012_07_04' '2012_07_11' '2012_07_18'};
        fileDates = {'2405_' '3005_' '2706_' '0407_' '1107_' '1807_'};
        load("data\BM\Tolkin\VSDI\datesConsCell.mat");
end


% parameters
full_length = 1:99; % frames for testing the model
fullTc = ((full_length)-27) * 10;
plt_frames = 20:60;
plt_time = fullTc(plt_frames);
conds_labels =  {'retino' 'noMaks' 'soa150' 'soa100' 'soa75' 'soa50'};
colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30"];

% Pre-allocation
all_con2Fig = cell(6,1);
all_con2Bg = cell(6,1);
all_con4Fig = cell(6,1);
all_con4Bg = cell(6,1);

for date_inx=1:numel(dates)
    DATE = dates{date_inx};
    disp(['Loading data from session ',DATE, '...'])

    vsd_folder = ['data\BM\', monkey, '\VSDI\', DATE];
    shaniRoi = importdata([vsd_folder, '\shani_ROIs.mat']);
    V1_pix = shaniRoi.rioPix_allV1;

    verOnly = intersect(V1_pix, setdiff(shaniRoi.roiPix_fit_ver, shaniRoi.roiPix_fit_hor));
    verOnly = setdiff(verOnly, shaniRoi.jointResidue);
    horOnly = intersect(V1_pix, setdiff(shaniRoi.roiPix_fit_hor, shaniRoi.roiPix_fit_ver));
    horOnly = setdiff(horOnly, shaniRoi.jointResidue);

    soaIn = {};

    figure; hold on
    for SOA_inx=2:6
        if datesConsCell{SOA_inx, date_inx} == "N"
            disp(['  condition ', conds_labels{SOA_inx}, ' - missing data'])
            continue
        end
        
        % load matrices
        try
            condsn2 = importdata(fullfile(vsd_folder, datesConsCell{SOA_inx, date_inx}, 'condsn2_new.mat'));
            condsn4 = importdata(fullfile(vsd_folder, datesConsCell{SOA_inx, date_inx}, 'condsn4_new.mat'));
            load(fullfile(vsd_folder, datesConsCell{SOA_inx, date_inx}, 'condsn.mat'))
            F02 = mean(condsn2(:, 25:27, :), 2);
            F04 = mean(condsn4(:, 25:27, :), 2);
            condsn2 = condsn2(:, full_length, :) ./ F02;
            condsn2 = condsn2 ./ mean(condsn(:,full_length,3), 3);
            condsn4 = condsn4(:, full_length, :) ./ F04;
            condsn4 = condsn4 ./ mean(condsn(:,full_length,3), 3);

        catch
            condsn2 = importdata(fullfile(vsd_folder, datesConsCell{SOA_inx, date_inx}, 'condsn2_bl.mat'));
            condsn4 = importdata(fullfile(vsd_folder, datesConsCell{SOA_inx, date_inx}, 'condsn4_bl.mat'));
            condsn2 = condsn2(:, full_length, :);
            condsn4 = condsn4(:, full_length, :);
        end

        % Balance data
        evenTrials = min(size(condsn2, 3), size(condsn4, 3));
        condsn2 = condsn2(:,:,1:evenTrials);
        condsn4 = condsn4(:,:,1:evenTrials);

        % Calculate mean(ROI) no norm.
        con2Fig = squeeze(mean(condsn2(verOnly,:,:), 1))-1;
        con2Bg = squeeze(mean(condsn2(horOnly,:,:), 1))-1;
        con4Fig = squeeze(mean(condsn4(horOnly,:,:), 1))-1;
        con4Bg = squeeze(mean(condsn4(verOnly,:,:), 1))-1;

        % Calculate FGm
        con2FGm = con2Fig - con2Bg;
        con4FGm = con4Fig - con4Bg;
        bothFGm = [con2FGm, con4FGm];
        
        % Plot Fig and Bg time-course
        soaIn = cat(1, soaIn, conds_labels{SOA_inx});
        limits = [-0.0003 0.0004];
        figTC = mean([con2Fig(plt_frames,:), con4Fig(plt_frames,:)], 2);
        bgTC = mean([con2Bg(plt_frames,:), con4Bg(plt_frames,:)], 2);
        if SOA_inx == 3
            condTitle = 'SOA150 fig vs. bg';
            subplot(2,3,1); hold on
            plot(plt_time, figTC, 'Color', colors(SOA_inx-1))
            plot(plt_time, bgTC, 'Color', colors(SOA_inx-1), 'LineStyle', '--', 'HandleVisibility','off')
            yline(0, '--', 'HandleVisibility','off')
            title(condTitle)
            ylim([-0.001 0.0025])
            grid on
            hold off
        elseif SOA_inx == 4
            condTitle = 'SOA100 fig vs. bg';
            subplot(2,3,2); hold on
            plot(plt_time, figTC, 'Color', colors(SOA_inx-1))
            plot(plt_time, bgTC, 'Color', colors(SOA_inx-1), 'LineStyle', '--', 'HandleVisibility','off')
            yline(0, '--', 'HandleVisibility','off')
            title(condTitle)
            ylim([-0.001 0.0025])
            grid on
            hold off
        elseif SOA_inx == 6
            condTitle = 'SOA50 fig vs. bg';
            subplot(2,3,3); hold on
            plot(plt_time, figTC, 'Color', colors(SOA_inx-1))
            plot(plt_time, bgTC, 'Color', colors(SOA_inx-1), 'LineStyle', '--', 'HandleVisibility','off')
            yline(0, '--', 'HandleVisibility','off')
            title(condTitle)
            ylim([-0.001 0.0025])
            grid on
            hold off
        end
 
        % Plot cond2 FGm
        subplot(2,3,4); hold on
        plot(fullTc(plt_frames), smoothdata(mean(con2FGm(plt_frames,:), 2), 'SmoothingFactor', 0.1), 'Color', colors(SOA_inx-1))
        yline(0, '--', 'HandleVisibility','off')
        title('ver target FGm')
        ylim(limits)
        grid on; hold off

        % Plot cond2 FGm
        subplot(2,3,5); hold on
        plot(fullTc(plt_frames), smoothdata(mean(con4FGm(plt_frames,:), 2), 'SmoothingFactor', 0.1), 'Color', colors(SOA_inx-1))
        yline(0, '--', 'HandleVisibility','off')
        title('hor target FGm')
        ylim(limits)
        grid on; hold off

        % Plot total FGm
        subplot(2,3,6); hold on
        plot(fullTc(plt_frames), smoothdata(mean(bothFGm(plt_frames,:), 2), 'SmoothingFactor', 0.1), 'Color', colors(SOA_inx-1))
        yline(0, '--', 'HandleVisibility','off')
        title('both target FGm')
        ylim(limits)
        grid on; hold off

        % Assigna to data struct
        n = size(con2Fig,2);
        rng(42)
        shuffledIdx = randperm(n);

        all_con2Fig{SOA_inx} = cat(2, all_con2Fig{SOA_inx}, con2Fig(:,shuffledIdx));
        all_con2Bg{SOA_inx} = cat(2, all_con2Bg{SOA_inx}, con2Bg(:,shuffledIdx));
        all_con4Fig{SOA_inx} = cat(2, all_con4Fig{SOA_inx}, con4Fig);
        all_con4Bg{SOA_inx} = cat(2, all_con4Bg{SOA_inx}, con4Bg);

        disp(['  condition ', conds_labels{SOA_inx}, ' - data loaded'])

    end

    legend(soaIn)
    sgtitle(DATE, 'Interpreter','none')
    hold off
    disp(['Session ',DATE, ' - data loaded'])

end

%% Calculate all fig_stats, bg_stats & fgm_stats

soa_labels = {'noMask'; 'SOA150'; 'SOA100'; 'SOA75'; 'SOA50'};
columnNames = {'soa_labels', 'bsln_value', 'bsln_SEM', 'bsln_rnksm', 'early_value', 'early_SEM', 'early_rnnksm', 'late_value', 'late_SEM', 'late_rnksm'};
time_points_ms = {-40:10:-20, 80:10:100, 180:10:200}; % Adjust as desired
time_points = cellfun(@(x) find(ismember(fullTc, x)), time_points_ms, 'UniformOutput', false);

% Calculate all FGms
allPairsFGm = cell(6,1);
shfldPairsFGm = cell(6,1);
fig_stats = table(soa_labels);  % Value of 'figure'
bg_stats = table(soa_labels);   % Value of 'back ground'
fgm_stats = table(soa_labels);  % Value of FGm

% Loop over releavant SOA
for SOA_inx=2:6

    pairedFig = mean(cat(3, all_con2Fig{SOA_inx}, all_con4Fig{SOA_inx}), 3);
    pairedBg = mean(cat(3, all_con2Bg{SOA_inx}, all_con4Bg{SOA_inx}), 3);
    all_con2FGm = all_con2Fig{SOA_inx} - all_con2Bg{SOA_inx};
    all_con4FGm = all_con4Fig{SOA_inx} - all_con4Bg{SOA_inx};
    pairedFGm = mean(cat(3, all_con2FGm, all_con4FGm), 3);
    allPairsFGm{SOA_inx} = pairedFGm;
    n = size(pairedFGm, 2);

    % Create shuffled ROIs, 100-iter shuffle
    shfldFig = nan(99,100);
    shfldBg = nan(99,100);
    shfldPairsFGm{SOA_inx} = nan(99,100);
    rng(42)
    random_mat = randi([0, 1], 100, n) * 2 - 1;
    for i=1:100
        random_vector = random_mat(i,:);
        singleShfldFig = cat(3, random_vector.*all_con2Fig{SOA_inx},...
                          random_vector.*all_con4Fig{SOA_inx});
        shfldFig(:, i) = mean(singleShfldFig, [2 3]);

        singleShfldBg = cat(3, random_vector.*all_con2Bg{SOA_inx},...
                               random_vector.*all_con4Bg{SOA_inx});
        shfldBg(:, i) = mean(singleShfldBg, [2 3]);

        singleShfldFGm = cat(3, random_vector.*(all_con2Fig{SOA_inx} - all_con2Bg{SOA_inx}),...
                             random_vector.*(all_con4Fig{SOA_inx} - all_con4Bg{SOA_inx}));
        shfldPairsFGm{SOA_inx}(:,i) = nanmean(singleShfldFGm, [2 3]);
    end

    % Assign values - bsln time-point
    for timeInx = 0:2
        tp = time_points{timeInx+1};

        % Value
        fig_stats{SOA_inx-1, 2 + 3*(timeInx)} = mean(pairedFig(tp,:), 'all');
        bg_stats{SOA_inx-1, 2 + 3*(timeInx)}  = mean(pairedBg(tp,:), 'all');
        fgm_stats{SOA_inx-1, 2 + 3*(timeInx)} = mean(pairedFGm(tp,:), 'all');

        % SEM
        fig_stats{SOA_inx-1, 3 + 3*(timeInx)} = std(mean(pairedFig(tp,:), 1)) / sqrt(n);
        bg_stats{SOA_inx-1, 3 + 3*(timeInx)}  = std(mean(pairedBg(tp,:), 1)) / sqrt(n);
        fgm_stats{SOA_inx-1, 3 + 3*(timeInx)} = std(mean(pairedFGm(tp,:), 1)) / sqrt(n);

        % RANKSM shuffled
        fig_stats{SOA_inx-1, 4 + 3*(timeInx)} = ranksum(mean(pairedFig(tp,:), 1), mean(shfldFig(tp,:), 1));
        bg_stats{SOA_inx-1, 4 + 3*(timeInx)}  = ranksum(mean(pairedBg(tp,:), 1), mean(shfldBg(tp,:), 1));
        fgm_stats{SOA_inx-1, 4 + 3*(timeInx)} = ranksum(mean(pairedFGm(tp,:), 1), mean(shfldPairsFGm{SOA_inx}(tp,:), 1));
    end
end

% Label all columns
fig_stats.Properties.VariableNames = columnNames;
bg_stats.Properties.VariableNames = columnNames;
fgm_stats.Properties.VariableNames = columnNames;



%% Plot all FGm curves one plot

fullTc = ((1:80)-27) * 10;
plt_frames = 20:60;

figure; hold on
for SOA_inx=2:6
%     y = smoothdata(mean(allPairsFGm{SOA_inx}(plt_frames,:), 2), 'SmoothingFactor', 0.1);
    y = mean(allPairsFGm{SOA_inx}(plt_frames,:), 2);

%       plot(fullTc(plt_frames), y, 'Color', colors(SOA_inx-1))

    sem = std(allPairsFGm{SOA_inx}(plt_frames,:), [], 2) ./ sqrt(size(allPairsFGm{SOA_inx}, 2));
    shadedErrorBar(fullTc(plt_frames), y, sem,...
        {'Color', colors(SOA_inx-1)}, 1)
    plot(fullTc(plt_frames), mean(shfldPairsFGm{SOA_inx}(plt_frames,:), 2),...
        'LineStyle','--', 'Color', colors(SOA_inx-1), 'HandleVisibility','off')

end
title('FG roiA_condA - roiB_condB', 'Interpreter','none')
legend({'' 'noMask' '' 'SOA150' '' 'SOA100' '' 'SOA75' '' 'SOA50'})
yline(0, '--', 'HandleVisibility','off')
xline(0', 'r', 'HandleVisibility','off')
ylabel('\Delta F/F')
xlabel('time (ms)')
grid on
% ylim([-0.05 0.08])

%% Calculate longSoa & shortSoa fig_stats, bg_stats & fgm_stats

soa_labels = {'longSoa'; 'shortSoa'};
columnNames = {'soa_labels', 'bsln_value', 'bsln_SEM', 'bsln_rnksm', 'early_value', 'early_SEM', 'early_rnnksm', 'late_value', 'late_SEM', 'late_rnksm'};
time_points_ms = {-40:10:-20
                  80:10:100
                  180:10:200}; % Adjust as desired
time_points = cellfun(@(x) find(ismember(fullTc, x)), time_points_ms, 'UniformOutput', false);

% Calculate all FGms
shfldPairsFGm = cell(2,1);
fgm_stats = table(soa_labels);  % Value of FGm

longSoa = cat(2, allPairsFGm{3}, allPairsFGm{4});
shortSoa = cat(2, allPairsFGm{5}, allPairsFGm{6});

for soaLength=1:2
    if soaLength==1
        combined2Fig = cat(2, all_con2Fig{3}, all_con2Fig{4});
        combined4Fig = cat(2, all_con4Fig{3}, all_con4Fig{4});
        combined2Bg = cat(2, all_con2Bg{3}, all_con2Bg{4});
        combined4Bg = cat(2, all_con4Bg{3}, all_con4Bg{4});

        combinedFGm = longSoa;
    else
        combined2Fig = cat(2, all_con2Fig{5}, all_con2Fig{6});
        combined4Fig = cat(2, all_con4Fig{5}, all_con4Fig{6});
        combined2Bg = cat(2, all_con2Bg{5}, all_con2Bg{6});
        combined4Bg = cat(2, all_con4Bg{5}, all_con4Bg{6});

        combinedFGm = shortSoa;
    end

    n = size(combinedFGm, 2);

    % Create shuffled ROIs
    shfldPairsFGm{soaLength} = nan(99,100);
    rng(4)
    random_mat = randi([0, 1], 100, n) * 2 - 1;
    for i=1:100
        random_vector = random_mat(i,:);
        singleShfldFGm = cat(3, random_vector.*(combined2Fig - combined2Bg),...
                                random_vector.*(combined4Fig - combined4Bg));
        shfldPairsFGm{soaLength}(:,i) = nanmean(singleShfldFGm, [2 3]);
        if soaLength == 1 && i == 1
            singleLongShfld = mean(singleShfldFGm, 3);
        end
        if soaLength == 2 && i == 1
            singleShortShfld = mean(singleShfldFGm, 3);
        end
    end
    
    % Assign values - bsln time-point
    for timeInx = 0:2
        tp = time_points{timeInx+1};
    
        % Value
        fgm_stats{soaLength, 2 + 3*(timeInx)} = mean(combinedFGm(tp,:), 'all');
    
        % SEM
        fgm_stats{soaLength, 3 + 3*(timeInx)} = std(mean(combinedFGm(tp,:), 1)) / sqrt(n);
    
        % RANKSM shuffled
        fgm_stats{soaLength, 4 + 3*(timeInx)} = ranksum(mean(combinedFGm(tp,:), 1), mean(shfldPairsFGm{soaLength}(tp,:), 1));
    end
end

% Label all columns
fgm_stats.Properties.VariableNames = columnNames;

%% Plot long vs. short FGm

fullTc = ((1:80)-27) * 10;
plt_frames = 20:60;
plt_time = fullTc(plt_frames);

switch monkey
    case 'Frodo'
        ylimits = [-4e-5 10e-5];
    case 'Tolkin'
        ylimits = [-2e-5 10e-5];
end


% Calculate timecourse sig.
rnksmLongShort = false(60,1);
for i = 1:60
    rnksmLong(i) = ranksum(longSoa(i,:), shfldPairsFGm{1}(i,:)) < 0.05;
    rnksmShort(i) = ranksum(shortSoa(i,:), shfldPairsFGm{2}(i,:)) < 0.05;
    rnksmLongShort(i) = ranksum(longSoa(i,:), shortSoa(i,:)) < 0.05;
end

% pltLong = smoothdata(mean(longSoa(plt_frames,:), 2), "movmean",3);
% pltShort = smoothdata(mean(shortSoa(plt_frames,:), 2), "movmean",3);
sem = @(x) std(x, [], 2)/sqrt(size(x, 2));

pltLong = mean(longSoa(plt_frames,:), 2);
pltShort = mean(shortSoa(plt_frames,:), 2);
semLong = sem(longSoa(plt_frames,:));
semShort = sem(shortSoa(plt_frames,:));

pltLongShfld = mean(singleLongShfld(plt_frames,:), 2);
pltShortShfld = mean(singleShortShfld(plt_frames,:), 2);
semLongShfld = sem(singleLongShfld(plt_frames,:));
semShortShfld = sem(singleShortShfld(plt_frames,:));

% Plot Shaded Error bars of long & short SOAs
figure; hold on
shadedErrorBar(plt_time, pltLong, semLong, {'Color', colors(1)}, 0)
shadedErrorBar(plt_time, pltShort, semLong, {'Color', colors(2)}, 0)
shadedErrorBar(plt_time, pltLongShfld, semLongShfld, {'Color', colors(1), 'LineStyle', '-'}, 0)
shadedErrorBar(plt_time, pltShortShfld, semShortShfld, {'Color', colors(2), 'LineStyle', '-'}, 0)


% Add sigstars to plot
try
    plot(fullTc(rnksmLong), 0.00012, 'b', 'Marker', '*')
    plot(fullTc(rnksmShort), 0.00011, 'r', 'Marker', '*')
    plot(fullTc(rnksmLongShort), 0.0001, 'k', 'Marker', '*')
end
title('FG roiA_condA - roiB_condB', 'Interpreter','none')
legend({'' 'longSoa' '' 'shortSoa'})
yline(0, 'HandleVisibility','off')
xline(0', 'r', 'HandleVisibility','off')
ylabel('\Delta F/F')
xlabel('time (ms)')
xlim([-30 350])
ylim(ylimits)
% grid on


% Bar plot first peak
time_points_ms = 80:10:100; % Adjust as desired
time_points = find(ismember(fullTc,time_points_ms));

figure; hold on
bar([1 2], fgm_stats.early_value)
errorbar([1 2], fgm_stats.early_value, fgm_stats.early_SEM, 'k', 'linestyle', 'none', 'linewidth', 2);
P = ranksum(mean(longSoa(time_points,:), 1), mean(shortSoa(time_points,:), 1));
sigstar([1 2], P)

% Bar plot 2nd peak
time_points_ms = 180:10:200; % Adjust as desired
time_points = find(ismember(fullTc,time_points_ms));

bar([4 5], fgm_stats.late_value)
errorbar([4 5], fgm_stats.late_value, fgm_stats.late_SEM, 'k', 'linestyle', 'none', 'linewidth', 2);

P = ranksum(mean(longSoa(time_points,:), 1), mean(shortSoa(time_points,:), 1));
sigstar([4 5], P)

set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'soa150&100' 'soa75&50' 'soa150&100' 'soa75&50'});
ylim(ylimits)
xlim([0 6])
hold off

%% Combine Frodo and Tolkin FGm

frodoFGm = importdata('data\frodo_paired_FGm.mat');
tolkinFGm = importdata('data\tolkin_paired_FGm.mat');
fullTc = ((1:80)-27) * 10;
plt_frames = 20:60;

figure; hold on
for SOA_inx=2:6
    rng(42)
    frodoNum = size(frodoFGm{SOA_inx}, 2);
    tolkinNum = size(tolkinFGm{SOA_inx}, 2);
    evenNum = min(frodoNum, tolkinNum);
    frodoSoaFgm = frodoFGm{SOA_inx}(:, randsample(frodoNum, evenNum));
    tolkinSoaFgm = tolkinFGm{SOA_inx}(:, randsample(tolkinNum, evenNum));
    FG = mean([frodoSoaFgm, tolkinSoaFgm], 2);

%     FG = mean([frodoFGm{SOA_inx}, tolkinFGm{SOA_inx}], 2);
%     y = smoothdata(mean(bothFGm{SOA_inx}(plt_frames,:), 2), 'SmoothingFactor', 0.1);
    y = FG(plt_frames,:);
    plot(fullTc(plt_frames), y)
end
title('FG roiA_condA - roiB_condB', 'Interpreter','none')
legend({'noMask' 'SOA150' 'SOA100' 'SOA75' 'SOA50'})
yline(0, '--', 'HandleVisibility','off')
xline(0', 'r', 'HandleVisibility','off')
ylabel('\Delta F/F')
xlabel('time (ms)')

%% SOA correlation with FGm - 2 peaks

bothFGearly = [6.00E-05
               6.61E-05
               5.84E-05
               5.14E-05];

bothFGlate = [9.21E-05
          4.72E-05
          1.61E-05
          1.28E-06];

soa = [150 100 75 50];

% Fit a polynomial of degree 1 (linear fit)
mainSeqEq = polyfit(soa, bothFGlate, 1);
mainSeqLine = polyval(mainSeqEq, soa);
shamSeqEq = polyfit(soa, bothFGearly, 1);
shamSeqLine = polyval(shamSeqEq, soa);

% Calculate R and P value
[r_early,P_early] = corr(bothFGearly, soa');
[r_late,P_late] = corr(bothFGlate, soa');

% Colors
colors = [
    [80 0 80]/255;
    [120 20 120]/255;
    [160 60 160]/255;
    [200 100 200]/255;
];

% Plot figure
figure; hold on
for k = 1:4
    plot(soa(k), bothFGlate(k), '.', 'MarkerSize', 30, 'Color', colors(k, :));
    plot(soa(k), bothFGearly(k), '.', 'MarkerSize', 30, 'Color', [128 128 128]/256);

end
plot(soa, mainSeqLine, 'r');
plot(soa, shamSeqLine, 'g');
xlim([35 165])
xticks(flip(soa))

% Titles and legends
title('FGm 1st & 2nd phases')
text(50, 0.00009, sprintf('early phase: r = %.3f, P = %.3f', r_early,P_early))
text(50, 0.00008, sprintf('late phase: r = %.3f, P = %.3f', r_late,P_late))
ylabel('FGm (\Delta F/F)')
xlabel('Performance (%)')
legend({'SOA150' 'SOA100' 'SOA75' 'SOA50'})
