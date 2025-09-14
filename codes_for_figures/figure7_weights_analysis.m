%%% Weights Dynamics Analysis for Registrated Data

% SET PROJECT DIRECTORY
project_root = 'D:\Backup\Itiel_BM\BM_article_project';
cd(project_root)

% Add all subfolders to path so MATLAB can find your functions
addpath(genpath(project_root))

% Change monkey
monkey = 'frodo';

% Load data
switch monkey
    case 'frodo'
        % Frodo early sessions
        weightsTC = importdata('data\normalized2_regist\frodo_early_sessions\weightsTC.mat');
        accuracyCell = importdata('data\normalized2_regist\frodo_early_sessions\accuracyCell.mat');
        roiHor = importdata('data\normalized2_regist\frodo_early_sessions\ellipseHor.mat');
        roiVer = importdata('data\normalized2_regist\frodo_early_sessions\ellipseVer.mat');
        V1_pix = importdata('data\normalized2_regist\frodo_early_sessions\V1_pix.mat');
        residue = importdata('data\normalized2_regist\frodo_early_sessions\residue.mat');
    case 'tolkin_early'
        % Tolkin early sessions
        weightsTC = importdata('data\normalized2_regist\tolkin_early_sessions\weightsTC.mat');
        roiHor = importdata('data\normalized2_regist\tolkin_early_sessions\roiHor.mat');
        roiVer = importdata('datanormalized2_regist\tolkin_early_sessions\roiVer.mat');
        V1_pix = importdata('data\normalized2_regist\tolkin_early_sessions\V1_pix.mat');
        residue = [];
    case 'tolkin_late'
        % Tolkin late sessions
        weightsTC = importdata('data\normalized2_regist\tolkin_late_sessions\weightsTC_456.mat');
        roiHor = importdata('data\normalized2_regist\tolkin_late_sessions\roiHor.mat');
        roiVer = importdata('data\normalized2_regist\tolkin_late_sessions\roiVer.mat');
        V1_pix = importdata('data\normalized2_regist\tolkin_late_sessions\V1_pix.mat');
        residue = [];
end


% ROI assign
[~,Iver,~] = intersect(V1_pix, roiVer);
[~,Ihor,~] = intersect(V1_pix, roiHor);
[~,Ires,~] = intersect(V1_pix, residue);
[bgHor,IbgHor] = setdiff(V1_pix, roiHor);
[bgVer,IbgVer] = setdiff(V1_pix, roiVer);

verOnly = setdiff(roiVer, union(roiHor, residue));
[~,IverOnly,~] = intersect(V1_pix, verOnly);
horOnly = setdiff(roiHor, union(roiVer, residue));
[~,IhorOnly,~] = intersect(V1_pix, horOnly);

sharedRegion = intersect(roiHor, roiVer);
[~,IsharedRegion,~] = intersect(V1_pix, sharedRegion);

bothTargetsOnly = union(verOnly, horOnly);
IbothTargetsOnly = union(IverOnly, IhorOnly);

bothTargets = union(roiVer, roiHor);
IbothTargets = union(Iver, Ihor);

[noTargets,InoTargets] = setdiff(V1_pix, bothTargets);

allTheRest = 1:numel(V1_pix);
allTheRest = setdiff(allTheRest, union(IverOnly, IhorOnly));

full_time_course = ((1:80)-27) * 10;
plt_frames = 23:62;
plt_time = full_time_course(plt_frames);

%% Weights Dynamics plot

weightsMeanTC = squeeze(nanmean(weightsTC, 2));
weigthsSemCell = cell(6, 3);

for ii=1:6
    weigthsSemCell{ii,1} = nanstd(squeeze(weightsTC(ii,IhorOnly,:)), 1) / sqrt(size(IhorOnly, 1));
    weigthsSemCell{ii,2} = nanstd(squeeze(weightsTC(ii,IverOnly,:)), 1) / sqrt(size(IverOnly, 1));
    weigthsSemCell{ii,3} = nanstd(squeeze(weightsTC(ii,IsharedRegion,:)), 1) / sqrt(size(IsharedRegion, 1));
end

% plot 4 SOAs graphs
soa_labels = {'SOA150' 'SOA100' 'SOA75' 'SOA50'};
for i=1:4
    figure(); hold on
    ylim([-1.5e-3 1.5e-3])
    h1 = shadedErrorBar(plt_time, squeeze(nanmean(weightsTC(i+2,IverOnly,:), 2)),...
        weigthsSemCell{i+2,1}, {'color', [114 147 203]./255, 'linewidth', 1}, 1);
    h2 = shadedErrorBar(plt_time, squeeze(nanmean(weightsTC(i+2,IhorOnly,:), 2)),...
        weigthsSemCell{i+2,2}, {'color', [211 94 96]./255, 'linewidth', 1}, 1);
    xline(0, 'r', 'HandleVisibility','off')
    yline(0, '--', 'HandleVisibility','off')
    xlabel('time (ms)'); ylabel('Mean Weights')
    title(soa_labels{i})
    legend([h1.mainLine h2.mainLine], {'vertical only', 'horizontal only'})
end

diff_tc = squeeze(nanmean(weightsTC(:,IverOnly,:), 2)) - squeeze(nanmean(weightsTC(:,IhorOnly,:), 2));

figure(); hold on
plot(plt_time, diff_tc')
xline(0, 'r', 'HandleVisibility','off')
xlabel('frames'); ylabel('Mean Weights Difference')
yline(0, '--', 'HandleVisibility','off')
conds_labels = {'retino' 'noMask' 'SOA150' 'SOA100' 'SOA75' 'SOA50'};
legend(conds_labels,...
    'Location', 'northwest')
title('Mean Weights Difference')

%% Combine Tolkin's early and late DS

% earlyDiffTc = diff_tc;
% lateDiffTc = diff_tc;

diff_tc = cat(3, earlyDiffTc, lateDiffTc);
diff_tc = mean(diff_tc, 3);

figure; hold on
plot(plt_time, diff_tc)
xline(0)
yline(0)
legend(conds_labels)
ylabel('/delta')
title('Mean Weights Difference')

%% Combine Frodo & Tolkin

% frdodDiffTc = diff_tc;
% tolkinDiffTc = diff_tc;

diff_tc = cat(3, frdodDiffTc, tolkinDiffTc);
diff_tc = mean(diff_tc, 3);

figure; hold on
plot(plt_time, diff_tc)
xline(0)
yline(0)
legend(conds_labels)
ylabel('/delta')
title('Mean Weights Difference')

%% First peak & Second peak correlation with Difference Weights ROI

tp1 = 80;
tp2 = 190;
pk1 = find(plt_time == tp1);
pk2 = find(plt_time == tp2);

x = [50 75 100 150];
firstPeakDiff = flip(diff_tc(3:end,pk1));
secondPeakDiff = flip(diff_tc(3:end,pk2));


% Fit a polynomial of degree 1 (linear fit)
firstSeqEq = polyfit(x, firstPeakDiff, 1);
firstSeqLine = polyval(firstSeqEq, x);
secondSeqEq = polyfit(x, secondPeakDiff, 1);
secondSeqLine = polyval(secondSeqEq, x);

% Calculate R and P value
[firstR,firstP] = corr(x', firstPeakDiff);
[secondR,secondP] = corr(x', secondPeakDiff);


figure; hold on
plot(x, firstPeakDiff, '.', 'MarkerSize', 25, 'color', 	"#4DBEEE")
plot(x, secondPeakDiff, '.', 'MarkerSize', 25, 'color', "#0072BD")
plot(x, firstSeqLine, 'color', 	"#4DBEEE")
plot(x, secondSeqLine, 'color', "#0072BD")
xticks(x)
xlim([25 175])
title('Mean Weights Difference Correlation with SOA')
yloc = max(secondPeakDiff);
text(35, yloc, sprintf('1st peak: R=%.3f P=%.3f', firstR, firstP))
text(35, yloc-0.0005, sprintf('2nd peak: R=%.3f P=%.3f', secondR, secondP))
legend({'1st peak' '2nd peak'})

%% weights map visualization & weights histogram - 2X2 grid

timePoints = {70, 180};
conds = [3 6];
FIRST_FRAME = 21;
NUM_TRAIN_FRAMES = 5;
S = 0.5;
sbpltIdx = 1;

% Define the minimum and maximum values for clipping the colormap
minValue = -0.0035; % Adjust as needed based on your data
maxValue = 0.0035; % Adjust as needed based on your data

% create contour matrices
horCont = false(10000,1); horCont(roiHor) = 1;
horCont = reshape(horCont, 100, 100)';
verCont = false(10000,1); verCont(roiVer) = 1;
verCont = reshape(verCont, 100, 100)';
contourV1 = false(100, 100); contourV1(V1_pix) = 1;
contourV1 = contourV1';

%  Reset stats
N_hor = numel(horOnly);
N_ver = numel(verOnly);
hor_early = nan(2,2); % soa150/50, mean/std
ver_early = nan(2,2); % soa150/50, mean/std
hor_late = nan(2,2); % soa150/50, mean/std
ver_late = nan(2,2); % soa150/50, mean/std
rnksm_early = nan(2,1); % soa150/50
rnksm_late = nan(2,1); % soa150/50
AUC_early = nan(2,1); % soa150/50
AUC_late = nan(2,1); % soa150/50
weights_eary_cell = cell(2,2); % soa150/50, hor/ver
weights_late_cell = cell(2,2); % soa150/50, hor/ver

h=figure(100); hold on; set(h,'WindowStyle','docked')
h=figure(101); hold on; set(h,'WindowStyle','docked')
h=figure(102); hold on; set(h,'WindowStyle','docked')

for gridRow = 1:2 % soa150/50
    SOA_inx = conds(gridRow);
    for gridCol = 1:2  % early/late
        % Adjust as desire
        time_point = timePoints{gridCol};
        frameIdx = find(plt_time == time_point(1)):find(plt_time == time_point(end));
        map = zeros(10000, 1);
        % map(V1_pix) = squeeze(mean(all_data{SOA_inx, 2}(:,frame,:), [1 2]));
        map(V1_pix) = squeeze(mean(weightsTC(SOA_inx,:,frameIdx), [1 3]));


        % select ROI pixels
        horValues = map(horOnly);
        verValues = map(verOnly);
        sharedValues = map(sharedRegion);

        % map(map < 1e-3 & map > -2e-3) = 0;
        % map(map > prctile(map, 10) & map < prctile(map, 90)) = 0;
        % [min(map), max(map)]

        % filter and prepare to plot
        map = mfilt2(map, 100, 100, S, 'lm');
        map(setdiff(1:10000, V1_pix)) = 1000;
        map = reshape(map, 100, 100);

        % plot vsd map
        figure(100);
        subplot(2, 2, sbpltIdx);
        fg = imagesc(reshape(map, 100, 100)', [minValue maxValue]); colormap mapgeog
        hold on
        contour(horCont, 'Color', '#808080')
        contour(verCont, 'Color', '#808080')
        contour(contourV1, 'k');
        title([conds_labels{SOA_inx}, ', time:', num2str(time_point)])
        if sbpltIdx==1
            text(10, 5, sprintf('frames: %d:%d', frameIdx + FIRST_FRAME-1, frameIdx + FIRST_FRAME + NUM_TRAIN_FRAMES -2))
            text(10, 10, sprintf('smoothing: S=%g', S))
            text(10,15, sprintf('cliping: [%g, %g]', minValue, maxValue))
        end
        hold off

        % plot histogram
        figure(101);
        subplot(2, 2, sbpltIdx);        hold on
        h1 = histogram(horValues, 'Normalization', 'probability');
        h2 = histogram(verValues, 'Normalization', 'probability');
        h1.BinWidth = 0.0005;
        h2.BinWidth = 0.0005;
        ylim([0 0.15]); xlim([-10e-3 10e-3])
        xline(mean(horValues), 'Color', '#0343DF')
        xline(mean(verValues), 'Color', '#8C000F')
        title(sprintf('Weights in %s, Hor & Ver', conds_labels{SOA_inx}))
        subtitle(ranksum(horValues, verValues))
        if sbpltIdx==1
            legend({'horOnly weights' 'verOnly weights'})
        end

        % plot ROC and AUC
        figure(102);
        subplot(2, 2, sbpltIdx);        hold on
        plotROC(verValues, horValues)

        sbpltIdx = sbpltIdx + 1;

        % enter stats
        if gridCol == 1
            hor_early(gridRow, 1) = mean(horValues);
            hor_early(gridRow, 2) = std(horValues);
            ver_early(gridRow, 1) = mean(verValues);
            ver_early(gridRow, 2) = std(verValues);
            rnksm_early(gridRow) = ranksum(horValues, verValues);
            AUC_early(gridRow) = calculateAUC(horValues, verValues);
            weights_eary_cell{gridRow, 1} = horValues;
            weights_eary_cell{gridRow, 2} = verValues;
        else
            hor_late(gridRow, 1) = mean(horValues);
            hor_late(gridRow, 2) = std(horValues);
            ver_late(gridRow, 1) = mean(verValues);
            ver_late(gridRow, 2) = std(verValues);
            rnksm_late(gridRow) = ranksum(horValues, verValues);
            AUC_late(gridRow) = calculateAUC(horValues, verValues);
            weights_late_cell{gridRow, 1} = horValues;
            weights_late_cell{gridRow, 2} = verValues;
        end
    end
end

figure(100); hold off
figure(101); hold off
figure(102); hold off


%% Calculate correlation with TC-retino-maps

CORRELATION_MODE = 'single'; % switch between single frame and TC ('single', 'TC')
RETINO_MAP_MS = 80;

full_time_course = ((1:80)-27) * 10;
plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
retinoMapsIndex = find(plt_time == RETINO_MAP_MS);
retinoMaps = squeeze(weightsTC(1,:,:));

% retinoMaps = zeros(10000,40); retinoMaps(V1_pix,:) = squeeze(retinoWeights);
% retinoMaps = mfilt2(retinoMaps, 100, 100, 1, 'lm');
% figure; mimg2(retinoMaps, 100, 100,-0.002,0.002,-40:10:400); colormap mapgeog

soa_names = ["SOA150"; "SOA100"; "SOA75"; "SOA50"];


loc = 1;
allCorls = nan(40,4);

figure; hold on
sgtitle("Normalized Cross-Correlation of Weight-maps Pairs")
for soaIdx=3:6
%     soaMap = zeros(10000,40); soaMap(V1_pix,:) = squeeze(weightsTC(soaIdx,:,:));
%     soaMap = mfilt2(soaMap, 100, 100, 1, 'lm');

    soaMap = squeeze(weightsTC(soaIdx,:,:));

    corrTc = nan(40,1);
    for i=1:40
        switch CORRELATION_MODE
            case 'single'
                referenceMap = retinoMaps(:,retinoMapsIndex);
            case 'TC'
                referenceMap = retinoMaps(:,i);
        end
        corl = corrcoef(soaMap(:,i), referenceMap);
        corrTc(i) = corl(2);

    end

    subplot(2,2,loc)
    plot(plt_time, corrTc)
    yline(0, '--')
    xline(0, 'r')
    ylim([-0.1, 0.7])
    title(sprintf("NCC of %s and retinoMap", soa_names(soaIdx-2)))

%     earlyCorrelation{inx_b-3, inx_a-1} = corrTc(time_points(1));
%     lateCorrelation{inx_b-3, inx_a-1} = corrTc(time_points(2));
    allCorls(:,soaIdx-2) = corrTc;
    loc = loc + 1;
end

figure; hold on
plot(plt_time, allCorls)
xline(0)
yline(0)
legend(soa_names)
ylabel('correlation')
title(sprintf('TC Correlation of SOA-weights maps with %s retino-weight map', CORRELATION_MODE))

%% Combine Tolkin's early and late DS

% earlyCorls = allCorls;
% lateCorls = allCorls;

allCorls = cat(3, earlyCorls, lateCorls);
allCorls = mean(allCorls, 3);

figure; hold on
plot(plt_time, allCorls)
xline(0)
yline(0)
legend(soa_names)
ylabel('correlatio')
title('TC Correlation of SOA-weights maps with single retino-weight map')

%% Combine Frodo & Tolkin

% frodoCorls = allCorls;
% tolkinCorls = allCorls;

allCorls = cat(3, frodoCorls, tolkinCorls);
allCorls = mean(allCorls, 3);

figure; hold on
plot(plt_time, allCorls)
xline(0)
yline(0)
legend(soa_names)
ylabel('correlatio')
title('TC Correlation of SOA-weights maps with single retino-weight map')

%% First peak & Second peak correlation with single retino-map

tp1 = 80;
tp2 = 190;
pk1 = find(plt_time == tp1);
pk2 = find(plt_time == tp2);

x = [50 75 100 150];
firstPeakCorr = flip(allCorls(pk1,:));
secondPeakCorr = flip(allCorls(pk2,:));


% Fit a polynomial of degree 1 (linear fit)
firstSeqEq = polyfit(x, firstPeakCorr, 1);
firstSeqLine = polyval(firstSeqEq, x);
secondSeqEq = polyfit(x, secondPeakCorr, 1);
secondSeqLine = polyval(secondSeqEq, x);

% Calculate R and P value
[firstR,firstP] = corr(x', firstPeakCorr');
[secondR,secondP] = corr(x', secondPeakCorr');


figure; hold on
plot(x, firstPeakCorr, '.', 'MarkerSize', 25, 'color', 	"#4DBEEE")
plot(x, secondPeakCorr, '.', 'MarkerSize', 25, 'color', "#0072BD")
plot(x, firstSeqLine, 'color', 	"#4DBEEE")
plot(x, secondSeqLine, 'color', "#0072BD")
xticks(x)
xlim([25 175])
title('Weight Map Correlation with Retino 190ms Map')
yloc = max(secondPeakCorr);
text(35, yloc, sprintf('1st peak: R=%.3f P=%.3f', firstR, firstP))
text(35, yloc-0.05, sprintf('2nd peak: R=%.3f P=%.3f', secondR, secondP))
legend({'1st peak' '2nd peak'})



%% functions

calculateAUC(horValues, verValues)

function auc = calculateAUC(dist1, dist2)
    labels = [zeros(size(dist1)); ones(size(dist2))]; % 1 for dist1, 0 for dist2
    scores = [dist1; dist2]; % Combine distributions
    [~,~,~,auc] = perfcurve(labels, scores, 1); % Compute AUC
end

function plotROC(d1, d2)
    % Compute true positive rate (TPR) and false positive rate (FPR)
    scores = [d1(:); d2(:)];
    labels = [ones(size(d1(:))); zeros(size(d2(:)))]; % 1 for d1, 0 for d2
    
    [X, Y, ~, AUC] = perfcurve(labels, scores, 1);
    
    % Plot ROC Curve
    plot(X, Y, 'b-', 'LineWidth', 2);
    hold on;
    plot([0 1], [0 1], 'k--'); % Diagonal reference line
    xlabel('False Positive Rate (FPR)');
    ylabel('True Positive Rate (TPR)');
    title(sprintf('ROC Curve (AUC = %.3f)', AUC));
    grid on;
    hold off;
    
    fprintf('AUC: %.3f\n', AUC);
end


function avg_vec = average_adjusted_vectors(vec1, vec2)
    % Averages two sorted vectors after adjusting them to the same length.
    % vec1, vec2: Input sorted vectors
    % avg_vec: Output vector that is the average of the adjusted vectors
    
    % Determine lengths of input vectors
    len1 = numel(vec1);
    len2 = numel(vec2);

    if len1 > len2
        % Downsample vec1 to match vec2's length
        vec1_adj = interp1(linspace(1, len1, len1), vec1, linspace(1, len1, len2));
        vec2_adj = vec2;
    elseif len2 > len1
        % Downsample vec2 to match vec1's length
        vec2_adj = interp1(linspace(1, len2, len2), vec2, linspace(1, len2, len1));
        vec1_adj = vec1;
    else
        % If already same length, just use them
        vec1_adj = vec1;
        vec2_adj = vec2;
    end

    % Make sure both vectors are columns
    vec1_adj = vec1_adj(:);
    vec2_adj = vec2_adj(:);

    % Compute the average vector
    avg_vec = (vec1_adj + vec2_adj) / 2;
end