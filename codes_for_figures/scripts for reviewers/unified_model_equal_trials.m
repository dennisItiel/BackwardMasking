%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Answer for rev2 - model of all SOAs together                    %%%
%%%           equal no. of trials for each SOA                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moving Window of Registrated Data - LOOCV
% Data created from register_all_data
% All SOAs combined together for training, but are tested separately

% hyper parameters
TEST_FRAMES = 1:80; % frames for testing the model
N_FRAMES = 40;

BASELINE_FRAMES = 5:20;
REGULIZER = 'ridge'; % L1 = 'lasso', L2 = 'ridge'
LAMBDA = 1;
NUM_TRAIN_FRAMES = 5;
FIRST_FRAME = 21;

% plotting options
SHUFFLE_LABELS = 0;

% useful variables
conds_labels = {'retino', 'noMask', '150', '100', '75', '50'};
all_frames = 1:256;
num_test_frames = size(TEST_FRAMES, 2);
num_bl_frames = size(BASELINE_FRAMES, 2);
num_SOA = size(conds_labels, 2);

%% Main loop - training and testing

% % Frodo early sessions
% load('E:\Itiel\BM_ID_analysis\normalized2_regist\frodo_early_sessions\all_data.mat')
% roiHor = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\frodo_early_sessions\roiHor.mat');
% roiVer = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\frodo_early_sessions\roiVer.mat');
% V1_pix = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\\frodo_early_sessions\V1_pix.mat');

% % Tolkin early sessions
% load('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_early_sessions\all_data.mat')
% roiHor = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_early_sessions\roiHor.mat');
% roiVer = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_early_sessions\roiVer.mat');
% V1_pix = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_early_sessions\V1_pix.mat');

% Tolkin late sessions
load('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_late_sessions\all_data_456.mat')
roiHor = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_late_sessions\roiHor.mat');
roiVer = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_late_sessions\roiVer.mat');
V1_pix = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\tolkin_late_sessions\V1_pix.mat');

disp('Initiating model training...')
fprintf('Number of ranges: %d\n', N_FRAMES)


SOAs_to_use = 3:6;
rng(45)
h = @(x) size(x,1);
min_of_all_trials = min(min(cellfun(h, all_data(3:6,:))));

remaining_data = cell(6,2);
% balance trials
for SOA_inx = SOAs_to_use
    vertical = all_data{SOA_inx, 1};
    horizontal = all_data{SOA_inx, 2};

    % Balance data
    vertical = vertical(randperm(size(vertical, 1)),:,:);
    horizontal = horizontal(randperm(size(horizontal, 1)),:,:);
    balanced_n = min(size(vertical, 1), size(horizontal, 1));

    % Put in structs
    all_data{SOA_inx, 1} = vertical(1:min_of_all_trials,:,:);
    all_data{SOA_inx, 2} = horizontal(1:min_of_all_trials,:,:);
    remaining_data{SOA_inx, 1} = vertical(min_of_all_trials+1:balanced_n,:,:);
    remaining_data{SOA_inx, 2} = horizontal(min_of_all_trials+1:balanced_n,:,:);
end

% split across cons to train/test sets
weights = nan(10, size(V1_pix,1), N_FRAMES);
test_accuracy = nan(10, N_FRAMES, 4);


fprintf('\n')
fprintf('Trainig of all sessions...\n')
reverseStr = '';
for fold=1:10
    % preallocate train/test data
    all_train_vertcial = [];
    all_train_horizontal = [];
    all_test_vertical = remaining_data(:, 1);
    all_test_horizontal = remaining_data(:, 2);

    for SOA_inx = SOAs_to_use
        vertical = all_data{SOA_inx, 1};
        horizontal = all_data{SOA_inx, 2};

        m = size(vertical, 1);
        cv = cvpartition(m, 'KFold', 10);
        train_idx = training(cv, fold);  % logical index vector for training
        test_idx = test(cv, fold);       % logical index vector for testing

        all_train_vertcial = cat(1, all_train_vertcial, vertical(train_idx,:,:));
        all_train_horizontal = cat(1, all_train_horizontal, horizontal(train_idx,:,:));
        all_test_vertical{SOA_inx} = cat(1, all_test_vertical{SOA_inx}, vertical(test_idx,:,:));
        all_test_horizontal{SOA_inx} = cat(1, all_test_horizontal{SOA_inx}, horizontal(test_idx,:,:));
    end

    rng(42)
    if SHUFFLE_LABELS
        num_all_train_trials = size(all_train_vertcial, 1);
        idx_vector = [zeros(num_all_train_trials,1); ones(num_all_train_trials, 1)];
        idx_vector = logical(idx_vector(randperm(num_all_train_trials*2)));
        all_train_vertcial = cat(1, all_train_vertcial(idx_vector(1:num_all_train_trials),:,:),...
                                    all_train_horizontal(idx_vector(num_all_train_trials+1:end),:,:));
        all_train_horizontal = cat(1, all_train_vertcial(~idx_vector(1:num_all_train_trials),:,:),...
                                    all_train_horizontal(~idx_vector(num_all_train_trials+1:end),:,:));
    end

    for window_inx=1:N_FRAMES
        train_frames = FIRST_FRAME-1+window_inx : FIRST_FRAME+NUM_TRAIN_FRAMES+window_inx-2;

        % assign data to trainable matrices
        [X_in, Y_in] = ID_assignXY(all_train_vertcial(:,train_frames,:),...
                                   all_train_horizontal(:,train_frames,:));

        % Train SVM model
        SVM_Model = fitclinear(X_in, Y_in, 'Learner', 'svm', ...
            'Regularization', REGULIZER, 'Lambda', LAMBDA);

        % assign weights & bias
        weights(fold, :, window_inx) = SVM_Model.Beta;

        for SOA_inx = SOAs_to_use

            % Arange test-set data
            test_set = [all_test_vertical{SOA_inx}; all_test_horizontal{SOA_inx}];
            test_label = [ones(1, size(all_test_vertical{SOA_inx}, 1)), -ones(1, size(all_test_horizontal{SOA_inx}, 1))];
            accuracy = nan(size(test_label, 2), 1);
    
            % Predict test data, trial by trial, and measure accuracy
            for single_trial=1:(size(test_label, 2))
                X_out = squeeze(test_set(single_trial,train_frames,:));
                Y_pred = predict(SVM_Model, X_out);
                %                 Y_pred = sign(mean(Y_pred));
                Y_pred = Y_pred == test_label(single_trial);
                accuracy(single_trial) = mean(Y_pred);
            end
            test_accuracy(fold, window_inx, SOA_inx-2) = nanmean(accuracy);
        end
    end
   
    % Print message
    msg = sprintf('Processed %d/%d folds', fold, 10);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

%% plot shadedErrorBar

plt_time = ((23:62) - 27) * 10;

figure(); hold on

ylim([0.35 1]) %([0.3 0.7])
% h1 = shadedErrorBar(time, nanmean(accuracyCell{1}, 1), nanstd(accuracyCell{1}, 1)/size(accuracyCell{1}, 1),...
%     {'color', [102 0 102]./255, 'linewidth', 1},1);
% h2 = shadedErrorBar(time, nanmean(accuracyCell{2}, 1), nanstd(accuracyCell{2}, 1)/size(accuracyCell{2}, 1),...
%     {'color', [132 186 91]./255,  'linewidth', 1},1);
h3 = shadedErrorBar(plt_time, nanmean(test_accuracy(:,:,1), 1), nanstd(test_accuracy(:,:,1), 1)/sqrt(10),...
    {'color', [114 147 203]./255, 'linewidth', 1},1);
h4 = shadedErrorBar(plt_time, nanmean(test_accuracy(:,:,2), 1), nanstd(test_accuracy(:,:,2), 1)/sqrt(10),...
    {'color', [255 150 0]./255, 'linewidth', 1},1);
h5 = shadedErrorBar(plt_time, nanmean(test_accuracy(:,:,3), 1), nanstd(test_accuracy(:,:,3), 1)/sqrt(10),...
    {'color', [171 104 87]./255, 'linewidth', 1},1);
h6 = shadedErrorBar(plt_time, nanmean(test_accuracy(:,:,4), 1), nanstd(test_accuracy(:,:,4), 1)/sqrt(10),...
    {'color', [255 150 193]./255, 'linewidth', 1},1);

yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')

conds_labels = {'150', '100', '75', '50'}; 
legend([h3.mainLine h4.mainLine h5.mainLine h6.mainLine], conds_labels,...
    'Location', 'northwest')

xlabel('frames'); ylabel('accuracy')
title('Grand Analysis of Models, 5-frames Moving Window')

%% Combine short and long SOAs

longSoaAccuracy = cat(3, test_accuracy(:,:,1), test_accuracy(:,:,2));
shortSoaAccuracy = cat(3, test_accuracy(:,:,3), test_accuracy(:,:,4));

plt_time = ((23:62) - 27) * 10;

figure(); hold on

ylim([0.35 1]) %([0.3 0.7])
h1 = shadedErrorBar(plt_time, nanmean(longSoaAccuracy(:,:,1), 1), nanstd(longSoaAccuracy(:,:,1), 1)/sqrt(10),...
    {'color', [102 0 102]./255, 'linewidth', 1},1);
h2 = shadedErrorBar(plt_time, nanmean(shortSoaAccuracy(:,:,2), 1), nanstd(shortSoaAccuracy(:,:,2), 1)/sqrt(10),...
    {'color', [132 186 91]./255, 'linewidth', 1},1);

yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')

conds_labels = {'long soa', 'short soa'}; 
legend([h1.mainLine h2.mainLine], conds_labels,...
    'Location', 'northwest')

xlabel('frames'); ylabel('accuracy')
title('Grand Analysis of Models, 5-frames Moving Window')


%% weights map visualization & weights histogram - 2X1 grid

residue = importdata('E:\Itiel\BM_ID_analysis\normalized2_regist\frodo_early_sessions\residue.mat');

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

full_time_course = ((1:80)-27) * 10;
plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
weightsTC = squeeze(nanmean(weights, 1));
timePoints = {70, 180};
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

for gridRow = 1 % soa150/50
    for gridCol = 1:2  % early/late
        % Adjust as desire
        time_point = timePoints{gridCol};
        frameIdx = find(plt_time == time_point(1)):find(plt_time == time_point(end));
        map = zeros(10000, 1);
        % map(V1_pix) = squeeze(mean(all_data{SOA_inx, 2}(:,frame,:), [1 2]));
        map(V1_pix) = squeeze(mean(weightsTC(:,frameIdx), 2));


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
        subplot(1, 2, sbpltIdx);
        fg = imagesc(reshape(map, 100, 100)', [minValue maxValue]); colormap mapgeog
        hold on
        contour(horCont, 'Color', '#808080')
        contour(verCont, 'Color', '#808080')
        contour(contourV1, 'k');
        title(['time:', num2str(time_point)])
        if sbpltIdx==1
            text(10, 5, sprintf('frames: %d:%d', frameIdx + FIRST_FRAME-1, frameIdx + FIRST_FRAME + NUM_TRAIN_FRAMES -2))
            text(10, 10, sprintf('smoothing: S=%g', S))
            text(10,15, sprintf('cliping: [%g, %g]', minValue, maxValue))
        end
        hold off

        % plot histogram
        figure(101);
        subplot(1, 2, sbpltIdx);        hold on
        h1 = histogram(horValues, 'Normalization', 'probability');
        h2 = histogram(verValues, 'Normalization', 'probability');
        h1.BinWidth = 0.0005;
        h2.BinWidth = 0.0005;
        ylim([0 0.15]); xlim([-10e-3 10e-3])
        xline(mean(horValues), 'Color', '#0343DF')
        xline(mean(verValues), 'Color', '#8C000F')
        title(sprintf('Weights, Hor & Ver'))
        subtitle(ranksum(horValues, verValues))
        if sbpltIdx==1
            legend({'horOnly weights' 'verOnly weights'})
        end

        % plot ROC and AUC
        figure(102);
        subplot(1, 2, sbpltIdx);        hold on
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


% distr1 = weights_eary_cell{1,1}-weights_eary_cell{2,1};
% distr2 = weights_eary_cell{1,2}-weights_eary_cell{2,2};
distr1 = weights_late_cell{1,1};
distr2 = weights_late_cell{2,1};

figure; hold on
h1 = histogram(distr1, 'Normalization', 'probability');
h2 = histogram(distr2, 'Normalization', 'probability');
h1.BinWidth = 0.0005;
h2.BinWidth = 0.0005;
ylim([0 0.15]); xlim([-10e-3 10e-3])
xline(mean(distr1), 'Color', '#0343DF')
xline(mean(distr2), 'Color', '#8C000F')

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