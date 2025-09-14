%% Figure 4 - Retinotopic
% in this script I produce all of figure 4

% SET PROJECT DIRECTORY
project_root = 'D:\Backup\Itiel_BM\BM_article_project';
cd(project_root)

% Add all subfolders to path so MATLAB can find your functions
addpath(genpath(project_root))

% Load data
ANALYSIS_WINDOW = 34:38;
SOA_inx = 1;

% Change monkey (frodo, tolkin)
monkey = 'tolkin';

% Load data
switch monkey
    case 'frodo'
        % Frodo early sessions
        all_data = importdata('data\normalized2_regist\frodo_early_sessions\all_data.mat');
        roiHor = importdata('data\normalized2_regist\frodo_early_sessions\ellipseHor.mat');
        roiVer = importdata('data\normalized2_regist\frodo_early_sessions\ellipseVer.mat');
        V1_pix = importdata('data\normalized2_regist\\frodo_early_sessions\V1_pix.mat');
        egSession = '23/09';
        egIdx = 2;
        sessionsIdx = {1:7, 1:10;
                   8:15, 11:16;
                   16:20, 17:25;
                   21:31, 26:33};
    case 'tolkin'
        all_data = importdata('data\normalized2_regist\tolkin_late_sessions\all_data.mat');
        roiHor = importdata('data\normalized2_regist\tolkin_late_sessions\ellipseHor.mat');
        roiVer = importdata('data\normalized2_regist\tolkin_late_sessions\ellipseVer.mat');
        V1_pix = importdata('data\normalized2_regist\tolkin_late_sessions\V1_pix.mat');
        egSession = '04/07';
        egIdx = 1;
        sessionsIdx = {1:8, 1:8;
                   9:17, 9:13;
                   18:28, 14:22;
                   29:38, 23:33};
end

% weights assign
[~,Iver,~] = intersect(V1_pix, roiVer);
[~,Ihor,~] = intersect(V1_pix, roiHor);
IverOnly = setdiff(Iver, Ihor);
IhorOnly = setdiff(Ihor, Iver);
IsharedRegion = intersect(Ihor, Iver);

% allocate data
vertical = all_data{SOA_inx, 1};
eg_vertical = vertical(sessionsIdx{egIdx,1}, :, :);
horizontal = all_data{SOA_inx, 2};
eg_horizontal = horizontal(sessionsIdx{egIdx,2}, :, :);
full_time_course = ((1:80)-27) * 10;

%% Show example session Maps for Vertical and Horizontal Targets
% Adjust as desire
time_point = full_time_course(ANALYSIS_WINDOW);
S = 1;

% Define the minimum and maximum values for clipping the colormap
minValue = 0; % Adjust as needed based on your data
maxValue = 7; % Adjust as needed based on your data

% Create contours
contourV1 = zeros(100, 100);
contourV1(V1_pix) = 1; contourV1 = contourV1';

contourVer = zeros(100, 100);
contourVer(roiVer) = 1; contourVer = contourVer';

contourHor = zeros(100, 100);
contourHor(roiHor) = 1; contourHor = contourHor';

% Create vertical map
verMap = zeros(10000, 1);
verMap(V1_pix) = squeeze(mean(eg_vertical(:,ANALYSIS_WINDOW,:), [1 2]));
verMap = mfilt2(verMap, 100, 100, S, 'lm');
verMap(setdiff(1:10000, V1_pix)) = 1000;
verMap = reshape(verMap, 100, 100)';

% Create horizontal map
horMap = zeros(10000, 1);
horMap(V1_pix) = squeeze(mean(eg_horizontal(:,ANALYSIS_WINDOW,:), [1 2]));
horMap = mfilt2(horMap, 100, 100, S, 'lm');
horMap(setdiff(1:10000, V1_pix)) = 1000;
horMap = reshape(horMap, 100, 100)';


figure;
imagesc(horMap, [minValue maxValue]); colormap mapgeog
hold on
contour(contourV1, 'k');
contour(contourHor, 'cyan');
title(sprintf('Horizontal Retino, time: %g:%g ms', time_point(1), time_point(end)))
text(10, 5, sprintf('frames: %d:%d', ANALYSIS_WINDOW(1), ANALYSIS_WINDOW(end)))
text(10, 10, sprintf('smoothing: S=%g', S))
text(10,15, sprintf('cliping: [%g, %g]', minValue, maxValue))
hold off

figure;
imagesc(verMap, [minValue maxValue]); colormap mapgeog
hold on
contour(contourV1, 'k');
contour(contourVer, 'cyan');
title(sprintf('Vertical Retino, time: %g:%g ms', time_point(1), time_point(end)))
text(10, 5, sprintf('frames: %d:%d', ANALYSIS_WINDOW(1), ANALYSIS_WINDOW(end)))
text(10, 10, sprintf('smoothing: S=%g', S))
text(10,15, sprintf('cliping: [%g, %g]', minValue, maxValue))
hold off

%% Show example session Signal TC

plt_frames = 23:80;
plt_time = full_time_course(plt_frames);
meanTC = nanmean([nanmean(eg_vertical(:,plt_frames,Iver), 3);...
            nanmean(horizontal(:,plt_frames,Ihor), 3)]);
semTC = nanstd([nanmean(eg_vertical(:,plt_frames,Iver), 3);...
                (nanmean(horizontal(:,plt_frames,Ihor), 3))]) / sqrt(size(eg_vertical, 1)+size(horizontal, 1));

figure;
shadedErrorBar(plt_time, meanTC, semTC,...
    {'color', [102 0 102]./255, 'linewidth', 1},1);
title('Signal TC of Both Targets in Retino Sessions')
xlabel('time (ms)')
ylabel('Zscore')


%% Train Model with e.g. Session
% Model setting
% hyper parameters
TEST_FRAMES = 1:80; % frames for testing the model
TRAIN_FRAMES = ANALYSIS_WINDOW;
BASELINE_FRAMES = 5:20;
PRE_STIM_FRAMES = 23:27;

REGULIZER = 'ridge'; % L1 = 'lasso', L2 = 'ridge'
LAMBDA = 1;
NUM_TRAIN_FRAMES = 5;
FIRST_FRAME = 21;

% useful variables
SOA_inx = 1;
all_frames = 1:256;
num_test_frames = size(TEST_FRAMES, 2);
num_train_frames = size(TRAIN_FRAMES, 2);
num_bsln_frames = size(BASELINE_FRAMES, 2);


%% Model training and testing:
% Real data training
[eg_acuracy, eg_weights] = trainAndTest(eg_vertical, eg_horizontal, V1_pix, TRAIN_FRAMES, TEST_FRAMES,...
    REGULIZER, LAMBDA, 0);
% Shuffled data training
[shfld_accuracy, shfld_weights] = trainAndTest(eg_vertical, eg_horizontal, V1_pix, TRAIN_FRAMES, TEST_FRAMES,...
    REGULIZER, LAMBDA, 1);

%% Visualise Outputs of Model
% Plot Accuracy Time-Course of e.g. Session

% plot shaded errorbars
plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
patch_frame = full_time_course(TRAIN_FRAMES(1));

figure(); hold on
rectangle('Position', [patch_frame 0.5 40 0.5], 'FaceColor', [204 155 253 10]./255,'EdgeColor','none')
shadedErrorBar(plt_time, nanmean(eg_acuracy(:,plt_frames), 1), nanstd(eg_acuracy(:,plt_frames), 1)/size(eg_acuracy(:,plt_frames), 1),...
    {'color', [102 0 102]./255, 'linewidth', 1}, 0);
shadedErrorBar(plt_time, nanmean(shfld_accuracy(:,plt_frames), 1), nanstd(shfld_accuracy(:,plt_frames), 1)/size(shfld_accuracy(:,plt_frames), 1),...
    {'color', [150 0 150]./255, 'linewidth', 1}, 0);

% grid on
ylim([0.2 1]) 
yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')
legend({'', 'real data', '', 'shuffled data'},...
    'Location', 'northwest')
xlabel('Time (ms)'); ylabel('Accuracy')
title('Accuracy Time-Course of e.g. session, 5-frames')

% Values and pValues
accuracyBsln = nanmean(eg_acuracy(:,PRE_STIM_FRAMES), 'all')
semBsln = std(nanmean(eg_acuracy(:,PRE_STIM_FRAMES), 2)) / sqrt(size(eg_acuracy, 1))
sgnrnkBsln = signrank(nanmean(eg_acuracy(:,PRE_STIM_FRAMES), 2), 0.5)
rnksmBsln = ranksum(nanmean(eg_acuracy(:,PRE_STIM_FRAMES), 2), nanmean(shfld_accuracy(:,PRE_STIM_FRAMES), 2))

accuracyWindow = nanmean(eg_acuracy(:,TRAIN_FRAMES), 'all')
semWindow = std(nanmean(eg_acuracy(:,TRAIN_FRAMES), 2)) / sqrt(size(eg_acuracy, 1))
sgnrnkWindow = signrank(nanmean(eg_acuracy(:,TRAIN_FRAMES), 2), 0.5)
rnksmWindow = ranksum(nanmean(eg_acuracy(:,TRAIN_FRAMES), 2), nanmean(shfld_accuracy(:,TRAIN_FRAMES), 2))

%% single map visualization & weights histogram
% Adjust as desire
time_point = full_time_course(TRAIN_FRAMES);
S = 1;
FIRST_FRAME = 21;
NUM_TRAIN_FRAMES = 5;

weightsMap = zeros(10000, 1);
weightsMap(V1_pix) = squeeze(mean(eg_weights, 1));

% Define the minimum and maximum values for clipping the colormap
minValue = -0.001; % Adjust as needed based on your data
maxValue = 0.001; % Adjust as needed based on your data

% create contour matrices
horCont = false(10000,1); horCont(roiHor) = 1;
horCont = reshape(horCont, 100, 100)';
verCont = false(10000,1); verCont(roiVer) = 1;
verCont = reshape(verCont, 100, 100)';
contourV1 = false(100, 100); contourV1(V1_pix) = 1;
contourV1 = contourV1';

% select ROI pixels
horPixelValues = mean(eg_weights(:, IhorOnly));
verPixelValues = mean(eg_weights(:, IverOnly));
sharedPixelValues =  mean(eg_weights(:, IsharedRegion));

% filter and prepare to plot
weightsMap = mfilt2(weightsMap, 100, 100, S, 'lm');
weightsMap(setdiff(1:10000, V1_pix)) = 1000;
weightsMap = reshape(weightsMap, 100, 100);

% plot vsd map
figure;
fg = imagesc(reshape(weightsMap, 100, 100)', [minValue maxValue]); colormap mapgeog
hold on
contour(horCont, 'c')
contour(verCont, 'c')
contour(contourV1, 'k');
title(sprintf('Retino, time: %g:%g ms', time_point(1), time_point(end)))
text(10, 5, sprintf('frames: %d:%d', TRAIN_FRAMES(1), TRAIN_FRAMES(end)))
text(10, 10, sprintf('smoothing: S=%g', S))
text(10,15, sprintf('cliping: [%g, %g]', minValue, maxValue))

% Values in ROIs
sz = [3 5];
varTypes = ["string","double","double","double","double"];
varNames = ["ROI" "value" "ROI std" "sem" "pValue"];
egStats = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames', varNames);
egStats.ROI = ["horizontal"; "vertical"; "shared"];

egStats{1,2} = mean(eg_weights(:,IhorOnly), 'all');
egStats{2,2} = mean(eg_weights(:,IverOnly), 'all');
egStats{3,2} = mean(eg_weights(:,IsharedRegion), 'all');

egStats{1,3} = std(mean(eg_weights(:,IhorOnly), 1));
egStats{2,3} = std(mean(eg_weights(:,IverOnly), 1));
egStats{3,3} = std(mean(eg_weights(:,IsharedRegion), 1));

egStats{1,4} = std(mean(eg_weights(:,IhorOnly), 2)) / sqrt(size(eg_weights, 1));
egStats{2,4} = std(mean(eg_weights(:,IverOnly), 2)) / sqrt(size(eg_weights, 1));
egStats{3,4} = std(mean(eg_weights(:,IsharedRegion), 2)) / sqrt(size(eg_weights, 1));

egStats{1,5} = ranksum(nanmean(eg_weights(:,IhorOnly), 2), nanmean(shfld_weights(:,IhorOnly), 2));
egStats{2,5} = ranksum(nanmean(eg_weights(:,IverOnly), 2), nanmean(shfld_weights(:,IverOnly), 2));
egStats{3,5} = ranksum(nanmean(eg_weights(:,IsharedRegion), 2), nanmean(shfld_weights(:,IsharedRegion), 2));

disp(egStats)
rnksmVerHor = ranksum(nanmean(eg_weights(:,IverOnly), 2), nanmean(eg_weights(:,IhorOnly), 2))

% plot histogram
figure;
h1 = histogram(horPixelValues, 'Normalization', 'probability');
hold on
h2 = histogram(verPixelValues, 'Normalization', 'probability');
h1.BinWidth = 0.00015;
h2.BinWidth = 0.00015;
sigstar([egStats{1,2}, egStats{2,2}], rnksmVerHor)
ylim([0 0.5]); xlim([-4e-3 4e-3])
title(sprintf('Weights Distribution in Retino Weight map, Horizontal & Vertical targets'))
legend({'Horizontal weights' 'Vertical weights'})

%% Train and Test Model with All Trials
% Real data training
[test_accuracy, weights] = trainAndTest(vertical, horizontal, V1_pix, TRAIN_FRAMES, TEST_FRAMES,...
    REGULIZER, LAMBDA, 0);
% Shuffled data training
[shfld_accuracy, shfld_weights] = trainAndTest(vertical, horizontal, V1_pix, TRAIN_FRAMES, TEST_FRAMES,...
    REGULIZER, LAMBDA, 1);
% Alocate data
verPixelValues = mean(weights(:, IverOnly), 1);
horPixelValues = mean(weights(:, IhorOnly), 1);
sharedPixelValues = mean(weights(:, IsharedRegion), 1);

% case Tolkin
if monkey == "tolkin"
    all_data_early = importdata('data\normalized2_regist\tolkin_early_sessions\all_data.mat');
    roiHor_early = importdata('data\normalized2_regist\tolkin_early_sessions\ellipseHor.mat');
    roiVer_early = importdata('data\normalized2_regist\tolkin_early_sessions\ellipseVer.mat');   
    V1_pix_early = importdata('data\normalized2_regist\tolkin_early_sessions\V1_pix.mat');

    [~,Iver_early,~] = intersect(V1_pix_early, roiVer_early);
    [~,Ihor_early,~] = intersect(V1_pix_early, roiHor_early);
    IverOnly_early = setdiff(Iver_early, Ihor_early);
    IhorOnly_early = setdiff(Ihor_early, Iver_early);
    IsharedRegion_early = intersect(Ihor_early, Iver_early);

    vertical_early = all_data_early{SOA_inx, 1};
    horizontal_early = all_data_early{SOA_inx, 2};

    % train and test model
    [test_accuracy_early, weights_early] = trainAndTest(vertical_early, horizontal_early, V1_pix_early, TRAIN_FRAMES, TEST_FRAMES,...
    REGULIZER, LAMBDA, 0);
    [shfld_accuracy_early, ~] = trainAndTest(vertical_early, horizontal_early, V1_pix_early, TRAIN_FRAMES, TEST_FRAMES,...
    REGULIZER, LAMBDA, 1);
    test_accuracy = [test_accuracy; test_accuracy_early];
    shfld_accuracy = [shfld_accuracy; shfld_accuracy_early];

    % Alocate data
    verPixelValues_early = mean(weights_early(:, IverOnly_early), 1);
    horPixelValues_early = mean(weights_early(:, IhorOnly_early), 1);
    sharedPixelValues_early = mean(weights_early(:, IsharedRegion_early), 1);
end

%% Visualise output
% Time-Course Accuracy
% plot shaded errorbars
full_time_course = ((1:80)-27) * 10;
plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
patch_frame = full_time_course(TRAIN_FRAMES(1));

figure(); hold on
rectangle('Position', [patch_frame 0.5 40 0.5], 'FaceColor', [204 155 253 10]./255,'EdgeColor','none')
shadedErrorBar(plt_time, nanmean(test_accuracy(:,plt_frames), 1), nanstd(test_accuracy(:,plt_frames), 1)/size(test_accuracy(:,plt_frames), 1),...
    {'color', [102 0 102]./255, 'linewidth', 1},0);
shadedErrorBar(plt_time, nanmean(shfld_accuracy(:,plt_frames), 1), nanstd(shfld_accuracy(:,plt_frames), 1)/size(shfld_accuracy(:,plt_frames), 1),...
    {'color', [150 0 150]./255, 'linewidth', 1},0);

% grid on
ylim([0.2 1]) 
yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')
legend({'', 'real data', '', 'shuffled data'},...
    'Location', 'northwest')
xlabel('Time (ms)'); ylabel('Accuracy')
title('Accuracy Time-Course GA, 5-frames')

accuracyBsln = nanmean(test_accuracy(:,PRE_STIM_FRAMES), 'all')
semBsln = std(nanmean(test_accuracy(:,PRE_STIM_FRAMES), 2)) / sqrt(size(test_accuracy, 1))
sgnrnkBsln = signrank(nanmean(test_accuracy(:,PRE_STIM_FRAMES), 2), 0.5)
rnksmBsln = ranksum(nanmean(test_accuracy(:,PRE_STIM_FRAMES), 2), nanmean(shfld_accuracy(:,PRE_STIM_FRAMES), 2))

accuracyWindow = nanmean(test_accuracy(:,TRAIN_FRAMES), 'all')
semWindow = std(nanmean(test_accuracy(:,TRAIN_FRAMES), 2)) / sqrt(size(test_accuracy, 1))
sgnrnkWindow = signrank(nanmean(test_accuracy(:,TRAIN_FRAMES), 2), 0.5)
rnksmWindow = ranksum(nanmean(test_accuracy(:,TRAIN_FRAMES), 2), nanmean(shfld_accuracy(:,TRAIN_FRAMES), 2))


%% Bar Plot of Weights
% Arrange bar data
X = categorical({'Vertical target','Horizontal target','Shared region'});
X = reordercats(X, {'Vertical target','Horizontal target','Shared region'});
Y = [mean(verPixelValues), mean(horPixelValues), mean(sharedPixelValues)];
sem_vec = [std(verPixelValues)/sqrt(numel(verPixelValues))
           std(horPixelValues)/sqrt(numel(horPixelValues))
           std(sharedPixelValues)/sqrt(numel(sharedPixelValues))];

% Get pValues
p_ver_hor = ranksum(verPixelValues, horPixelValues);
p_hor_shr = ranksum(horPixelValues, sharedPixelValues);
p_ver_shr = ranksum(verPixelValues, sharedPixelValues);

if monkey == "tolkin"
    Y_early = [mean(verPixelValues_early), mean(horPixelValues_early), mean(sharedPixelValues_early)];
    sem_vec_early = [std(verPixelValues_early)/sqrt(numel(verPixelValues_early))
                       std(horPixelValues_early)/sqrt(numel(horPixelValues_early))
                       std(sharedPixelValues_early)/sqrt(numel(sharedPixelValues_early))];
    Y = mean([Y; Y_early], 1);
    sem_vec = mean([sem_vec, sem_vec_early], 2);

    p_ver_hor_early = ranksum(verPixelValues_early, horPixelValues_early);
    p_hor_shr_early = ranksum(horPixelValues_early, sharedPixelValues_early);
    p_ver_shr_early = ranksum(verPixelValues_early, sharedPixelValues_early);

    p_ver_hor = max(p_ver_hor, p_ver_hor_early);
    p_hor_shr = max(p_hor_shr, p_hor_shr_early);
    p_ver_shr = max(p_ver_shr, p_ver_shr_early);
end

% Colors
colors = [
    [80 0 80]/255;
    [120 20 120]/255;
    [160 60 160]/255;
];

% Plot figure
figure; hold on
for k = 1:size(Y,2)
    bar(k, Y(k), 'FaceColor', colors(k, :));
end
errorbar(1:3,Y,sem_vec,'k','linestyle','none');
set(gca, 'XTick', 1:3);
set(gca,'XTickLabel',X)

% Plot the stars
groups = {[1 2] [1 3] [2 3]};
sigstar(groups,[(p_ver_hor), (p_ver_shr), (p_hor_shr)]);

ylim([-0.001 0.001])
ylabel('Mean Value (a.u)')
title(sprintf('Weights Mean Value, Horizontal & Vertical targets, %d ms', time_point))
hold off

% Stats
sz = [3 5];
varTypes = ["string","double","double","double","double"];
varNames = ["ROI" "value" "ROI std" "sem" "pValue"];
gaStats = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames', varNames);
gaStats.ROI = ["vertical"; "horizontal"; "shared"];

gaStats{:,2} = Y';
gaStats{:,4} = sem_vec;

gaStats{1,5} = ranksum(nanmean(weights(:,IhorOnly), 1), nanmean(shfld_weights(:,IhorOnly), 1));
gaStats{2,5} = ranksum(nanmean(weights(:,IverOnly), 1), nanmean(shfld_weights(:,IverOnly), 1));
gaStats{3,5} = ranksum(nanmean(weights(:,IsharedRegion), 1), nanmean(shfld_weights(:,IsharedRegion), 1));

disp(gaStats)

%% SVM training and testing function
function [test_accuracy, weights] = trainAndTest(vertical,horizontal,V1_pix,TRAIN_FRAMES,TEST_FRAMES,...
    REGULIZER, LAMBDA,shuffle)
rng(42)
% balance data
vertical = vertical(randperm(size(vertical, 1)),:,:);
horizontal = horizontal(randperm(size(horizontal, 1)),:,:);

balanced_n = min(size(vertical, 1), size(horizontal, 1));
vertical = vertical(1:balanced_n,:,:);
horizontal = horizontal(1:balanced_n,:,:);

% label shuffling
if shuffle
    pool = [vertical;horizontal];
    pool = pool(randperm(size(pool, 1)), :, :);
    vertical = pool(1:balanced_n, :, :);
    horizontal = pool(balanced_n+1:end, :, :);
end

num_test_frames = numel(TEST_FRAMES);
weights = nan(balanced_n, size(V1_pix,1));
test_accuracy = nan(balanced_n, num_test_frames);

reverseStr = '';
for perm_inx=1:balanced_n
    msg = sprintf('Processed %d/%d permutaions', perm_inx, balanced_n);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    % preallocate train/test data
    out_inx = perm_inx;
    in_inx = setdiff(1:balanced_n, out_inx);
    % split vertical data
    train_vertical = vertical(in_inx,:,:);
    test_vertical = vertical(out_inx,:,:);

    % split horizontal data
    train_horizontal = horizontal(in_inx,:,:);
    test_horizontal = horizontal(out_inx,:,:);

    num_trials_in = numel(in_inx);
    num_trials_out = numel(out_inx);



    % assign data to trainable matrices
    [X_in, Y_in] = ID_assignXY(train_vertical(:,TRAIN_FRAMES,:),...
        train_horizontal(:,TRAIN_FRAMES,:));

    % Train SVM model
    SVM_Model = fitclinear(X_in, Y_in, 'Learner', 'svm', ...
        'Regularization', REGULIZER, 'Lambda', LAMBDA);

    % assign weights & bias
    weights(perm_inx, :) = SVM_Model.Beta;

    % Arange test-set data
    test_set = [test_vertical; test_horizontal];
    test_label = [ones(1, size(test_vertical, 1)), -ones(1, size(test_horizontal, 1))];
    accuracy = nan(size(test_label, 2), num_test_frames);

    % Predict test data, trial by trial, and measure accuracy
    for single_trial=1:(size(test_label, 2))
        X_out = squeeze(test_set(single_trial, TEST_FRAMES, :));
        Y_pred = predict(SVM_Model, X_out);
        accuracy(single_trial, :) = Y_pred == test_label(single_trial);
    end
    test_accuracy(perm_inx, :) = nanmean(accuracy, 1);
end
disp(' ')
end

