%%% Moving Window of Registrated Data - LOOCV
% Data created from register_all_data
% Each SOA trained & tested seperately

% SET PROJECT DIRECTORY
project_root = 'D:\Backup\Itiel_BM\BM_article_project';
cd(project_root)

% Add all subfolders to path so MATLAB can find your functions
addpath(genpath(project_root))

% hyper parameters
TEST_FRAMES = 1:80; % frames for testing the model
N_FRAMES = 40;

BASELINE_FRAMES = 5:20;
REGULIZER = 'ridge'; % L1 = 'lasso', L2 = 'ridge'
LAMBDA = 1;
NUM_TRAIN_FRAMES = 5;
FIRST_FRAME = 21;

% testing options
TEST_ERROR_DATA = 1;
SHUFFLE_LABELS = 0;

% useful variables
conds_labels = {'retino', 'noMask', '150', '100', '75', '50'};
all_frames = 1:256;
num_test_frames = size(TEST_FRAMES, 2);
num_bl_frames = size(BASELINE_FRAMES, 2);
num_SOA = size(conds_labels, 2);


%% Main loop - training and testing

% Change monkey
monkey = 'frodo';

% Load data
switch monkey
    case 'frodo'
        load('data\normalized2_regist\frodo_early_sessions\all_data.mat')
        load('data\normalized2_regist\frodo_early_sessions\error_data.mat')
        roiHor = importdata('data\normalized2_regist\frodo_early_sessions\roiHor.mat');
        roiVer = importdata('data\normalized2_regist\frodo_early_sessions\roiVer.mat');
        V1_pix = importdata('data\normalized2_regist\\frodo_early_sessions\V1_pix.mat');
    case 'tolkin_early'
        load('data\normalized2_regist\tolkin_early_sessions\all_data.mat')
        load('data\normalized2_regist\tolkin_early_sessions\error_data.mat')
        roiHor = importdata('data\normalized2_regist\tolkin_early_sessions\roiHor.mat');
        roiVer = importdata('data\normalized2_regist\tolkin_early_sessions\roiVer.mat');
        V1_pix = importdata('data\normalized2_regist\tolkin_early_sessions\V1_pix.mat');
    case 'tolkin_late'
        load('data\normalized2_regist\tolkin_late_sessions\all_data_456.mat')
        load('data\normalized2_regist\tolkin_late_sessions\error_data.mat')
        roiHor = importdata('data\normalized2_regist\tolkin_late_sessions\roiHor.mat');
        roiVer = importdata('data\normalized2_regist\tolkin_late_sessions\roiVer.mat');
        V1_pix = importdata('data\normalized2_regist\tolkin_late_sessions\V1_pix.mat');
end


% preallocate
weightsTC = nan(num_SOA, size(V1_pix,1), N_FRAMES);

% weights assign
[~,Iver,~] = intersect(V1_pix, roiVer);
[~,Ihor,~] = intersect(V1_pix, roiHor);
verOnly = setdiff(Iver, Ihor);
horOnly = setdiff(Ihor, Iver);
sharedRegion = intersect(Ihor, Iver);

accuracyCell = cell(num_SOA, 1);
errAccuracyCell = cell(num_SOA, 1);

disp('Initiating model training...')
fprintf('Number of ranges: %d\n', N_FRAMES)

% split across cons to train/test sets
for SOA_inx=3:6
    % allocate data
    vertical = all_data{SOA_inx, 1};
    horizontal = all_data{SOA_inx, 2};

    % balance data
    rng(45)
    vertical = vertical(randperm(size(vertical, 1)),:,:);
    horizontal = horizontal(randperm(size(horizontal, 1)),:,:);

    balanced_n = min(size(vertical, 1), size(horizontal, 1));
    vertical = vertical(1:balanced_n,:,:);
    horizontal = horizontal(1:balanced_n,:,:);
    
    % label shuffling
    if SHUFFLE_LABELS
        pool = [vertical;horizontal];
        pool = pool(randperm(size(pool, 1)), :, :);
        vertical = pool(1:balanced_n, :, :);
        horizontal = pool(balanced_n+1:end, :, :);
    end

    weights = nan(balanced_n, size(V1_pix,1), N_FRAMES);
    test_accuracy = nan(balanced_n, N_FRAMES);

    % Arange error data
    error_set = [error_data{SOA_inx, 1}; error_data{SOA_inx, 2}];
    error_label = [ones(1, size(error_data{SOA_inx, 1}, 1)), -ones(1, size(error_data{SOA_inx, 2}, 1))];
    n_err = size(error_data{SOA_inx,1}, 1) + size(error_data{SOA_inx,2}, 1);
    err_accuracy = nan(size(error_label, 2), N_FRAMES);
    
    fprintf('\n')
    fprintf('Trainig of %s session...\n', conds_labels{SOA_inx})
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

        num_trials_in(SOA_inx) = numel(in_inx);
        num_trials_out(SOA_inx) = numel(out_inx);

        for window_inx=1:N_FRAMES
            train_frames = FIRST_FRAME-1+window_inx : FIRST_FRAME+NUM_TRAIN_FRAMES+window_inx-2;

            % assign data to trainable matrices
            [X_in, Y_in] = ID_assignXY(train_vertical(:,train_frames,:),...
                                       train_horizontal(:,train_frames,:));
        
            % Train SVM model
            SVM_Model = fitclinear(X_in, Y_in, 'Learner', 'svm', ...
                'Regularization', REGULIZER, 'Lambda', LAMBDA);

            % assign weights & bias
            weights(perm_inx, :, window_inx) = SVM_Model.Beta;

            % Arange test-set data
            test_set = [test_vertical; test_horizontal];
            test_label = [ones(1, size(test_vertical, 1)), -ones(1, size(test_horizontal, 1))];
            accuracy = nan(size(test_label, 2), 1);
            
            % Predict test data, trial by trial, and measure accuracy
            for single_trial=1:(size(test_label, 2))
                X_out = squeeze(test_set(single_trial,train_frames,:)); 
                Y_pred = predict(SVM_Model, X_out);
                Y_pred = Y_pred == test_label(single_trial);
                accuracy(single_trial) = mean(Y_pred);
            end
            test_accuracy(perm_inx, window_inx) = nanmean(accuracy);

            if TEST_ERROR_DATA
                accuracy = nan(size(error_label, 2), 1);
        
                % Predict test data, trial by trial, and measure accuracy
                for single_trial=1:(size(error_label, 2))
                    X_out = squeeze(error_set(single_trial,train_frames,:)); 
                    Y_pred = predict(SVM_Model, X_out);
                    Y_pred = Y_pred == error_label(single_trial);
                    accuracy(single_trial) = mean(Y_pred);
                end
                err_accuracy(:, perm_inx, window_inx) = accuracy;
            end
        end
        weightsTC(SOA_inx,:,:) = squeeze(nanmean(weights, 1));
    end
    errAccuracyCell{SOA_inx} = squeeze(mean(err_accuracy, 2));
    accuracyCell{SOA_inx} = test_accuracy;
end


%% plot shadedErrorBar

colors = {[102 0 102]./255; [132 186 91]./255; [114 147 203]./255;
          [255 150 0]./255; [171 104 87]./255; [255 150 193]./255};
time = -40:10:350;

figure(); hold on
for soa_inx=[3 4 5 6]
    shadedErrorBar(time, nanmean(accuracyCell{soa_inx}, 1), nanstd(accuracyCell{soa_inx}, 1)/size(accuracyCell{soa_inx}, 1),...
        {'color', colors{soa_inx}, 'linewidth', 1},1);
end
ylim([0.35 1]) 
yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')
legend({'', '150', '', '100', '', '75', '', '50'},...
    'Location', 'northwest')
xlabel('Time (ms)'); ylabel('Accuracy')
title(sprintf('Grand Analysis of Models, %g-frames Moving Window', NUM_TRAIN_FRAMES))

