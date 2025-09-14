%% get values of key points TABLE & Bars

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
        trueAccuracy = importdata('data\normalized2_regist\frodo_early_sessions\accuracyCell.mat');
        shfldAccuracy = importdata('data\normalized2_regist\frodo_early_sessions\accuracyCellShfld.mat');

    case 'tolkin'
        accuracy_early = importdata('data\normalized2_regist\tolkin_early_sessions\accuracyCell.mat');
        accuracy_early_shlfd = importdata('data\normalized2_regist\tolkin_early_sessions\accuracyCellShfld.mat');
        
        accuracy_late = importdata('data\normalized2_regist\tolkin_late_sessions\accuracyCell_456.mat');
        accuracy_late_shfld = importdata('data\normalized2_regist\tolkin_late_sessions\accuracyCellShfld_456.mat');
        
        for i=1:6
            trueAccuracy{i,1} = cat(1, accuracy_early{i,1}, accuracy_late{i,1});
            shfldAccuracy{i,1} = cat(1, accuracy_early_shlfd{i,1}, accuracy_late_shfld{i,1});
        end
end

full_time_course = ((1:80)-27) * 10;


%% plot accuracy time-course

plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
colors = {[102 0 102]./255; [132 186 91]./255; [114 147 203]./255;
          [255 150 0]./255; [171 104 87]./255; [255 150 193]./255};
sem = @(x) std(x, [], 1)/sqrt(size(x, 1));

figure(); hold on
ylim([0.35 1])

for soa_inx=[3 4 5 6]
    shadedErrorBar(plt_time, nanmean(trueAccuracy{soa_inx}, 1), sem(trueAccuracy{soa_inx}),...
        {'color', colors{soa_inx}, 'linewidth', 1},1);
end
        
yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')
legend({'', '150', '', '100', '', '75', '', '50'},...
    'Location', 'northwest')
xlabel('Time (ms)'); ylabel('Accuracy')
title('Grand Analysis of Models, 5-frames Moving Window')

%% Table of Values

% initiate table of statistics
soa_labels = {'retino'; 'noMask'; 'SOA150'; 'SOA100'; 'SOA75'; 'SOA50'};
accuracyStats = table(soa_labels);
shfldStats = table(soa_labels);

% time points to analyze
time_points_ms = {-40:10:-20, 100:10:120, 170:10:190};
time_points = cellfun(@(x) find(ismember(plt_time, x)), time_points_ms, 'UniformOutput', false);

rng(42);
for soa_inx=[3 4 5 6]
    for tp=1:3
        frames = time_points{tp};
        trueAcc = trueAccuracy{soa_inx}(:, frames);
        shfldAcc = shfldAccuracy{soa_inx}(:, frames);
        
        switch tp
            case 1
                accuracyStats.value_0(soa_inx) = mean(trueAcc, 'all');
                shfldStats.value_0(soa_inx) = mean(shfldAcc, 'all');

                accuracyStats.sem_0(soa_inx) = std(mean(trueAcc, 2)) / sqrt(size(trueAcc, 1));
                shfldStats.sem_0(soa_inx) = std(mean(shfldAcc, 2)) / sqrt(size(shfldAcc, 1));

                accuracyStats.rnksm_0(soa_inx) = ranksum(mean(trueAcc, 2), mean(shfldAcc, 2));

            case 2
                accuracyStats.value_1(soa_inx) = mean(trueAcc, 'all');
                shfldStats.value_1(soa_inx) = mean(shfldAcc, 'all');

                accuracyStats.sem_1(soa_inx) = std(mean(trueAcc, 2)) / sqrt(size(trueAcc, 1));
                shfldStats.sem_1(soa_inx) = std(mean(shfldAcc, 2)) / sqrt(size(shfldAcc, 1));

                accuracyStats.rnksm_1(soa_inx) = ranksum(mean(trueAcc, 2), mean(shfldAcc, 2));

            case 3
                accuracyStats.value_2(soa_inx) = mean(trueAcc, 'all');
                shfldStats.value_2(soa_inx) = mean(shfldAcc, 'all');

                accuracyStats.sem_2(soa_inx) = std(mean(trueAcc, 2)) / sqrt(size(trueAcc, 1));
                shfldStats.sem_2(soa_inx) = std(mean(shfldAcc, 2)) / sqrt(size(shfldAcc, 1));

                accuracyStats.rnksm_2(soa_inx) = ranksum(mean(trueAcc, 2), mean(shfldAcc, 2));
        end
    end
end

%% Combine 150&100, 75&50, show plot

switch monkey
    case 'frodo'
        ylimits = [0.4 1];
        ySig = 0.95;
    case 'tolkin'
        ylimits = [0.45 0.7];
        ySig = 0.65;
end

plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
sem = @(x) std(x, [], 1)/sqrt(size(x, 1));

trueAccuracy{1} = cat(1, trueAccuracy{3}, trueAccuracy{4});
trueAccuracy{2} = cat(1, trueAccuracy{5}, trueAccuracy{6});
trueAccuracy = trueAccuracy(1:2, 1);

shfldAccuracy{1} = cat(1, shfldAccuracy{3}, shfldAccuracy{4});
shfldAccuracy{2} = cat(1, shfldAccuracy{5}, shfldAccuracy{6});
shfldAccuracy = shfldAccuracy(1:2, 1);


colors = {[102 0 102]./255; [132 186 91]./255; [114 147 203]./255;
          [255 150 0]./255; [171 104 87]./255; [255 150 193]./255};

rnksmLong = false(40,1);
rnksmShort = false(40,1);
rnksmLongShort = false(40,1);
for i = 1:40
    rnksmLong(i) = ranksum(trueAccuracy{1}(:,i), shfldAccuracy{1}(:,i)) < 0.05;
    rnksmShort(i) = ranksum(trueAccuracy{2}(:,i), shfldAccuracy{2}(:,i)) < 0.05;
    rnksmLongShort(i) = ranksum(trueAccuracy{1}(:,i), trueAccuracy{2}(:,i)) < 0.05;
end

figure(); hold on
ylim(ylimits)

for soa_inx=[1 2]
    shadedErrorBar(plt_time, nanmean(trueAccuracy{soa_inx}, 1), sem(trueAccuracy{soa_inx}),...
        {'color', colors{soa_inx}, 'linewidth', 1},0);
end
plot(plt_time(rnksmLong), ySig, '*', 'color', [102 0 102]./255)
plot(plt_time(rnksmShort), ySig+0.02, '*', 'color', [132 186 91]./255)
plot(plt_time(rnksmLongShort), ySig+0.04, '*', 'color', 'k')

yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')
legend({'', '150 & 100', '', '75 & 50'},...
    'Location', 'northwest')
xlabel('Time (ms)'); ylabel('Accuracy')
title('Grand Analysis of Models, 5-frames Moving Window')

%% early vs. late peaks, short & long SOAs

% Data points
switch monkey
    case 'frodo'
        tp1 = 70:10:90;
        tp2 = 190:10:210;
        yLim = 1;
    case 'tolkin'
        tp1 = 90:10:110;
        tp2 = 150:10:170;
        yLim = 0.7;
end


pk1 = find(ismember(plt_time, tp1));
pk2 = find(ismember(plt_time, tp2));


x = [1, 2, 4, 5];  % X-coordinates for the points
y = [mean(trueAccuracy{2}(:,pk1), 'all');
     mean(trueAccuracy{1}(:,pk1), 'all');
     mean(trueAccuracy{2}(:,pk2), 'all');
     mean(trueAccuracy{1}(:,pk2), 'all')];

shuffledY = [mean(shfldAccuracy{2}(:,pk1), 'all');
             mean(shfldAccuracy{1}(:,pk1), 'all');
             mean(shfldAccuracy{2}(:,pk2), 'all');
             mean(shfldAccuracy{1}(:,pk2), 'all')]; 

sem = @(x) std(mean(x,2))/sqrt(numel(mean(x,2)));

errors = [sem(trueAccuracy{2}(:,pk1));
          sem(trueAccuracy{1}(:,pk1));
          sem(trueAccuracy{2}(:,pk2));
          sem(trueAccuracy{1}(:,pk2))];

shuffledErrors = [sem(shfldAccuracy{2}(:,pk1));
                  sem(shfldAccuracy{1}(:,pk1));
                  sem(shfldAccuracy{2}(:,pk2));
                  sem(shfldAccuracy{1}(:,pk2))];

rnksm = [ranksum(mean(trueAccuracy{2}(:,pk1), 2), mean(shfldAccuracy{2}(:,pk1), 2));
        ranksum(mean(trueAccuracy{1}(:,pk1), 2), mean(shfldAccuracy{1}(:,pk1), 2));
        ranksum(mean(trueAccuracy{2}(:,pk2), 2), mean(shfldAccuracy{2}(:,pk2), 2));
        ranksum(mean(trueAccuracy{1}(:,pk2), 2), mean(shfldAccuracy{1}(:,pk2), 2))];

% Get pValues
p1 = ranksum(mean(trueAccuracy{1}(:,pk1), 2), mean(trueAccuracy{2}(:,pk1), 2));
p2 = ranksum(mean(trueAccuracy{1}(:,pk2), 2), mean(trueAccuracy{2}(:,pk2), 2));

% SOA labels
soaLabels = {'Short', 'Long', 'Short', 'Long'};

% Create a figure
figure; hold on;

% Plot the data points with error bars
errorbar(x, y, errors, 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
errorbar(x, shuffledY, shuffledErrors, 'o', 'MarkerFaceColor', '#929591', 'Color', '#929591', 'LineWidth', 1.5);
yline(0.5, '--')

% Connect each pair with a line
line([x(1), x(2)], [y(1), y(2)], 'Color', 'k', 'LineWidth', 1.5); % First pair
line([x(3), x(4)], [y(3), y(4)], 'Color', 'k', 'LineWidth', 1.5); % Second pair

line([x(1), x(2)], [shuffledY(1), shuffledY(2)], 'Color', '#929591', 'LineWidth', 1.5); % First pair
line([x(3), x(4)], [shuffledY(3), shuffledY(4)], 'Color', '#929591', 'LineWidth', 1.5); % Second pair

% Add significance markers above each pair
sigstar({[1, 2], [4, 5]}, [p1, p2]);  % Replace p-values with your actual values

% Set axis labels and title
% xlabel('SOA Conditions');
ylabel('Peak Values');
title('Short vs. Long SOAs Accuracy');

% Set X-axis limits and labels
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', {sprintf('%g ms', tp1), sprintf('%g ms', tp2)});
xlim([0.5, 5.5]);
ylim([0.45 yLim])

% Add labels underneath each data point
y_loc = ones(4,1) * min(shuffledY) - 0.02;
text(x, y_loc, soaLabels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hold off;


%% Plot Signal in a Colorcode vs. Performance

% Frodo
load('data\normalized2_regist\frodo_early_sessions\accuracyCell.mat');
load('data\normalized2_regist\frodo_early_sessions\all_data.mat')

conInx = 3;
frames = 23:57;
time_course = ((frames) - 27) * 10;
cond_labels = {'retino' 'noMoask' 'SOA150' 'SOA100' 'SOA75' 'SOA50'};

siganl = nanmean([nanmean(all_data{conInx,1}(:,frames,:), [1 3]);...
                  nanmean(all_data{conInx,2}(:,frames,:), [1 3])]);
accuracy = nanmean(accuracyCell{conInx}(:,1:35), 1);
accuracy(end) = NaN;

% signal time-course
figure; hold on
subplot(2,1,1)
patch(time_course, siganl, accuracy, 'EdgeColor', 'interp', 'LineWidth', 5);
caxis([0.5 0.85])
ylim([-1 8])
xline(0', 'r', 'HandleVisibility','off')
colorbar('westoutside');
title(sprintf('Signal TC in zScore All V1, %s', cond_labels{conInx}))
ylabel('\Delta ZS')

% accuracy time-course
subplot(2,1,2)
plot(time_course, accuracy)
title('Accuracy of Models Moving Window')
ylim([0.4, 1])
yline(0.5, '--')
ylabel('Accuracy')
xline(0', 'r', 'HandleVisibility','off')
xlabel('Time (ms)')


%% Load correct/error data

full_time_course = ((1:80)-27) * 10;

plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
colors = {[102 0 102]./255; [132 186 91]./255; [114 147 203]./255;
          [255 150 0]./255; [171 104 87]./255; [255 150 193]./255};

% Load data
switch monkey
    case 'frodo'
        accuracyCell = importdata('data\normalized2_regist\frodo_early_sessions\correctError\accuracyCell_correct_only.mat');
        accuracyCellShfld = importdata('data\normalized2_regist\frodo_early_sessions\correctError\accuracyCell_correct_only_shfld.mat');
        errAccuracyCell = importdata('data\normalized2_regist\frodo_early_sessions\correctError\accuracyCell_error_only.mat');
        errAccuracyCellShfld = importdata('data\normalized2_regist\frodo_early_sessions\correctError\accuracyCell_error_only_shfld.mat');
    case 'tolkin'
        accuracyCell = importdata('data\normalized2_regist\tolkin_early_sessions\correctError\accuracyCell_correct_only.mat');
        accuracyCellShfld = importdata('data\normalized2_regist\tolkin_early_sessions\correctError\accuracyCell_correct_only_shfld.mat');
        errAccuracyCell = importdata('data\normalized2_regist\tolkin_early_sessions\correctError\accuracyCell_error_only.mat');
        errAccuracyCellShfld = importdata('data\normalized2_regist\tolkin_early_sessions\correctError\accuracyCell_error_only_shfld.mat');

        accuracy_late = importdata('data\normalized2_regist\tolkin_late_sessions\correctError\accuracyCell_correct_only.mat');
        accuracy_late_shfld = importdata('data\normalized2_regist\tolkin_late_sessions\correctError\accuracyCell_correct_only_shfld.mat');
        accuracy_late_error = importdata('data\normalized2_regist\tolkin_late_sessions\correctError\accuracyCell_error_only.mat');
        accuracy_late_error_shfld = importdata('data\normalized2_regist\tolkin_late_sessions\correctError\accuracyCell_error_only_shfld.mat');

        for i=3:6
            accuracyCell{i,1} = cat(1, accuracyCell{i,1}, accuracy_late{i,1});
            accuracyCellShfld{i,1} = cat(1, accuracyCellShfld{i,1}, accuracy_late_shfld{i,1});
            errAccuracyCell{i,1} = cat(1, errAccuracyCell{i,1}, accuracy_late_error{i,1});
            errAccuracyCellShfld{i,1} = cat(1, errAccuracyCellShfld{i,1}, accuracy_late_error_shfld{i,1});
        end
end



sem = @(x) std(x, [], 1)/sqrt(size(x, 1));

soa_inx = 5;
length = 'long'; 

switch length
    case 'single'
        correctTC = accuracyCell{soa_inx,1};
        errTC = errAccuracyCell{soa_inx,1};
    case 'long'
        correctTC = cat(1, accuracyCell{3,1}, accuracyCell{4,1});
        errTC = cat(1, errAccuracyCell{3,1}, errAccuracyCell{4,1});
    case 'short'
        correctTC = cat(1, accuracyCell{5,1}, accuracyCell{6,1});
        errTC = cat(1, errAccuracyCell{5,1}, errAccuracyCell{6,1});
end


rnksm = false(40,1);
for i=1:40
    rnksm(i) = ranksum(correctTC(:,i), errTC(:,i)) < 0.05;
end

figure; hold on
shadedErrorBar(plt_time, nanmean(correctTC, 1), sem(correctTC),...
    {'color', colors{3}, 'linewidth', 1},1);
shadedErrorBar(plt_time, nanmean(errTC, 1), sem(errTC),...
    {'color', colors{4}, 'linewidth', 1},1);
plot(plt_time(rnksm), 0.95, '*', 'color', 'k')

yline(0.5)
legend({'', 'correct', '', 'error'})
title('Tolkin longSOA, correct vs. error (trained on correct only)')
ylim([0.3 1])


%% early vs. late peaks, short & long SOAs

% Data points
switch monkey
    case 'frodo'
        tp1 = 70:10:90;
        tp2 = 190:10:210;
        yLims = [0.45 1];

        % Combine 150&100, 75&50
        accuracyCell{1} = cat(1, accuracyCell{3}, accuracyCell{4});
        accuracyCell{2} = cat(1, accuracyCell{5}, accuracyCell{6});
        accuracyCell = accuracyCell(1:2, 1);

        errAccuracyCell{1} = cat(1, errAccuracyCell{3}, errAccuracyCell{4});
        errAccuracyCell{2} = cat(1, errAccuracyCell{5}, errAccuracyCell{6});
        errAccuracyCell = errAccuracyCell(1:2, 1);

    case 'tolkin'
        tp1 = 90:10:110;
        tp2 = 150:10:170;
        yLims = [0.4 0.7];

        % Combine 150&100, 75&50
        accuracyCell{1} = cat(1, accuracyCell{3}, accuracyCell{4});
        accuracyCell{2} = accuracyCell{5};
        accuracyCell = accuracyCell(1:2, 1);

        errAccuracyCell{1} = cat(1, errAccuracyCell{3}, errAccuracyCell{4});
        errAccuracyCell{2} = errAccuracyCell{5};
        errAccuracyCell = errAccuracyCell(1:2, 1);
end


pk1 = find(ismember(plt_time, tp1));
pk2 = find(ismember(plt_time, tp2));



x = [1, 2, 4, 5];  % X-coordinates for the points
y_correct = [mean(accuracyCell{1}(:,pk1), 'all');
     mean(accuracyCell{2}(:,pk1), 'all');
     mean(accuracyCell{1}(:,pk2), 'all');
     mean(accuracyCell{2}(:,pk2), 'all')];

y_error = [mean(errAccuracyCell{1}(:,pk1), 'all');
         mean(errAccuracyCell{2}(:,pk1), 'all');
         mean(errAccuracyCell{1}(:,pk2), 'all');
         mean(errAccuracyCell{2}(:,pk2), 'all')];

sem = @(x) std(mean(x,2))/sqrt(numel(mean(x,2)));

sem_correct = [sem(accuracyCell{1}(:,pk1));
              sem(accuracyCell{2}(:,pk1));
              sem(accuracyCell{1}(:,pk2));
              sem(accuracyCell{2}(:,pk2))];

sem_error = [sem(errAccuracyCell{1}(:,pk1));
              sem(errAccuracyCell{2}(:,pk1));
              sem(errAccuracyCell{1}(:,pk2));
              sem(errAccuracyCell{2}(:,pk2))];

rnksm = [ranksum(mean(accuracyCell{1}(:,pk1), 2), mean(errAccuracyCell{1}(:,pk1), 2));
        ranksum(mean(accuracyCell{2}(:,pk1), 2), mean(errAccuracyCell{2}(:,pk1), 2));
        ranksum(mean(accuracyCell{1}(:,pk2), 2), mean(errAccuracyCell{1}(:,pk2), 2));
        ranksum(mean(accuracyCell{2}(:,pk2), 2), mean(errAccuracyCell{2}(:,pk2), 2))];

% SOA labels
predictorLabels = {'Correct', 'Error', 'Correct', 'Error'};

% Create a figure
figure; hold on;

% Plot the data points with error bars
errorbar(1, y_correct(1), sem_correct(1), 'ko', 'MarkerFaceColor', "#800000", 'LineWidth', 1.5);
errorbar(1, y_correct(2), sem_correct(2), 'ko', 'MarkerFaceColor', "#FFA500", 'LineWidth', 1.5);
errorbar(4, y_correct(3), sem_correct(3), 'ko', 'MarkerFaceColor', "#800000", 'LineWidth', 1.5);
errorbar(4, y_correct(4), sem_correct(4), 'ko', 'MarkerFaceColor', "#FFA500", 'LineWidth', 1.5);

errorbar(2, y_error(1), sem_error(1), 'ko', 'MarkerFaceColor', "#800000", 'LineWidth', 1.5);
errorbar(2, y_error(2), sem_error(2), 'ko', 'MarkerFaceColor', "#FFA500", 'LineWidth', 1.5);
errorbar(5, y_error(3), sem_error(3), 'ko', 'MarkerFaceColor', "#800000", 'LineWidth', 1.5);
errorbar(5, y_error(4), sem_error(4), 'ko', 'MarkerFaceColor', "#FFA500", 'LineWidth', 1.5);

yline(0.5, '--')

% Connect each pair with a line
line([x(1), x(2)], [y_correct(1), y_error(1)], 'Color', "#800000", 'LineWidth', 1.5); % First pair
line([x(1), x(2)], [y_correct(2), y_error(2)], 'Color', "#FFA500", 'LineWidth', 1.5); % First pair

line([x(3), x(4)], [y_correct(3), y_error(3)], 'Color', "#800000", 'LineWidth', 1.5); % Second pair
line([x(3), x(4)], [y_correct(4), y_error(4)], 'Color', "#FFA500", 'LineWidth', 1.5); % Second pair

% Add significance markers above each pair
H1 = sigstar({[1, 2], [4, 5]}, [rnksm(1) rnksm(3)]);  % Replace p-values with your actual values
set(H1,'color',"#800000")
H2 = sigstar({[1, 2], [4, 5]}, [rnksm(2) rnksm(4)]);  % Replace p-values with your actual values
set(H2,'color',"#FFA500")

% Set axis labels and title
% xlabel('SOA Conditions');
ylabel('Peak Values');
title('Correct vs. Error SOAs Accuracy');

% Set X-axis limits and labels
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', {sprintf('%g ms', tp1), sprintf('%g ms', tp2)});
xlim([0.5, 5.5]);
ylim(yLims)

% Add labels underneath each data point
y_loc = ones(4,1) * yLims(1) + 0.02;
text(x, y_loc, predictorLabels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
legend({'long SOA' 'short SOA'})

hold off;


%% early vs. late peaks, ALL SOAs

% Combine all SOAs
correctTC = [];
correctShfldTC = [];
errorTC = [];
errorShfldTC = [];

switch monkey
    case 'frodo'
        tp1 = 70:10:90;
        tp2 = 190:10:210;
        yLims = [0.45 1];

        for i=3:6
            correctTC = cat(1, correctTC, accuracyCellShfld{i});
            correctShfldTC = cat(1, correctShfldTC, accuracyCellShfld{i});
            errorTC = cat(1, errorTC, errAccuracyCellShfld{i});
            errorShfldTC = cat(1, errorShfldTC, errAccuracyCellShfld{i});
        end

    case 'tolkin'
        tp1 = 90:10:110;
        tp2 = 150:10:170;
        yLims = [0.4 0.7];

        for i=3:5
            correctTC = cat(1, correctTC, accuracyCell{i});
            correctShfldTC = cat(1, correctShfldTC, accuracyCellShfld{i});
            errorTC = cat(1, errorTC, errAccuracyCell{i});
            errorShfldTC = cat(1, errorShfldTC, errAccuracyCellShfld{i});
        end
end


pk1 = find(ismember(plt_time, tp1));
pk2 = find(ismember(plt_time, tp2));



x = [1, 2, 4, 5];  % X-coordinates for the points
sem = @(x) std(mean(x,2))/sqrt(numel(mean(x,2)));
rnksm =  @(x,y) ranksum(mean(x,2), mean(y,2));

% find timepoint values
y_correct = [mean(correctTC(:,pk1), 'all');
            mean(correctTC(:,pk2), 'all')];
sem_correct = [sem(correctTC(:,pk1));
              sem(correctTC(:,pk2))];
rnskm_correct = [rnksm(correctTC(:,pk1), correctShfldTC(:,pk1));
                 rnksm(correctTC(:,pk2),  correctShfldTC(:,pk2))];


y_error = [mean(errorTC(:,pk1), 'all');
         mean(errorTC(:,pk2), 'all')];
sem_error = [sem(errorTC(:,pk1));
              sem(errorTC(:,pk2))];
rnskm_error = [rnksm(errorTC(:,pk1), errorShfldTC(:,pk1));
                 rnksm(errorTC(:,pk2),  errorShfldTC(:,pk2))];

rnksm_correct_error = [rnksm(correctTC(:,pk1), errorTC(:,pk1));
                       rnksm(correctTC(:,pk2), errorTC(:,pk2))];

% Create a figure
figure; hold on;

% Plot the data points with error bars
% bar(1, y_correct(1), 'FaceColor', "#069AF3");
errorbar(1, y_correct(1), sem_correct(1), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', "#069AF3", 'LineWidth', 1.5);
% bar(4, y_correct(2), 'FaceColor', "#069AF3");
errorbar(4, y_correct(2), sem_correct(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', "#069AF3", 'LineWidth', 1.5);

% bar(2, y_error(1), 'FaceColor', "#FF6347");
errorbar(2, y_error(1), sem_error(1), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', "#FF6347", 'LineWidth', 1.5);
% bar(5, y_error(2), 'FaceColor', "#FF6347");
errorbar(5, y_error(2), sem_error(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', "#FF6347", 'LineWidth', 1.5);

yline(0.5, '--')

% Connect each pair with a line
line([x(1), x(2)], [y_correct(1), y_error(1)], 'Color', "#929591", 'LineWidth', 1.5); % First pair
line([x(3), x(4)], [y_correct(2), y_error(2)], 'Color', "#929591", 'LineWidth', 1.5); % First pair


% Add significance markers above each pair
sigstar({[1, 2], [4, 5]}, [rnksm_correct_error(1) rnksm_correct_error(2)]);  % Replace p-values with your actual values

% Set axis labels and title
% xlabel('SOA Conditions');
ylabel('Peak Values');
title('Correct vs. Error SOAs Accuracy');

% Set X-axis limits and labels
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', {sprintf('%g ms', tp1), sprintf('%g ms', tp2)});
xlim([0.5, 5.5]);
ylim(yLims)

% Add labels underneath each data point
legend({'correct trials' '' 'error trials'})

hold off;

% (y_correct - y_error) ./ y_correct