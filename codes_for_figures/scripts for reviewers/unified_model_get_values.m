% Change monkey
monkey = 'tolkin';

% Load data
switch monkey
    case 'frodo'
        trueAccuracy = importdata('E:\Itiel\BM_ARTICLE\figures\Revision answers\unified_model\frodo\test_accuracy.mat');
        shfldAccuracy = importdata('E:\Itiel\BM_ARTICLE\figures\Revision answers\unified_model\frodo\shfld_test_accuracy.mat');

    case 'tolkin'
        accuracy_early = importdata('E:\Itiel\BM_ARTICLE\figures\Revision answers\unified_model\tolkin_early_sessions\test_accuracy.mat');
        accuracy_early_shlfd = importdata('E:\Itiel\BM_ARTICLE\figures\Revision answers\unified_model\tolkin_early_sessions\shfld_test_accuracy.mat');
        
        accuracy_late = importdata('E:\Itiel\BM_ARTICLE\figures\Revision answers\unified_model\tolkin_late_sessions\test_accuracy.mat');
        accuracy_late_shfld = importdata('E:\Itiel\BM_ARTICLE\figures\Revision answers\unified_model\tolkin_late_sessions\shfld_test_accuracy.mat');
        
        for i=1:4
            trueAccuracy(:,:,i) = cat(1, accuracy_early(:,:,i), accuracy_late(:,:,i));
            shfldAccuracy(:,:,i) = cat(1, accuracy_early_shlfd(:,:,i), accuracy_late_shfld(:,:,i));
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

for soa_inx=[1 2 3 4]
    shadedErrorBar(plt_time, nanmean(trueAccuracy(:,:,soa_inx), 1), sem(trueAccuracy(:,:,soa_inx)),...
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

%% pValues & Bar plot

% initiate pValue table

time_point = time_points_ms{3};     %[180 200];  % Change to desired values
frames = time_points{3};

sz = [3 4];
varTypes = ["string","double","double","double"];
varNames = ["soa_labels" "SOA150" "SOA100" "SOA75"];
pValues = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames', varNames);
pValues.soa_labels = ["SOA100"; "SOA75"; "SOA50"];
% frames = find(plt_time == time_point(1)):find(plt_time == time_point(2));

for inx_a=3:5
    for inx_b=inx_a+1:6
        pValues{inx_b-3, inx_a-1} = ranksum(mean(trueAccuracy{inx_a}(:,frames), 2), mean(trueAccuracy{inx_b}(:,frames), 2));
    end
end
disp(pValues)


% Arrange data for bars
X = categorical({'SOA50','SOA75','SOA100', 'SOA150'});
X = reordercats(X, {'SOA50','SOA75','SOA100', 'SOA150'});
Y = flipud(accuracyStats.value_2(3:6));
sem_vec = flipud(accuracyStats.sem_2(3:6));

% Get pValues
p50_75 = pValues{3, 4};
p50_100 = pValues{3, 3};
p50_150 = pValues{3, 2};
p75_100 = pValues{2, 3};
p75_150 = pValues{2, 2};
p100_150 = pValues{1, 2};

% Colors
colors = [
    [80 0 80]/255;
    [120 20 120]/255;
    [160 60 160]/255;
    [200 100 200]/255;
];

% Plot figure
figure; hold on
for k = 1:size(Y,1)
    bar(k, Y(k), 'FaceColor', colors(k, :));
end
errorbar(1:4,Y,sem_vec,'k','linestyle','none');
set(gca, 'XTick', 1:4);
set(gca,'XTickLabel',X)

% Plot the stars
groups = {[1 2] [1 3] [1 4] [2 3] [2 4] [3 4]};
sigstar(groups,[(p50_75), (p50_100), (p50_150), (p75_100), (p75_150), (p100_150)]);

ylim([0.5 1])
ylabel('Figure-Ground modulation')
title(sprintf('SVM Accuracy of Frodo 2nd peak, %d-%d ms', time_point(1), time_point(end)))
hold off




%% Plot 1st & 2nd peaks for 50 & 150 SOAs


% % Data points
% soa_idx = [3 6];
% x = [1, 2, 4, 5];  % X-coordinates for the points
% y = [flipud(accuracyStats.value_1(soa_idx));
%     flipud(accuracyStats.value_2(soa_idx))]; 
% shuffledY = [flipud(shfldStats.value_1(soa_idx));
%             flipud(shfldStats.value_2(soa_idx))]; 
% errors = [flipud(accuracyStats.sem_1(soa_idx));
%            flipud(accuracyStats.sem_1(soa_idx))];
% shuffledErrors = [flipud(shfldStats.sem_1(soa_idx));
%            flipud(shfldStats.sem_1(soa_idx))];

% Data points
shortSoaIdx = 6;
longSoaIdx = 3;
soa_labels = {'retino' 'noMask' '150' '100' '75' '50'};

colors = {[102 0 102]./255; [132 186 91]./255; [114 147 203]./255;
          [255 150 0]./255; [171 104 87]./255; [255 150 193]./255};

rnksm = false(40,1);
for i=1:40
    rnksm(i) = ranksum(trueAccuracy{longSoaIdx}(:,i), trueAccuracy{shortSoaIdx}(:,i)) < 0.05;
end

figure(); hold on
ylim([0.35 1])

for soa_inx=[longSoaIdx shortSoaIdx]
    shadedErrorBar(plt_time, nanmean(trueAccuracy{soa_inx}, 1), sem(trueAccuracy{soa_inx}),...
        {'color', colors{soa_inx}, 'linewidth', 1},1);
end
% plot(plt_time(rnksm), 0.7, '*', 'color', 'k')
plot(plt_time(rnksm), 0.95, '*', 'color', 'k')

yline(0.5, 'k', 'linestyle', '--', 'HandleVisibility','off')
xline(0, 'r', 'HandleVisibility','off')
legend({'', soa_labels{longSoaIdx}, '', soa_labels{shortSoaIdx}},...
    'Location', 'northwest')
xlabel('Time (ms)'); ylabel('Accuracy')
title('Grand Analysis of Models, 5-frames Moving Window')


tp1 = 80;
tp2 = 190;
pk1 = find(plt_time == tp1);
pk2 = find(plt_time == tp2);

x = [1, 2, 4, 5];  % X-coordinates for the points
y = [mean(trueAccuracy{shortSoaIdx}(:,pk1));
     mean(trueAccuracy{longSoaIdx}(:,pk1));
     mean(trueAccuracy{shortSoaIdx}(:,pk2));
     mean(trueAccuracy{longSoaIdx}(:,pk2))];

shuffledY = [mean(shfldAccuracy{shortSoaIdx}(:,pk1));
             mean(shfldAccuracy{longSoaIdx}(:,pk1));
             mean(shfldAccuracy{shortSoaIdx}(:,pk2));
             mean(shfldAccuracy{longSoaIdx}(:,pk2))]; 

sem = @(x) std(x)/sqrt(numel(x));

errors = [sem(trueAccuracy{shortSoaIdx}(:,pk1));
          sem(trueAccuracy{shortSoaIdx}(:,pk1));
          sem(trueAccuracy{longSoaIdx}(:,pk2));
          sem(trueAccuracy{longSoaIdx}(:,pk2))];

shuffledErrors = [sem(shfldAccuracy{shortSoaIdx}(:,pk1));
                  sem(shfldAccuracy{longSoaIdx}(:,pk1));
                  sem(shfldAccuracy{shortSoaIdx}(:,pk2));
                  sem(shfldAccuracy{longSoaIdx}(:,pk2))];
% Get pValues
p1 = ranksum(trueAccuracy{shortSoaIdx}(:,pk1), trueAccuracy{longSoaIdx}(:,pk1));
p2 = ranksum(trueAccuracy{shortSoaIdx}(:,pk2), trueAccuracy{longSoaIdx}(:,pk2));
% % Get pValues
% p1 = ranksum(trueAccuracy{soa_idx(1)}(:,12), trueAccuracy{soa_idx(2)}(:,12));
% p2 = ranksum(trueAccuracy{soa_idx(1)}(:,24), trueAccuracy{soa_idx(2)}(:,24));

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
xlabel('SOA Conditions');
ylabel('Peak Values');
title('Short vs. Long SOAs Accuracy');

% Set X-axis limits and labels
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', {sprintf('%g ms', tp1), sprintf('%g ms', tp2)});
xlim([0.5, 5.5]);

% Add labels underneath each data point
y_loc = ones(4,1) * min(shuffledY) - 0.02;
text(x, y_loc, soaLabels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hold off;

%% Combine 150&100, 75&50, show plot

switch monkey
    case 'frodo'
        ylimits = [0.4 1];
        ySig = 0.95;
    case 'tolkin'
        ylimits = [0.45 0.8];
        ySig = 0.65;
end

plt_frames = 23:62;
plt_time = full_time_course(plt_frames);
sem = @(x) std(x, [], 1)/sqrt(size(x, 1));

trueLong = cat(1, trueAccuracy(:,:,1), trueAccuracy(:,:,2));
trueShort = cat(1, trueAccuracy(:,:,3), trueAccuracy(:,:,4));
shfldLong = cat(1, shfldAccuracy(:,:,1), shfldAccuracy(:,:,2));
shfldShort = cat(1, shfldAccuracy(:,:,3), shfldAccuracy(:,:,4));


colors = {[102 0 102]./255; [132 186 91]./255; [114 147 203]./255;
          [255 150 0]./255; [171 104 87]./255; [255 150 193]./255};

% Calculate sig.
rnksmLong = false(40,1);
rnksmShort = false(40,1);
rnksmLongShort = false(40,1);
for i = 1:40
    rnksmLong(i) = ranksum(trueLong(:,i), shfldLong(:,i)) < 0.05;
    rnksmShort(i) = ranksum(trueShort(:,i), shfldShort(:,i)) < 0.05;
    rnksmLongShort(i) = ranksum(trueLong(:,i), trueShort(:,i)) < 0.05;
end

% Plot figure
figure(); hold on
ylim(ylimits)

% Curves
shadedErrorBar(plt_time, nanmean(trueLong, 1), sem(trueLong),...
    {'color', colors{1}, 'linewidth', 1},0);
shadedErrorBar(plt_time, nanmean(trueShort, 1), sem(trueShort),...
    {'color', colors{2}, 'linewidth', 1},0);

% Sigstars
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
        yLim = [0.4 1];
    case 'tolkin'
        tp1 = 90:10:110;
        tp2 = 150:10:170;
        yLim = [0.45 0.8];
end


pk1 = find(ismember(plt_time, tp1));
pk2 = find(ismember(plt_time, tp2));


x = [1, 2, 4, 5];  % X-coordinates for the points
y = [mean(trueShort(:,pk1), 'all');   % short
     mean(trueLong(:,pk1), 'all');   % long
     mean(trueShort(:,pk2), 'all');
     mean(trueLong(:,pk2), 'all')];

shuffledY = [mean(shfldShort(:,pk1), 'all');
             mean(shfldLong(:,pk1), 'all');
             mean(shfldShort(:,pk2), 'all');
             mean(shfldLong(:,pk2), 'all')]; 

sem = @(x) std(mean(x,2))/sqrt(numel(mean(x,2)));

errors = [sem(trueShort(:,pk1));   % short
          sem(trueLong(:,pk1));   % long
          sem(trueShort(:,pk2));
          sem(trueLong(:,pk2))];

shuffledErrors = [sem(shfldShort(:,pk1));
                  sem(shfldLong(:,pk1));
                  sem(shfldShort(:,pk2));
                  sem(shfldLong(:,pk2))];

rnksm = [ranksum(mean(trueShort(:,pk1), 2), mean(shfldShort(:,pk1), 2));   % short
        ranksum(mean(trueLong(:,pk1), 2), mean(shfldLong(:,pk1), 2));   % long
        ranksum(mean(trueShort(:,pk2), 2), mean(shfldShort(:,pk2), 2));
        ranksum(mean(trueLong(:,pk2), 2), mean(shfldLong(:,pk2), 2))];

% Get pValues
p1 = ranksum(mean(trueShort(:,pk1), 2), mean(trueLong(:,pk1), 2));
p2 = ranksum(mean(trueShort(:,pk2), 2), mean(trueLong(:,pk2), 2));

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
ylim(yLim)

% Add labels underneath each data point
y_loc = ones(4,1) * min(shuffledY) - 0.02;
text(x, y_loc, soaLabels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hold off;