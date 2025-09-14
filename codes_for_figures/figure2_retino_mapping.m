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
fileDates = {'1609_' '2309_' '0610_' '1410_' '2511_' '0212_' '1703_'};
load("data\BM\Frodo\VSDI\datesConsCell.mat");
load("data\BM\Frodo\VSDI\maskDatesCons.mat");


TEST_FRAMES = 1:100;
MAP_WINDOW = 34:38;
S = 1;

% Load matrices
vertical = importdata([vsd_folder, '\a\condsn2_bl.mat']);
horizontal = importdata([vsd_folder, '\a\condsn4_bl.mat']);
load([vsd_folder, '\shani_ROIs'])
chamber = importdata([vsd_folder, '\ID_chamber.mat']);

vertical = mfilt2(vertical, 100, 100, S, 'lm');
horizontal = mfilt2(horizontal, 100, 100, S, 'lm');

% Adjust as desire
full_time_course = ((1:80)-27) * 10;
time_point = full_time_course(MAP_WINDOW);

verMap = squeeze(mean(vertical(:,MAP_WINDOW,:), [2 3]));
verMap = reshape(verMap, 100, 100)';
% circle_contour_gui(verMap) % save a variable 'chamber'

%% Draw Retino Maps
% Create contours
contourChamber = zeros(100, 100);
contourChamber(chamber) = 1; contourChamber = contourChamber';

contourVer = zeros(100, 100);
contourVer(roiPix_fit_ver) = 1; contourVer = contourVer';

contourHor = zeros(100, 100);
contourHor(roiPix_fit_hor) = 1; contourHor = contourHor';

% Create vertical map
verMap = ones(10000, 1)*1000;
verMap(chamber) = squeeze(mean(vertical(chamber,MAP_WINDOW,:), [2 3])) -1;
verMap = reshape(verMap, 100, 100)';

% Create horizontal map
horMap = ones(10000, 1)*1000;
horMap(chamber) = squeeze(mean(horizontal(chamber,MAP_WINDOW,:), [2 3])) -1;
horMap = reshape(horMap, 100, 100)';


% Define the minimum and maximum values for clipping the colormap
minValue = 0.0001; % Adjust as needed based on your data
maxValue = 0.0025; % Adjust as needed based on your data


figure;
imagesc(verMap, [minValue maxValue]); colormap mapgeog
hold on
contour(contourChamber, 'k');
contour(contourVer, 'cyan');
title(sprintf('Vertical Retino, time: %g:%g ms', time_point(1), time_point(end)))
text(10, 5, sprintf('frames: %d:%d', MAP_WINDOW(1), MAP_WINDOW(end)))
text(10, 10, sprintf('smoothing: S=%g', S))
text(10,15, sprintf('cliping: [%g, %g]', minValue, maxValue))
hold off


figure;
imagesc(horMap, [minValue maxValue]); colormap mapgeog
hold on
contour(contourChamber, 'k');
contour(contourHor, 'cyan');
title(sprintf('Horizontal Retino, time: %g:%g ms', time_point(1), time_point(end)))
text(10, 5, sprintf('frames: %d:%d', MAP_WINDOW(1), MAP_WINDOW(end)))
text(10, 10, sprintf('smoothing: S=%g', S))
text(10,15, sprintf('cliping: [%g, %g]', minValue, maxValue))
hold off

%% Draw mimg of Retino e.g. session
plotFrames = 31:2:49;
timeLabel = (plotFrames - 27)*10;

pltMaps = ones(10000, 10) * 1000;
smoothMaps = mfilt2(nanmean(vertical(:,plotFrames,:), 3), 100, 100)-1;
pltMaps(chamber, :) = smoothMaps(chamber, :);

figure;
mimg2(pltMaps, 100, 100, 0.00005, 0.0025, timeLabel, 'k'); colormap mapgeog

smoothMaps = mfilt2(nanmean(horizontal(:,plotFrames,:), 3), 100, 100)-1;
pltMaps(chamber, :) = smoothMaps(chamber, :);

figure;
mimg2(pltMaps, 100, 100, 0.00005, 0.0025, timeLabel, 'k'); colormap mapgeog