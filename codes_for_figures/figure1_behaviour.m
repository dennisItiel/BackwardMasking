%%% Calculate statiscs of animal's behaviour, all sessions and SOA's
% namely performance and reaction time, with their relationship

monkey = 'frodo'; %'tolkin'
cond_labels = {'retino' 'noMaks' 'soa150' 'soa100' 'soa75' 'soa50'};
pth = 'C:\Users\itiel\Desktop\BM_article_project\behaviour_files\';

switch monkey
    case 'frodo'
        fileDates = {'1609' '2309' '710' '1410'};
        load("E:\Itiel\BM\Frodo\VSDI\datesConsCell.mat")
    case 'tolkin'
        fileDates = {'2405' '3005' '2706' '0407' '1107' '1807'};
        load("E:\Itiel\BM\Tolkin\VSDI\datesConsCell.mat");
end

corrRT = cell(4,1);
errRT = cell(4,1);
allPerf = cell(4,1);

for SOA_inx=6:-1:3
    SOA = cond_labels{SOA_inx};

    soaCorrRT = [];
    soaErrRT = [];
    for date_inx=1:numel(fileDates)
        letter = datesConsCell{SOA_inx, date_inx};
        if letter == "N"
            continue
        end
    
        tbl = importdata(strcat(pth, monkey,'\', SOA, '\', fileDates{date_inx}, letter, '.xlsx'));
        tbl = tbl.data;
        tbl = tbl(~isnan(tbl(:,4)) & ~isnan(tbl(:,5)), :);
    
        singleCorrRT = tbl(tbl(:,5) == 1, 7);
        singleErrRT = tbl(tbl(:,5) == 0, 7);
        
        soaCorrRT = [soaCorrRT; singleCorrRT];
        soaErrRT = [soaErrRT; singleErrRT];

        correctsVer = sum(tbl(:,2) == 2 & tbl(:,5) == 1);
        correctsHor = sum(tbl(:,2) == 4 & tbl(:,5) == 1);
        errorsVer = sum(tbl(:,2) == 2 & tbl(:,5) == 0);
        errorsHor = sum(tbl(:,2) == 4 & tbl(:,5) == 0);
        
        allPerf{SOA_inx-2} = cat(1, allPerf{SOA_inx-2}, correctsVer/(correctsVer+errorsVer));
        allPerf{SOA_inx-2} = cat(1, allPerf{SOA_inx-2}, correctsHor/(correctsHor+errorsHor));
    end
    corrRT{SOA_inx-2} = soaCorrRT;
    errRT{SOA_inx-2} = soaErrRT;
end

%% plot RT

x = [50 75 100 150];
rts = flipud(cellfun(@mean, corrRT));
rtSEM = flipud(cellfun(@(x) std(x)/sqrt(numel(x)), corrRT));


% Colors palette
colors = [
    [200 100 200]/255;
    [160 60 160]/255;
    [120 20 120]/255;
    [80 0 80]/255;
];

figure; hold on
errorbar(x, rts, rtSEM, 'k')
plot(x, rts, 'LineWidth', 2, 'color', [0.7 0.7 0.7])
for i=1:4
    plot(x(i), rts(i), '.', 'MarkerSize', 30, 'Color', colors(i, :))
end

% ylim([0.6, 1])
xlim([40, 160])
xticks([50 75 100 150])
ylabel('RT + SEM (ms)')
xlabel('SOA')
title(sprintf('%s RT, 4-sessions cross-trials average', monkey))

% Rt ranksums
sz = [3 4];
varTypes = ["string","double","double","double"];
varNames = ["soa_labels" "SOA150" "SOA100" "SOA75"];
RT_PV = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames', varNames);
RT_PV.soa_labels = ["SOA100"; "SOA75"; "SOA50"];

for inx_a=1:3
    for inx_b=inx_a+1:4
        RT_PV{inx_b-1, inx_a+1} = ranksum(corrRT{inx_a}, corrRT{inx_b});
    end
end


% Get pValues
p50_75 = RT_PV{3, 4};
p50_100 = RT_PV{3, 3};
p50_150 = RT_PV{3, 2};
p75_100 = RT_PV{2, 3};
p75_150 = RT_PV{2, 2};
p100_150 = RT_PV{1, 2};

% Plot the stars
groups = {[50 75] [50 100] [50 150] [75 100] [75 150] [100 150]};
sigstar(groups,[(p50_75), (p50_100), (p50_150), (p75_100), (p75_150), (p100_150)], 1);

%% plot Performance

x = [50 75 100 150];
perfs = flipud(cellfun(@mean, allPerf));
perfSEM = flipud(cellfun(@(x) std(x)/sqrt(numel(x)), allPerf));

% Colors palette
colors = [
    [200 100 200]/255;
    [160 60 160]/255;
    [120 20 120]/255;
    [80 0 80]/255;
];

figure; hold on
errorbar(x, perfs, perfSEM, 'k')
plot(x, perfs, 'LineWidth', 2, 'color', [0.7 0.7 0.7])
for i=1:4
    plot(x(i), perfs(i), '.', 'MarkerSize', 30, 'Color', colors(i, :))
end
ylim([0.6, 1])
xlim([40, 160])
xticks([50 75 100 150])
ylabel('performance + SEM (ms)')
xlabel('SOA')
title(sprintf('%s performance, 4-sessions cross-sessions average', monkey))


% Performance ranksums
sz = [3 4];
varTypes = ["string","double","double","double"];
varNames = ["soa_labels" "SOA150" "SOA100" "SOA75"];
perfs_PV = table('Size',sz, 'VariableTypes',varTypes, 'VariableNames', varNames);
perfs_PV.soa_labels = ["SOA100"; "SOA75"; "SOA50"];

for inx_a=1:3
    for inx_b=inx_a+1:4
        perfs_PV{inx_b-1, inx_a+1} = ranksum(allPerf{inx_a}, allPerf{inx_b});
    end
end

% disp(perfs_PV)

% Get pValues
p50_75 = perfs_PV{3, 4};
p50_100 = perfs_PV{3, 3};
p50_150 = perfs_PV{3, 2};
p75_100 = perfs_PV{2, 3};
p75_150 = perfs_PV{2, 2};
p100_150 = perfs_PV{1, 2};

% Plot the stars
groups = {[50 75] [50 100] [50 150] [75 100] [75 150] [100 150]};
sigstar(groups,[(p50_75), (p50_100), (p50_150), (p75_100), (p75_150), (p100_150)], 1);



x = [50 75 100 150];
% perfsBoth = [];
% perfSEM = std(perfsBoth, 0, 2)/sqrt(2);

figure; hold on
errorbar(x, mean(perfs, 2), perfSEM)
plot(x, perfs, 'x', 'LineWidth', 3, 'Color', 'b')
ylim([0.6, 1])
xlim([40, 160])
xticks([50 75 100 150])
ylabel('performance + SEM (ms)')
xlabel('SOA')
title('Frodo performance, 4-sessions cross-sessions average')

%% Error vs Correct RT

meanCorrectRT = cellfun(@nanmean, corrRT);
semCorrectRT = cellfun(@(x) nanstd(x)/sqrt(numel(x)), corrRT);
meanErrorRT = cellfun(@nanmean, errRT)';
semErrorRT = cellfun(@(x) nanstd(x)/sqrt(numel(x)), errRT)';

% Colors palette
colors = [
    [80 0 80]/255;
    [120 20 120]/255;
    [160 60 160]/255;
    [200 100 200]/255;
];

pValues = [];

figure; hold on
switch monkey
    case 'frodo'
        plot(80:260, 80:260, 'color', [0.5 0.5 0.5], 'LineStyle', '-')
    case 'tolkin'
        plot(180:300, 180:300, 'color', [0.5 0.5 0.5], 'LineStyle', '-')
end


for k = 1:4
    semX = semCorrectRT(k);
    semY = semErrorRT(k);
    
    errorbar(meanCorrectRT(k), meanErrorRT(k), semY, semY, semX, semX, '.', 'MarkerSize', 30, 'Color', 'k')
    plot(meanCorrectRT(k), meanErrorRT(k), '.', 'MarkerSize', 30, 'Color', colors(k, :))

    pValues(k) = ranksum(corrRT{k}, errRT{k});
end


% Plot the stars
groups = {[meanCorrectRT(1)-5 meanCorrectRT(1)+5],...
          [meanCorrectRT(2)-5 meanCorrectRT(2)+5],...
          [meanCorrectRT(3)-5 meanCorrectRT(3)+5],...
          [meanCorrectRT(4)-5 meanCorrectRT(4)+5]};
sigstar(groups,[(pValues(1)), (pValues(2)), (pValues(3)), (pValues(4))]);

xlabel('correct trials RT')
ylabel('error trials RT')
legend({'' '' 'SOA150' '' 'SOA100' '' 'SOA75' '' 'SOA50'})
title('Reaction Time of Correct vs. Error trials')

% Pool all SOAs together, correct vs. error
pooledCorr = [];
pooledErr = [];
for i=1:4
    pooledCorr = [pooledCorr; corrRT{i}];
    pooledErr = [pooledErr; errRT{i}];
end
diffCorrErr = nanmean(pooledErr) - nanmean(pooledCorr);
P = ranksum(pooledErr, pooledCorr);

%% Accuracy rate vs. Reaction Time

% Colors palette
colors = [
    [80 0 80]/255;
    [120 20 120]/255;
    [160 60 160]/255;
    [200 100 200]/255;
];

minThreshold = 100;
maxThreshold = 300;
binSize = 40;
minTrials = 5;

edges = minThreshold:binSize:1000;
numDataPoints = numel(edges)-1;

% Preallocate
allBinnedPerf = nan(4, numDataPoints);
allBinnedSEM = nan(4, numDataPoints);
pooledPerf = cell(numDataPoints,1);
trialsPerBin = nan(4, numDataPoints);

for soaInx=1:4
    soa_RT = cat(1, corrRT{soaInx}, errRT{soaInx});
    soa_resp = cat(1, ones(numel(corrRT{soaInx}), 1), zeros(numel(errRT{soaInx}), 1));
    
    % Sort RT and respones by ascending order
    [sortedRT,I] = sort(soa_RT);
    sortedResp = soa_resp(I);

    % Threshold by min RT
    minInx = find(sortedRT < minThreshold, 1, 'last') +1;
    maxInx = find(sortedRT >= maxThreshold, 1, 'first') -1;
    sortedRT = sortedRT(minInx:maxInx);
    sortedResp = sortedResp(minInx:maxInx);
    
    % Count into bins
    [N,edges,bin] = histcounts(sortedRT, edges);
    trialsPerBin(soaInx,:) = N;
    
    % Mean and SEM of %correct
    startIdx = 1;
    for binInx = 1:numDataPoints
        if N(binInx) < minTrials
            continue
        end
        responses = sortedResp(startIdx:startIdx + N(binInx) - 1);
        allBinnedPerf(soaInx, binInx) = nanmean(responses);
        allBinnedSEM(soaInx, binInx) = nanstd(responses) ./ sqrt(numel(responses));
        pooledPerf{binInx} = [pooledPerf{binInx}; responses];
        startIdx = startIdx + N(binInx);
    end
end

% important stats
nanIdx = find(cellfun(@isempty, pooledPerf), 1, 'first');

meanPooledPerf = cellfun(@mean, pooledPerf, 'UniformOutput', 1);
sem = @(x) std(x)./sqrt(numel(x));
semPooledPerf = cellfun(sem, pooledPerf, 'UniformOutput', 1);
P_pooled = ranksum(pooledPerf{1}, pooledPerf{nanIdx-1});
numTrialsTotal = sum(trialsPerBin, 2)';

pltBins = edges(1:numDataPoints) + 0.5*binSize;
figure; hold on
for soaInx=1:4
    shadedErrorBar(pltBins, allBinnedPerf(soaInx, :), allBinnedSEM(soaInx, :),...
        {'color', colors(soaInx, :), 'linewidth', 2}, 0);
end
plot(pltBins, meanPooledPerf)
yline(0.5, 'k', 'Linestyle', '--')
legend({'' '150' '' '100' '' '75' '' '50'})
xlabel('Binned RT length (ms)')
ylabel('Percent correct')
ylim([0.4 1])
xlim([minThreshold maxThreshold])
xticks(edges(1:numDataPoints))
title('Fraction correct versus reaction time')
