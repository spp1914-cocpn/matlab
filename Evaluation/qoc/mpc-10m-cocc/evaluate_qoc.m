numRepetitions = 200;
simTime = 300; % in seconds
numNcs = 2;

prefix = 'mpc-10m-cocc';
%scenarioName = 'mpc-10m-cocc';
scenarioName = 'mpc-10m-cocc-decrease-increase';
filePrefix = ['matlab/Evaluation/qoc/' prefix '/' scenarioName];
translatorPath = 'libncs_matlab/matlab/config/translators/mpc/hamming/translator_inverted_pendulum_short_mpc_hamming_10sec.mat';
translatorData = load(translatorPath, 'translator');

% read the Omnet results
targetQocs = readtable([filePrefix '-target-qoc.csv'], 'ReadVariableNames', false);
ncs1TargetQocTimes = cell2mat(cellfun(@str2double, table2array(targetQocs(2:end, 1)), 'UniformOutput', false));
ncs1TargetQocTimes = ncs1TargetQocTimes(~isnan(ncs1TargetQocTimes));
ncs1TargetQocs = table2array(targetQocs(2:numel(ncs1TargetQocTimes)+1, 2));

ncs2TargetQocTimes = cell2mat(cellfun(@str2double, table2array(targetQocs(2:end, 3)), 'UniformOutput', false));
ncs2TargetQocTimes = ncs2TargetQocTimes(~isnan(ncs2TargetQocTimes));
ncs2TargetQocs = table2array(targetQocs(2:numel(ncs2TargetQocTimes)+1, 4));

clear targetQocs

states = datastore([filePrefix '-controller-state-1.csv'], 'ReadVariableNames', false);
allVariableNames = states.VariableNames;

% stored run-wise
% first read the times (index is odd number, starting add 1)
states.SelectedVariableNames = allVariableNames(1);    
data = readall(states);    
ncs1Times = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
ncs1Times = ncs1Times(~isnan(ncs1Times));
states.SelectedVariableNames = allVariableNames(3);    
data = readall(states);    
ncs2Times = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
ncs2Times = ncs2Times(~isnan(ncs2Times));

%now the controller states
states.SelectedVariableNames = allVariableNames([2:4:length(states.VariableNames)]);
data = readall(states);
ncs1States =  table2array(data(2:numel(ncs1Times)+1, :));
states.SelectedVariableNames = allVariableNames([4:4:length(states.VariableNames)]);
data = readall(states);
ncs2States =  table2array(data(2:numel(ncs2Times)+1, :));

states = datastore([filePrefix '-controller-state-2.csv'], 'ReadVariableNames', false);
states.SelectedVariableNames = allVariableNames([2:4:length(states.VariableNames)]);
data = readall(states);
ncs1States =  cat(3, ncs1States, table2array(data(2:numel(ncs1Times)+1, :)));
states.SelectedVariableNames = allVariableNames([4:4:length(states.VariableNames)]);
data = readall(states);
ncs2States =  cat(3, ncs2States, table2array(data(2:numel(ncs2Times)+1, :)));

states = datastore([filePrefix '-controller-state-3.csv'], 'ReadVariableNames', false);
states.SelectedVariableNames = allVariableNames([2:4:length(states.VariableNames)]);
data = readall(states);
ncs1States =  cat(3, ncs1States, table2array(data(2:numel(ncs1Times)+1, :)) -pi);
states.SelectedVariableNames = allVariableNames([4:4:length(states.VariableNames)]);
data = readall(states);
ncs2States =  cat(3, ncs2States, table2array(data(2:numel(ncs2Times)+1, :)) -pi);

states = datastore([filePrefix '-controller-state-4.csv'], 'ReadVariableNames', false);
states.SelectedVariableNames = allVariableNames([2:4:length(states.VariableNames)]);
data = readall(states);
ncs1States =  cat(3, ncs1States, table2array(data(2:numel(ncs1Times)+1, :)));
states.SelectedVariableNames = allVariableNames([4:4:length(states.VariableNames)]);
data = readall(states);
ncs2States =  cat(3, ncs2States, table2array(data(2:numel(ncs2Times)+1, :)));

assert(size(ncs1States, 2) == numRepetitions);
assert(size(ncs2States, 2) == numRepetitions);

ncs1ErrorNorms = zeros(numel(ncs1Times), numRepetitions);
ncs2ErrorNorms = zeros(numel(ncs2Times), numRepetitions);
ncs1Qocs = zeros(numel(ncs1Times), numRepetitions);
ncs2Qocs = zeros(numel(ncs2Times), numRepetitions);
for r=1:numRepetitions
    ncs1ErrorNorms(:, r) = vecnorm(squeeze(ncs1States(:, r, :)), 2, 2);
    ncs2ErrorNorms(:, r) = vecnorm(squeeze(ncs2States(:, r, :)), 2, 2);
    
%     ncs1Qocs(:, r) = translatorData.translator.translateControlError(ncs1ErrorNorms(:, r));
%     ncs2Qocs(:, r) = translatorData.translator.translateControlError(ncs2ErrorNorms(:, r));
end
% 
ncs1MeanErrors = mean(ncs1ErrorNorms, 2, 'omitnan');
ncs1MedianErrors = median(ncs1ErrorNorms, 2, 'omitnan');
ncs1UpperQuantiles = quantile(ncs1ErrorNorms, 0.9, 2);
ncs1LowerQuantiles = quantile(ncs1ErrorNorms, 0.1, 2);
ncs1ErrorVars = var(ncs1ErrorNorms, 0, 2, 'omitnan');
ncs2MeanErrors = mean(ncs2ErrorNorms, 2, 'omitnan');
ncs2MedianErrors = median(ncs2ErrorNorms, 2, 'omitnan');
ncs2UpperQuantiles = quantile(ncs2ErrorNorms, 0.9, 2);
ncs2LowerQuantiles = quantile(ncs2ErrorNorms, 0.1, 2);
ncs2ErrorVars = var(ncs2ErrorNorms, 0, 2, 'omitnan');

% translate into qoc (could be improved -> first compute QoC, then take the
% mean)
ncs1MeanQocs = translatorData.translator.translateControlError(ncs1MeanErrors);
ncs1MedianQocs = translatorData.translator.translateControlError(ncs1MedianErrors);
ncs1UpperQuantilesQoc = translatorData.translator.translateControlError(ncs1UpperQuantiles);
ncs1LowerQuantilesQoc = translatorData.translator.translateControlError(ncs1LowerQuantiles);
ncs2MeanQocs = translatorData.translator.translateControlError(ncs2MeanErrors);
ncs2MedianQocs = translatorData.translator.translateControlError(ncs2MedianErrors);
ncs2UpperQuantilesQoc = translatorData.translator.translateControlError(ncs2UpperQuantiles);
ncs2LowerQuantilesQoc = translatorData.translator.translateControlError(ncs2LowerQuantiles);

figure('Name', scenarioName, 'NumberTitle','off');
t = tiledlayout(2, 1);
title(t, 'Ncs 1');
xlabel(t, 'Simulation Time (in seconds)');

nexttile;
hold on;
fill([ncs1Times', fliplr(ncs1Times')], [ncs1UpperQuantiles', fliplr(ncs1LowerQuantiles')], 'g', 'FaceAlpha', 0.5);
plot(ncs1Times, ncs1MeanErrors, 'LineWidth', 2);
plot(ncs1Times, ncs1MedianErrors, 'LineWidth', 2);

xlim([0 simTime]);
xticks([0:10:simTime]);
%ylim([0 0.005]);    
legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex', 'Location', 'southwest'); 
ylabel('Control Error');    
hold off;

nexttile;
hold on;
%fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantilesQoc', fliplr(ncsLowerQuantilesQoc')], 'g', 'FaceAlpha', 0.5);
plot(ncs1Times, ncs1MeanQocs, 'LineWidth', 2);
plot(ncs1Times, ncs1MedianQocs, 'LineWidth', 2);
plot(ncs1TargetQocTimes, ncs1TargetQocs, 'LineWidth', 2);

xlim([0 simTime]);
xticks([0:10:simTime]);
ylim([0 1]);
ylabel('QoC'); 
%legend('Mean', 'Target', 'interpreter', 'latex', 'Location', 'southwest'); 
legend('Mean', 'Median', 'Target', 'interpreter', 'latex'); 
%legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Target', 'interpreter', 'latex');
hold off; 

%
%
figure('Name', scenarioName, 'NumberTitle','off');
t = tiledlayout(2, 1);
title(t, 'Ncs 2');
xlabel(t, 'Simulation Time (in seconds)');

nexttile;
hold on;
fill([ncs2Times', fliplr(ncs2Times')], [ncs2UpperQuantiles', fliplr(ncs2LowerQuantiles')], 'g', 'FaceAlpha', 0.5);
plot(ncs2Times, ncs2MeanErrors, 'LineWidth', 2);
plot(ncs2Times, ncs2MedianErrors, 'LineWidth', 2);

xlim([0 simTime]);
xticks([0:10:simTime]);
%ylim([0 0.005]);    
legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex', 'Location', 'southwest'); 
ylabel('Control Error');    
hold off;

nexttile;
hold on;
%fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantilesQoc', fliplr(ncsLowerQuantilesQoc')], 'g', 'FaceAlpha', 0.5);
plot(ncs2Times, ncs2MeanQocs, 'LineWidth', 2);
plot(ncs2Times, ncs2MedianQocs, 'LineWidth', 2);
plot(ncs2TargetQocTimes, ncs2TargetQocs, 'LineWidth', 2);

xlim([0 simTime]);
xticks([0:10:simTime]);
ylim([0 1]);
ylabel('QoC'); 
%legend('Mean', 'Target', 'interpreter', 'latex', 'Location', 'southwest'); 
legend('Mean', 'Median', 'Target', 'interpreter', 'latex'); 
%legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Target', 'interpreter', 'latex');
hold off;   

% also plot the error variance
% also plot the error variance over time
figure('Name', 'Error Variances', 'NumberTitle','off');
t = tiledlayout(2,1);
title(t, 'Error Variance over Time');
xlabel(t, 'Simulation time (in seconds)');
ylabel(t, 'Error Variance');

nexttile;
hold on;
plot(ncs1Times, ncs1ErrorVars, 'LineWidth', 2);
xticks([0:10:simTime]);
xlim([0 simTime]);
ylim([0 0.00002]);
title('Ncs 1');
hold off;

nexttile;
hold on;
plot(ncs2Times, ncs2ErrorVars, 'LineWidth', 2);
xticks([0:10:simTime]);
xlim([0 simTime]);
ylim([0 0.00002]);
title('Ncs 2');
hold off;



