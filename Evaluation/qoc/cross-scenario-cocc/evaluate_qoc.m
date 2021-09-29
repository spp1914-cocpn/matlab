prefix = 'cross-scenario-cocc';
numRepetitions = 100;
numNcs = 1; % one real ncs and a mock
simTime = 200; % in seconds

name = 'constrate'; %constflows % %constrate-mpc

scenarioName = [prefix '-' name];
filePrefix = ['matlab/Evaluation/qoc/' prefix '/' scenarioName];
if ~contains(scenarioName, 'mpc')
    translatorPath = 'libncs_matlab/matlab/config/translators/lqr/hamming/translator_inverted_pendulum_short_lqr_hamming_10sec.mat';
else   
   translatorPath = 'libncs_matlab/matlab/config/translators/mpc/hamming/translator_inverted_pendulum_short_mpc_hamming_10sec.mat';
end
translatorData = load(translatorPath, 'translator');

% read the Omnet results
targetQocs = readtable([filePrefix '-target-qoc.csv'], 'ReadVariableNames', false);
ncsTargetQocTimes = cell2mat(cellfun(@str2double, table2array(targetQocs(2:end, 1)), 'UniformOutput', false));
ncsTargetQocTimes = ncsTargetQocTimes(~isnan(ncsTargetQocTimes));
ncsTargetQocs = table2array(targetQocs(2:numel(ncsTargetQocTimes)+1, 2));

clear targetQocs

states = datastore([filePrefix '-controller-state-1.csv'], 'ReadVariableNames', false);
allVariableNames = states.VariableNames;

% stored run-wise
% first read the times (index is odd number, starting add 1)
states.SelectedVariableNames = allVariableNames(1);    
data = readall(states);    
ncsTimes = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
ncsTimes = ncsTimes(~isnan(ncsTimes));
states.SelectedVariableNames = allVariableNames([2:2:length(states.VariableNames)]);
data = readall(states);

ncsStates =  table2array(data(2:numel(ncsTimes)+1, :));

states = datastore([filePrefix '-controller-state-2.csv'], 'ReadVariableNames', false);
states.SelectedVariableNames = allVariableNames([2:2:length(states.VariableNames)]);
data = readall(states);

ncsStates =  cat(3, ncsStates, table2array(data(2:numel(ncsTimes)+1, :)));

states = datastore([filePrefix '-controller-state-3.csv'], 'ReadVariableNames', false);
states.SelectedVariableNames = allVariableNames([2:2:length(states.VariableNames)]);
data = readall(states);

ncsStates =  cat(3, ncsStates, table2array(data(2:numel(ncsTimes)+1, :)) - pi);

states = datastore([filePrefix '-controller-state-4.csv'], 'ReadVariableNames', false);
states.SelectedVariableNames = allVariableNames([2:2:length(states.VariableNames)]);
data = readall(states);

ncsStates =  cat(3, ncsStates, table2array(data(2:numel(ncsTimes)+1, :)));

clear data states
assert(size(ncsStates, 2) == numRepetitions);

% compute the error norm directly, per time step and run

ncsErrorNorms = zeros(numel(ncsTimes), numRepetitions);
for r=1:numRepetitions
    ncsErrorNorms(:, r) = vecnorm(squeeze(ncsStates(:, r, :)), 2, 2);
end

ncsMeanErrors = mean(ncsErrorNorms, 2, 'omitnan');
ncsMedianErrors = median(ncsErrorNorms, 2, 'omitnan');
ncsUpperQuantiles = quantile(ncsErrorNorms, 0.9, 2);
ncsLowerQuantiles = quantile(ncsErrorNorms, 0.1, 2);
ncsErrorVars = var(ncsErrorNorms, 0, 2, 'omitnan');

% translate into qoc (could be improved -> first compute QoC, then take the
% mean)
ncsMeanQocs = translatorData.translator.translateControlError(ncsMeanErrors);
ncsMedianQocs = translatorData.translator.translateControlError(ncsMedianErrors);
ncsUpperQuantilesQoc = translatorData.translator.translateControlError(ncsUpperQuantiles);
ncsLowerQuantilesQoc = translatorData.translator.translateControlError(ncsLowerQuantiles);

figure('Name', scenarioName, 'NumberTitle','off');
t = tiledlayout(2, 1);
title(t, 'Ncs 1');
xlabel(t, 'Simulation Time (in seconds)');

nexttile;
hold on;
fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantiles', fliplr(ncsLowerQuantiles')], 'g', 'FaceAlpha', 0.5);
plot(ncsTimes, ncsMeanErrors, 'LineWidth', 2);
plot(ncsTimes, ncsMedianErrors, 'LineWidth', 2);

xlim([0 simTime]);
xticks([0:10:simTime]);
%ylim([0 0.002]);    
legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex', 'Location', 'southwest'); 
ylabel('Control Error');    
hold off;

nexttile;
hold on;
%fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantilesQoc', fliplr(ncsLowerQuantilesQoc')], 'g', 'FaceAlpha', 0.5);
plot(ncsTimes, ncsMeanQocs, 'LineWidth', 2);
%plot(ncsTimes, ncsMedianQocs, 'LineWidth', 2);
plot(ncsTargetQocTimes, ncsTargetQocs, 'LineWidth', 2);

xlim([0 simTime]);
xticks([0:10:simTime]);
ylim([0 1]);
ylabel('QoC'); 
legend('Mean', 'Target', 'interpreter', 'latex', 'Location', 'southwest'); 
%legend('Mean', 'Median', 'Target', 'interpreter', 'latex'); 
%legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Target', 'interpreter', 'latex');
hold off;   

% also plot the error variance
figure('Name', scenarioName, 'NumberTitle','off');
title('Ncs 1');

hold on;
plot(ncsTimes, ncsErrorVars, 'LineWidth', 2);
xlim([0 simTime]);
xticks([0:10:simTime]);
%ylim([0 0.0000002]);
ylabel('Control Error Variance');
xlabel('Simulation Time (in seconds)');
hold off;
