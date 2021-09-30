class ='cloud-scenario';
prefix = 'with-cross';
simTime = 100; % in seconds
numRepetitions = 20%100;
numNcs = 1; % one real ncs and plenty of mocks

names = {'seq2', 'seq3', 'seq4', 'seq6'}; % indicate the sequence length in use
types = {'-mpc', '-lqr'}; % indicate the controller in use
name = names{2};
type = types{2};
%type='';
%name = 'control-increase-mpc';
scenarioName = [class '-' prefix type '-' name];
%scenarioName = [prefix type];
filePrefix = ['matlab/Evaluation/qoc/' class '/' prefix '/' scenarioName];
if contains(scenarioName, 'mpc')    
    translatorPath = 'libncs_matlab/matlab/config/translators/mpc/hamming/translator_inverted_pendulum_short_mpc_hamming_10sec.mat';
elseif contains(scenarioName, 'robust')
   translatorPath = 'libncs_matlab/matlab/config/translators/robust/hamming/translator_inverted_pendulum_short_robust_hamming_10sec.mat';
else
    translatorPath = 'libncs_matlab/matlab/config/translators/lqr/hamming/translator_inverted_pendulum_short_lqr_hamming_10sec.mat';
end
translatorData = load(translatorPath, 'translator');

% read the Omnet results
targetQocs = readtable([filePrefix '-target-qoc.csv'], 'ReadVariableNames', false);
ncsTargetQocTimes = cell2mat(cellfun(@str2double, table2array(targetQocs(2:end, 1)), 'UniformOutput', false));
ncsTargetQocTimes = ncsTargetQocTimes(~isnan(ncsTargetQocTimes));
ncsTargetQocs = table2array(targetQocs(2:numel(ncsTargetQocTimes)+1, 2));

clear targetQocs

% read the actual delays in a single run
actDelays = datastore([filePrefix, '-actual-delay.csv'], 'ReadVariableNames', false);
allVariableNames = actDelays.VariableNames;
actDelays.SelectedVariableNames = allVariableNames(1);
data = readall(actDelays);
actDelayTimes = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
actDelayTimes = actDelayTimes(~isnan(actDelayTimes));
actDelays.SelectedVariableNames = allVariableNames(2);
data = readall(actDelays);
actPacketDelays =  table2array(data(2:numel(actDelayTimes)+1, :))* 1000; % in milliseconds

% read all delays
actDelays.SelectedVariableNames = allVariableNames(2:4:end);
data = readall(actDelays);
allDelays = cell2mat(table2cell(data));
allDelays = allDelays(~isnan(allDelays));

clear actDelays data

% read the control period
controlPeriods = readtable([filePrefix '-control-period.csv'], 'ReadVariableNames', false);
controlPeriodTimes = cell2mat(cellfun(@str2double, table2array(controlPeriods(2:end, 1)), 'UniformOutput', false));
controlPeriodTimes = controlPeriodTimes(~isnan(controlPeriodTimes));
ncsPeriods = table2array(controlPeriods(2:numel(controlPeriodTimes)+1, 2));

clear controlPeriods

% read the observed packet delays
% first read the times (index is odd number, starting add 1)
delays = datastore([filePrefix, '-observed-delay.csv'], 'ReadVariableNames', false);
allVariableNames = delays.VariableNames;
delays.SelectedVariableNames = allVariableNames(1);   
data = readall(delays);
delayTimes = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
delayTimes = delayTimes(~isnan(delayTimes));

delays.SelectedVariableNames = allVariableNames([2:2:length(delays.VariableNames)]);
data = readall(delays);
packetDelays =  table2array(data(2:numel(delayTimes)+1, :))* 1000; % in milliseconds

packetDelaysTimesteps = zeros(size(packetDelays, 1), numRepetitions);
for j=1:size(packetDelays, 1)
    idx = find(controlPeriodTimes(controlPeriodTimes <= delayTimes(j)), 1, 'last');
    packetDelaysTimesteps(j, :) = packetDelays(j, :) / 1000 / ncsPeriods(idx);
end

clear delays

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
ylim([0 0.002]);    
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

% and plot the actual packet delays during single run
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays (Single Run)');

hold on;
plot(actDelayTimes, actPacketDelays, 'LineWidth', 2);
xlim([0 simTime]);
xticks([0:10:simTime]);    
ylim([0 20]);
yticks([0:2:20]);
ylabel('Packet Delay (in milliseconds)');
xlabel('Simulation Time (in seconds)');
hold off;

% and plot the actual packet delays during single run, starting from 50s
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays (Single Run, 50-100s)');

hold on;
scatter(actDelayTimes((actDelayTimes >= 50)), actPacketDelays(actDelayTimes >= 50), 'LineWidth', 2);
xlim([50 simTime]);
xticks([50:10:simTime]);
ylim([0 20]);
yticks([0:2:20]);
ylabel('Packet Delay (in milliseconds)');
xlabel('Simulation Time (in seconds)');
hold off;
% 
% plot a histogram of the actual packet delays in that run until 50s
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (Single Run, 0-50s)');
hold on;
histogram(actPacketDelays(actDelayTimes < 50), 0:2:20);
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

 % plot a histogram of the actual packet delays in that after 50s
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (Single Run, 50-100s)');
hold on;
histogram(actPacketDelays(actDelayTimes >= 50), 0:2:20);
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

% plot a histogram of the actual packet delays in that run
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (Single Run)');
hold on;
histogram(actPacketDelays, 0:2:20);
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

% plot a histogram of all actual packet delays in all runs
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (All runs)');
hold on;
histogram(allDelays * 1000, 0:2:20);
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

% plot a histogram of all observed packet delays in all runs
figure('Name', scenarioName, 'NumberTitle','off');
title('Observed Packet Delays - Histogram (Single run)');
hold on;
histogram(packetDelaysTimesteps(:, 1), 0:1:6, 'Normalization', 'probability');
ylabel('Probability');
xlabel('Observed Delay (in time steps)');
hold off;


% plot a histogram of all observed packet delays in all runs
figure('Name', scenarioName, 'NumberTitle','off');
title('Observed Packet Delays - Histogram (All runs)');
hold on;
histogram(packetDelaysTimesteps, 0:1:6, 'Normalization', 'probability');
ylabel('Probability');
xlabel('Observed Delay (in time steps)');
hold off;
