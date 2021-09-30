class ='cloud-scenario';
prefix = 'with-cross3';
simTime = 100; % in seconds
numRepetitions = 200;
numNcs = 1; % one real ncs and plenty of mocks
% starting from seq5, we additionally consider a bottleck line delay of 10ms
% hence, the maxMeasDelay=4
% cross-traffic starts 20s after warm-up period
types = { ...
    '-seq3-mpc', '-seq3-lqr', ...
    '-seq4-mpc', '-seq4-lqr', ...
    '-seq5-mpc', '-seq5-lqr', ...
    '-seq15-mpc', '-seq15-lqr'
    };
type = types{5};

scenarioName = [class '-' prefix type];
filePrefix = ['matlab/Evaluation/qoc/' class '/' prefix '/' type(2:end) '/' scenarioName];

if contains(scenarioName, 'mpc')    
    translatorPath = 'libncs_matlab/matlab/config/translators/mpc/hamming/translator_double_pendulum_mpc_hamming_10sec.mat';
elseif contains(scenarioName, 'robust')
    translatorPath = 'libncs_matlab/matlab/config/translators/robust/hamming/translator_inverted_pendulum_short_robust_hamming_10sec.mat';
else
    translatorPath = 'libncs_matlab/matlab/config/translators/lqr/hamming/translator_double_pendulum_lqr_hamming_10sec.mat';
end
translatorData = load(translatorPath, 'translator');

% read the Omnet results
allIdx = 1:numRepetitions;

% first get the crashed runs
plantAdmissible = tabularTextDatastore([filePrefix '-plant-admissible.csv'], 'ReadVariableNames', false);
allVariableNames = plantAdmissible.VariableNames;
plantAdmissible.SelectedVariableNames = allVariableNames(end);
data = plantAdmissible.readall();
plantStateAdmissible = table2array(data(2:end, 1));

crashedRuns = find(~plantStateAdmissible); % flag plant state admissible =0 -> run terminated prematurely

controllerAdmissible = tabularTextDatastore([filePrefix '-controller-admissible.csv'], 'ReadVariableNames', false);
allVariableNames = controllerAdmissible.VariableNames;
controllerAdmissible.SelectedVariableNames = allVariableNames(end);
data = controllerAdmissible.readall();
controllerStateAdmissible = table2array(data(2:end, 1));

crashedRuns = union(crashedRuns, find(~controllerStateAdmissible));
properRuns = setdiff(allIdx, crashedRuns); %the runs that did not crash

clear data plantAdmissible controllerAdmissible

states = tabularTextDatastore([filePrefix '-true-state-1.csv'], 'ReadVariableNames', false);
allVariableNames = states.VariableNames;

allIdx = 1:numRepetitions;

% stored run-wise
% first read the times (index is odd number, starting add 1)
% we read them from the first run that did not crash
states.SelectedVariableNames = allVariableNames(2 * properRuns(1)- 1);    
data = states.readall();     
ncsTimes = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
ncsTimes = ncsTimes(~isnan(ncsTimes));
states.SelectedVariableNames = allVariableNames([2:2:length(states.VariableNames)]);
data = states.readall();

data = table2array(data(2:numel(ncsTimes)+1, :));
ncsStates =  data(:, properRuns);

for j=2:6    
    states = tabularTextDatastore([filePrefix sprintf('-true-state-%d.csv', j)], 'ReadVariableNames', false);
    states.SelectedVariableNames = allVariableNames([2:2:length(states.VariableNames)]);
    data = states.readall();    
    data = table2array(data(2:numel(ncsTimes)+1, :));        
    
    ncsStates =  cat(3, ncsStates, data(:, properRuns));
end

clear data states
assert(size(ncsStates, 2) == numRepetitions - numel(crashedRuns));

% target QoC is a random variable in this scenario since cross traffic is random
% likewise the corresponding times are jittered
targetQocs = tabularTextDatastore([filePrefix '-target-qoc.csv'], 'ReadVariableNames', false);
allVariableNames = targetQocs.VariableNames;
targetQocs.SelectedVariableNames = allVariableNames(1:2:end); % use the times from the first run that did not crash
data = targetQocs.readall();

ncsTargetQocTimes = cell2mat(cellfun(@str2double, table2array(data(2:end, properRuns(1))), 'UniformOutput', false));
ncsTargetQocTimes = ncsTargetQocTimes(~isnan(ncsTargetQocTimes));
targetQocs.SelectedVariableNames = allVariableNames(2:2:end);
data = targetQocs.readall();
ncsTargetQocsAll = table2array(data(2:numel(ncsTargetQocTimes)+1, properRuns));

assert(size(ncsTargetQocTimes, 2) == 1);
assert(size(ncsTargetQocsAll, 2) == numRepetitions - numel(crashedRuns));

% we use the mean of the target QoCs as an indicator
ncsTargetQocs = mean(ncsTargetQocsAll, 2);

clear targetQocs data

% read the actual delays in a single run
actDelays = tabularTextDatastore([filePrefix, '-actual-delay.csv'], 'ReadVariableNames', false);
allVariableNames = actDelays.VariableNames;
actDelays.SelectedVariableNames = allVariableNames(2*properRuns(1)-1); % take the first proper run
data = actDelays.readall();
actDelayTimes = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
actDelayTimes = actDelayTimes(~isnan(actDelayTimes));
actDelays.SelectedVariableNames = allVariableNames(2 * properRuns(1));
data = actDelays.readall();
actPacketDelays = table2array(data(2:numel(actDelayTimes)+1, :)) * 1000; % in milliseconds

% read all delays
actDelays.SelectedVariableNames = allVariableNames(2:4:end);
data = actDelays.readall();
allDelays = cell2mat(table2cell(data));
allDelays = allDelays(~isnan(allDelays));

clear actDelays data

% read the control period (should be jittered)
% use a single run (the first proper run)
controlPeriods = tabularTextDatastore([filePrefix '-control-period.csv'], 'ReadVariableNames', false);
allVariableNames = controlPeriods.VariableNames;
controlPeriods.SelectedVariableNames = allVariableNames(2*properRuns(1)-1);
data = controlPeriods.readall();
controlPeriodTimes = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
controlPeriodTimes = controlPeriodTimes(~isnan(controlPeriodTimes));
controlPeriods.SelectedVariableNames = allVariableNames(2*properRuns(1));
data = controlPeriods.readall();
ncsPeriods = table2array(data(2:numel(controlPeriodTimes)+1, :));

clear controlPeriods data

% % read the observed packet delays (in time steps)
% % first read the times (index is odd number, starting add 1)
% delays = tabularTextDatastore([filePrefix, '-observed-delay.csv'], 'ReadVariableNames', false);
% allVariableNames = delays.VariableNames;
% delays.SelectedVariableNames = allVariableNames(2*properRuns(1)-1); % take the first proper run
% data = delays.readall();
% delayTimes = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
% delayTimes = delayTimes(~isnan(delayTimes));
% 
% delays.SelectedVariableNames = allVariableNames([2:2:length(delays.VariableNames)]);
% data = delays.readall();
% packetDelaysTimesteps =  table2array(data(2:numel(delayTimes)+1, :)); % in timesteps
% 
% clear delays data delayTimes

% compute the error norm directly, per time step and run
dim=3;
ncsErrorNorms = vecnorm(ncsStates, 2, dim);
ncsQocs = translatorData.translator.translateControlError(ncsErrorNorms);

ncsMeanErrors = mean(ncsErrorNorms, 2, 'omitnan');
ncsMedianErrors = median(ncsErrorNorms, 2, 'omitnan');
ncsUpperQuantiles = quantile(ncsErrorNorms, 0.9, 2);
ncsLowerQuantiles = quantile(ncsErrorNorms, 0.1, 2);
ncsErrorVars = var(ncsErrorNorms, 0, 2, 'omitnan');

% translate into qoc
ncsMeanQocs = mean(ncsQocs, 2, 'omitnan');
ncsMedianQocs = median(ncsQocs, 2, 'omitnan');
ncsUpperQuantilesQoc = quantile(ncsQocs, 0.9, 2);
ncsLowerQuantilesQoc = quantile(ncsQocs, 0.1, 2);

% norms = tabularTextDatastore([filePrefix '-state-norms.csv'], 'ReadVariableNames', false);
% norms.SelectedVariableNames = allVariableNames([2:2:length(norms.VariableNames)]);
% data = norms.readall();    
% data = table2array(data(2:numel(ncsTimes)+1, :));        
% 
% ncsNorms =  data(:, properRuns);
% clear norms data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % figure('Name', scenarioName, 'NumberTitle','off');
% % t = tiledlayout(1, 1);
% % title(t, 'Ncs 1');
% % xlabel(t, 'Simulation Time (in seconds)');
% % 
% % nexttile;
% % hold on;
% % fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantiles', fliplr(ncsLowerQuantiles')], 'g', 'FaceAlpha', 0.5);
% % plot(ncsTimes, ncsMeanErrors, 'LineWidth', 2);
% % plot(ncsTimes, ncsMedianErrors, 'LineWidth', 2);
% % 
% % xlim([0 simTime]);
% % xticks([0:10:simTime]);
% % ylim([0 0.06]);    
% % legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex', 'Location', 'southwest'); 
% % ylabel('Control Error');
% % hold off;
% % 
% % 
% % figure('Name', scenarioName, 'NumberTitle','off');
% % t = tiledlayout(1, 1);
% % title(t, 'Ncs 1');
% % xlabel(t, 'Simulation Time (in seconds)');
% % 
% % nexttile;
% % hold on;
% % %fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantiles', fliplr(ncsLowerQuantiles')], 'g', 'FaceAlpha', 0.5);
% % %plot(ncsTimes, ncsMeanErrors, 'LineWidth', 2);
% % %plot(ncsTimes, ncsMedianErrors, 'LineWidth', 2);
% % plot(ncsTargetQocTimes, ncsTargetQocsAll(:, 1), 'LineWidth', 2);
% % 
% % xlim([0 simTime]);
% % xticks([0:10:simTime]);
% % %ylim([0 0.06]);    
% % %legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex', 'Location', 'southwest'); 
% % ylabel('Target QoC');
% % hold off;
% % 
% % 
% % figure('Name', scenarioName, 'NumberTitle','off');
% % t = tiledlayout(1, 1);
% % title(t, 'Ncs 1');
% % xlabel(t, 'Simulation Time (in seconds)');
% % 
% % nexttile;
% % hold on;
% % %fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantiles', fliplr(ncsLowerQuantiles')], 'g', 'FaceAlpha', 0.5);
% % %plot(ncsTimes, ncsMeanErrors, 'LineWidth', 2);
% % %plot(ncsTimes, ncsMedianErrors, 'LineWidth', 2);
% % plot(controlPeriodTimes, ncsPeriods * 1000, 'LineWidth', 2);
% % 
% % xlim([0 simTime]);
% % xticks([0:10:simTime]);
% % %ylim([0 0.06]);    
% % %legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex', 'Location', 'southwest'); 
% % ylabel('Control Period (in milli-seconds)');
% % hold off;
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', scenarioName, 'NumberTitle','off');
t = tiledlayout(2, 1);
title(t, 'Ncs 1');
xlabel(t, 'Simulation Time (in seconds)');

nexttile;
hold on;
fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantiles', fliplr(ncsLowerQuantiles')], 'g', 'FaceAlpha', 0.5, ...
    'DisplayName', '$Q_{0.9}$ - $Q_{0.1}$');
hold on;
plot(ncsTimes, ncsMeanErrors, 'LineWidth', 2, 'DisplayName', 'Mean');
hold on;
plot(ncsTimes, ncsMedianErrors, 'LineWidth', 2, 'DisplayName', 'Median');

xlim([0 simTime]);
xticks([0:10:simTime]);
%ylim([0 0.002]);    
%legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex', 'Location', 'southwest'); 
legend('interpreter', 'latex', 'Location', 'southwest');
ylabel('Control Error');
hold off;

nexttile;
hold on;
fill([ncsTimes', fliplr(ncsTimes')], [ncsUpperQuantilesQoc', fliplr(ncsLowerQuantilesQoc')], 'g', 'FaceAlpha', 0.5, ...
    'DisplayName', '$Q_{0.9}$ - $Q_{0.1}$');
hold on
plot(ncsTimes, ncsMeanQocs, 'LineWidth', 2, 'DisplayName', 'Mean');
hold on
plot(ncsTimes, ncsMedianQocs, 'LineWidth', 2, 'DisplayName', 'Median');
hold on
plot(ncsTargetQocTimes, ncsTargetQocs, 'LineWidth', 2, 'DisplayName', 'Target');

xlim([0 simTime]);
xticks([0:10:simTime]);
ylim([0 1]);
ylabel('QoC'); 
legend('interpreter', 'latex', 'Location', 'southwest');
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
%ylim([0 20]);
%yticks([0:2:20]);
ylabel('Packet Delay (in milliseconds)');
xlabel('Simulation Time (in seconds)');
hold off;

% and plot the actual packet delays during single run, starting from 50s
figure('Name', scenarioName, 'NumberTitle', 'off');
title('Actual Packet Delays (Single Run, 50-100s)');

hold on;
scatter(actDelayTimes((actDelayTimes >= 50)), actPacketDelays(actDelayTimes >= 50), 'LineWidth', 2);
xlim([50 simTime]);
xticks([50:10:simTime]);
%ylim([0 20]);
%yticks([0:2:20]);
ylabel('Packet Delay (in milliseconds)');
xlabel('Simulation Time (in seconds)');
hold off;
% 
% plot a histogram of the actual packet delays in that run until 50s
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (Single Run, 0-50s)');
hold on;
histogram(actPacketDelays(actDelayTimes < 50));
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

% plot a histogram of the actual packet delays in that after 50s
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (Single Run, 50-100s)');
hold on;
histogram(actPacketDelays(actDelayTimes >= 50));
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

% plot a histogram of the actual packet delays in that run
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (Single Run)');
hold on;
histogram(actPacketDelays);
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

% plot a histogram of all actual packet delays in all runs
figure('Name', scenarioName, 'NumberTitle','off');
title('Actual Packet Delays - Histogram (All runs)');
hold on;
histogram(allDelays * 1000);
ylabel('Count');
xlabel('Packet Delay (in milliseconds)');
hold off;

% % % plot a histogram of all observed packet delays in all runs
% % figure('Name', scenarioName, 'NumberTitle','off');
% % title('Observed Packet Delays - Histogram (Single run)');
% % hold on;
% % histogram(packetDelaysTimesteps(:, 1), 'Normalization', 'probability');
% % ylabel('Probability');
% % xlabel('Observed Delay (in time steps)');
% % hold off;
% % 
% % 
% % % plot a histogram of all observed packet delays in all runs
% % figure('Name', scenarioName, 'NumberTitle','off');
% % title('Observed Packet Delays - Histogram (All runs)');
% % hold on;
% % histogram(packetDelaysTimesteps, 'Normalization', 'probability');
% % ylabel('Probability');
% % xlabel('Observed Delay (in time steps)');
% % hold off;
