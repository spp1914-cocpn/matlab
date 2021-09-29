numRepetitions = 100;
simTime = 300; % in seconds
startOffset = 30;

prefix = 'matlab/Evaluation/qoc/lqr-10m-cocc/';
translatorPath = 'libncs_matlab/matlab/config/translators/lqr/hamming/translator_inverted_pendulum_short_lqr_hamming_10sec.mat';

% % TODO: translate into QoC (all the curves)
translatorData = load(translatorPath, 'translator');

% target Qoc is the same for all runs
targetQocs =  readtable([prefix 'lqr-10m-cocc-target-qoc.csv'], 'ReadVariableNames', false);
ncs1TargetQocTimes = cell2mat(cellfun(@str2double, table2array(targetQocs(2:end, 1)), 'UniformOutput', false));
ncs2TargetQocTimes = cell2mat(cellfun(@str2double, table2array(targetQocs(2:end, 3)), 'UniformOutput', false));
ncs2TargetQocTimes = ncs2TargetQocTimes(~isnan(ncs2TargetQocTimes));

ncs1TargetQoc = table2array(targetQocs(2:end, 2));
ncs2TargetQoc = table2array(targetQocs(2:numel(ncs2TargetQocTimes)+1, 4));

% % read the Omnet results
errors = readtable([prefix 'lqr-10m-cocc-actual-error.csv'], 'ReadVariableNames', false);

ncs1Times = cell2mat(cellfun(@str2double, table2array(errors(2:end, 1)), 'UniformOutput', false));
ncs2Times = cell2mat(cellfun(@str2double, table2array(errors(2:end, 3)), 'UniformOutput', false));
ncs2Times = ncs2Times(~isnan(ncs2Times));
% 
% ncs1Errors = table2array(errors(2:end, 2:4:end));
% ncs2Errors = table2array(errors(2:numel(ncs2Times)+1, 4:4:end));
% 
% % compute the mean error and the variances, and quantiles
% ncs1MeanError = mean(ncs1Errors, 2);
% ncs1MedianError = median(ncs1Errors, 2);
% ncs2MeanError = mean(ncs2Errors, 2);
% ncs2MedianError = median(ncs2Errors, 2);
% 
% ncs1ErrorVar = var(ncs1Errors, 0, 2);
% ncs2ErrorVar = var(ncs2Errors, 0, 2);
% 
% ncs1UpperQuantile = quantile(ncs1Errors, 0.9, 2);
% ncs1LowerQuantile = quantile(ncs1Errors, 0.1, 2);
% ncs2UpperQuantile = quantile(ncs2Errors, 0.9, 2);
% ncs2LowerQuantile = quantile(ncs2Errors, 0.1, 2);
% 

% 
% ncs1MeanQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1MeanError);
% ncs1MedianQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1MedianError);
% ncs1UpperQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1UpperQuantile);
% ncs1LowerQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1LowerQuantile);
% 
% ncs2MeanQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2MeanError);
% ncs2MedianQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2MedianError);
% ncs2UpperQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2UpperQuantile);
% ncs2LowerQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2LowerQuantile);


% controller estimates of the true plant state
states = readtable([prefix 'lqr-10m-cocc-controller-state-1.csv'], 'ReadVariableNames', false);
ncs1States = table2array(states(2:end, 2:4:end));
ncs2States = table2array(states(2:numel(ncs2Times)+1, 4:4:end));

states = readtable([prefix 'lqr-10m-cocc-controller-state-2.csv'], 'ReadVariableNames', false);
ncs1States = cat(3, ncs1States, table2array(states(2:end, 2:4:end)));
ncs2States = cat(3, ncs2States, table2array(states(2:numel(ncs2Times)+1, 4:4:end)));

states = readtable([prefix 'lqr-10m-cocc-controller-state-3.csv'], 'ReadVariableNames', false);
ncs1States = cat(3, ncs1States, table2array(states(2:end, 2:4:end))-pi);
ncs2States = cat(3, ncs2States, table2array(states(2:numel(ncs2Times)+1, 4:4:end))-pi);

states = readtable([prefix 'lqr-10m-cocc-controller-state-4.csv'], 'ReadVariableNames', false);
ncs1States = cat(3, ncs1States, table2array(states(2:end, 2:4:end)));
ncs2States = cat(3, ncs2States, table2array(states(2:numel(ncs2Times)+1, 4:4:end)));

clear states errors targetQocs

assert(size(ncs1States, 2) == numRepetitions);
assert(size(ncs2States, 2) == numRepetitions);


% ncs1Times = cell2mat(cellfun(@str2double, table2array(controllerStates{1}(2:end, 1)), 'UniformOutput', false));
% ncs2Times = cell2mat(cellfun(@str2double, table2array(controllerStates{1}(2:end, 3)), 'UniformOutput', false));
% ncs2Times = ncs2Times(~isnan(ncs2Times));


% compute the error norm directly, per time step and run
ncs1ErrorNorms = zeros(size(ncs1States, 1),  numRepetitions);
ncs2ErrorNorms = zeros(size(ncs2States, 1),  numRepetitions);
ncs1Qocs = zeros(size(ncs1States, 1),  numRepetitions);
ncs2Qocs = zeros(size(ncs2States, 1),  numRepetitions);

for j=1:numRepetitions    
    ncs1ErrorNorms(:, j) = vecnorm(squeeze(ncs1States(:, j, :)), 2, 2);
    ncs2ErrorNorms(:, j) = vecnorm(squeeze(ncs2States(:, j, :)), 2, 2);
    
    ncs1Qocs(:, j) = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1ErrorNorms(:, j));
    ncs2Qocs(:, j) = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2ErrorNorms(:, j));
end

ncs1MeanError = mean(ncs1ErrorNorms, 2);
ncs1MedianError = median(ncs1ErrorNorms, 2);
ncs2MeanError = mean(ncs2ErrorNorms, 2);
ncs2MedianError = median(ncs2ErrorNorms, 2);

ncs1ErrorVar = var(ncs1ErrorNorms, 0, 2);
ncs2ErrorVar = var(ncs2ErrorNorms, 0, 2);

ncs1UpperQuantile = quantile(ncs1ErrorNorms, 0.9, 2);
ncs1LowerQuantile = quantile(ncs1ErrorNorms, 0.1, 2);
ncs2UpperQuantile = quantile(ncs2ErrorNorms, 0.9, 2);
ncs2LowerQuantile = quantile(ncs2ErrorNorms, 0.1, 2);

% TODO: translate into QoC (all the curves)
translatorData = load(translatorPath, 'translator');

ncs1MeanQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1MeanError);
%ncs1MeanQoc = mean(ncs1Qocs, 2);
ncs1MedianQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1MedianError);
%ncs1MedianQoc = median(ncs1Qocs, 2);
%ncs1UpperQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1UpperQuantile);
ncs1UpperQuantileQoc = quantile(ncs1Qocs, 0.9, 2);
%ncs1LowerQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs1LowerQuantile);
ncs1LowerQuantileQoc = quantile(ncs1Qocs, 0.1, 2);

ncs2MeanQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2MeanError);
%ncs2MeanQoc = mean(ncs2Qocs, 2);
ncs2MedianQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2MedianError);
%ncs2MedianQoc = median(ncs2Qocs, 2);
%ncs2UpperQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2UpperQuantile);
ncs2UpperQuantileQoc = quantile(ncs2Qocs, 0.9, 2);
%ncs2LowerQuantileQoc = arrayfun(@(e) translatorData.translator.translateControlError(e), ncs2LowerQuantile);
ncs2LowerQuantileQoc = quantile(ncs2Qocs, 0.1, 2);

% plot the results
% ncs 1
figure('Name','Ncs 1','NumberTitle','off');
t = tiledlayout(2,1);
title(t, 'Ncs 1');
xlabel(t, 'Simulation time (in seconds)');

nexttile;
hold on;
n = numel(ncs1Times);
l = ncs1MeanError(1:n) - sqrt(ncs1ErrorVar(1:n));
u = ncs1MeanError(1:n) + sqrt(ncs1ErrorVar(1:n));
%fill([ncs1Times', fliplr(ncs1Times')], [l', fliplr(u')], 'g', 'FaceAlpha', 0.5);
fill([ncs1Times', fliplr(ncs1Times')], [ncs1UpperQuantile(1:n)', fliplr(ncs1LowerQuantile(1:n)')], 'g', 'FaceAlpha', 0.5);
plot(ncs1Times, ncs1MeanError(1:n), 'LineWidth', 2);
plot(ncs1Times, ncs1MedianError(1:n), 'LineWidth', 2);

ylim([0 0.016]);
xlim([0 300]);
xticks([0:10:300]);
ylabel('Control Error');
legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex'); 
hold off;

nexttile;
hold on;
%fill([ncs1Times', fliplr(ncs1Times')], [ncs1UpperQuantileQoc(1:n)', fliplr(ncs1LowerQuantileQoc(1:n)')], 'g', 'FaceAlpha', 0.5);
plot(ncs1Times, ncs1MeanQoc(1:n), 'LineWidth', 2);
plot(ncs1Times, ncs1MedianQoc(1:n), 'LineWidth', 2);
plot(ncs1TargetQocTimes, ncs1TargetQoc, 'LineWidth', 2);

ylim([0 1]);
xlim([0 300]);
xticks([0:10:300]);
ylabel('QoC');
legend('Mean', 'Median', 'Target', 'interpreter', 'latex'); 
%legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex'); 
hold off;

% ncs 2
figure('Name','Ncs 2','NumberTitle','off');
t = tiledlayout(2,1);
title(t, 'Ncs 2');
xlabel(t, 'Simulation time (in seconds)');

nexttile;
hold on;
fill([ncs2Times', fliplr(ncs2Times')], [ncs2UpperQuantile', fliplr(ncs2LowerQuantile')], 'g', 'FaceAlpha', 0.5);
plot(ncs2Times, ncs2MeanError, 'LineWidth', 2);
plot(ncs2Times, ncs2MedianError, 'LineWidth', 2);

ylim([0 0.016]);
xlim([0 300]);
xticks([0:10:300]);
ylabel('Control Error');
legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex'); 
hold off;

nexttile;
hold on;
%fill([ncs2Times', fliplr(ncs2Times')], [ncs2UpperQuantileQoc', fliplr(ncs2LowerQuantileQoc')], 'g', 'FaceAlpha', 0.5);
plot(ncs2Times, ncs2MeanQoc, 'LineWidth', 2);
plot(ncs2Times, ncs2MedianQoc, 'LineWidth', 2);
plot(ncs2TargetQocTimes, ncs2TargetQoc, 'LineWidth', 2);

ylim([0 1]);
xlim([0 300]);
xticks([0:10:300]);
ylabel('QoC');
legend('Mean', 'Median', 'Target', 'interpreter', 'latex'); 
%legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex'); 
hold off;


% also plot the error variance over time
figure('Name','Error Variance','NumberTitle','off');
t = tiledlayout(2,1);
title(t, 'Error Variance over Time');
xlabel(t, 'Simulation time (in seconds)');
ylabel(t, 'Error Variance');

nexttile;
hold on;
n = numel(ncs1Times);
plot(ncs1Times, ncs1ErrorVar(1:n), 'LineWidth', 2);
xticks([0:10:300]);
xlim([0 300]);
ylim([0 0.00002]);
title('Ncs 1');
hold off;

nexttile;
hold on;
plot(ncs2Times, ncs2ErrorVar, 'LineWidth', 2);
xticks([0:10:300]);
xlim([0 300]);
ylim([0 0.00002]);
title('Ncs 2');
hold off;


