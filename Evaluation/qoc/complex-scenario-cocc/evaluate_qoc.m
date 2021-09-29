numRepetitions = 100;
simTime = 480; % in seconds
numNcs = 8;

prefix = 'matlab/Evaluation/qoc/complex-scenario-cocc/';
%scenarioName = 'complex-scenario-cocc'; % 'complex-scenario-cocc-mpc'
%scenarioName = 'complex-scenario-cocc-mpc';
scenarioName = 'complex-scenario-cocc-mpc-alt';
filePrefix = [prefix scenarioName];
if strcmp(scenarioName, 'complex-scenario-cocc')
    translatorPath = 'libncs_matlab/matlab/config/translators/lqr/hamming/translator_inverted_pendulum_short_lqr_hamming_10sec.mat';
elseif strcmp(scenarioName, 'complex-scenario-cocc-mpc-alt')
    translatorPath = 'libncs_matlab/matlab/config/translators/mpc_alt/hamming/translator_inverted_pendulum_short_mpc_alt_hamming_10sec.mat';
else
    translatorPath = 'libncs_matlab/matlab/config/translators/mpc/hamming/translator_inverted_pendulum_short_mpc_hamming_10sec.mat';
end
translatorData = load(translatorPath, 'translator');

% read the Omnet results
targetQocs = readtable([filePrefix '-target-qoc.csv'], 'ReadVariableNames', false);

ncsTargetQocTimes = cell(1, numNcs);
ncsTargetQocs = cell(1, numNcs);

idx = 1;
for j=1:numNcs   
    ncsTargetQocTimes{j} = cell2mat(cellfun(@str2double, table2array(targetQocs(2:end, idx)), 'UniformOutput', false));
    ncsTargetQocTimes{j} = ncsTargetQocTimes{j}(~isnan(ncsTargetQocTimes{j}));
    ncsTargetQocs{j} = table2array(targetQocs(2:numel(ncsTargetQocTimes{j})+1, idx+1));

    idx = idx + 2;
end

clear targetQocs

ncsTimes = cell(1, numNcs);
ncsStates = cell(1, numNcs);
states = datastore([filePrefix '-controller-state-1.csv'], 'ReadVariableNames', false);
allVariableNames = states.SelectedVariableNames;
idx = 1;
for j=1:numNcs    
    % read the times first
    states.SelectedVariableNames = allVariableNames(idx);    
    data = readall(states);    
    ncsTimes{j} = cell2mat(cellfun(@str2double, table2array(data(~any(ismissing(data), 2), :)), 'UniformOutput', false));
    ncsTimes{j} = ncsTimes{j}(~isnan(ncsTimes{j}));
    
    states.SelectedVariableNames = allVariableNames([idx+1:2*numNcs:length(states.VariableNames)]);
    data = readall(states);    
    ncsStates{j} = table2array(data(~any(ismissing(data), 2), :));    
    idx = idx + 2;
end

states = datastore([filePrefix '-controller-state-2.csv'], 'ReadVariableNames', false);
idx = 1;
for j=1:numNcs
    states.SelectedVariableNames = allVariableNames([idx+1:2*numNcs:length(states.VariableNames)]);
    data = readall(states);    
    ncsStates{j} = cat(3, ncsStates{j}, ...
        table2array(data(~any(ismissing(data), 2), :)));
    idx = idx + 2;
end

states = datastore([filePrefix '-controller-state-3.csv'], 'ReadVariableNames', false);
idx = 1;
for j=1:numNcs
    states.SelectedVariableNames = allVariableNames([idx+1:2*numNcs:length(states.VariableNames)]);
    data = readall(states);    
    ncsStates{j} = cat(3, ncsStates{j}, ...
        table2array(data(~any(ismissing(data), 2), :)) -pi);
    idx = idx + 2;
end

states = datastore([filePrefix '-controller-state-4.csv'], 'ReadVariableNames', false);
idx = 1;
for j=1:numNcs
    states.SelectedVariableNames = allVariableNames([idx+1:2*numNcs:length(states.VariableNames)]);
    data = readall(states);    
    ncsStates{j} = cat(3, ncsStates{j}, ...
        table2array(data(~any(ismissing(data), 2), :)));
    idx = idx + 2;
end

assert(sum(cellfun(@(c) size(c, 2), ncsStates)) == numRepetitions * numNcs);

% compute the error norm directly, per time step and run
ncsErrorNorms = cell(1, numNcs);
ncsMeanErrors = cell(1, numNcs);
ncsMedianErrors = cell(1, numNcs);
ncsUpperQuantiles = cell(1, numNcs);
ncsLowerQuantiles = cell(1, numNcs);
%ncsQocs = cell(1, numNcs);
ncsErrorVars = cell(1, numNcs);


for j=1:numNcs
    ncsErrorNorms{j} = zeros(size(ncsStates{j}, 1), numRepetitions);
    %ncsQocs{j} = zeros(size(ncsStates{j}, 1), numRepetitions);
    for r=1:numRepetitions
        ncsErrorNorms{j}(:, r) = vecnorm(squeeze(ncsStates{j}(:, r, :)), 2, 2);
        %ncsQocs{j}(:, r) = translatorData.translator.translateControlError(ncsErrorNorms{j}(:, r));
    end
    ncsMeanErrors{j} = mean(ncsErrorNorms{j}, 2, 'omitnan');
    ncsMedianErrors{j} = median(ncsErrorNorms{j}, 2, 'omitnan');
    ncsUpperQuantiles{j} = quantile(ncsErrorNorms{j}, 0.9, 2);
    ncsLowerQuantiles{j} = quantile(ncsErrorNorms{j}, 0.1, 2);
    ncsErrorVars{j} = var(ncsErrorNorms{j}, 0, 2, 'omitnan');
end

% translate into qoc
ncsMeanQocs = cell(1, numNcs);
ncsMedianQocs = cell(1, numNcs);
ncsUpperQuantileQocs = cell(1, numNcs);
ncsLowerQuantileQocs = cell(1, numNcs);
for j=1:numNcs
    ncsMeanQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsMeanErrors{j});
    ncsMedianQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsMedianErrors{j}); 
    %ncsUpperQuantileQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsUpperQuantiles{j});
    %ncsLowerQuantileQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsLowerQuantiles{j});
end

% % compute the mean error and the variances, and quantiles
% ncsMeanErrors = cellfun(@(errors) mean(errors, 2, 'omitnan'), ncsErrors, 'UniformOutput', false);
% ncsMedianErrors = cellfun(@(errors) median(errors, 2, 'omitnan'), ncsErrors, 'UniformOutput', false);
% ncsUpperQuantiles = cellfun(@(errors) quantile(errors, 0.9, 2), ncsErrors, 'UniformOutput', false);
% ncsLowerQuantiles = cellfun(@(errors) quantile(errors, 0.1, 2), ncsErrors, 'UniformOutput', false);
% 
% ncsMeanQocs = cell(1, numNcs);
% ncsMedianQocs = cell(1, numNcs);
% ncsUpperQuantileQocs = cell(1, numNcs);
% ncsLowerQuantileQocs = cell(1, numNcs);
% for j=1:numNcs
%     ncsMeanQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsMeanErrors{j});
%     ncsMedianQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsMedianErrors{j}); 
%     ncsUpperQuantileQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsUpperQuantiles{j});
%     ncsLowerQuantileQocs{j} = arrayfun(@(e) translatorData.translator.translateControlError(e), ncsLowerQuantiles{j});
% end

for j=1:numNcs
    figure('Name', sprintf('Ncs %d', j), 'NumberTitle','off');
    t = tiledlayout(2, 1);
    title(t, sprintf('Ncs %d', j));
    xlabel(t, 'Simulation Time (in seconds)');
    
    nexttile;
    hold on;
    n = min(numel(ncsTimes{j}), numel(ncsMeanErrors{j}));
    fill([ncsTimes{j}(1:n)', fliplr(ncsTimes{j}(1:n)')], [ncsUpperQuantiles{j}(1:n)', fliplr(ncsLowerQuantiles{j}(1:n)')], 'g', 'FaceAlpha', 0.5);
    plot(ncsTimes{j}(1:n), ncsMeanErrors{j}(1:n), 'LineWidth', 2);
    plot(ncsTimes{j}(1:n), ncsMedianErrors{j}(1:n), 'LineWidth', 2);
    
    xlim([0 simTime]);
    xticks([0:20:simTime]);
    ylim([0 0.005]);    
    legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'interpreter', 'latex'); 
    ylabel('Control Error');    
    hold off;
    
    nexttile;
    hold on;
    %fill([ncsTimes{j}', fliplr(ncsTimes{j}')], [ncsUpperQuantileQocs{j}', fliplr(ncsLowerQuantileQocs{j}')], 'g', 'FaceAlpha', 0.5);
    plot(ncsTimes{j}(1:n), ncsMeanQocs{j}(1:n), 'LineWidth', 2);
    plot(ncsTimes{j}(1:n), ncsMedianQocs{j}(1:n), 'LineWidth', 2);
    plot(ncsTargetQocTimes{j}, ncsTargetQocs{j}, 'LineWidth', 2);
    
    xlim([0 simTime]);
    xticks([0:20:simTime]);
    ylim([0 1]);
    ylabel('QoC'); 
    %legend('Mean', 'Target', 'interpreter', 'latex'); 
    legend('Mean', 'Median', 'Target', 'interpreter', 'latex'); 
    %legend('$Q_{0.9}$ - $Q_{0.1}$', 'Mean', 'Median', 'Target', 'interpreter', 'latex');
    hold off;     
end


for j = 1:numNcs     
    % also plot the error variance
    figure('Name', sprintf('Ncs %d', j), 'NumberTitle','off');
    title(sprintf('Ncs %d', j));
    
    n = min(numel(ncsTimes{j}), numel(ncsErrorVars{j}));
    
    hold on;
    plot(ncsTimes{j}(1:n), ncsErrorVars{j}(1:n), 'LineWidth', 2);
    xlim([0 simTime]);
    ylim([0 0.00002]);
    xticks([0:20:simTime]);
    ylabel('Control Error Variance');
    xlabel('Simulation Time (in seconds)');
    hold off;
end
