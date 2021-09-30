% parameters as defined in Omnet config file ncs-testbench/simulations/controller/characterization.ini
% results are from config SamplingRateSeriesDoublePendulumMPC
% corresponding to stabilization of of a double inverted pendulum
packetRates = 50:5:300;%50:5:350; 
numRepetitions = 50;
numPacketRates = numel(packetRates);
simTime = 30; % in seconds
numEls = simTime * packetRates;
prefix = 'matlab/Evaluation/controller-characteristics/characteristics/double-pendulum/robust/';

% read the Omnet results
dataRes = tabularTextDatastore([prefix 'double-pendulum-robust-10sec.csv'], 'ReadVariableNames',false);
allVariableNames = dataRes.VariableNames;
dataRes.SelectedVariableNames = allVariableNames(2:2:numel(allVariableNames));
allData = readall(dataRes);

controlErrors = cell(1, numPacketRates);
mappingRateRun = zeros(1, numPacketRates*numRepetitions);

for counter=1:size(allData, 2)  
    data = table2array(allData(:, counter));
    data = data(~isnan(data));     % is a column vector
    rate = find(numEls == size(data, 1));
    if ~isempty(rate)
        controlErrors{rate} = [controlErrors{rate}; data'];
        mappingRateRun(counter) = rate;
    end
end
clear allData dataRes

assert(all(cellfun(@(c) size(c, 1), controlErrors) == numRepetitions));

% % dataRes = tabularTextDatastore([prefix 'double-pendulum-robust-10sec-state-norms.csv'], 'ReadVariableNames',false);
% % allVariableNames = dataRes.VariableNames;
% % dataRes.SelectedVariableNames = allVariableNames(2:2:numel(allVariableNames));
% % allData = readall(dataRes);
% % 
% % numPlantSteps = simTime * 1000;
% % errorNorms = zeros(numPacketRates, numRepetitions, numPlantSteps + 1);
% % idx = ones(1, numPacketRates);
% % 
% % for counter=1:size(allData, 2)  
% %     data = table2array(allData(2:end, counter));
% %     data = data(~isnan(data));     % is a column vector
% %     rate = mappingRateRun(counter);
% %     if rate ~= 0
% %         errorNorms(rate, idx(rate), :) = data;
% %         idx(rate) = idx(rate) + 1;
% %     end
% % end
% % clear allData dataRes


% get rid of first 10 seconds and consider only steady state
skip = 10; % wait for 10 seconds
avgErrors = zeros(numPacketRates, 1);

avgErrorsOverTime = cell(numPacketRates, 1);
% %avgErrorNorms = zeros(numPacketRates, 1);

for j=1:numPacketRates    
    timeInstants = (0:numEls(j)-1) / packetRates(j);    

    avgErrors(j) = mean(controlErrors{j}(:, timeInstants > skip), 'all');  
    avgErrorsOverTime{j} = mean(controlErrors{j});
    
% %     avgErrorNorms(j) = mean(errorNorms(j, :, (0:numPlantSteps) / 1000 > skip), 'all');
end

% % avgErrorNormsOverTime = squeeze(mean(errorNorms, 2));


qoc = avgErrors(end) ./ avgErrors;
qoc = rescale(qoc, 0, 1); % so that 50Hz corresponds to QoC = 0 and 350 Hz is QoC=1


% fit the error functions: rate -> control error
[xData, yData] = prepareCurveData(packetRates, avgErrors);
% Set up fittype and options.
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');
opts.StartPoint = [0.0704686609696171 -0.0338799221308433 0.0119180702561441 -0.00438121598062712];
% Fit model to data.
[fitresult, ~] = fit(xData, yData, fittype('exp2'), opts);

% fit the QoC functions: rate -> QoC
[xData, yData] = prepareCurveData(packetRates, qoc);
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');
opts.StartPoint = [0.222111923535457 0.00542877168428098 -0.372075377609209 -0.01760065657094];
% Fit model to data.
[fitresult_qoc, ~] = fit(xData, yData, fittype('exp2'), opts);

% but we need the inverse mapping QoC-> rate
[xData, yData] = prepareCurveData(qoc, packetRates);
opts = fitoptions('Method', 'SmoothingSpline');
opts.SmoothingParam = 0.999908842787108;
% Fit model to data.
[fitresult_qocRate, gof] = fit(xData, yData, fittype('smoothingspline'), opts);

% finally, relate error to QoC
[xData, yData] = prepareCurveData(avgErrors, qoc);
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');    
opts.StartPoint = [0.521135830804001 0.231594386708524 0.488897743920167];
[fitresult_qocError, ~] = fit(xData, yData, fittype('rat11'), opts);

translatorRobustDoublePendulum = NcsTranslator(fitresult_qocRate, fitresult_qocError, packetRates(end));

% plot something
figure();
plot(packetRates, avgErrors, 'Marker', 'o', 'LineWidth', 2');
hold on;
plot(fitresult);
legend('true', 'fitted', 'Location', 'NorthEast', 'Interpreter', 'none');
ylabel('Avg. Control Error');
xlabel('Controller Sampling Rate (in Hz)');
title('Double Pendulum (MPC)');


legends = cell(1, numPacketRates);
figure();
for j=1:numPacketRates
    timeInstants = (0:numEls(j)-1) / packetRates(j);
    hold on;
    legends{j} = sprintf('%d Hz', packetRates(j));    
    plot(timeInstants(timeInstants > 0), avgErrorsOverTime{j}(:, timeInstants > 0), 'LineWidth', 2);
    grid on;
end
ylabel('Average Actual Control Error');
xlabel('Simulation time (in seconds)');
legend(legends, 'Location', 'westoutside', 'FontSize', 16);
title('Double Pendulum (Robust)');
