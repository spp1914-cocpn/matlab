% parameters as defined in Omnet config file ncs-testbench/simulations/controller/characterization.ini
% results are from config SamplingRateSeriesDoublePendulumMPCAlt
% corresponding to stabilization of a double inverted pendulum
packetRates = 50:5:300;
numRepetitions = 50;
numPacketRates = numel(packetRates);
simTime = 30; % in seconds
numEls = simTime * packetRates;
prefix = 'matlab/Evaluation/controller-characteristics/characteristics/double-pendulum/mpc-alt/';

% read the Omnet results
dataRes = tabularTextDatastore([prefix 'double-pendulum-mpc-alt-10sec.csv'], 'ReadVariableNames',false);
allVariableNames = dataRes.VariableNames;
dataRes.SelectedVariableNames = allVariableNames(2:2:numel(allVariableNames));
allData = readall(dataRes);

controlErrors = cell(1, numPacketRates);

for counter=1:size(allData, 2)  
    data = table2array(allData(:, counter));
    data = data(~isnan(data));     % is a column vector
    rate = find(numEls == size(data, 1));
    controlErrors{rate} = [controlErrors{rate}; data'];
end
clear allData dataRes

assert(all(cellfun(@(c) size(c, 1), controlErrors) == numRepetitions));

% get rid of first 10 seconds and consider only steady state
skip = 10; % wait for 10 seconds
avgErrors = zeros(numPacketRates, 1);

avgErrorsOverTime = cell(numPacketRates, 1);

for j=1:numPacketRates    
    timeInstants = (0:numEls(j)-1) / packetRates(j);    

    avgErrors(j) = mean(controlErrors{j}(:, timeInstants > skip), 'all');    
    avgErrorsOverTime{j} = mean(controlErrors{j});
end

qoc = avgErrors(end) ./ avgErrors;
qoc = rescale(qoc, 0, 1); % so that 50Hz corresponds to QoC = 0 and 300 Hz is QoC=1

% fit the error functions: rate -> control error
[xData, yData] = prepareCurveData(packetRates, avgErrors);
% Set up fittype and options.
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');
opts.StartPoint = [0.40800394072395 0.103601841487092 0.32199319271624];
% Fit model to data.
[fitresult, ~] = fit(xData, yData, fittype('rat11'), opts);

% fit the QoC functions: rate -> QoC
[xData, yData] = prepareCurveData(packetRates, qoc);
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');
opts.StartPoint = [0.475073084688092 0.221807301771909 0.564419272401814];
% Fit model to data.
[fitresult_qoc, ~] = fit(xData, yData, fittype('rat11'), opts);

% but we need the inverse mapping QoC-> rate
[xData, yData] = prepareCurveData(qoc, packetRates);
opts = fitoptions('Method', 'SmoothingSpline');
opts.SmoothingParam = 0.998478167166173;
% Fit model to data.
[fitresult_qocRate, gof] = fit(xData, yData, fittype('smoothingspline'), opts);

% finally, relate error to QoC
[xData, yData] = prepareCurveData(avgErrors, qoc);
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');    
opts.StartPoint = [0.000178934722006612 -1.42473305826902 0.0276083733077218];
[fitresult_qocError, ~] = fit(xData, yData, fittype('power2'), opts);

translatorMpcAltDoublePendulum = NcsTranslator(fitresult_qocRate, fitresult_qocError, packetRates(end));

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
title('Double Pendulum (MPC)');
