% parameters as defined in Omnet config file ncs-testbench/simulations/controller/characterization.ini
% results are from config SamplingRateSeriesShortPendulumRobust
% corresponding to stabilization of short pendulum with the robust
% controller using delta = 0.1 and sequences of length 2
packetRates = 20:5:200;
numRepetitions = 50;
numPacketRates = numel(packetRates);
simTime = 30; % in seconds
numEls = simTime * packetRates;
prefix = 'matlab/Evaluation/controller-characteristics/characteristics/short-pendulum/robust/';

% read the Omnet results
b = readtable([prefix 'pendulum-short-robust-10sec.csv'], 'ReadVariableNames',false);

controlErrorsShort = cell(1, numPacketRates);

for counter=1:numPacketRates*numRepetitions-1
    data = table2array(b(2:end, 2*counter));
    data = data(~isnan(data));     % is a column vector
    rate = find(numEls == size(data, 1));
    controlErrorsShort{rate} = [controlErrorsShort{rate}; data'];
end

% get rid of first 10 seconds and consider only steady state
skip = 10; % wait for 10 seconds
avgErrorsShort = zeros(numPacketRates, 1);

avgErrorsOverTimeShort = cell(numPacketRates, 1);

for j=1:numPacketRates    
    timeInstants = (0:numEls(j)-1) / packetRates(j);    

    avgErrorsShort(j) = mean(controlErrorsShort{j}(:, timeInstants > skip), 'all');
    
    avgErrorsOverTimeShort{j} = mean(controlErrorsShort{j});
end

qocShort = avgErrorsShort(end) ./ avgErrorsShort;
qocShort = rescale(qocShort, 0, 1); % so that 20Hz corresponds to QoC = 0 and 200 Hz is QoC=1

% fit the error functions: rate -> control error
[xData, yData] = prepareCurveData(packetRates, avgErrorsShort);
% Set up fittype and options.
% opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');
% opts.StartPoint = [0.890903252535798 0.959291425205444 0.547215529963803];

opts = fitoptions('Normalize', 'on');
% Fit model to data.
[fitresult_short, ~] = fit(xData, yData, fittype('pchipinterp'), opts);

% fit the QoC functions: rate -> QoC
[xData, yData] = prepareCurveData(packetRates, qocShort);
excludedPoints = excludedata( xData, yData, 'Indices', [29 30 31 32 33 34 36] );
opts = fitoptions( 'Method', 'PchipInterpolant' );
opts.Normalize = 'on';
opts.Exclude = excludedPoints;
%opts.StartPoint = [0.00476025855939998 1.07562754713306 -0.0806990088369622];
% Fit model to data.
[fitresult_qocShort, ~] = fit(xData, yData, fittype('pchipinterp'), opts);

% but we need the inverse mapping QoC-> rate
[xData, yData] = prepareCurveData(qocShort, packetRates);
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.999992506973253;

[fitresult_qocRateShort, gof] = fit(xData, yData, fittype('smoothingspline'), opts);

% finally, relate error to QoC
[xData, yData] = prepareCurveData(avgErrorsShort, qocShort);
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');    
opts.StartPoint = [6.88518449391243e-05 -1.40165016902307 0.0163853735906044];
[fitresult_qocErrorShort, ~] = fit(xData, yData, fittype('power2'), opts);

translatorRobustShort = NcsTranslator(fitresult_qocRateShort, fitresult_qocErrorShort, packetRates(end));

% plot something
figure();
plot(packetRates, avgErrorsShort, 'Marker', 'o', 'LineWidth', 2');
hold on;
plot(fitresult_short);
legend('true', 'fitted', 'Location', 'NorthEast', 'Interpreter', 'none');
ylabel('Avg. Control Error');
xlabel('Controller Sampling Rate (in Hz)');
title('Short Pendulum');


legends = cell(1, numPacketRates);
figure();
for j=1:numPacketRates
    timeInstants = (0:numEls(j)-1) / packetRates(j);
    hold on;
    legends{j} = sprintf('%d Hz', packetRates(j)); 
    %plot(timeInstants(timeInstants > skip), avgErrorsOverTimeShort{j}(:, timeInstants > skip), 'LineWidth', 2);
    plot(timeInstants(timeInstants > 0), avgErrorsOverTimeShort{j}(:, timeInstants > 0), 'LineWidth', 2);
    grid on;
end
ylabel('Average Actual Control Error');
xlabel('Simulation time (in seconds)');
legend(legends, 'Location', 'westoutside', 'FontSize', 16);
title('Short Pendulum');
