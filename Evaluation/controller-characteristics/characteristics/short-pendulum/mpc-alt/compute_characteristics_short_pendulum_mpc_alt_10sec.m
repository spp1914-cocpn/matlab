% parameters as defined in Omnet config file ncs-testbench/simulations/controller/characterization.ini
% results are from config SamplingRateSeriesShortPendulumMPCAlternative
% corresponding to stabilization of short pendulum
packetRates = 20:5:200;
numRepetitions = 50;
numPacketRates = numel(packetRates);
simTime = 30; % in seconds
numEls = simTime * packetRates;
prefix = 'matlab/Evaluation/controller-characteristics/characteristics/short-pendulum/mpc-alt/';

% read the Omnet results
b = readtable([prefix 'pendulum-short-mpc-alt-10sec.csv'], 'ReadVariableNames',false);

controlErrorsShort = cell(1, numPacketRates);

for counter=1:numPacketRates*numRepetitions-1
    % read the first time
    %t = str2double(table2array(b(3, 2*counter-1)))-str2double(table2array(b(2, 2*counter-1)))
    data = table2array(b(2:end, 2*counter));
    data = data(~isnan(data));     % is a column vector
    rate = find(numEls == size(data, 1));   
    %rate = find(abs(packetRates-(1/t)) <= 1e-5);
    try
        controlErrorsShort{rate} = [controlErrorsShort{rate}; data'];
    catch
    end
end

skip = 12; % wait for 12 seconds
avgErrorsShort = zeros(numPacketRates, 1);

avgErrorsOverTimeShort = cell(numPacketRates, 1);

for j=1:numPacketRates
    nums = size(controlErrorsShort{j}, 2);
    timeInstants = (0:nums-1) / packetRates(j);    

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
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.StartPoint = [5.97083360317071e-05 2.00197405462814 -0.311078507804225];
% Fit model to data.
[fitresult_qocShort, ~] = fit(xData, yData, fittype('power2'), opts);

% but we need the inverse mapping QoC-> rate
[xData, yData] = prepareCurveData(qocShort, packetRates);
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.StartPoint = [97.2211706927963 1.17908435943609 -75.0191575587005 0.432661691237592];
% Fit model to data.
[fitresult_qocRateShort, gof] = fit(xData, yData, fittype('exp2'), opts);

% finally, relate error to QoC
[xData, yData] = prepareCurveData(avgErrorsShort, qocShort);
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');    
opts.StartPoint = [0.000105761311645691 -1.21189100684468 0.0200698131102982];
[fitresult_qocErrorShort, ~] = fit(xData, yData, fittype('power2'), opts);

translatorMpcShort = NcsTranslator(fitresult_qocRateShort, fitresult_qocErrorShort, packetRates(end));

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
