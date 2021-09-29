% parameters as defined in Omnet config file ncs-testbench/simulations/controller/characterization.ini
% results are from config SamplingRateSeriesShortPendulumLQR
% corresponding to stabilization of short
packetRates = 20:5:200;
numRepetitions = 50;
numPacketRates = numel(packetRates);
simTime = 30; % in seconds
numEls = simTime * packetRates;
prefix = 'matlab/Evaluation/controller-characteristics/characteristics/short-pendulum/lqr/';

% read the Omnet results
b = readtable([prefix 'pendulum-short-lqr-10sec.csv'], 'ReadVariableNames',false);

controlErrorsShort = cell(1, numPacketRates);

clamp100= false;
clamp100Above = false;

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

if clamp100
    r = find(packetRates == 100);
    packetRates = packetRates(1:r);
    avgErrorsShort = avgErrorsShort(1:r);
    numPacketRates = r;
    numEls = numEls(1:r);
elseif clamp100Above
    r = find(packetRates > 100);
    packetRates = packetRates(r);
    avgErrorsShort = avgErrorsShort(r);
    numPacketRates = numel(r);
    numEls = numEls(r);     
end
qocShort = avgErrorsShort(end) ./ avgErrorsShort;
qocShort = rescale(qocShort, 0, 1); % so that 20Hz corresponds to QoC = 0 and 200 Hz is QoC=1

if clamp100
    % fit the error functions: rate -> control error
    [xData, yData] = prepareCurveData(packetRates, avgErrorsShort);
    % Set up fittype and options.

    opts = fitoptions('Normalize', 'on');
    % Fit model to data.
    [fitresult_short, ~] = fit(xData, yData, fittype('pchipinterp'), opts);
    
    % fit the QoC functions: rate -> QoC
    [xData, yData] = prepareCurveData(packetRates, qocShort);
    opts = fitoptions('Normalize', 'on');
    %opts.StartPoint = [0.00476025855939998 1.07562754713306 -0.0806990088369622];
    % Fit model to data.
    [fitresult_qocShort, ~] = fit(xData, yData, fittype('pchipinterp'), opts);
    
    % but we need the inverse mapping QoC-> rate
    [xData, yData] = prepareCurveData(qocShort, packetRates);
    opts = fitoptions('Normalize', 'on');
    %opts.SmoothingParam = 0.999999662184475;
    % Fit model to data.
    [fitresult_qocRateShort, gof] = fit(xData, yData, fittype('pchipinterp'), opts);
    
    % finally, relate error to QoC
    [xData, yData] = prepareCurveData(avgErrorsShort, qocShort);
    opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');    
    opts.StartPoint = [6.23185104133919e-05 -1.34624215393511 0.0159302656232381];
    [fitresult_qocErrorShort, ~] = fit(xData, yData, fittype('power2'), opts);

    translatorLqrShort = NcsTranslator(fitresult_qocRateShort, fitresult_qocErrorShort, packetRates(end));
else
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
    opts = fitoptions('Normalize', 'on');
    %opts.StartPoint = [0.00476025855939998 1.07562754713306 -0.0806990088369622];
    % Fit model to data.
    [fitresult_qocShort, ~] = fit(xData, yData, fittype('pchipinterp'), opts);

    % but we need the inverse mapping QoC-> rate
    [xData, yData] = prepareCurveData(qocShort, packetRates);
    opts = fitoptions('Normalize', 'on');
    %opts.SmoothingParam = 0.999999662184475;
    % Fit model to data.
    [fitresult_qocRateShort, gof] = fit(xData, yData, fittype('pchipinterp'), opts);

    % finally, relate error to QoC
    [xData, yData] = prepareCurveData(avgErrorsShort, qocShort);
    opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');    
    opts.StartPoint = [9.94972816084908e-05 -1.21653802397541 0.014228410827563];
    [fitresult_qocErrorShort, ~] = fit(xData, yData, fittype('power2'), opts);

    translatorLqrShort = NcsTranslator(fitresult_qocRateShort, fitresult_qocErrorShort, packetRates(end));

end

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
