% example for an event based approach to stabilize an inverted pendulum over 3 seconds, 
% where the impact of different threshold values (called deadband) is investigated

% requires libncs_matlab to be on the Matlab path, so check and add 
p = genpath('libncs_matlab/matlab');
if ~contains(path, split(strip(p, 'right', ':'), ':'))
    addpath(p);
end
% System properties: discrete time inverted pendulum with default
% parameters 
config = CreateInvertedPendulumScenario(); 

caDelayProb = config.caDelayProbs;
scDelayProb = caDelayProb;
acDelayProb = caDelayProb;

simTime = 3; % seconds to be simulated
config.numTimeSteps = simTime / config.samplingInterval;
config.plantSamplingInterval = config.samplingInterval;
config.ignoreControllerCache = true;

config.controllerClassName = 'InfiniteHorizonController';
config.controllerEventBased = true;
config.controllerEventTrigger = EventBasedControllerTriggerCriterion.Sequence;
controllerDeadbands = [0.0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4];

% seed ftw
rng(42);

horizonlength = 10; % time steps for the control error
runsPerDeadband = 10;

% draw the packet delays in advance
caPacketDelays = DiscreteSample(caDelayProb,config.numTimeSteps * runsPerDeadband)-1;
caPacketDelays = reshape(caPacketDelays, config.numTimeSteps, runsPerDeadband);
scPacketDelays = DiscreteSample(scDelayProb,config.numTimeSteps * runsPerDeadband)-1;
scPacketDelays = reshape(scPacketDelays, config.numTimeSteps, runsPerDeadband);
acPacketDelays = DiscreteSample(acDelayProb,config.numTimeSteps * runsPerDeadband)-1;
acPacketDelays = reshape(acPacketDelays, config.numTimeSteps, runsPerDeadband);

packetRates = zeros(runsPerDeadband, numel(controllerDeadbands));
avgCosts = zeros(runsPerDeadband, numel(controllerDeadbands));
qualityOfControl = zeros(runsPerDeadband, numel(controllerDeadbands));
controlError = zeros(runsPerDeadband, numel(controllerDeadbands));

printer = SimulationInfoPrinter('Ncs Example: Event-based stabilization of inverted pendulum', ...
    numel(controllerDeadbands), runsPerDeadband);
printer.printSimulationStart();

for i=1:numel(controllerDeadbands)
    config.controllerDeadband = controllerDeadbands(i);
    for r=1:runsPerDeadband        
        printer.printProgress(i, r);
        example = NcsPendulumExample(config, simTime);
        [statsPerStep, stats, costs] = example.simulate(caPacketDelays(:, r), scPacketDelays(:, r), acPacketDelays(:, r));

        packetRates(r, i) = sum([statsPerStep(:).ca_sent]) / simTime;
        avgCosts(r, i) = costs;  

        % compute the value of the cost function (finite-horizon, 10 time steps)
        sums = movsum([statsPerStep(:).actual_stagecosts], [horizonlength - 1 0], 'Endpoints', 'discard');
        qualityOfControl(r, i) = mean(sums);
        
        errors = movsum([statsPerStep(:).actual_control_error], [horizonlength - 1 0], 'Endpoints', 'discard');
        controlError(r, i) = mean(errors);
    end
end

printer.printSimulationEnd();

figure();
hold on;

ax = gca;
%ax.FontSize = fontSize;
ax.Box = 'On';
%set(gca, 'DefaultLineLineWidth', lineWidth);

xlabel('Packets / $s$', 'Interpreter', 'latex');
ylabel('$J_{avg}$', 'Interpreter', 'latex');

plot(mean(packetRates), mean(avgCosts), '-o');

figure();
hold on;

ax = gca;
%ax.FontSize = fontSize;
ax.Box = 'On';
%set(gca, 'DefaultLineLineWidth', lineWidth);

xlabel('Packets / $s$', 'Interpreter', 'latex');
ylabel('$J_{mean}$', 'Interpreter', 'latex');

plot(mean(packetRates), mean(qualityOfControl), '-x');

figure();
hold on;

ax = gca;
%ax.FontSize = fontSize;
ax.Box = 'On';
%set(gca, 'DefaultLineLineWidth', lineWidth);

xlabel('Packets / $s$', 'Interpreter', 'latex');
ylabel('Norm of Cumulated Control Error', 'Interpreter', 'latex');

plot(mean(packetRates), mean(controlError), '-x');
