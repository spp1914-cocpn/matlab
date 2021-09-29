% function to highlight how API functions ins libncs_matlab/matlab/api can
% be called from Matlab to carry out simulations without Omnet
% this can for instance be used to show that, if same realizations of
% packet delays and losses are used, simulation results with and without Omnet
% are exactly the same

simTime = 600;
maxSimTime = ConvertToPicoseconds(simTime);
configStruct.maxMeasDelay = int64(1);
configStruct.controllerClassName='MSSController';
%configStruct.filterClassName='DelayedModeIMMF';
%configStruct.maxMeasDelay = 1;
configStruct.networkType = 2;
configStruct.controlSequenceLength = int64(2);
%configStruct.plantSamplingInterval = 0.001;
configStruct.caDelayProbs = zeros(1, 10) + 1/10;
configStruct.scDelayProbs = [0 1 0 0 0 0 0];

configStruct.mpcHorizon = 100;
%configStruct.Q = eye(4);%10 * blkdiag(100, 0, 100, 0);
%configStruct.R = 1;%10;

% arbitrarily place the poles
configStruct.polesCont = [-1 -1.1 -1.2 -1.3]%[-5.61 -0.143 -6 -5];%-5+1.5j -5-1.5j];
% configStruct.polesCont = [-1+5i; 
%             -1-5i
%             -6+2*1i;
%             -6-2*1i;            
%          ];

packetRates = [100]%[50 100 200];%[20:5:200]%[20 25 50 100 200];%[40:5:200]%50:5:200; %200 [20:5:200]
numPacketRates = numel(packetRates);
samplingIntervals = 1 ./ packetRates
% numControllerTimesteps = simTime ./ samplingIntervals;
% ratios = numTimesteps ./ numControllerTimesteps;
numRepetitions = 1;

caDelays = [0 1 0 0 0 0 0];
trueStates = cell(numRepetitions, numPacketRates);
costs = cell(numRepetitions, numPacketRates);
controllerStats = cell(numRepetitions, numPacketRates);
plantStats = cell(numRepetitions, numPacketRates);
controlErrors = cell(numRepetitions, numPacketRates); % true vs estimated
for r =1:numRepetitions
    for j=1:numPacketRates    
        configStruct.plantSamplingInterval = samplingIntervals(j);
        configStruct.plantSamplingInterval = 0.001;
        configStruct.samplingInterval = samplingIntervals(j);
        
        %configStruct.caDelayProbs = [0 .5 0.3 0.2];% 0 0 0 0];
        configStruct.caDelayProbs = [1e-20 1-1e-20 0 0 0];
        %configStruct.caDelayProbs = [1e-20 0.85 1-0.85-1e-20 0 0];        
        %configStruct.caDelayProbs = [0 0.25 0.5 0.25 0];
        
        ncs_seedRng(r-1); % same as in omnet
        ncs = ncs_initialize(maxSimTime, 1, configStruct, 'libncs_matlab/matlab/config/inverted_pendulum_short.mat');
        %ncs = ncs_initialize(maxSimTime, 1, configStruct, 'libncs_matlab/matlab/config/double_integrator.mat');
        
        numControllerTimesteps = maxSimTime / uint64(ncs_getTickerInterval(ncs));
        numTimesteps = maxSimTime / uint64(ncs_getPlantTickerInterval(ncs));
        ratio = numTimesteps / numControllerTimesteps

        caNetwork = CommunicationNetwork(numTimesteps, 12);
        scNetwork = CommunicationNetwork(numTimesteps, 1); % constant delay, 1 time step
        
        controlErrors{r, j} = zeros(2, numControllerTimesteps); % true vs estimated
        packetList = cell(1,0);
        nextControllerInvocation = ratio;
        controllerK = 1;
        for k=1:uint64(numTimesteps)            
            % first a plant step
            %timestamp = ConvertToPicoseconds(k * configStruct.plantSamplingInterval)
            timestamp = int64(ConvertToPicoseconds(double(k) * configStruct.plantSamplingInterval));
            if~ncs_doPlantStep(ncs, timestamp)
                % plant state is inadmissible, e.g. pendulum has fallen over
                k
                break
            end
            if k == nextControllerInvocation
                % receive all the packets first
                [numRcvdCa, pktsInCa] = caNetwork.receivePackets(controllerK);
                for ii=1:numRcvdCa
                    % construct the packets to be received        
                    ncs_doHandlePacket(ncs, timestamp, pktsInCa(ii));
                end
                [numRcvdSc, pktsInSc] = scNetwork.receivePackets(controllerK);
                for ii=1:numRcvdSc
                    % construct the packets to be received
                    ncs_doHandlePacket(ncs, timestamp, pktsInSc(ii));
                end
                params = [];
                %pktsIn = cell(1, numel(packetList));
%                 for ii=1:numel(packetList)
%                     % construct the packets to be received
%                     packet = ncs_pktCreate(packetList{ii}.src, packetList{ii}.dst, packetList{ii}.payload);
%                     pktsIn{ii} = packet;                
%                     ncs_doHandlePacket(ncs, timestamp, pktsOut{ii});
%                 end            
                params = [];
%                 if k == uint64(1 / configStruct.plantSamplingInterval) % 1 seconds
% %                     params.samplingInterval = 1/100;
% %                     params
% %                     ratio = ratio / 2;
% 
%                 elseif k == uint64(20 / configStruct.plantSamplingInterval) % 2 seconds
%                     params.samplingInterval = 1/50;
% %                     params
%                     ratio = ratio * 2;
%                 elseif k == uint64(5 / configStruct.plantSamplingInterval) % 5 seconds
%                     params.samplingInterval = 1/25;
% %                     params
%                     ratio = ratio * 2;
% %                       res = ncs_doHandleQocTarget(ncs, 0);
% %                       res
% %                       ratio = ratio * 5
%                 elseif k == uint64(15 / configStruct.plantSamplingInterval) % 15 seconds
%                     params.samplingInterval = 1/100;
% %                     params
%                     ratio = ratio / 4;
%                       % factor 10, expected a sampling rate of 200 Hz
% %                       res = ncs_doHandleQocTarget(ncs, 1);
% %                       res.samplingInterval
% %                       ratio = ratio / 10
%                 end
%                 %params.samplingInterval = samplingIntervals(j);
%                 if controllerK == 1
%                     params.caDelayProbs = zeros(1, 10) + 1/10;
%                 else
%                     params.caDelayProbs = [0 1 0 0 0 0 0 0 0 0];
%                 end
                seconds = ConvertToSeconds(double(timestamp));
                [pktsOut, ~, stats] = ncs_doLoopStep(ncs, timestamp, params);
                controlErrors{r, j}(1, controllerK) = stats.actual_control_error;
                controlErrors{r, j}(2, controllerK) = stats.estimated_control_error;
                                
                % do directly emulate what's happening in Omnet, we must do a
                % bit more
                % extract all the information from the packets
%                 packetList = cell(1, numel(pktsOut));
                for ii=1:numel(pktsOut)
                    if ncs_pktGetDstAddr(pktsOut{ii}) == 1 % to the actuator
                        % send the packet
                        caNetwork.sendPacket(pktsOut{ii}, DiscreteSample(configStruct.caDelayProbs, 1)-1, controllerK);
                    else
                       scNetwork.sendPacket(pktsOut{ii}, 1, controllerK); 
                    end
%                     packetData.src = ncs_pktGetSrcAddr(pktsOut{ii});
%                     packetData.dst = ncs_pktGetDstAddr(pktsOut{ii});
%                     packetData.payload = ncs_pktGetPayload(pktsOut{ii});                
%                     packetData.id = ncs_pktGetId(pktsOut{ii});
%                     packetData.isAck = ncs_pktIsAck(pktsOut{ii});
% 
%                     packetList{ii} = packetData;
                end

                controllerK = controllerK + 1
                nextControllerInvocation = nextControllerInvocation + ratio;        
             end
        end
        [costs{r, j}, controllerStats{r, j}, plantStats{r, j}] = ncs_finalize(ncs);
        trueStates{r, j} = plantStats{r, j}.trueStates;    
    end
end
angleDeviationsDeg = cell(numRepetitions, numPacketRates);
maxAngleDeviations = zeros(numRepetitions, numPacketRates);
posDeviations = cell(numRepetitions, numPacketRates);
maxPosDeviations = zeros(1, numPacketRates);

avgControlErrors = zeros(numRepetitions, numPacketRates);
maxControlError = zeros(numRepetitions, 1);
minControlError = inf(numRepetitions, 1); % min possible control error is zero, upper bound, yet unrealistic
%qocs = zeros(1, numPacketRates);
skip = 5; % wait for 5 seconds

for r=1:numRepetitions
    for j=1:numPacketRates
        angleDeviationsDeg{r, j} = rad2deg(trueStates{r, j}(3, :))-180;
        maxAngleDeviations(r, j) = max(abs(angleDeviationsDeg{r, j}));
        posDeviations{r, j} = trueStates{r, j}(1,:);
        maxPosDeviations(r, j) = max(abs(posDeviations{r, j}));

        % get the actual control error but discard first five seconds
        numEl = numel(controlErrors{r, j}(1, :));
        timeInstants = samplingIntervals(j)*(0:numEl-1);

        errors = controlErrors{r, j}(1, timeInstants > skip);
        avgControlErrors(r, j) = mean(errors);
        minControlError(r) = min(minControlError(r), avgControlErrors(r, j));
        maxControlError(r) = max(maxControlError(r), avgControlErrors(r, j));        
    end
end

% TODO: compute average per time step and sampling rate, so far only
% average per run and sampling rate computed


% translate into qoc
qocs = arrayfun(@(avg) abs((avg - maxControlError) / (maxControlError - minControlError)), avgControlErrors);

legends = cell(1, numPacketRates);
figure();
for j=1:numPacketRates
    hold on;
    legends{j} = sprintf('%d Hz', packetRates(j));
    numEl = size(trueStates{j}(3, :), 2);
    stairs(configStruct.plantSamplingInterval*(0:numEl-1), rad2deg(trueStates{j}(3, :))-180, 'LineWidth', 2);
    grid on;
end
ylabel('Angle Deviation (in degrees)');
xlabel('Simulation time (in seconds)');
legend(legends, 'Location', 'westoutside', 'FontSize', 14);

figure();
for j=1:numPacketRates
    hold on;
    numEl = size(trueStates{j}(1, :), 2);
    stairs(configStruct.plantSamplingInterval*(0:numEl-1), trueStates{j}(1, :), 'LineWidth', 2);
    grid on;
end
ylabel('Position Deviation (in meters)');
xlabel('Simulation time (in seconds)');
legend(legends, 'Location', 'westoutside', 'FontSize', 14);

figure();
for j=1:numPacketRates
    hold on;
    numEl = numel(controlErrors{j}(1, :));
    stairs(samplingIntervals(j)*(0:numEl-1), controlErrors{j}(1, :), 'LineWidth', 2);
    grid on;
end
ylabel('Actual Control Error');
xlabel('Simulation time (in seconds)');
legend(legends, 'Location', 'westoutside', 'FontSize', 14);

figure();
for j=1:numPacketRates
    hold on;
    numEl = numel(controlErrors{j}(2, :));
    stairs(samplingIntervals(j)*(0:numEl-1), controlErrors{j}(2, :), 'LineWidth', 2);
    grid on;
end
ylabel('Estimated Control Error');
xlabel('Simulation time (in seconds)');
legend(legends, 'Location', 'westoutside', 'FontSize', 14);

figure();
plot(packetRates, maxPosDeviations, 'Marker', 'o', 'LineWidth', 2');
ylabel('Max Absolute Position Deviation (in meters)');
xlabel('Controller Sampling Rate (in Hz)');

figure();
plot(packetRates, maxAngleDeviations, 'Marker', 'o', 'LineWidth', 2');
ylabel('Max Absolute Angle Deviation (in degrees)');
xlabel('Controller Sampling Rate (in Hz)');

figure();
plot(packetRates, avgControlErrors, 'Marker', 'o', 'LineWidth', 2');
ylabel('Avg. Control Error');
xlabel('Controller Sampling Rate (in Hz)');

figure();
plot(packetRates, qocs, 'Marker', 'o', 'LineWidth', 2');
ylabel('Quality of Control');
xlabel('Controller Sampling Rate (in Hz)');

