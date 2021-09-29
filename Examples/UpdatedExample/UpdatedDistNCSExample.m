classdef UpdatedDistNCSExample < handle
    
    methods (Static, Access = public)
        function [groundTruthResults, fixedResults, variableResults] = simulate(...
                simTime, flipIntervals, distributionUpdateIntervals)
            if nargin == 2
                distributionUpdateIntervals = 10;
            end
            
            setUpRoutine = SetUpScenarios(simTime, "Histogram");
            [controller, groundTruthEstimator, fixedEstimator, variableEstimator] = ...
                setUpRoutine.setUp();
            [plant, augmentedPlantModel, sensor, actuator] = setUpRoutine.setUpPlant();
            
            [dimX, dimY, dimU] = setUpRoutine.getDimensions();
            
            [maxMeasurementDelay, controlSequenceLength] = setUpRoutine.getDelay();
            
            caNetwork = CommunicationNetwork(simTime, controlSequenceLength + 1 );
            scNetwork = CommunicationNetwork(simTime, maxMeasurementDelay);
            acNetwork = CommunicationNetwork(simTime, maxMeasurementDelay + 1 );
            
            groundTruthResults = UpdatedDistNCSExample.initializeResults(dimX, simTime);
            groundTruthResults.appliedU = zeros(dimU, simTime);
            groundTruthResults.trueModes = zeros(1, simTime + 1);
            groundTruthResults.measurements = zeros(dimY, simTime);
            groundTruthResults.realTrajectory = zeros(dimX, simTime + 1);
            
            fixedResults = UpdatedDistNCSExample.initializeResults(dimX, simTime);
            variableResults = UpdatedDistNCSExample.initializeResults(dimX, simTime);
            
            [plantState, ~] = setUpRoutine.getInitialPlantStateAndEstimate();
            xReal = groundTruthResults.realTrajectory(:,1); % initial plant state
            groundTruthResults.realTrajectory(:,1) = plantState;

            defaultInput = zeros(dimU, 1); % the input if the actuator ran empty
            bufferedInputSequences = repmat(defaultInput, ...
                [1 controlSequenceLength controlSequenceLength]);

            printer = SimulationInfoPrinter('Updated (Distribution) Ncs Example', 1, simTime);
            %printer.turnOff;
            printer.printSimulationStart();
            
            [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                setUpRoutine.getPackageDelays();

            save 'caDelays.mat' caPacketDelays;
            save 'acDelays.mat' acPacketDelays;
            
            setUpUpdate = ChangeUpScenarios(simTime, 6, 'Histogram');
            
            delayHistory = cell(1, uint8(distributionUpdateIntervals / 4));
            modeHistogram = zeros(1, controlSequenceLength + 1);
            transitionMatrix = [];
            
            % Array of zeros and ones.
            % Displays which control packets of the last few steps was sent
            % but not successfully acknowledged (1)
            % lengthcaPacketHistory is calculated by maxDelay of caNetwork
            % + maxDelay of acNetwork + 1 per network because to be lost, the packets
            % needs be longer away than the sum of the max delays of both
            % networks.
            lengthcaPacketHistory = controlSequenceLength + maxMeasurementDelay + 4;
            sentcaPacketHistory = zeros(1, lengthcaPacketHistory);
            for k = 1 : simTime
                printer.printProgress(1, k);
                
                % flip network characteristics
                if mod(k, (2 * flipIntervals)) == flipIntervals
                    [transitionMatrix, delayProbs] = ...
                        setUpUpdate.getTransitionMatrix();
                    groundTruthEstimator.setModeTransitionMatrix(transitionMatrix);
                    [controller, ~, ~, ~] = setUpUpdate.setUp();
                    %variableEstimator.updateDefaultDelays(delayProbs);
                    
                    [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                        setUpUpdate.getPackageDelays();
                elseif mod(k, (2 * flipIntervals)) == 0
                    [transitionMatrix, delayProbs] = ...
                        setUpRoutine.getTransitionMatrix();
                    groundTruthEstimator.setModeTransitionMatrix(transitionMatrix);
                    [controller, ~, ~, ~] = setUpRoutine.setUp();
                    %variableEstimator.updateDefaultDelays(delayProbs);  

                    [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                       setUpRoutine.getPackageDelays();
                end
                
                % take a measurement and transmit it
                measurement = sensor.simulate(xReal);  
                scPacket = DataPacket(measurement, k);
                % add source and destination (addresses are arbitraty)
                scPacket.sourceAddress = 3;
                scPacket.destinationAddress = 2;
                
                scNetwork.sendPacket(scPacket, scPacketDelays(k), k);                
                groundTruthResults.measurements(:, k) = measurement;

                % lost ACK and lost control packet cannot be distinguished.
                % There is some penalty needed! TODO: What is a acceptable
                % penalty.
                if sentcaPacketHistory(1) == 1
                    modeHistogram(1, controlSequenceLength + 1) = ...
                            modeHistogram(1, controlSequenceLength + 1) + 1;
                    sentcaPacketHistory(1) = 0;
                end
                
                sentcaPacketHistory = circshift(sentcaPacketHistory, lengthcaPacketHistory - 1);
                
                % receive measurements and mode observations
                [~, arrivedAcPackets] = acNetwork.receivePackets(k);
                [~, arrivedScPackets] = scNetwork.receivePackets(k); 

                % use this to update the controller's state estimate
                if k > 1
                    [modeObservations, modeDelays] = ...
                        UpdatedDistNCSExample.processAcPackets(k, arrivedAcPackets);
                    [measurements, measurementDelays] = ...
                        UpdatedDistNCSExample.processScPackets(arrivedScPackets);

                    timer = tic;
                    % fill delays into histogram
                    for observation = 1:length(modeObservations)
                        modeHistogram(1, modeObservations(observation)) = ...
                            modeHistogram(1, modeObservations(observation)) + 1;
                        
                        index = lengthcaPacketHistory - (modeObservations(observation) ...
                            + modeDelays(observation));
                        sentcaPacketHistory(index) = 0;
                    end
                    % modeHistogram
                    % fill shifting window
                    UpdatedDistNCSExample.updateDelayHistory(delayHistory, ...
                        modeObservations);
                    
                    if mod(k, distributionUpdateIntervals) == 0
                        modeHistogram = (1/sum(modeHistogram)) ...
                            * modeHistogram;
                        
                        idx = find(modeHistogram <= 1e-12);
                        normalizedModeProbs = modeHistogram;
                        if ~isempty(idx)
                            %this.warning('NormalizingModeProbabilities', ...
                            %    '** %d of %d mode probablities too small. Perform normalization **', numel(idx), this.numModes); 
                            normalizedModeProbs(idx) = 1e-12;
                            normalizedModeProbs = normalizedModeProbs / sum(normalizedModeProbs);
                        end
                        modeHistogram = normalizedModeProbs;
                        %fprintf('\n histogramm  %d \n', modeHistogram);
                        transitionMatrix = ...
                            Utility.calculateDelayTransitionMatrix(modeHistogram);
                        modeHistogram = zeros(1, controlSequenceLength + 1);
                        for delays = 1:length(delayHistory)
                            for delay = 1 : length(delays)
                                modeHistogram(1, delay) = ...
                            modeHistogram(1, delay) + 1;
                            end
                        end
                        
                        variableEstimator.setModeTransitionMatrix(...
                            transitionMatrix);
                    end
                    variableResults.time = variableResults.time + toc(timer); 
                    

                    % distribute the possible inputs to all modes
                    modeSpecificInputs = arrayfun(@(mode) bufferedInputSequences(:, mode, mode), ...
                            1:controlSequenceLength, 'UniformOutput', false);
                    % include the default input for the last mode 
                    augmentedPlantModel.setSystemInput([cell2mat(modeSpecificInputs) defaultInput]);
                    timer = tic;
                    groundTruthEstimator.step(augmentedPlantModel, sensor, ...
                        measurements, measurementDelays, modeObservations, modeDelays);
                    groundTruthResults.time = groundTruthResults.time + toc(timer);
                    
                    timer = tic;
                    fixedEstimator.step(augmentedPlantModel, sensor, ...
                        measurements, measurementDelays, modeObservations, modeDelays);
                    fixedResults.time = fixedResults.time + toc(timer);
                    
                    timer = tic;
                    variableEstimator.step(augmentedPlantModel, sensor, ...
                        measurements, measurementDelays, modeObservations, modeDelays);                    
                    variableResults.time = variableResults.time + toc(timer); 
                    
                    previousModeEstimate = groundTruthEstimator.getPreviousModeEstimate();
                else
                    % initial timestep, no update required
                     previousModeEstimate = groundTruthEstimator.getPreviousModeEstimate(false);
                end
                [groundTruthResults.estimatedTrajectory(:,k), groundTruthResults.covs(:, :, k)] = ...
                    groundTruthEstimator.getStateMeanAndCov();
                
                [fixedResults.estimatedTrajectory(:,k), fixedResults.covs(:, :, k)] = ...
                    fixedEstimator.getStateMeanAndCov();
                
                [variableResults.estimatedTrajectory(:,k), variableResults.covs(:, :, k)] = ...
                    variableEstimator.getStateMeanAndCov();

                % now compute the input sequence based on state estimate
                inputs = controller.computeControlSequence(...
                    groundTruthEstimator.getState(), previousModeEstimate, k);
                % and store as matrix
                inputSequence = reshape(inputs, [], controlSequenceLength);    
                caPacket = DataPacket(inputSequence, k);
                caPacket.sourceAddress = 2;
                caPacket.destinationAddress = 1;
                % send the packet to the actuator
                caNetwork.sendPacket(caPacket, caPacketDelays(k), k);

                sentcaPacketHistory(lengthcaPacketHistory) = 1;
                
                bufferedInputSequences = circshift(bufferedInputSequences, 1, 3);            
                bufferedInputSequences(:,:, 1) = inputSequence;    

                % receive packets that arrive at the actuator and process them to
                % obtain the actual control input u_k
                % and send ACK if required
                [~, arrivedCaPackets] = caNetwork.receivePackets(k);    
                controllerAck = actuator.processControllerPackets(arrivedCaPackets);
                if ~isempty(controllerAck)     
                    acNetwork.sendPacket(controllerAck, acPacketDelays(k), k);
                end
                [actualInput, actualMode] = actuator.getCurrentInput(k);

                groundTruthResults.appliedInputs(:, k) = actualInput;
                groundTruthResults.trueModes(:, k) = actualMode;

                % apply the input and obtain new plant state (i.e. proceed to k+1)
                plant.setSystemInput(actualInput);
                xReal = plant.simulate(xReal);

                groundTruthResults.realTrajectory(:, k + 1) = xReal; 
            end % main simulation loop
            % compute and return costs
            groundTruthResults.costs = setUpRoutine.ComputeCosts(...
              groundTruthResults.realTrajectory, groundTruthResults.appliedU);

            groundTruthResults.lastTransitionMatrix = ...
                groundTruthEstimator.getModeTransitionMatrix();
           
            fixedResults.lastTransitionMatrix = ...
                fixedEstimator.getModeTransitionMatrix();
            
            variableResults.lastTransitionMatrix = ...
                variableEstimator.getModeTransitionMatrix();

            printer.printSimulationEnd();
        end
        
        function results = initializeResults(dimX, simTime)
            results = [];
            results.estimatedTrajectory = zeros(dimX, simTime+1);
            results.covs = zeros(dimX, dimX, simTime + 1);
            results.time = 0;
            results.lastTransitionMatrix = [];
        end
        
        function newHistory = updateDelayHistory(delayHistory, delays)
            newHistory = circshift(delayHistory, 1, 2);
            newHistory{1} = delays;
        end
    end
    
    methods (Access = private, Static)
        %% processAcPackets
        function [modeObservations, modeDelays] = processAcPackets(timestep, acPackets)
            modeObservations = [];
            modeDelays = [];
            if numel(acPackets) ~= 0
                ackPayloads = cell(1, numel(acPackets));
                ackDelays = cell(1, numel(acPackets));
                [ackPayloads{:}] = acPackets(:).payload;
                [ackDelays{:}] = acPackets(:).packetDelay;
                % get the observed modes from the ACK packets
                ackedPacketTimes = cell2mat(ackPayloads);
                modeDelays = cell2mat(ackDelays);
                % get the time stamps of the ACK packets
                ackTimeStamps = timestep - modeDelays;

                modeObservations = ackTimeStamps - ackedPacketTimes + 1;
            end
        end  

        %% processScPackets
        function [measurements, measDelays] = processScPackets(scPackets)
            % Extract the measurements and corresponding delays from the
            % data packets received from the sensor. 
            %
            % Parameters:
            %   >> scPackets (Array of DataPackets, might be empty)
            %      An array of DataPackets containing measurements taken and transmitted from the sensor.
            %
            % Returns:
            %   << measurements (Matrix, might be empty)
            %      The measurements extracted from the data packets,
            %      column-wise arranged.
            %      Empty matrix is returned if scPackets is empty.
            %
            %   << measDelays (Vector, might be empty)
            %      A vector containing the delays of the received measurements (in timesteps), 
            %      where the i-th element denotes the delay of the i-th measurement.
            %      Empty matrix is returned if scPackets is empty.

            measurements = [];
            measDelays = [];
            if numel(scPackets) ~= 0
                delays = cell(1, numel(scPackets));
                meas = cell(1, numel(scPackets));
                [delays{:}] = scPackets(:).packetDelay;
                [meas{:}] = scPackets(:).payload;
                measurements = cell2mat(meas);
                measDelays = cell2mat(delays);
            end
        end
    end
end