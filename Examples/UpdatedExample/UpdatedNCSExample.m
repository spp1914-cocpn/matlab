classdef UpdatedNCSExample < handle
    % Example code which demonstrates the setup of an networked control
    % system with a updatable transition probability matrix. The class is
    % used to initialize a network containing several filters to be compared.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
    %                        Joanna Mueller <joanna.mueller@student.kit.edu>
    %
    %                        Institute for Anthropomatics and Robotics
    %                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
    %                        Karlsruhe Institute of Technology (KIT), Germany
    %
    %                        https://isas.iar.kit.edu
    %
    %    This program is free software: you can redistribute it and/or modify
    %    it under the terms of the GNU General Public License as published by
    %    the Free Software Foundation, either version 3 of the License, or
    %    (at your option) any later version.
    %
    %    This program is distributed in the hope that it will be useful,
    %    but WITHOUT ANY WARRANTY; without even the implied warranty of
    %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %    GNU General Public License for more details.
    %
    %    You should have received a copy of the GNU General Public License
    %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    methods (Static, Access = public)
        function [groundTruthResults, fixedResults, variableResults] = simulate(...
                simTime, flipIntervals, estimatorType)
            % Static method to run the simulation of an updated networked 
            % control system.
            %
            % parameters:
            %   >> simTime: duration of the simulation in number of time 
            %   steps 
            %   >> flipIntervals: length after which the properties of the
            %   network are changing
            %   >> estimatorType: type determining if and how the
            %   transition probability matrix is updated
            % Returns:
            %   << results: array containing the estimation state, the
            %   unscertainties, the estimation time,
            %   as well as the last transition matrix
            setUpRoutine = SetUpScenarios(simTime, estimatorType);
            [controller, groundTruthEstimator, fixedEstimator, variableEstimator] = ...
                setUpRoutine.setUp();
            [plant, augmentedPlantModel, sensor, actuator] = setUpRoutine.setUpPlant();
            
            [dimX, dimY, dimU] = setUpRoutine.getDimensions();
            
            [maxMeasurementDelay, controlSequenceLength] = setUpRoutine.getDelay();
            
            caNetwork = CommunicationNetwork(simTime, controlSequenceLength + 1 );
            scNetwork = CommunicationNetwork(simTime, maxMeasurementDelay);
            acNetwork = CommunicationNetwork(simTime, maxMeasurementDelay + 1 );
            
            groundTruthResults = UpdatedNCSExample.initializeResults(dimX, simTime);
            groundTruthResults.appliedU = zeros(dimU, simTime);
            groundTruthResults.trueModes = zeros(1, simTime + 1);
            groundTruthResults.measurements = zeros(dimY, simTime);
            groundTruthResults.realTrajectory = zeros(dimX, simTime + 1);
            
            fixedResults = UpdatedNCSExample.initializeResults(dimX, simTime);
            variableResults = UpdatedNCSExample.initializeResults(dimX, simTime);
            
            [plantState, ~] = setUpRoutine.getInitialPlantStateAndEstimate();
            xReal = groundTruthResults.realTrajectory(:,1); % initial plant state
            groundTruthResults.realTrajectory(:,1) = plantState;

            defaultInput = zeros(dimU, 1); % the input if the actuator ran empty
            bufferedInputSequences = repmat(defaultInput, ...
                [1 controlSequenceLength controlSequenceLength]);

            printer = SimulationInfoPrinter('Updated Ncs Example', 1, simTime);
            printer.printSimulationStart();
            
            [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                setUpRoutine.getPackageDelays();
            
            setUpUpdate = ChangeUpScenarios(simTime, 1, estimatorType); 
            setUpUpdate2 = ChangeUpScenarios(simTime, 6, estimatorType);   
            for k = 1 : simTime
                printer.printProgress(1, k);
                
                % change delays  
                if mod(k, (3 * flipIntervals)) == flipIntervals
                    [transitionMatrix, ~] = ...
                        setUpUpdate.getTransitionMatrix();
                    groundTruthEstimator.setModeTransitionMatrix(transitionMatrix);
                    %controller = setUpUpdate.updateTransitionMatrix(controller);
                    
                    [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                        setUpUpdate.getPackageDelays();
                elseif mod(k, (3 * flipIntervals)) == (2 * flipIntervals)
                    [transitionMatrix, ~] = ...
                        setUpUpdate2.getTransitionMatrix();
                    groundTruthEstimator.setModeTransitionMatrix(transitionMatrix);
                    %controller = setUpUpdate2.updateTransitionMatrix(controller);
                    
                    [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                        setUpUpdate2.getPackageDelays();
                elseif mod(k, (3 * flipIntervals)) == 0
                    [transitionMatrix, ~] = ...
                        setUpRoutine.getTransitionMatrix();
                    groundTruthEstimator.setModeTransitionMatrix(transitionMatrix);
                    %controller = setUpRoutine.updateTransitionMatrix(controller);
                    
                    [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                       setUpRoutine.getPackageDelays();
                end
                
                % take a measurement and transmit it
                measurement = sensor.simulate(xReal);  
                scPacket = DataPacket(measurement, k);
                % add source and destination (addresses are arbitraty)
                scPacket.sourceAddress = 3;
                scPacket.destinationAddress = 2;
                % add source and destination (addresses are arbitraty)
                scNetwork.sendPacket(scPacket, scPacketDelays(k), k);                
                groundTruthResults.measurements(:, k) = measurement;

                % receive measurements and mode observations
                [~, arrivedAcPackets] = acNetwork.receivePackets(k);
                [~, arrivedScPackets] = scNetwork.receivePackets(k); 

                % use this to update the controller's state estimate
                if k > 1
                    [modeObservations, modeDelays] = ...
                        UpdatedNCSExample.processAcPackets(k, arrivedAcPackets);
                    [measurements, measurementDelays] = ...
                        UpdatedNCSExample.processScPackets(arrivedScPackets);

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
            % Method to initialize results array.
            %
            % Parameters:
            %   >> dimX: length of state vector
            %   >> simTime: duration of the simulation in discrete time
            %   steps
            % Returns:
            %   << results: array containing the estimation state, the
            %   unscertainties, the estimation time,
            %   as well as the last transition matrix
            results = [];
            results.estimatedTrajectory = zeros(dimX, simTime+1);
            results.covs = zeros(dimX, dimX, simTime + 1);
            results.time = 0;
            results.lastTransitionMatrix = [];
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
