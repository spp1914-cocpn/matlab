classdef NCSTemplate < handle
    % Example code which demonstrates the setup of an networked control
    % system with a updatable transition probability matrix. The class is
    % used to initialize a network containing one filter.
    
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
        function results = simulate(...
                simTime, flipIntervals, estimatorType, setUpRoutine, setUpUpdate)
            % Static method to run the simulation of an updated networked 
            % control system.
            %
            % Parameters:
            %   >> simTime: duration of the simulation in number of time 
            %   steps 
            %   >> flipIntervals: length after which the properties of the
            %   network are changing
            %   >> estimatorType: type determining if and how the
            %   transition probability matrix is updated
            %   >> setUpRoutine: initialized components of the networked
            %   control system
            %   >> setUpUpdate: predefined changes of the network
            % Returns:
            %   << results: array containing the estimation state, the
            %   unscertainties, the real trajectory, the estimation time,
            %   the measurements, the inputs, and the true modes of each 
            %   time step as well as the last transition matrix
            
            % The filter in this control loop is chosen beforehand and
            % initialized here.            
            if strcmp(estimatorType, 'groundtruth')
                [controller, estimator, ~, ~] = ...
                    setUpRoutine.setUp();
            elseif strcmp(estimatorType, 'fixed')
                [controller, ~, estimator, ~] = ...
                    setUpRoutine.setUp();
            elseif strcmp(estimatorType, 'variable')
                [controller, ~, ~, estimator] = ...
                    setUpRoutine.setUp();
            else
                error('** Not a valid estimator type **')
            end
            
            [plant, augmentedPlantModel, sensor, actuator] = setUpRoutine.setUpPlant();
            
            [dimX, dimY, dimU] = setUpRoutine.getDimensions();
            
            [maxMeasurementDelay, controlSequenceLength] = setUpRoutine.getDelay();
            
            caNetwork = CommunicationNetwork(simTime, controlSequenceLength + 1 );
            scNetwork = CommunicationNetwork(simTime, maxMeasurementDelay);
            acNetwork = CommunicationNetwork(simTime, maxMeasurementDelay + 1 );
            
            results = NCSTemplate.initializeResults(dimX, dimU, dimY, simTime);
            
            [plantState, ~] = setUpRoutine.getInitialPlantStateAndEstimate();
            xReal = results.realTrajectory(:,1); % initial plant state
            results.realTrajectory(:,1) = plantState;

            defaultInput = zeros(dimU, 1); % the input if the actuator ran empty
            bufferedInputSequences = repmat(defaultInput, ...
                [1 controlSequenceLength controlSequenceLength]);

            printer = SimulationInfoPrinter('Updated Ncs Example', 1, simTime);
            printer.printSimulationStart();
            
            [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                setUpRoutine.getPackageDelays();
            
            divisor = 2 * flipIntervals;
            for k = 1 : simTime
                printer.printProgress(1, k);
                
                % change delays  
                if mod(k, divisor) == flipIntervals
                    [transitionMatrix, ~] = ...
                        setUpUpdate.getTransitionMatrix();
                    if strcmp(estimatorType, 'groundtruth')
                        estimator.setModeTransitionMatrix(transitionMatrix);
                    end
                    controller = setUpUpdate.updateTransitionMatrix(controller);
                    
                    [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                        setUpUpdate.getPackageDelays();
                elseif mod(k, divisor) == 0
                    [transitionMatrix, ~] = ...
                        setUpRoutine.getTransitionMatrix();
                    if strcmp(estimatorType, 'groundtruth')
                        estimator.setModeTransitionMatrix(transitionMatrix);
                    end
                    controller = setUpRoutine.updateTransitionMatrix(controller);
                    
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
                results.measurements(:, k) = measurement;

                % receive measurements and mode observations
                [~, arrivedAcPackets] = acNetwork.receivePackets(k);
                [~, arrivedScPackets] = scNetwork.receivePackets(k); 

                % use this to update the controller's state estimate
                if k > 1
                    [modeObservations, modeDelays] = ...
                        NCSTemplate.processAcPackets(k, arrivedAcPackets);
                    [measurements, measurementDelays] = ...
                        NCSTemplate.processScPackets(arrivedScPackets);

                    % distribute the possible inputs to all modes
                    modeSpecificInputs = arrayfun(@(mode) bufferedInputSequences(:, mode, mode), ...
                            1:controlSequenceLength, 'UniformOutput', false);
                    % include the default input for the last mode 
                    augmentedPlantModel.setSystemInput([cell2mat(modeSpecificInputs) defaultInput]);
                    timer = tic;
                    estimator.step(augmentedPlantModel, sensor, ...
                        measurements, measurementDelays, modeObservations, modeDelays);
                    results.time(1,k) = toc(timer);
                    
                    previousModeEstimate = estimator.getPreviousModeEstimate();
                else
                    % initial timestep, no update required
                     previousModeEstimate = estimator.getPreviousModeEstimate(false);
                end
                [results.estimatedTrajectory(:,k), results.covs(:, :, k)] = ...
                    estimator.getStateMeanAndCov();
                
                % now compute the input sequence based on state estimate
                inputs = controller.computeControlSequence(...
                    estimator.getState(), previousModeEstimate, k);
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

                results.appliedInputs(:, k) = actualInput;
                results.trueModes(:, k) = actualMode;

                % apply the input and obtain new plant state (i.e. proceed to k+1)
                plant.setSystemInput(actualInput);
                xReal = plant.simulate(xReal);

                results.realTrajectory(:, k + 1) = xReal; 
            end % main simulation loop
            % compute and return costs
            
            results.costs = setUpRoutine.ComputeCosts(...
              results.realTrajectory, results.appliedU);

            results.lastTransitionMatrix = ...
                estimator.getModeTransitionMatrix();
         
            printer.printSimulationEnd();
        end
        
        function results = initializeResults(dimX, dimU, dimY, simTime)
            % Method to initialize results array.
            %
            % Parameters:
            %   >> dimX: length of state vector
            %   >> dimU: length of input vector
            %   >> dimY: length of measurement vector
            %   >> simTime: duration of the simulation in discrete time
            %   steps
            % Returns:
            %   << results: array containing the estimation state, the
            %   unscertainties, the real trajectory, the estimation time,
            %   the measurements, the inputs, and the true modes of each 
            %   time step as well as the last transition matrix
            results = [];
            results.estimatedTrajectory = zeros(dimX, simTime+1);
            results.covs = zeros(dimX, dimX, simTime + 1);
            results.time = zeros(1, simTime+1);
            results.lastTransitionMatrix = [];
            results.appliedU = zeros(dimU, simTime);
            results.trueModes = zeros(1, simTime + 1);
            results.measurements = zeros(dimY, simTime);
            results.realTrajectory = zeros(dimX, simTime + 1);
        end
    end
    
    methods (Access = private, Static)
        %% see other Examples for more information.
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