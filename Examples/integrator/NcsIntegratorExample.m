classdef NcsIntegratorExample < handle
    % This class is an example how to simulate a complete networked control
    % system (NCS) (in this case a double integrator plant) in Matlab.
     properties (Access = private)
      %% input properties
        % state cost matrix; (matrix of dimension <dimX> x <dimX>)
        Q = [];
        % control cost matrix; (matrix of dimension <dimU> x <dimU>)
        R = [];
        % system noise covariance; (matrix of dimension <dimX> x <dimX>)
        W = []; 
        % measurement noise covariance; (matrix of dimension <dimX> x <dimX>)
        V = [];
        % simulation time; (positive integer)
        simTime = [];
        % delay probability distribution of the controller-actuator network;
        % (vector of dimension >= <caPacketLength>)
        caDelayProb = [];
        % control sequence length; (positive integer)
        controlSeqLength = [];
        % actual delays of the controller-actuator network; 
        % (vector of dimension <simTime>)
        caPacketDelays = [];
        % delay probability distribution of the sensor-controller network;
        % (vector)
        scDelayProb = [];
        % actual delays of the sensor-controller network; 
        % (vector of dimension <simTime>)
        scPacketDelays = [];
        % maximum allowed measurement delay; (positive integer or zero)
        scMaxDelay = [];
      %% derived properties and components
        % dimension of the system state; (positive integer)
        dimX = -1;
        % dimension of the control input; (positive integer)
        dimU = -1;
        % dimension of the measurement; (positive integer)
        dimY = -1;
        % handle to the controller;
        controller = [];
        % handle to the controller-actuator network;
        caNetwork = [];
        % handle to the actuator; 
        actuator = [];
        plant = [];
        augmentedPlantModel = []; % the MJLS model used by the stimator
        sensor = [];
        scNetwork = [];
        acNetwork = [];
        estimator = [];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods 
        %% NcsIntegratorExample
        function this = NcsIntegratorExample(A,B,C,Q,R,W,V, simTime, caDelayProb,controlSeqLength,scMaxDelay)

            this.dimX = size(A,1);
            this.dimU = size(B,2);
            this.dimY = size(C,1);
            this.Q = Q;
            this.R = R;
            this.W = W;
            this.V = V;
            this.simTime = simTime;
            this.controlSeqLength = controlSeqLength;
            this.scMaxDelay = scMaxDelay;
    
            
            Validator.validateDiscreteProbabilityDistribution(caDelayProb);
            this.caDelayProb = Utility.truncateDiscreteProbabilityDistribution(caDelayProb, controlSeqLength + 1);
                      
            transitionMatrix = Utility.calculateDelayTransitionMatrix(this.caDelayProb);      
            
            % define the network for measurements, control inputs and ACKs
            % from the actuator, respectively
            this.caNetwork = CommunicationNetwork(this.simTime, this.controlSeqLength + 1);
            this.scNetwork = CommunicationNetwork(this.simTime,this.scMaxDelay);
            this.acNetwork = CommunicationNetwork(this.simTime,this.scMaxDelay + 1);

            % initialize plant and sensor
            this.actuator = BufferingActuator(this.controlSeqLength, Inf, zeros(this.dimU, 1));
            this.plant = LinearPlant(A, B, W);
            this.sensor = LinearMeasurementModel(C);
            % zero mean measurement noise
            this.sensor.setNoise(Gaussian(zeros(size(C, 1), 1), V));    

            % now define the controller and the state estimator used by it
            numModes = this.controlSeqLength + 1;
            modeFilters = arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:numModes, 'UniformOutput', false);
            this.estimator = DelayedModeIMMF(modeFilters, transitionMatrix, this.scMaxDelay);

            this.augmentedPlantModel = JumpLinearSystemModel(numModes, ...
                arrayfun(@(~) LinearPlant(A, B, W), ...
                        1:numModes, 'UniformOutput', false));
            this.augmentedPlantModel.setNoise(Gaussian(zeros(size(C, 1), 1), V))

            this.controller = InfiniteHorizonController(A, B, this.Q, this.R, ...
                this.caDelayProb,this.controlSeqLength);
            switch this.controller.status
              case -2
                error(['** Stabilizing infinite time horizon controller', ...
                  'does not exist. **']);
              case -1
                disp(['** Warning: Numerical Problems; Controller gain might',...
                  'be corrupted. Working with best gain found. **']);
              case 0
                disp(['** Warning: Controller diverged or convergence was not ',...
                  'detected. Working with best gain found. **']);
              case 1
                % everything fine
            end
        end

        %% simulate
        function results = simulate(this, initialPlantState, initialEstimate, caPacketDelays,scPacketDelays, acPacketDelays)
            % reset reusable components
            this.controller.reset();
            this.actuator.reset();
            this.caNetwork.reset();
            this.scNetwork.reset();

            % set the initial estimate
            this.estimator.setState(initialEstimate);

            % construct output structure
            results.realTrajectory = zeros(this.dimX, this.simTime+1);

            results.estimatedTrajectory = zeros(this.dimX, this.simTime+1);
            results.covs = zeros(this.dimX, this.dimX, this.simTime +1);
            results.appliedInputs = zeros(this.dimU, this.simTime);
            results.trueModes = zeros(1, this.simTime+1);
            results.measurements = zeros(this.dimY, this.simTime);

            xReal = results.realTrajectory(:,1); % initial plant state

            results.realTrajectory(:,1) = initialPlantState;

            defaultInput = zeros(this.dimU, 1); % the input if the actuator ran empty
            bufferedInputSequences = repmat(defaultInput, [1 this.controlSeqLength this.controlSeqLength]);

            printer = SimulationInfoPrinter('Ncs Double Integrator Example', 1, this.simTime);
            printer.printSimulationStart();

          % main simulation loop
          for k = 1 : this.simTime
            printer.printProgress(1, k);

            % take a measurement and transmit it
            measurement = this.sensor.simulate(xReal);  
            scPacket = DataPacket(measurement, k);
            % add source and destination (addresses are arbitraty)
            scPacket.sourceAddress = 3;
            scPacket.destinationAddress = 2;

            this.scNetwork.sendPacket(scPacket, scPacketDelays(k), k);                
            results.measurements(:, k) = measurement;

            % receive measurements and mode observations
            [~, arrivedAcPackets] = this.acNetwork.receivePackets(k);
            [~, arrivedScPackets] = this.scNetwork.receivePackets(k); 

            % use this to update the controller's state estimate
            if k > 1
                [modeObservations, modeDelays] = NcsIntegratorExample.processAcPackets(k, arrivedAcPackets);
                [measurements, measDelays] = NcsIntegratorExample.processScPackets(arrivedScPackets);

                % distribute the possible inputs to all modes
                modeSpecificInputs = arrayfun(@(mode) bufferedInputSequences(:, mode, mode), ...
                        1:this.controlSeqLength, 'UniformOutput', false);
                % include the default input for the last mode 
                this.augmentedPlantModel.setSystemInput([cell2mat(modeSpecificInputs) defaultInput]);

                this.estimator.step(this.augmentedPlantModel, this.sensor, ...
                    measurements, measDelays, modeObservations, modeDelays);

                previousModeEstimate = this.estimator.getPreviousModeEstimate();
            else
                % initial timestep, no update required
                 previousModeEstimate = this.estimator.getPreviousModeEstimate(false);
            end
            [results.estimatedTrajectory(:,k), results.covs(:, :, k)] = this.estimator.getStateMeanAndCov();

            % now compute the input sequence based on state estimate
            inputs = this.controller.computeControlSequence(this.estimator.getState(), previousModeEstimate, k);
            % and store as matrix
            inputSequence = reshape(inputs, [], this.controlSeqLength);    
            caPacket = DataPacket(inputSequence, k);
            caPacket.sourceAddress = 2;
            caPacket.destinationAddress = 1;
            % send the packet to the actuator
            this.caNetwork.sendPacket(caPacket, caPacketDelays(k), k);

            bufferedInputSequences = circshift(bufferedInputSequences, 1, 3);            
            bufferedInputSequences(:,:, 1) = inputSequence;    

            % receive packets that arrive at the actuator and process them to
            % obtain the actual control input u_k
            % and send ACK if required
            [~, arrivedCaPackets] = this.caNetwork.receivePackets(k);    
            controllerAck = this.actuator.processControllerPackets(arrivedCaPackets);
            if ~isempty(controllerAck)     
                this.acNetwork.sendPacket(controllerAck, acPacketDelays(k), k);
            end
            [actualInput, actualMode] = this.actuator.getCurrentInput(k);

            results.appliedInputs(:, k) = actualInput;
            results.trueModes(:, k) = actualMode;

            % apply the input and obtain new plant state (i.e. proceed to k+1)
            this.plant.setSystemInput(actualInput);
            xReal = this.plant.simulate(xReal);

            results.realTrajectory(:, k + 1) = xReal; 
          end % main simulation loop
          % compute and return costs
          results.costs = this.ComputeCosts(results.realTrajectory, results.appliedInputs);

          printer.printSimulationEnd();
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

    methods (Access = private)


      %% ComputeCosts
      function costs = ComputeCosts(this,xTrajectory,uSeq)
      %========================================================================
      % ComputeCosts: computes average costs per time step according to an
      % LQG cost function
      %
      % Synopsis:
      % costs = ComputeCosts(xTrajectory,uSeq);
      % 
      % Input Parameters: 
      % - xTrajectory: state trajectory; 
      %                (matrix of dimension <dimX> x <simTime+1>)
      % - uSeq: sequence of applied control inputs;
      %         (matrix of dimension <dimU> x <simTime>)
      %
      % Output Parameters::
      % - costs: average costs per time step; (positive integer)
      %========================================================================
         costs = Utility.computeLQGCosts(this.simTime, xTrajectory, uSeq, this.Q, this.R) / this.simTime;     
      end
    end
end