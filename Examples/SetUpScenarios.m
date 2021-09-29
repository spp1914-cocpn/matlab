classdef SetUpScenarios
    % Class for initializing controllers, estimators, networks, and plant.
    % In this class, there are default setups for each component which can
    % be used in example code.
    
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
    properties(Access = protected, Constant)
        % Default values for size of the control packets and considered
        % measurement delays.
        controlSequenceLength = 4;
        maxMeasurementDelay = 6;
    end
    
    properties(Access = protected)
        % Probability distributions for the three communication channels
        % considered in this framework.
        caDelayProb = [];
        scDelayProb = [];
        acDelayProb = [];
        
        % Drawn packet delays for each channel and every time step.
        caPacketDelays = [];
        scPacketDelays = [];
        acPacketDelays = [];
        
        % Configuration describing the behavior of the plant (e.g., system 
        % matrix, input matrix, and measurement matrix)
        config = [];
        
        % True transition probabilities in the beginning of the simulation.
        transitionMatrix = [];
        
        % Initial plant state and estimate. Estimate is usually described
        % by a Gaussian distribution i.e., a mean and a covariance.
        initialPlantState = [];
        initialEstimate = [];
        
        % Duration time of the simulation in steps.
        simTime;
        % Initial guess of the delay distribution of the
        % controller-actuator link for variable and fixed estimators.
        startingNetworkDelayDistribution = [];
        
        % Number of modes of the Markov jump system.
        numModes;
        % Set of filters with a filter for each mode used in the IMMF.
        modeFilters;
        
        % String describing if either the Random updater, the RKL updater,
        % the ConvexOptimization updater, or the general updater (super
        % class which does not update the transition probability matrix) is
        % used.
        estimatorType = '';
    end
    
    methods(Access = public)
        function this = SetUpScenarios(simTime, estimatorType)
            % Constructor of class which is used to initialize the
            % components of the updated networked control system.
            %
            % parameters:
            %   >> simTime: duration of the simulation in number of time 
            %   steps 
            %   >> estimatorType: type determining if and how the
            %   transition probability matrix is updated
            % Returns:
            %   << this: object of type SetUpScenarios
            this.simTime = simTime;
            this.estimatorType = estimatorType;
                
            this.numModes = this.controlSequenceLength + 1;
            % The mode filters are initializes as extended Kalman filters. 
            this.modeFilters = arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), ...
            1:this.numModes, 'UniformOutput', false);
            
            % In the SimulatedNetworks.mat is a set of different example
            % distributions which are used in this example setup.
            NetworkDelayDistribution = load('SimulatedNetworks.mat');
            NetworkDelayDistribution = NetworkDelayDistribution.NetworkDelayDistribution;
            this.startingNetworkDelayDistribution = NetworkDelayDistribution(8,:);
            % Ground truth distribution.
            NetworkDelayDistribution = NetworkDelayDistribution(1,:);
            
            this.caDelayProb = NetworkDelayDistribution;
            this.scDelayProb = NetworkDelayDistribution;
            this.acDelayProb = NetworkDelayDistribution;
            
            this.caPacketDelays = DiscreteSample(this.caDelayProb,this.simTime)-1;
            this.scPacketDelays = DiscreteSample(this.scDelayProb,this.simTime)-1;
            this.acPacketDelays = DiscreteSample(this.acDelayProb,this.simTime)-1;

            % System properties: discrete time double integrator plant (+noise) with default
            % parameters 
            this.config = [];
            this.config.mass = 0.1; % kg
            this.config.samplingInterval = 0.01; % sec, 100 Hz
            this.config.initialPlantState = [100; 0];
            this.config.plantNoiseCov = 0.5 * eye(2); % taken from ACC 2013 paper by Jörg and Maxim
            this.config.measNoiseCov = 0.2^2; % taken from ACC 2013 paper by Jörg and Maxim
            
            this.config.Q = eye(2); 
            this.config.R = 1;
            this.config.A = [0 1; 0 0];
            this.config.B = [0; 1 / this.config.mass];
            this.config.C = [1 0];
            W_cont = 0.1 * eye(2); % process noise, continuous time
            
            this.config.contSys = ss(this.config.A, this.config.B, this.config.C, 0); % no D matrix
            this.config.discreteSys = c2d(this.config.contSys, this.config.samplingInterval);

            this.config.W = integral(@(x) expm(this.config.A*x) * W_cont * ...
                expm(this.config.A'*x), 0, this.config.samplingInterval, ...
                'ArrayValued', true);

            this.config.V = 0.2^2; % variance of the measurement noise
            
            this.caDelayProb = Utility.truncateDiscreteProbabilityDistribution(...
                this.caDelayProb, this.numModes );
            this.transitionMatrix = Utility.calculateDelayTransitionMatrix(this.caDelayProb);
            
            this.initialPlantState = this.config.initialPlantState;
            this.initialEstimate = Gaussian(this.initialPlantState, 0.5 * eye(2));
        end
        
        function [controller, groundtruthFilter, fixedFilter, ...
                variableFilter] = setUp(this)
            % This method sets up the controller and filters used in the
            % example code.
            %
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << controller: infinite horizon controller object
            %   << groundtruthFilter: delayed mode IMMF with true
            %   transition matrix
            %   << fixedFilter: delayed mode IMMF with fixed transition
            %   matrix which is not the ground truth
            %   << variableFilter: delayed mode IMMF with a updated
            %   transition matrix which is generally not the true
            %   transition matrix
            controller = this.setUpController();
            % Filter with the true transition matrix which ca change over
            % time.
            groundtruthFilter = this.setUpGroundtruthEstimator();
            % Filter with a fixed transition matrix.
            fixedFilter = this.setUpFixedEstimator();
            % Filter with a transition matrix updated over time.
            variableFilter = this.setUpVariableEstimator();
        end
        
        function estimator = setUpGroundtruthEstimator(this) 
            % This method sets up the ground truth filter used in the
            % example code.
            %
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << estimator: delayed mode IMMF with true
            %   transition matrix
            estimator = DelayedModeIMMF(this.modeFilters, this.transitionMatrix, ...
                this.maxMeasurementDelay);
            
            estimator.setState(this.initialEstimate);
        end
        
        function estimator = setUpFixedEstimator(this)
            % This method sets up the fixed filter used in the example 
            % code.
            %
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << estimator: delayed mode IMMF with fixed transition
            %   matrix which is not the ground truth
            delayProbs = Utility.truncateDiscreteProbabilityDistribution( ...
                this.startingNetworkDelayDistribution,  this.numModes );
            localTransitionMatrix = Utility.calculateDelayTransitionMatrix(delayProbs);
                
            estimator = DelayedModeIMMF(this.modeFilters, localTransitionMatrix, ...
                this.maxMeasurementDelay);
            
            estimator.setState(this.initialEstimate);
        end
        
        function estimator = setUpVariableEstimator(this) 
            % This method sets up the variable filter used in the
            % example code. The filter's transition matrix is either
            % updater by the RKL updater, the ConvexOptimization updater,
            % the Random updater, or the Super class updater.
            %
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << estimator: delayed mode IMMF with a updated
            %   transition matrix which is generally not the true
            %   transition matrix
            initialWeights = (1/ this.numModes) * ones(1, this.numModes);

            delayProbs = Utility.truncateDiscreteProbabilityDistribution( ...
                this.startingNetworkDelayDistribution,  this.numModes );
            
            localTransitionMatrix = Utility.calculateDelayTransitionMatrix(delayProbs);
            if strcmp(this.estimatorType, 'RKL')
                tpmUpdater = RKLUpdater(initialWeights, 0.02, ...
                    this.config.C, this.config.A, this.config.B, this.config.W, this.config.V );
            elseif strcmp(this.estimatorType, 'ConvexOptimization')             
                tpmUpdater = ConvexOptimizationUpdater(40, ...
                    this.numModes , this.config.C, ...
                    this.config.A, this.config.B, this.config.W, this.config.V);
            elseif strcmp(this.estimatorType, 'Random')
                tpmUpdater = RandomUpdater();
            elseif strcmp(this.estimatorType, 'Super')
                tpmUpdater = Updater();
            elseif strcmp(this.estimatorType, 'Histogram')
                tpmUpdater = Updater();
            else
                error(['** Please choose an existing updater like ', ...
                    'ConvexOptimization, RKL, or Random **']);
            end 
            estimator = DelayedModeIMMF(this.modeFilters, localTransitionMatrix, ...
                    this.maxMeasurementDelay, tpmUpdater);
            estimator.setState(this.initialEstimate);
        end
        
        function controller = setUpController(this)
            % This method sets up the infinite horizon controller used in
            % the example code.
            %
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << controller: infinite horizon controller object
            controller = InfiniteHorizonController(this.config.A, this.config.B, ...
                this.config.Q, this.config.R, this.caDelayProb, this.controlSequenceLength);
            
            switch controller.status
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
        
        function controller = updateTransitionMatrix(this, controller)
            % This method returns an updated controller object if the true
            % transition matrix is changed.
            %
            % Parameters:
            %   >> this: SetUp object which uses this method
            %   >> controller: to be updated
            % Returns:
            %   << controller: updated infinite horizon controller object
            controller.changeTransitionMatrix(this.transitionMatrix);
        end
        
        function [plant, augmentedPlant, sensor, actuator] = setUpPlant(this)
            % This method initializes the plant as well as the sensor and
            % the actuator.
            %
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << plant: linear plant object with behavior specified by
            %   this.config
            %   << augmentedPlant: Markov jump system
            %   << sensor: measurementModel object describing the structure
            %   of the measurements
            %   << actuator: actuator object with a buffer for control
            %   sequences
            actuator = BufferingActuator(this.controlSequenceLength, Inf, zeros(...
                size(this.config.B, 2), 1));
            plant = LinearPlant(this.config.A, this.config.B, this.config.W);
            sensor = LinearMeasurementModel(this.config.C);
            % zero mean measurement noise
            sensor.setNoise(Gaussian(zeros(size(this.config.C, 1), 1), this.config.V));  
            
            augmentedPlant = JumpLinearSystemModel(this.numModes, ...
                arrayfun(@(~) LinearPlant(this.config.A, this.config.B, this.config.W), ...
                        1:this.numModes, 'UniformOutput', false));
            augmentedPlant.setNoise(Gaussian(zeros(...
                size(this.config.C, 1), 1), this.config.V))
        end
        
        function [dimX, dimY, dimU] = getDimensions(this)
            % returns dimension of the state (X),
            %  the measurement vector (Y), and the input vector (U).
            % 
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << dimX: length of the state vector
            %   << dimY: length of the measurement vector
            %   << dimU: length of the input vector
            dimX = size(this.config.A,1);
            dimU = size(this.config.B,2);
            dimY = size(this.config.C,1);
        end
        
        function [delay, sequenceLength] = getDelay(this)
            % returns the predetermined maximum for the considered delays 
            % of measurement packets as well as the the maximum delays (-1)
            % of the control packets.
            % 
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << delay: maximum of considered measurement delays
            %   << sequenceLength: number of control inputs sent per
            %   control packet
            delay = this.maxMeasurementDelay;
            sequenceLength = this.controlSequenceLength;
        end
        
        function [initialPlantState, initialEstimate] = getInitialPlantStateAndEstimate(this)
            % returns the initial plant state and the corresponding
            % estimate.
            % 
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << initialPlantState: state vector initializing the plant
            %   << initialEstimate: initial estimate for filters
            initialPlantState = this.initialPlantState;
            initialEstimate = this.initialEstimate;
        end
        
        function [caPacketDelays, scPacketDelays, acPacketDelays] = ...
                getPackageDelays(this)
            % returns initially drawn delays for control packets,
            % measurement packets, and acknowledgment of control packets.
            % 
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << caPacketDelays: drawn delays of controller-actuator link
            %   << scPacketDelays: drawn delays of sensor-controller link
            %   << acPacketDelays: drawn delays of actuator-controller link
            caPacketDelays = this.caPacketDelays;
            scPacketDelays = this.scPacketDelays;
            acPacketDelays = this.acPacketDelays;
        end
        
        function [transitionMatrix, delayDistribution] = getTransitionMatrix(this)
            % returns the ground truth transition probability matrix and the 
            % true delay distribution of the controller-actuator link.
            % 
            % Parameters:
            %   >> this: SetUp object which uses this method
            % Returns:
            %   << transitionMatrix: ground truth matrix
            %   << delayDistribution: true delay distribution
            transitionMatrix = this.transitionMatrix;
            delayDistribution = this.caDelayProb;
        end
        
        %% ComputeCosts
        function costs = ComputeCosts(this, xTrajectory, uSeq)
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
         costs = Utility.computeLQGCosts(this.simTime, xTrajectory, uSeq, ...
             this.config.Q, this.config.R) / this.simTime;     
        end
    end
end