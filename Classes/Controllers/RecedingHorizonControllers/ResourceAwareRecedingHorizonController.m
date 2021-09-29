classdef ResourceAwareRecedingHorizonController < SequenceBasedTrackingController & CaDelayProbsChangeable
    % Implementation of a resource-aware sequence-based model predictive controller for
    % linear NCS with networks connecting the controller and actuator, 
    % and sensor and controller, respectively, where application layer ACKs
    % are sent out by the actuator upon reception of applicable control
    % inputs.
    %
    % Literature: 
    %   Florian Rosenthal, Fabio Broghammer, Benjamin Noack, and Uwe D. Hanebeck,
    %   Resource-Aware Sequence-Based MPC using Simulated Annealing,
    %   submitted to IEEE Control System Letters
    %
    %   Fabio Broghammer, 
    %   Ereignisgesteuerte Folgeregelung über verlustbehaftete Netzwerke mit Anwendung auf holonome Roboter (in German), 
    %   Master thesis, Karlsruhe Institute of Technology (KIT), 
    %   2021.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
    %                             Fabio Broghammer <fabio.broghammer@student.kit.edu>
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
    
    properties (GetAccess = ?ResourceAwareRecedingHorizonControllerTest, SetAccess = immutable)
        horizonLength double {Validator.validateHorizonLength(horizonLength)} = 1;        
        
        A;
        B;        
        C;
        W;
        V;
        % for fixed sequence length, F, G, H, and J are independent of (A,B)
        F;
        G;
        H; % mode-dependent, one for each mode (in total sequence length + 1, first and last one are zero)
        J; % we only need the one corresponding to the first mode (all others are zero)
        
        % Precomputed unrolled network system matrix.
        rollF;
        % Precomputed unrolled network input matrix.
        rollG;
 
    end
       
    properties (SetAccess = immutable, GetAccess = public)
        maxMeasDelay;
    end
    
    properties (Access = private)
        initialized = true; % flag to indicate that controller state has been freshly set
        % set to false once a control sequence has been computed (and the
        % state has been updated)        
        
        % the input related part of the state, evolves according to %eta_k+1=F*eta_k + G*U_k
        etaState; % holds eta_{k-1}
         
        delayTransitionMatrix; % for the delays of the controller packets (if delays are Markovian)
               
        conditionalInputProbs; % corresponds to mode theta in case of uncorrelated delays
        
        augQ; % Q_tilde in the paper, time-invariant diagonal of augQ can be prefilled (Q and terminalQ)
        %Q;
        %terminalQ; % use unique stabilizing solution of DARE 
        R;
        alpha(1,1) double {mustBeNonnegative, mustBeFinite} = 1; % Weight factor for trace of covariance C_{k+Horizon}^p.
        augA; % time-invariant portions A and F on diagonal can be prefilled
        augB; %time-invariant portion G (lower part of augB) can be prefilled
        
        lastNumUsedMeas = 0;
        lastNumDiscardedMeas = 0;          
    end
    
    properties (SetAccess = private, GetAccess = ?ResourceAwareRecedingHorizonControllerTest)        
        K; % precomputed matrices K, defined in eq. (3.22)
        S; % precomputed part of matrices Shat_k, defined in Appendix A.2.2
        T; % precomputed part of matrices That_k, defined in Appendix A.2.2
        
        % state history to store previous states (required for filter)
        % we store them in the form of Gaussians
        % (cell array of Gaussians)
        % i-th element indicates the prior(!) at time k-i+1, where k is
        % the current time
        % stateHistory{1} is thus x_k^p, the state estimate used by the controller at time k
        % stateHistory{2} is thus x_{k-1}^p
        % stateHistory{3} is thus x_{k-2}^p
        % ...
        % stateHistory{maxMeasDelay+1} is thus x_{k-maxMeasDelay}^p
        stateHistory;
        % measurement history: store the received measurements
        % cell array where each entry is a matrix of measurements (might be empty)
        % i-th element indicates the measurements taken at time k-i for update at k-i, where k is
        % the current time
        % measurementHistory{1} is thus y_{k-1} or empty
        % measurementHistory{2} is thus y_{k-2} or empty
        % ...
        % measurementHistory{meaxMeasDelay} is thus y_{k-maxMeasDelay}
        measurementHistory;        
        
        % inputHistory to store previous possible inputs (sequence length) column-wise (required for
        % time update of the filter, which uses expected input)
        % last possible input is always default input zero, which we thus need not store
        % inputHistory(:,:, 1) are possible inputs at time k-1, needed to compute u_{k-1}, needed for x_k^p
        % inputHistory(:,:, 2) are possible inputs at time k-2, needed to compute u_{k-2}, needed for x_{k-1}^p
        % ...
        % inputHistory(:,:, meaxMeasDelay) are possible inputs at time k-maxMeasDelay, needed to compute u_{k-maxMeasDelay}, needed for x_{k-maxMeasDelay+1}^p
        inputHistory;       
        % inputProbsHistory to store previous input probabilities (sequence length +1 per time step) as row vectors (required for
        % time update of the filter, which uses expected input) 
        % inputProbsHistory(1,:) are probabilities for inputs at time k-1, needed to compute u_{k-1}, needed for x_k^p
        % inputProbsHistory(2,:) are probabilities for inputs at time k-2, needed to compute u_{k-2}, needed for x_{k-1}^p
        % ...
        % inputProbsHistory(maxMeasDelay,:) are probabilities for inputs at time k-maxMeasDelay, needed to compute u_{k-maxMeasDelay}, needed for x_{k-maxMeasDelay+1}^p
        inputProbsHistory;
        
        delayProbs; % stores the current and the previous delay probs
        % delayProbs(:, 1) = current delay probs (time k)
        % delayProbs(:, 2) = previous delay probs(time k-1)
        % ...
        % delayProbs(:, N) = delay probs at time k-N, where N+1 is sequence length
        
        sendingCostFunction;
        neighborFunction;
        startScheduleFunction;
        % event history [e_{k-N}, ..., e_{k-1}], where N is sequence length -1
        % used to calculate the sending cost.
        eventHistory;
        % last computed schedule
        lastSchedule;        
                
        saOptions;
        quadprogOptions; 
    end
    
    methods (Access = public)
        %% ResourceAwareRecedingHorizonController
        function this = ResourceAwareRecedingHorizonController(A, B, C, Q, R, alpha, caDelayProbs, ...
                sequenceLength, horizonLength, maxMeasDelay, W, V, x0, x0Cov, ...
                sendingCostFunction, neighborFunction, initScheduleFunction, startScheduleFunction, refTrajectory)
            % Class constructor.
            %
            % Parameters:
            %   >> A (Square Matrix)
            %      The system matrix of the plant.
            %
            %   >> B (Matrix)
            %      The input matrix of the plant.
            %
            %   >> C (Matrix)
            %      The measurement matrix of the sensor.
            %
            %   >> Q (Positive semi-definite matrix)
            %      The state weighting matrix Q_k in the controller's underlying cost function.
            %      For the terminal state x_K, where K is the horizon
            %      length, the corresponding weighting matrix Q_K is chosen as the
            %      stabilizing solution of the DARE associated with A, B,
            %      Q, R. If this solution does not exist, Q is used instead.
            %
            %   >> R (Positive definite matrix)
            %      The input weighting matrix in the controller's underlying cost function.
            %
            %   >> alpha (Positive scalar)
            %      The penalty term for the predicted state covariance C_K^p at the terminal stage in the
            %      controller's underlying cost function.
            %
            %   >> caDelayProb (Nonnegative vector)
            %      The vector describing the delay distribution of the
            %      CA-network.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %   >> horizonLength (Positive integer)
            %      The horizon length K considered by the controller for optimization.
            %
            %   >> maxMeasDelay (Nonnegative integer)
            %      The maximum delay a measurement may experience before it
            %      is discarded by the controller.
            %
            %   >> W (Square matrix)
            %      The covariance matrix of the plant noise.
            %
            %   >> V (Square matrix)
            %      The covariance matrix of the measurement noise.
            %            
            %
            %   >> x0 (Vector)
            %      The controller's initial estimate of the plant state.
            %
            %   >> x0Cov (Positive definite matrix)
            %      The covariance matrix denoting the uncertainty of the controller's initial estimate.
            %
            %   >> sendingCostFunction (Function handle)
            %      Function used by the controller to determine the
            %      transmission cost S_k for a particular event schedule E_k,
            %      given a history of previous transmissions e_k, e_{k-1}, ..., e_{k-N}, where N is sequence length -1.
            %      This function must accept two arguments, <eventSchedule>
            %      and <eventHistory>, and return a nonnegative scalar <cost>.
            %            
            %   >> neighborFunction (Function handle)
            %      Function used by the controller during simulated
            %      annealing to obtain a neighboring solution candidate E_k' for
            %      a particular event schedule E_k.
            %      This function must accept one argument, <eventSchedule>,
            %      and return a new one of the same length.
            %            
            %   >> initScheduleFunction (Function handle)
            %      Function used by the controller during to generate an
            %      intial guess (i.e., an initial event schedule E_0) prior
            %      to the very first invocation of simulated annealing (i.e.,
            %      at the first time step).
            %      This function must accept one argument, <length>,
            %      and return a valid transmission schedule E_0 of length <length>.
            %            
            %   >> startScheduleFunction (Function handle)
            %      Function used by the controller during to, at time k, generate an
            %      intial solution candidate E_k for simulated annealing, 
            %      given the solution E_{k-1} from the previous time step k-1.
            %      This function must accept one argument, <eventSchedule>,
            %      and return a new (sligthly different) one of the same length.
            %
            %   >> referenceTrajectory (Matrix, optional)
            %      The state reference to track, given as a matrix
            %      with the reference states column-wise arranged, if a tracking task is performed.
            %      If the state shall be driven to the origin, this
            %      argument should be left out.      
            %
            % Returns:
            %   << this (ResourceAwareRecedingHorizonController)
            %      A new ResourceAwareRecedingHorizonController instance.
            
            % do not require changing setpoints or a reference to be
            % tracked            
            if nargin == 18
                refTrajectory = [];                
            elseif nargin ~= 19 
                error('ResourceAwareRecedingHorizonController:InvalidNumberOfArguments', ...
                    '** Constructor must be called with either 18 or 18 arguments **');
            end
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            
            Validator.validateMeasurementMatrix(C, dimX);
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            
            assert(Checks.isNonNegativeScalar(maxMeasDelay) && mod(maxMeasDelay, 1) == 0, ...
                'ResourceAwareRecedingHorizonController:InvalidMaxMeasDelay', ...
                ['** Input parameter <maxMeasDelay> (maximum measurement', ...
                'delay (M)) must be a nonnegative integer **']);
            
            % check the noise
            Validator.validateSysNoiseCovarianceMatrix(W, dimX);
            Validator.validateMeasNoiseCovarianceMatrix(V, size(C, 1));
                        
            % check the states
            assert(Checks.isVec(x0, dimX) && all(isfinite(x0)), ...
                'ResourceAwareRecedingHorizonController:InvalidX0', ...
                '** Input parameter <x0> (initial estimate of plant state) must be a %d-dimensional vector **', dimX);
            assert(Checks.isCov(x0Cov, dimX), ...
                'ResourceAwareRecedingHorizonController:InvalidX0Cov', ...
                '** Input parameter <x0Cov> (initial covariance) must be positive definite and %d-by-%d-dimensional **', dimX, dimX);
             
            assert(isa(sendingCostFunction, 'function_handle') && nargin(sendingCostFunction) == 2, ...
                'ResourceAwareRecedingHorizonController:InvalidSendingCostFunction', ...
                '** <sendingCostFunction> must be a function handle accepting two arguments **');
            
            assert(isa(neighborFunction, 'function_handle') && nargin(neighborFunction) == 1, ...
                'ResourceAwareRecedingHorizonController:InvalidNeighborFunction', ...
                '** <neighborFunction> must be a function handle accepting one argument **');
            
            assert(isa(initScheduleFunction, 'function_handle') && nargin(initScheduleFunction) == 1, ...
                'ResourceAwareRecedingHorizonController:InvalidInitScheduleFunction', ...
                '** <initScheduleFunction> must be a function handle accepting one argument **');
            
            assert(isa(startScheduleFunction, 'function_handle') && nargin(startScheduleFunction) == 1, ...
                'ResourceAwareRecedingHorizonController:InvalidStartScheduleFunction', ...
                '** <startScheduleFunction> must be a function handle accepting one argument **');           
             
            if isempty(refTrajectory)
                Z = [];
            else
                Z = eye(dimX);
            end
                  
            this = this@SequenceBasedTrackingController(dimX, dimU, sequenceLength, false, Z, refTrajectory, []);
            
            this.maxMeasDelay = maxMeasDelay;
            this.A = A;
            this.B = B;
            this.C = C;
            this.W = W;
            this.V = V;
            
            Validator.validateDiscreteProbabilityDistribution(caDelayProbs);
            truncatedProbs = Utility.truncateDiscreteProbabilityDistribution(caDelayProbs, this.sequenceLength + 1);
            this.delayProbs = repmat(truncatedProbs(:), 1, this.sequenceLength); % like a history
            this.conditionalInputProbs = [zeros(this.sequenceLength, 1); 1]; % default input at the beginning, last mode
            
            %Validator.validateTransitionMatrix(caDelayProbs, this.sequenceLength + 1);            
            %this.delayTransitionMatrix = caDelayProbs;
            %this.delayProbs = zeros(this.sequenceLength + 1, this.sequenceLength);
            %this.delayProbs(:, 1) = 1 / (this.sequenceLength + 1) + zeros(this.sequenceLength + 1, 1); % uniform is assumed
            
            this.horizonLength = horizonLength;            
            [this.F, this.G, this.H, allJ] = Utility.createActuatorMatrices(this.sequenceLength, dimU);            
            this.J = allJ(:, :, 1); % we only need the first of the J matrices
            
            this.etaState = zeros(size(this.F, 1), 1);
            % fill time-invariant portions of augA and augB
            % remaining blocks are computed each time a new sequence is
            % computed (dependent on transmission history/schedule)
            this.augA = repmat(blkdiag(this.A, this.F), 1, 1, this.horizonLength); % diagonal is time-invariant
            this.augB = repmat([zeros(dimX, size(this.G, 2)); this.G], 1, 1, this.horizonLength); % G is time-invariant 
            
            %this.Q = Q;
            % for stability: use stabilizing solution of DARE as terminal
            % weighting -> only exists in case (A,B) stabilizable + the
            % associated symplectic matrix has no eigenvalues on the unit
            % circle     
            [terminalQ,~,~,info] = idare(A, B, Q, R);
            if info.Report > 1
                % either report == 2: The solution is not finite
                % or report == 3: No solution found since the symplectic
                % matrix has eigenvalues on the unit circle
                terminalQ = Q;
                warning('ResourceAwareRecedingHorizonController:InvalidPlant', ...
                    '** (A,B) seems not stabilizable, cannot use stabilizing solution of associated DARE as terminal weighting Q_N **');
            end
            % fill time-invariant portions of augQ
            % remaining blocks are computed each time a new sequence is computed (dependent on schedule)
            this.augQ = repmat(blkdiag(Q, zeros(size(this.F))), 1, 1, this.horizonLength);
            % terminalQ is different
            this.augQ = cat(3, this.augQ, blkdiag(terminalQ, zeros(size(this.F))));            
            this.R = R;
            this.alpha = alpha;            
            
            this.rollF = this.rolloutSystemMatrixF();       
            this.rollG = this.rolloutInputMatrixG();            
            this.computeAndSetKMatrices();
            this.computeAndSetSAndTMatrices();
            
            this.sendingCostFunction = sendingCostFunction;
            this.neighborFunction = neighborFunction;
            this.startScheduleFunction = startScheduleFunction;
            
            this.stateHistory = cell(1, maxMeasDelay + 1);
            this.stateHistory{1} = Gaussian(x0, x0Cov); % all other entries of the cell array are empty
            this.measurementHistory = cell(1, this.maxMeasDelay); % no measurements available
            this.inputHistory = zeros(dimU, this.sequenceLength, this.maxMeasDelay); % no inputs available, default input zero
            % we do not store the last possible input (would be column
            % sequence length +1) as it is default input
            % no input applied, zero default input, corresponding probability is 1
            this.inputProbsHistory = [zeros(this.maxMeasDelay, this.sequenceLength), ones(this.maxMeasDelay, 1)];
                        
            this.eventHistory = zeros(1, this.sequenceLength - 1);
            
            % this function is only used in case of correlated delays
            %this.conditionalInputProbs = this.computeConditionalProbabilitiesForEvents();
            this.lastSchedule = initScheduleFunction(this.horizonLength);
                        
            this.saOptions = optimoptions(@simulannealbnd, ...
                'DataType', 'custom', ...
                'AnnealingFcn', @this.permuteSchedule, ...
                'AcceptanceFcn', @this.acceptanceFunction, ...
                'InitialTemperature', this.horizonLength, ...
                'MaxStallIterations', 10, ...
                'ReannealInterval', 800, ...
                'Display', 'Off' ...
            );
            this.quadprogOptions = optimoptions('quadprog', 'Display', 'off');            
        end
    
        %% setControllerPlantState
        function setControllerPlantState(this, state)
            % Set the estimate of the plant state.
            %
            % This function is mainly used to set an initial plant state or for testing purposes.            
            %
            % Parameters:
            %   >> state (Subclass of Distribution)
            %      The new system state.
            
            assert(Checks.isClass(state, 'Distribution') && state.getDim() == this.dimPlantState, ...
                'ResourceAwareRecedingHorizonController:SetControllerPlantState:InvalidState', ...
                '** <state> must be a %d-dimensional Distribution **', this.dimPlantState);
            
            [x0, x0Cov] = state.getMeanAndCov(); % moment matching
                       
            this.stateHistory = cell(1, this.maxMeasDelay + 1);
            this.stateHistory{1} = Gaussian(x0, x0Cov); % all other entries of the cell array are empty             
            this.initialized = true;
        end
        
        %% getControllerPlantState
        function plantState = getControllerPlantState(this)
            % Get the plant state as currently perceived by the controller, i.e., its internal estimate of the plant state.
            %
            % Returns:
            %   << plantState (Column vector)
            %      The controller's current estimate of the plant state.
            %
            if ~isempty(this.stateHistory{1})            
                % the state has been freshly set
                [plantState, ~] = this.stateHistory{1}.getMeanAndCov();
            else
                [plantState, ~] = this.stateHistory{2}.getMeanAndCov();
            end
            % should be called after a control sequence has been computed
        end
        
        %% getLastComputationMeasurementData
        function [numUsedMeas, numDiscardedMeas] = getLastComputationMeasurementData(this)
            % Get information on the last performed control sequence computation (due to a call of computeControlSequence).
            %
            % Returns:
            %   << numUsedMeas (Nonnegative integer)
            %      The number of measurements used to update the controller state.
            %
            %   << numDiscardedMeas (Nonnegative integer)
            %      The number of measurements discarded due to their delays.
            
            numUsedMeas = this.lastNumUsedMeas;
            numDiscardedMeas = this.lastNumDiscardedMeas;
        end
        
        %% changeCaDelayProbs
        function changeCaDelayProbs(this, newCaDelayProbs)
            % Change the distribution of the control packet delays to be
            % assumed by the controller from the current time step on until
            % this method is called again.
            %
            % Parameters:
            %  >> newCaDelayProbs (Nonnegative vector)
            %     Vector specifiying the new delay distribution.
            %
            Validator.validateDiscreteProbabilityDistribution(newCaDelayProbs);
            truncatedProbs = Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.sequenceLength + 1);
            this.delayProbs(:, 1) = truncatedProbs(:);
        end
        
     %% reset
        function reset(this)
            
        end
        
        %% computeControlSequence
        function inputSeq = computeControlSequence(this, measurements, measurementDelays, modeMeas, modeDelays, timestep)
            % time step is only required for picking the correct reference
            % trajectory (if present)
            assert(Checks.isPosScalar(timestep) && mod(timestep, 1) == 0, ...
                'ResourceAwareRecedingHorizonController:computeControlSequence:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be a positive integer **');
            
            if this.initialized
                % no update of the controller state required
                % only happens at the initial time step
                % or if state has been set using setControllerPlantState
                this.lastNumDiscardedMeas = 0;
                this.lastNumUsedMeas = 0;
                this.initialized = false;
                currentMeasurements = [];                
            else
                % update the controller state first
                [currentMeasurements, ~] = ...
                    updateControllerPlantState(this, measurements, measurementDelays, modeMeas, modeDelays);
            end
            % this.stateHistory{1} is x_k^p (the prior after prediction), used to compute the input sequence
            inputSeq = this.doControlSequenceComputation(this.getControllerPlantState(), timestep);
            
            oldEta = this.etaState;
            if isempty(inputSeq) % we do not send
                transmissionDecision = 0;
                this.etaState = this.F * this.etaState;  %eta_k+1
            else
                transmissionDecision = 1;
                this.etaState = this.F * this.etaState + this.G * inputSeq;
            end
            
            % prepate the histories for next time step
            this.measurementHistory = circshift(this.measurementHistory, 1, 2);
            this.measurementHistory{1} = currentMeasurements; % prepared for next time step k+1, this is y_k             
            
            this.stateHistory = circshift(this.stateHistory, 1, 2); % prepared for k+1
            this.stateHistory{1} = []; % not yet available, will hold x_{k+1}^p
            
            % prepare input history for next time step
            this.inputHistory = circshift(this.inputHistory, 1, 3); % shift each page to the right
            this.inputProbsHistory = circshift(this.inputProbsHistory, 1, 1); % shift the rows downwards            
            this.inputProbsHistory(1, :) = this.calculateInputProbabilitiesForSchedule(1, transmissionDecision, timestep)'; % input probs at time k
                        
            % H(:, :, end) is always zero and H(:, :, 1) is always zero
            possibleInputs = reshape(mtimesx(this.H(:, :, 2:end-1), oldEta), this.dimPlantInput, this.sequenceLength-1);
            if transmissionDecision == 0
                % do not append zero default input to possible inputs
                possibleInputs = [zeros(this.dimPlantInput, 1), possibleInputs];
            else
                possibleInputs = [inputSeq(1:this.dimPlantInput), possibleInputs];
            end
            this.inputHistory(:, :, 1) = possibleInputs; % without the default input, possible inputs at time k
                
            % predict the mode probs
            this.conditionalInputProbs = this.predictModeProbs(this.conditionalInputProbs, ...
                [this.eventHistory, transmissionDecision], timestep);
            
            % update event history and the ca delay probs
            this.eventHistory = [this.eventHistory(2:end), transmissionDecision];            
            this.delayProbs = [this.delayProbs(:, 1), this.delayProbs(:, 1:end-1)];            
            % predict the network delays            
            %this.delayProbs = this.predictNetworkDelays();       
        end        
    end
        
    methods(Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, timestep)            
            % augment the reference first
            % reference for eta portion is zero
            if ~isempty(this.refTrajectory)
                % the preview of the reference must cover the whole optimization horizon
                if timestep + this.horizonLength > size(this.refTrajectory, 2)
                    % too short, so append the last value
                    xRef = [this.refTrajectory(:, timestep:end), repmat(this.refTrajectory(:, end), 1, ...
                        timestep + this.horizonLength - size(this.refTrajectory, 2))];
                else
                    xRef = this.refTrajectory(:, timestep:timestep + this.horizonLength);
                end                
                augRef = [xRef; zeros(size(this.etaState, 1), this.horizonLength + 1)];
            else
                % no reference present, track zero
                augRef = zeros(this.dimPlantState + size(this.etaState, 1), this.horizonLength + 1);
            end           
            % prepare simulated annealing
            startSchedule = this.startScheduleFunction(this.lastSchedule);            
            
            % solve using SA
            this.lastSchedule = simulannealbnd(@(x) this.calculcateScheduleCost(x, state, augRef(:)), ...
                startSchedule, [], [], this.saOptions);
           
            transmissionDecision = this.lastSchedule(1);
            if transmissionDecision == 1
                [inputSequence, ~] = this.calculateControllerOutputForSchedule(state, this.lastSchedule, augRef(:));
            else
                % do not send
                inputSequence = [];
            end
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, timestep)
            if ~isempty(this.Z)
                % compute performance output and difference to reference
                % trajectory
                % directly track state reference, Z=I
                performance = state - this.refTrajectory(:, timestep);
            else
                performance = state;
            end            
            stageCosts = Utility.computeStageCosts(performance, input, this.augQ(1:this.dimPlantState, 1:this.dimPlantState, 1), this.R);
        end
        
        %% doCostsComputation
        function costs = doCostsComputation(this, stateTrajectory, appliedInputs)
            numInputs = size(appliedInputs, 2);
            assert(size(stateTrajectory, 2) == numInputs + 1, ...
                'ResourceAwareRecedingHorizonController:DoCostsComputation:InvalidStateTrajectory', ...
                '** <stateTrajectory> is expected to have %d columns ', numInputs + 1);
            
            if ~isempty(this.Z)
                % compute performance output and difference to reference
                % trajectory
                % directly track state reference, Z=I
                performance = stateTrajectory - this.refTrajectory(:, 1:numInputs + 1);
            else
                performance = stateTrajectory;
            end 
            Q = this.augQ(1:this.dimPlantState, 1:this.dimPlantState, 1);
            costs = Utility.computeLQGCosts(numInputs, performance, appliedInputs, Q, this.R);
        end
        
        %% doGetDeviationFromRefForState
        function deviation = doGetDeviationFromRefForState(this, state, timestep)
            if ~isempty(this.Z)
                % we track a reference trajectory
                assert(Checks.isScalarIn(timestep, 1, size(this.refTrajectory, 2)) && mod(timestep, 1) == 0, ...
                    'ResourceAwareRecedingHorizonController:GetDeviationFromRefForState:InvalidTimestep', ...
                    '** Input parameter <timestep> must be in {1, ... %d} **', ...
                    size(this.refTrajectory, 2));
                % directly track state reference, Z=I
                deviation = state - this.refTrajectory(:, timestep);
            else
                % we track the origin
                deviation = state;
            end
        end
    end
    
    methods (Access = private)        
        
        %% computeAndSetSAndTMatrices
        function computeAndSetSAndTMatrices(this)
            this.S = zeros(size(this.G, 2), size(this.G, 2), this.horizonLength);
            this.T = zeros(size(this.G, 1), size(this.G, 1), this.sequenceLength-2, this.horizonLength);
            for k=1:this.horizonLength
                this.S(:, :, k) = this.J' * this.K(:, :, k) * this.J;
                this.S(:, :, k) = (this.S(:, :, k) + this.S(:, :, k)') / 2; % ensure symmetry
                for j=2:this.sequenceLength
                    this.T(:, :, j, k) = this.H(:, :, j)' * this.K(:, :, k) * this.H(:, :, j);
                    this.T(:, :, j, k) = (this.T(:, :, j, k) + this.T(:, :, j, k)') / 2; % ensure symmetry
                end
            end            
        end
        
        %% computeAndSetKMatrices
        function computeAndSetKMatrices(this)
            % K_j, as defined in eq. (3.22)
            % K_j = B'*(A^(horizon-1-j))'*(A^(horizon-1-j))*B for j=0, 1,...,horizon-1
            powersA = zeros(this.dimPlantState, this.dimPlantState, this.horizonLength);
            powersA(:, :, 1) = eye(this.dimPlantState);
            powersA(:, :, 2) = this.A;
            for j=3:this.horizonLength
                powersA(:, :, j) = powersA(:, :, j-1) * this.A;
            end
            
            this.K = zeros(this.dimPlantInput, this.dimPlantInput, this.horizonLength);

            for i = 1:this.horizonLength
                temp = powersA(:, :, this.horizonLength - i + 1) * this.B;
                this.K(:, :, i) = temp' * temp;
            end
        end        
               
        %% rolloutInputMatrixG
        function rollG = rolloutInputMatrixG(this)
            % Creates network input matrix of the unrolled system.
            %
            % Literature:
            %   Fabio Broghammer, 
            %   Ereignisgesteuerte Folgeregelung über verlustbehaftete Netzwerke mit Anwendung auf holonome Roboter (in German), 
            %   Master thesis, Karlsruhe Institute of Technology (KIT), 
            %   2021.
            %
            % Returns:
            %   << rollG (Matrix, dimEta*horizon-by-dimInput*controlSequenceLength*horizon)
            %      The input matrix of the unrolled system.
            dimEta = size(this.F, 1);
            dimControlSequence = size(this.G, 2);

            rollG = zeros(dimEta * this.horizonLength, dimControlSequence * this.horizonLength);
            for column = 1:this.horizonLength -1 
                startColumnIdx = (column - 1) * dimControlSequence + 1;
                endColumnIdx = column * dimControlSequence;
                entry = this.G;

                for row = column + 1:this.horizonLength
                    if row ~= column + 1
                        entry = this.F * entry;
                    end                    
                    rollG((row - 1) * dimEta + 1:row * dimEta, startColumnIdx:endColumnIdx) = entry;
                    % short cut for large horizons: we can break the inner
                    % loop once entry is zero
                    if isequal(zeros(size(entry)), entry)
                        break
                    end
                end

            end          
        end
        
        %% rolloutSystemMatrixF
        function rollF = rolloutSystemMatrixF(this)
            % F_hat, as per eq. (3.23)
            % exploit that F as defined by (3.11) is a nilpotent matrix of
            % dimension <dimEta>x<dimEta>, thus F^b = 0 for some b <= dimEta
            % so we need only calculate the powers up to F^dimEta in case
            % dimEta < horizonLength 
            dimEta = size(this.F, 1);
            rollF = zeros(dimEta * this.horizonLength, dimEta);
            
            entry = eye(dimEta); 
            for k = 1:min(this.horizonLength, dimEta)             
                rollF((k - 1) * dimEta + 1:k * dimEta, :) = entry;
                entry = this.F * entry;
                if isequal(entry, zeros(dimEta))
                    % we can return, will remain zero
                    return
                end
            end
        end
                
        % for simulated annealing
        
        %% calculcateScheduleCost
        function cost = calculcateScheduleCost(this, schedule, plantState, augRef)
            [~, cost] = this.calculateControllerOutputForSchedule(plantState, schedule, augRef);
        end
        
        %% acceptanceFunction
        function acceptPoint = acceptanceFunction(this, optimValues, newx, newfval)
            % Defines a custom acceptanceFunction, as defined in line 15 of Algorithm 1 in the paper
            
            currT = optimValues.temperature;
            fval = optimValues.fval;
            temp = min(exp(-(newfval - fval) / currT), 1);
            acceptPoint = temp > rand();
        end

        %% permuteSchedule
        function schedule = permuteSchedule(this, optimValues, problemData)
            % Calls the neighbor function depending on the current temperature.
            % See literature for further details.
            schedule = optimValues.x;

            for i = 1:floor(optimValues.temperature) + 1
                schedule = this.neighborFunction(schedule);
            end
        end
        
        %% calculateControllerOutputForSchedule
        function [inputSeq, cost] = calculateControllerOutputForSchedule(this, plantState, eventSchedule, augRef)
            dimEta = numel(this.etaState);
            dimAugState  = dimEta + this.dimPlantState;
            dimInputSeq = size(this.G, 2);
            
            inputProbs = this.calculateInputProbabilitiesForSchedule(this.horizonLength, eventSchedule, this.sequenceLength);
            
            % augmented cost matrices (Q_tilde, R_tilde and O_tilde in the paper)            
            augR = zeros(dimInputSeq, dimInputSeq, this.horizonLength);
                        
            % covariance part C_{k+horizon}^p            
            rollAugS = zeros(dimInputSeq * this.horizonLength); % as defined in appendix A.2.2
            rollAugT = zeros(dimEta * this.horizonLength); % as defined in appendix A.2.2            
            rollAugJ = zeros(dimInputSeq * this.horizonLength); % as defined in appendix A.2.2
            rollAugH = zeros(dimEta * this.horizonLength); % as defined in appendix A.2.2
            rollAugX = zeros(dimEta * this.horizonLength, dimInputSeq * this.horizonLength); % as defined in appendix A.2.2
            
            % rollout matrix O (cross term O_hat in cost function, eq. (3.19) in thesis)
            % consists of matrices O_tilde penalizing the cross-terms
            rollAugO = zeros(dimAugState * (this.horizonLength + 1), dimInputSeq * this.horizonLength);
            for k = 1:this.horizonLength
                % H_k and J_k according to eq. (9) in the paper 
                expectedH = sum(reshape(inputProbs(:, k), 1, 1, []) .* this.H, 3);
                expectedJ = inputProbs(1, k) * this.J;
                
                % construct the nominal augmented dynamics (A_tilde and B_tilde in the paper, according to eq. (10)
                % only the upper-right corner of augA is time-varying
                this.augA(1:this.dimPlantState, this.dimPlantState + 1:end, k) = this.B * expectedH;
                
                this.augB(1:this.dimPlantState, :, k)  = this.B * expectedJ; %time-varying part
                          
                % lower diagonal block of augQ (Q_tilde in the thesis, eq. (3.14))
                HRH = expectedH' * this.R * expectedH;
                % ensure symmetry               
                this.augQ(this.dimPlantState + 1:end, this.dimPlantState + 1:end, k) ...
                    = (HRH + HRH') / 2;
                
                % augR (R_tilde in the thesis, eq. (3.14))
                augR(:, :, k) = expectedJ' * this.R * expectedJ;
                % ensure symmetry
                augR(:, :, k) = (augR(:, :, k) + augR(:, :, k)') / 2;
               
                startIdx = (k - 1) * dimInputSeq + 1;
                endIdx = k * dimInputSeq;
                
                % O_tilde in the thesis, eq. (3.14)
                %augO = [zeros(this.dimPlantState, dimInputSeq); expectedH' * this.R * expectedJ];
                rollAugO((k - 1) * dimAugState + 1:k * dimAugState, startIdx:endIdx) ...
                    = [zeros(this.dimPlantState, dimInputSeq); expectedH' * this.R * expectedJ];
                                                
                % augS is block-diagonal, as defined in appendix A.2.2                
                rollAugS(startIdx:endIdx, startIdx:endIdx) = inputProbs(1, k) * this.S(:, :, k);
                
                startIdxT = (k - 1) * dimEta + 1;
                endIdxT = k * dimEta;                
                % augT is block-diagonal, as defined in appendix A.2.2
                for j=2:this.sequenceLength                    
                    rollAugT(startIdxT:endIdxT, startIdxT:endIdxT) = rollAugT(startIdxT:endIdxT, startIdxT:endIdxT) ...
                        + inputProbs(j, k) * this.T(:, :, j, k); 
                end   
                
                % augJ is block-diagonal, as defined in appendix A.2.2
                JKJ = expectedJ' * this.K(:, :, k) * expectedJ;
                rollAugJ(startIdx:endIdx, startIdx:endIdx) = (JKJ + JKJ') / 2;
                
                % augH is block-diagonal, as defined in appendix A.2.2
                HKH = expectedH' * this.K(:, :, k) * expectedH;
                rollAugH(startIdxT:endIdxT, startIdxT:endIdxT) =  (HKH + HKH') / 2;
                
                % augX is block-diagonal, as defined in appendix A.2.2
                rollAugX(startIdxT:endIdxT, startIdx:endIdx) = expectedH' * this.K(:, :, k) * expectedJ;
            end
            
            % rollout over horizon            
            entry = eye(dimAugState);
            % these are the huge matrices A_tilde and B_tilde defined in Appendix A.1 in thesis
            rollAugA = zeros(dimAugState * (this.horizonLength + 1), dimAugState);
            rollAugB = zeros(dimAugState * (this.horizonLength + 1), dimInputSeq * this.horizonLength);
            for k = 1:this.horizonLength + 1                
                if k > 1
                    entry = this.augA(:, :, k - 1) * entry;
                end

                rollAugA((k - 1) * dimAugState + 1:k * dimAugState, :) = entry;
            end
                        
            for k=1:this.horizonLength 
                startColumnIdx = (k - 1) * dimInputSeq + 1;
                endColumnIdx = k * dimInputSeq;
                entry = this.augB(:, :, k);

                for row = k + 1:this.horizonLength + 1
                    if row ~= k + 1
                        entry = this.augA(:, :, row - 1) * entry;
                    end
                    rollAugB((row - 1) * dimAugState + 1:row * dimAugState, startColumnIdx:endColumnIdx) = entry;
                end
            end
            
            % rollout remaining cost matrices Q_hat and R_hat (eq. (3.19) in thesis)
            rollAugQ = Utility.createBlockdiagonalMatrix(this.augQ);
            rollAugR = Utility.createBlockdiagonalMatrix(augR);
                        
            Z_dyn = rollAugB' * rollAugQ * rollAugB + rollAugR + 2 * rollAugB' * rollAugO; % Z_k^Dyn in eq. (3.19)
            % Z_dyn is not symmetric due to last addend; but we can use its symmetric part w.l.o.g. (Lemma 3.1 in thesis) 
            Z_dyn = 0.5 * (Z_dyn + Z_dyn');
            a = rollAugA * [plantState; this.etaState];
            difference = a - augRef;

            f_dyn = difference' * rollAugQ * rollAugB + a' * rollAugO; % transpose of f_k^Dyn in eq. (3.19)
            c_dyn = a' * rollAugQ * a + augRef' * rollAugQ * augRef - 2 * augRef'*rollAugQ * a; % c_k^Dyn in eq. (3.19)
            
            % covariance part C_{k+horizon}^p of cost            
            rollAugW = rollAugT - rollAugH; % eq. (3.24)
            % Z_k^Cov in eq. (3.25)
            Z_cov = (rollAugS - rollAugJ) + this.rollG' * rollAugW * this.rollG - 2* this.rollG' * rollAugX;
            % Z_cov is not symmetric due to last addend; but we can use its symmetric part w.l.o.g. (Lemma 3.1 in thesis)
            Z_cov = (Z_cov + Z_cov') / 2;
            % transpose of f_k^Cov in eq. (3.25)
            f_cov = this.etaState' * (this.rollF' * (rollAugW * this.rollG - rollAugX)); 
            % c_k^Cov in eq. (3.25)
            c_cov = this.etaState' * this.rollF' * rollAugW * this.rollF * this.etaState;
            
            Z_k = Z_dyn + this.alpha * Z_cov; % eq. (3.26)
            f_k = f_dyn + this.alpha * f_cov; % eq. (3.26), row vector
            c_k = c_dyn + this.alpha * c_cov; % eq. (3.26)
            
            % solve using quadprog (in fact analytical solution is
            % possible, as there are no constraints)
            % resulting cost is half the cost as quadprog uses 0.5*x'*H*x +f'*x
            % but we want to solve x'*H*x+2*f'*x with H=Z_k and f'=f_k
            %[inputSequences, cost] = quadprog(Z_k, f_k, [], [], [], [], [], [], [], this.quadprogOptions);
            [inputSequences, cost] = quadprog(Z_k, f_k, [], [], [], [], [], [], []); % for mosek
            
            % analytical solution
%             sequences = -pinv(Z_k) * f_k';
%             cost2 =  sequences' * Z_k * sequences + 2*f_k*sequences;% + c_k
            
            cost = 2 * cost + c_k + this.sendingCostFunction(eventSchedule, this.eventHistory); % add c_k and transmission cost
            inputSeq = inputSequences(1:dimInputSeq);
        end       

        %% calculateInputProbabilitiesForSchedule
        function inputProbabilities = calculateInputProbabilitiesForSchedule(this, horizon, eventSchedule, timestep)
            % Calculates input probabilities for the given horizon K.
            %            
            %
            % Parameters:
            %   >> horizon (Positive Integer)
            %      The considered horizon K.
            %
            %   >> eventSchedule (Binary vector, horizon)
            %      The event schedule e_{k},..., e_{k + K - 1}.
            %
            %   >> timestep (Positive Integer)
            %      The current time step.
            %
            % Returns:
            %   << inputProbabilities (Matrix, controlSequenceLength + 1-by-K)
            %      The input probabilities for control inputs for k, ..., k
            %      + K, where K is the horizon length.
            %
            numModes = this.sequenceLength + 1;
            inputProbabilities = zeros(numModes, horizon);
            thetaState = this.conditionalInputProbs; % probabilities theta_{k-1}
            time = timestep;
            evtHistory = this.eventHistory;
            
            for k = 1:horizon
                events = [evtHistory, eventSchedule(k)];
                evtHistory = events(2:end);                
                time = time + 1;
                % propagate to current time                
                thetaState = this.predictModeProbs(thetaState, events, time);
                
                inputProbabilities(1:this.sequenceLength, k) = thetaState(1:this.sequenceLength);
                inputProbabilities(end, k) = 1 - sum(inputProbabilities(1:this.sequenceLength, k));
            end
        end
        
        %% predictModeProbs
        function newModeProbs = predictModeProbs(this, modeProbs, evtHistory, timestep)
            numModes = this.sequenceLength + 1;
            newTheta = this.calculateEventDependentModeTransitionMatrix(evtHistory)' * modeProbs;
            newModeProbs = zeros(size(modeProbs));
            % Probability is zero if simulation time is smaller than the age of the input.
            newModeProbs(1:min(timestep + 1, numModes)) = newTheta(1:min(timestep + 1, numModes));
            newModeProbs(numModes) = 1 - sum(newModeProbs(1:this.sequenceLength));
        end
        
        %% calculateEventDependentModeTransitionMatrix
        function eventTriggeredTransitionMatrix = calculateEventDependentModeTransitionMatrix(this, events)
            % Calculate the event triggered transition matrix when not all control sequences are send.
            % Transition from k-1 to k.
            % See chapter 4.
            %
            % Parameters:
            %
            %   >> events (Vector of size <sequenceLength>)
            %      [e_{k-N|k-N}, ..., e_{k|k}], where e_{i|i} is 0 or 1.
            %      If e{i|i} is 0, the control sequence at the time step i wasn't send.
            %
            % Returns:
            %   << eventTriggeredTransitionMatrix (Matrix of size <sequenceLength>+1 x <sequenceLength>+1)
            %      The event triggered transition matrix

            %delayMatrix = -1 * ones(this.sequenceLength + 1, this.sequenceLength);
            delayMatrix = zeros(this.sequenceLength + 1, this.sequenceLength);
            
            events = flip(events); % so first element is e_k, last element is e_{k-N}, where N is sequence length - 1
            % corresponding delay probs are stored in delayProbs(:, 1) to delayProbs(:, N+1)
            for i = 1:this.sequenceLength
                if events(i) == 1
                    delayMatrix(:, i) = this.delayProbs(:, i);
                else
                    % not sent, like a "loss"
                    delayMatrix(end, i) = 1;                    
                end

            end            
            eventTriggeredTransitionMatrix = Utility.calculateDelayTransitionMatrix(delayMatrix);
        end
        
%         %% predictNetworkDelays
%         function newDelays = predictNetworkDelays(this)
%             newDelays = [this.delayTransitionMatrix' * this.delayProbs(:, 1), this.delayProbs(:, 1:end-1)];
%         end
        
%         %% calculateInputProbabilitiesForSchedule
%         function inputProbabilities = calculateInputProbabilitiesForSchedule(this, horizon, eventSchedule, timestep)
%             % Calculate input probabilities for the given horizon.
%             %
%             % N + 1 must be bigger than controlSequenceLength.
%             % If you want to use a larger control sequence you have to extend the network.
%             % You can add states that are not reachable (see Evaluation/chapter4/eventSequenceLengthEvaluation).
%             %
%             % Parameters:
%             %   >> horizon (Positive Integer)
%             %      The horizon to consider.
%             %
%             %   >> eventSchedule (Binary vector, horizon)
%             %      The event schedule e_{k},..., e_{k + horizon - 1}.
%             %
%             %   >> timestep (Positive Integer)
%             %      The current time step.
%             %
%             % Returns:
%             %   << inputProbabilities (Matrix, N+2-by-horizon)
%             %      The input probabilities for k, ..., k + horizon - 1.
%             
%             inputProbabilities = zeros(this.sequenceLength + 1, horizon);           
%             time = timestep;
%             evtHistory = this.eventHistory;
%             for k = 1:horizon
%                 state = this.predictNetworkDelays(); % Update state to time step k;
%                 time = time + 1;
%                 events = [evtHistory, eventSchedule(k)];
%                 evtHistory = events(2:end);
% 
%                 eventsDecimal = bi2de(events);
%                 precomputedCondProbs = this.conditionalInputProbs(:, :, eventsDecimal + 1);
% 
%                 for i = 0:this.sequenceLength - 1
%                     if time < i
%                         % Probability is zero if simulation time is smaller than the age of the input.
%                         inputProbabilities(i + 1, k) = 0;
%                     else
%                         inputProbabilities(i + 1, k) = events(end - i) * precomputedCondProbs(:, i + 1)' * state(:, i + 1);
%                     end
% 
%                 end
%                 inputProbabilities(this.sequenceLength + 1, k) = 1 - sum(inputProbabilities(1:this.sequenceLength, k));
%             end
%         end
        
% %         %% computeConditionalProbabilitiesForEvents
% %         function condProbs = computeConditionalProbabilitiesForEvents(this)
% %             % Compute all conditional probabilities.
% %             %
% %             % Returns:
% %             %   << condProbs (Vector, N+2-by-N+1-by-2^(N+1))            
% %             %      The conditional probabilities for all input ages and
% %             %      delays and all possible events, where N+1 is the
% %             %      sequence length.
% %             %      [1, P(\tau_k >= e_k*1|\tau_{k-1} = 0), ..., P(\tau_k >= e_k*1, \tau_{k-1} >= e_{k-1}*2,... \tau_{k - age + 1} >= e_{k-N+1}*N | \tau_{k-N} = 0);
% %             %       1, P(\tau_k >= e_k*1|\tau_{k-1} = 1), ..., P(\tau_k >= e_k*1, \tau_{k-1} >= e_{k-1}*2,... \tau_{k - age + 1} >= e_{k-N+1}*N | \tau_{k-N} = 1);;
% %             %                                                       :
% %             %                                                       :
% %             %       1, P(\tau_k >= e_k*1|\tau_{k-1} = N+1), ..., P(\tau_k >= e_k*1, \tau_{k-1} >= e_{k-1}*2,... \tau_{k - age + 1} >= e_{k-N+1}*N | \tau_{k-N} = N+1);]
% %                         
% %             condProbs = zeros(this.sequenceLength + 1, this.sequenceLength, 2^this.sequenceLength);
% %             for eventDecimal = 0:(2^this.sequenceLength)-1
% %                 eventBinary = de2bi(eventDecimal, this.sequenceLength);
% %                 
% %                 for age=0:this.sequenceLength -1
% %                     % compute the possible bounds (depending on the transmissions) for the event
% %                     bounds = [1:age age];                    
% %                     for i=1:age
% %                         if eventBinary(end - i + 1) == 0
% %                             bounds(i) = 0;
% %                         end
% %                     end                   
% %                     condProbs(:, age + 1, eventDecimal + 1) = ...
% %                         this.computeConditionalProbabilitiesForInput(...
% %                             this.calculatePossibleCompounds(bounds)...
% %                             );
% %                 end
% %             end
% %         end
        
% %         %% calculatePossibleCompounds
% %         function possibleCompounds = calculatePossibleCompounds(this, bounds)
% %             % Calculate all possible compounds for
% %             % P(\tau_k >= 1, \tau_{k-1} >= 2, ... \tau_{k - age + 1} >= age | \tau_{k-age} <= age).
% %             %
% %             % Parameters:
% %             %   >> bounds (Vector)
% %             %      [e_k*1 e_{k-1}*2 e_{k-2}*3 ... {e_k-age}*age age]
% %             % Returns:
% %             %   << possibleCompounds (Matrix)
% %             %       Example: e_{k-1} = 0
% %             %       [1 0 ... age   0;
% %             %        1 0 ... age   1;
% %             %               :
% %             %               :
% %             %        1 0 ... age   age
% %             %        1 0 ... age+1 0
% %             %               :
% %             %               :
% %             %       N+1 N+1 ... N+1 age]
% %             %       (combinations of delays for given the bounds).
% %             if length(bounds) == 1
% %                 possibleCompounds = (0:bounds(1))';
% %                 return;
% %             end
% % 
% %             compounds = this.calculatePossibleCompounds(bounds(2:end));
% %             currentBound = bounds(1);
% %             possibleCompounds = [];
% % 
% %             for i = currentBound:this.sequenceLength
% %                 possibleCompounds = [possibleCompounds; ...
% %                                     [i * ones(size(compounds, 1), 1), compounds]]; 
% %             end
% %         end
        
% %         %% computeConditionalProbabilitiesForInput
% %         function condProbs = computeConditionalProbabilitiesForInput(this, possibleCompounds)
% %             % Precomputes the conditional probability for given the compounds.
% %             %
% %             % Parameters:
% %             %   >> possibleCompounds (Matrix, possibleCombination-by-'age')
% %             %      All possible compounds.
% %             %      Example: e_{k-1} = 0
% %             %       [1 0 ... age   0;
% %             %        1 0 ... age   1;
% %             %               :
% %             %               :
% %             %        1 0 ... age   age
% %             %        1 0 ... age+1 0
% %             %               :
% %             %               :
% %             %       N+1 N+1 ... N+1 age]
% %             %       (combinations of delays for given the bounds).
% %             %
% %             % Returns:
% %             %   << condProbs (Vector, N + 2-dimensional)
% %             %      Conditional probability
% %             %      [P(\tau_k >= e_{k}*1, \tau_{k-1} >= e_{k-1}*2, ... \tau_{k - age + 1} >= e_{k-age+1}*age | \tau_{k - age} = 0);
% %             %                                       ...
% %             %       P(\tau_k >= e_{k}*1, \tau_{k-1} >= e_{k-1}*2, ... \tau_{k - age + 1} >= age | e_{k-age+1}*\tau_{k - age} = age);
% %             %                                        0
% %             %                                       ...
% %             %                                        0
% %             numOfCombinations = size(possibleCompounds, 1);
% %             condProbs = zeros(this.sequenceLength + 1, 1);
% % 
% %             for i = 1:numOfCombinations                
% %                 condProbs(possibleCompounds(i, end) + 1, 1) = condProbs(possibleCompounds(i, end) + 1, 1) + this.calculateConditionalProbabilityForJoint(possibleCompounds(i, :));
% %             end
% %         end
% % 
% %         
% %         %% calculateConditionalProbabilityForJoint
% %         function probability = calculateConditionalProbabilityForJoint(this, joint)
% %             % Calculate the conditional probability
% %             % P(\tau_k = a, \tau_{k-1} = b, ..., \tau_{k-age + 1} = c | \tau_{age} = d)
% %             % = P(\tau_k = a | \tau_{k-1} = b) * ... * P(\tau_{k-age + 1} = c | \tau_{age} = d)
% %             %
% %             % Parameters:
% %             %   >> compound (Row-vector)
% %             %      The compound [a, b, ..., c, d] for \tau_k, ..., \tau_{k-age}.
% %             % Returns:
% %             %   << probability
% %             %      The conditional probability.
% %             
% %             probability = 1;
% %             for i=1:length(joint) -1
% %                 % P(\tau_{k - i} = a, \tau_{k - i - 1} = b)             %b               %  a
% %                 probability = probability * this.delayTransitionMatrix(joint(i + 1) + 1, joint(i) + 1);
% %             end
% %         end
        
        %% updateControllerPlantState
        function [currentMeasurements, currentMode] = updateControllerPlantState(this, measurements, measurementDelays, modes, modeDelays)
            % update the controller state first
            if isempty(measurements)
                currentMeasurements = [];
                numIterations = 1;
                % no measurement information available
                this.lastNumDiscardedMeas = 0;
                this.lastNumUsedMeas = 0;
            else
                this.checkMeasurementsAndDelays(measurements, measurementDelays);
                [applicableMeas, applicableDelays] = this.getApplicableMeasurements(measurements, measurementDelays);
                numIterations = max(applicableDelays);
                % store measurements with zero delay, if available
                % will not be processed, only stord for use at next time step
                currentMeasurements = this.incorporateDelayedMeasurements(applicableMeas, applicableDelays);

                this.lastNumUsedMeas = numel(applicableDelays);
                this.lastNumDiscardedMeas = numel(measurementDelays) - numel(applicableDelays);
            end    
            currentMode = [];
            if ~isempty(modes)
                this.checkModeMeasurementsAndDelays(modes, modeDelays);
                [applicableModes, applicableModeDelays] = this.getApplicableModeMeasurements(modes, modeDelays);
                if ~isempty(applicableModeDelays) % empty if delays were too large
                    numIterations = max(numIterations, max(applicableModeDelays));
                    currentMode = this.incorporateDelayedModeMeasurements(applicableModes, applicableModeDelays);
                end
            end               
            
            % if numIterations = i, 1<i<=maxMeasDelay, then y_{k-i} is the oldest received
            % measurement, it is stored in measurementHistory{i}
            % if maxDelay = 1, then either no new measurement information from the past is available
            % so that only a single prediction step is carried out, or
            % y_{k-1} is the oldest received measurement 
            delays = numIterations:-1:1;
            for i=delays
                % get the prior state at time k-i and the measurement y_k-i
                [stateMean, stateCov, priorCovSqrt] = this.stateHistory{i+1}.getMeanAndCov(); % x_{k-i}^p
                meas = this.measurementHistory{i};
                if ~isempty(meas)
                    % Kalman filter measurement update                    
                    Gmat = this.C * priorCovSqrt;                    
                    [stateMean, stateCov] = Utils.kalmanUpdate(stateMean, stateCov, meas(:), this.C * stateMean, ...
                        this.V + Gmat * Gmat', stateCov * this.C');
                end
                % predict to k-i+1 using the expected input
                % so compute it first
                [expectedInput, inputCov] = Utils.getMeanAndCov([this.inputHistory(:, :, i), zeros(this.dimPlantInput, 1)], ...
                    this.inputProbsHistory(i, :));
                predMean = this.A * stateMean + this.B * expectedInput; % noise is zero mean
                predCov = this.A * stateCov * this.A' + this.B * inputCov * this.B' + this.W;
                predCov = (predCov + predCov') / 2; % ensure symmetry
                if i == 1
                    % predicted to the current time step k, slot in
                    % cell array is empty      
                    this.stateHistory{i} = Gaussian(predMean, predCov);
                else
                    % update previous state x_{k-i+1} already present
                    this.stateHistory{i}.set(predMean, predCov);
                end
            end
        end
                
        %% checkModeMeasurementsAndDelays
        function checkModeMeasurementsAndDelays(this, trueModes, modeDelays)
            if ~isempty(trueModes)
                numCaModes = this.sequenceLength + 1;
                assert(Checks.isVec(trueModes) ...
                    && all(arrayfun(@(mode) Checks.isScalarIn(mode, 1, numCaModes) && mod(mode, 1) == 0, trueModes)), ...
                    'ResourceAwareRecedingHorizonController:CheckModeMeasurementsAndDelays:InvalidModeObservations', ...
                    '** Mode observations must be a given as a vector of integers from {1, ..., %d} **', ...
                        numCaModes);
                
                assert(Checks.isVec(modeDelays, numel(trueModes)) ...
                    && all(arrayfun(@(x) Checks.isNonNegativeScalar(x)&& mod(x, 1) == 0, modeDelays)), ...
                    'ResourceAwareRecedingHorizonController:CheckModeMeasurementsAndDelays:InvalidModeDelays', ...
                    '** Mode observation delays must be a given as a vector of %d nonnegative integers **', ...
                    numel(trueModes));
            end
        end
        
        %% getApplicableModeMeasurements
        function [applicableModes, applicableModeDelays] = getApplicableModeMeasurements(this, trueModes, modeDelays)
            % find the mode observations with valid delays
            idx = find(modeDelays <= this.maxMeasDelay);
            applicableModes = trueModes(idx);
            applicableModeDelays = modeDelays(idx);
            if numel(idx) ~= numel(modeDelays)
                warning('ResourceAwareRecedingHorizonController:IgnoringModeObservations', ...
                    '** %s\nIgnoring %d of %d mode observations. **', 'Delay too large', ...
                    numel(modeDelays) - numel(idx), numel(modeDelays));
            end
        end
        
        %% checkMeasurementsAndDelays
        function checkMeasurementsAndDelays(this, measurements, measDelays)
            if ~isempty(measurements)
                assert(Checks.isFixedRowMat(measurements, size(this.C, 1)), ...
                    'ResourceAwareRecedingHorizonController:CheckMeasurementsAndDelays:InvalidMeas', ...
                    '** Individual measurements must be %d-dimensional **', ...
                    size(this.C, 1));

                numMeas = size(measurements, 2);
                assert(Checks.isNonNegativeVec(measDelays, numMeas) ...
                        && all(arrayfun(@(delay) mod(delay, 1) == 0, measDelays)), ...
                    'ResourceAwareRecedingHorizonController:CheckMeasurementsAndDelays:InvalidMeasDelay', ...
                    '** Each measurement delay (%d in total) must be a nonnegative integer **', numMeas);
            end            
        end
        
        %% getApplicableMeasurements
        function [applicableMeas, applicableDelays] = getApplicableMeasurements(this, measurements, measDelays)
            numMeas = numel(measDelays);
            % the controller assumes only a single measurement per time step
            [~, uniqueIdx, ~] = unique(measDelays, 'stable');
            assert(numel(uniqueIdx) == numMeas, ...
                'ResourceAwareRecedingHorizonController:GetApplicableMeasurements:IgnoringMeasurementsNotUnique', ...
                '** %s\nIgnoring %d of %d measurements. **', ...
                'Controller assumes only one measurement per time step', ...
                numMeas - numel(uniqueIdx), numMeas);
  
            % find the measurements with valid delays
            idx = find(measDelays(uniqueIdx) <= this.maxMeasDelay);
            applicableMeas = measurements(:, idx);
            applicableDelays = measDelays(idx);
            if numel(idx) ~= numel(uniqueIdx)
                warning('ResourceAwareRecedingHorizonController:GetApplicableMeasurements:IgnoringMeasurementsDelayTooLarge', ...
                    '** %s\nIgnoring %d of %d measurements. **', ...
                    'Delay too large', ...
                    numel(uniqueIdx) - numel(idx), numel(uniqueIdx));
            end
        end
        
        %% incorporateDelayedModeMeasurements
        function currentMode = incorporateDelayedModeMeasurements(this, applicableModes, applicableModeDelays)
            currentMode = [];
            numModes = this.sequenceLength + 1;
            if ~isempty(applicableModeDelays)
                [contained, idx] = ismember(0, applicableModeDelays);
                if contained
                    currentMode = applicableModes(idx);
                    indices = [1:idx-1 idx+1:numel(applicableModeDelays)]; 
                else
                    indices = 1:numel(applicableModeDelays);
                end
                for i = indices
                    % delay of j time steps-> we now the input at time k-j
                    % so we can update the mode probs at that time to unit vector with 1 at position corresponding to that input, 
                    % which are stored in inputProbsHistory(j, :)
                    this.inputProbsHistory(applicableModeDelays(i), :) ...
                        = [zeros(1, applicableModes(i) - 1), 1 zeros(1, numModes - applicableModes(i))];
                end
            end
        end
        
        %% incorporateDelayedMeasurements
        function nonDelayedMeasurements = incorporateDelayedMeasurements(this, applicableMeasurements, applicableDelays)
            nonDelayedMeasurements = [];
            if ~isempty(applicableDelays)
                % get the available measurement delays
                occuringDelays = unique(applicableDelays); % sorted in ascending order (lowest delay first)
                % get the measurements per delay as a cell array of matrices
                % (lowest delay first); grouping by operation (column-wise)
                measurementsPerDelay = splitapply(@(group) {group}, applicableMeasurements, grp2idx(applicableDelays)');
                % update the history of measurements
                startIndex = 1;
                if occuringDelays(1) == 0
                    % measurements that are not delayed (i.e., from the current time step)
                    nonDelayedMeasurements = measurementsPerDelay{1};
                    startIndex = 2;
                end
                for i = startIndex:numel(occuringDelays)
                    % delay of j time steps-> store at j-th element of the
                    % history
                    this.measurementHistory{occuringDelays(i)} = measurementsPerDelay{i};
                end
            end
        end
        
    end   
end

