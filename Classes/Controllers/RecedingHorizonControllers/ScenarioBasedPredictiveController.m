classdef ScenarioBasedPredictiveController ...
        < SequenceBasedTrackingController & CaDelayProbsChangeable & ModelParamsChangeable
    % Implementation of a predictive controller (based on a quadratic cost function) for linear networked
    % control systems which computes control sequences by optimizing all possible scenarios.
    % Also, tracking a reference trajectory z_ref according to a
    % performance output z_k = Z*x_k is supported.
    %
    % Literature:    
    %   Isabel Jurado, Pablo Millán, Daniel Quevedo, and Franciso R. Rubio
    %   Stochastic MPC with Applications to Process Control,
    %   International Journal of Control,
    %   Volume 88, Issue 4, pp. 792-800, 2015
    %
    %   Achim Hekler, Jörg Fischer, and Uwe D. Hanebeck,
    %   Control over Unreliable Networks Based on Control Input Densities,
    %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
    %   Singapore, Singapore, July 2012.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2022  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (SetAccess = immutable, GetAccess = protected)        
        Q;
        R;                
    end
    
    properties (SetAccess = private, GetAccess = ?ScenarioBasedPredictiveControllerTest)
        A;
        B;        
        
        modeTransitionMatrix;
        
        alphas;
        scenarioProbabilities; % probabilities for each scenario/branch
                
        % store all relevant control sequences from time k to k-N, where N=sequence length
        % assuming zero default input
        bufferedControlInputSequences;
    end
    
    properties (SetAccess = private, GetAccess = protected)       
        recreateScenarioTree = false;
        solver; % yalmip optimizer, uses the scenario tree
    end
    
    properties (SetAccess = private, GetAccess = public)
        horizonLength;
    end   
    
    methods (Access = public)
        %% ScenarioBasedPredictiveController
        function this = ScenarioBasedPredictiveController(A, B, Q, R, sequenceLength, modeTransitionMatrix, ...
                lazyInitOptimProb, Z, refTrajectory)
            % Class constructor.
            %
            % Parameters:
            %   >> A (Square Matrix)
            %      The system matrix of the plant.
            %
            %   >> B (Matrix)
            %      The input matrix of the plant.
            %
            %   >> Q (Positive semi-definite matrix)
            %      The performance/plant output weighting matrix in the
            %      controller's underlying cost function, if a tracking
            %      task is performed, the state weighting matrix otherwise.
            %
            %   >> R (Positive definite matrix)
            %      The input weighting matrix in the controller's underlying cost function.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller, which
            %      (by default) equals the prediction horizon employed.
            %
            %   >> modeTransitionMatrix (Stochastic matrix, i.e. a square matrix with nonnegative entries whose rows sum to 1)
            %      The transition matrix of the mode theta_k of the augmented dynamics.
            %
            %   >> lazyInitOptimProb (Flag, i.e., a logical scalar)
            %      Flag to indicate whether the setup of the scenario tree used during the optimization shall be
            %      done during object instantiation, i.e., during this
            %      constructor call (lazyInitOptimProb = false), or shall be deferred until the first control sequence is to be computed (lazyInitOptimProb = true).
            %
            %   >> Z (Matrix, n-by-dimPlantState, optional)
            %      The time-invariant plant output (performance) matrix, 
            %      i.e., z_k = Z*x_k, if a tracking task is performed.
            %      If the state shall be driven to the origin, this
            %      argument should be left out.
            %
            %   >> referenceTrajectory (Matrix, optional)
            %      The reference trajectory to track, given as a matrix
            %      with the reference plant outputs column-wise arranged, if a tracking task is performed.
            %      If the state shall be driven to the origin, this
            %      argument should be left out.      
            %
            % Returns:
            %   << this (ScenarioBasedPredictiveController)
            %      A new ScenarioBasedPredictiveController instance.
            %
            
            % do not require changing setpoints or a reference to be
            % tracked
            if nargin == 7
                Z = [];
                refTrajectory = [];
            elseif nargin ~= 9 
                error('ScenarioBasedPredictiveController:InvalidNumberOfArguments', ...
                    '** Constructor must be called with either 7 or 9 arguments **');
            end
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
                        
            this = this@SequenceBasedTrackingController(dimX, dimU, sequenceLength, true, Z, refTrajectory, []);
                                                         
            Validator.validateCostMatrices(Q, R, this.dimRef, dimU);
            this.Q = Q;
            this.R = R;
            this.A = A;
            this.B = B;
            this.horizonLength = this.sequenceLength;
            
            % assume zero default input
            this.bufferedControlInputSequences = repmat(zeros(dimU, 1), [1 this.sequenceLength this.sequenceLength]);
            
            Validator.validateTransitionMatrix(modeTransitionMatrix, this.sequenceLength + 1);
            this.modeTransitionMatrix = modeTransitionMatrix;
            
            % initially, no inputs at actuator, so input is default input
            inputProbs = zeros(this.sequenceLength, 1);
            inputProbs(end + 1) = 1;
            this.updateInputProbs(inputProbs);
            this.updateScenarioProbs();
            
            if lazyInitOptimProb
                this.recreateScenarioTree = true;
            else
                this.solver = this.createConstrainedOptimizationProblem();           
            end
        end
        
        %% setBufferedControlSequences
        function setBufferedControlSequences(this, newBufferedSequences)
            % mainly for testing purposes

            assert(isequal(size(newBufferedSequences), size(this.bufferedControlInputSequences)) ...
                && all(isfinite(newBufferedSequences(:))), ...
                'ScenarioBasedPredictiveController:SetBufferedControlSequences:InvalidBufferedSequences', ...
                '** <newBufferedSequences> must be %d-by-%d-by-%d-dimensional **', ...
                size(newBufferedSequences, 1), size(newBufferedSequences, 2), size(newBufferedSequences, 3));
            
            this.bufferedControlInputSequences = newBufferedSequences;
        end
        
        %% reset
        function reset(this)
            % reset buffer
            % assume zero default input
            this.bufferedControlInputSequences = repmat(zeros(this.dimPlantInput, 1), [1 this.sequenceLength this.sequenceLength]);
        end
        
        %% changeModelParameters
        function changeModelParameters(this, newA, newB, ~)
            Validator.validateSystemMatrix(newA, this.dimPlantState);
            Validator.validateInputMatrix(newB, this.dimPlantState, this.dimPlantInput);                         
            
            this.A = newA;
            this.B = newB;
            
            % the optimization problem to solve changes
            this.recreateScenarioTree = true;
        end
        
        %% changeCaDelayProbs
        function changeCaDelayProbs(this, newCaDelayProbs)
            % Change the distribution of the control packet delays to be assumed by the controller.
            %
            % Parameters:
            %  >> newCaDelayProbs (Nonnegative vector)
            %     Vector specifiying the new delay distribution.
            %
            
            newMat = Utility.calculateDelayTransitionMatrix(...
                Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.sequenceLength + 1));
            if ~isequal(newMat, this.modeTransitionMatrix)
                this.modeTransitionMatrix = newMat;
            end
        end             
        
        %% changeHorizonLength
        function changeHorizonLength(this, newHorizonLength, lazyInitOptimProb)
            % Change the length of the optimization horizon to be considered in
            % the future.
            %
            % Parameters:
            %  >> newHorizonLength (Positive integer)
            %     The new length of the optimization horizon to be considered by the controller.
            %
            %   >> lazyInitOptimProb (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the new scenario tree shall be setup now, i.e., during this
            %      call (lazyInitOptimProb = false), or shall be deferred until the next control sequence is to be computed (lazyInitOptimProb = true).
            %      If left out, the default value true is used.
            %
            if nargin < 3
                lazyInitOptimProb = true;
            end
            
            this.horizonLength = newHorizonLength;
            % the optimization problem to solve changes
            if lazyInitOptimProb
                this.recreateScenarioTree = true;
            else
                this.solver = this.createConstrainedOptimizationProblem();
                this.recreateScenarioTree = false;
            end
        end 
        
        %% setScenarioTree
        function setScenarioTree(this, scenarioTree)
            % mainly for testing purposes
            
            this.solver = scenarioTree;
            this.recreateScenarioTree = false;
        end
        
        %% getScenarioTree
        function scenarioTree = getScenarioTree(this)
             % Get current scenario tree used for the compuation of U_k.
            %
            % Returns:
            %   << scenarioTree (Optimizer, i.e., a yalmip container for optimization problems)
            %      The current scenario tree, might be empty.
            %            
            scenarioTree = this.solver;
        end
    end
    
    methods (Access = protected)
         %% validateSequenceLengthOnSet
        function validateSequenceLengthOnSet(this, seqLength)
            % This function is called upon change of the sequence length to
            % validate the new value. 
            % It checks whether the value is a
            % positive integer and does not exceed the current optimization horizon, and errors if not.
            %
            % Parameters:
            %   >> seqLength (Positive integer)
            %      The desired new value of the sequence length.
            %
            if ~isempty(this.horizonLength)
                assert(seqLength <= this.horizonLength, ...
                    'ScenarioBasedPredictiveController:SetSequenceLength:InvalidSequenceLength', ...
                    '** <seqLength> must not exceed the current optimization horizon, i.e., {1, ... %d} **', ...
                    this.horizonLength);
            end
        end       
        
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, ~, timestep)
            % time step is only required for picking the correct reference
            % trajectory (if present)
            assert(Checks.isPosScalar(timestep) && mod(timestep, 1) == 0, ...
                'ScenarioBasedPredictiveController:DoControlSequenceComputation:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be a positive integer **');                      
            
            if this.recreateScenarioTree                 
                this.solver = this.createConstrainedOptimizationProblem();
                this.recreateScenarioTree = false;
            end
            
            if Checks.isClass(state, 'GaussianMixture')
                % filter provides Gaussian mixture
                % input probs are equal to mode probs
                [means, covs, modeProbs] = state.getComponents();
                [stateMean, ~] = Utils.getGMMeanAndCov(means, covs, modeProbs);
            else
                [stateMean, ~] = state.getMeanAndCov();
                % compute the mode probs directly
                modeProbs = this.modeTransitionMatrix' * this.alphas{1};
            end
            
            % update input probs for the scenario tree and, accordingly, the scenario pobs
            this.updateInputProbs(modeProbs);
            this.updateScenarioProbs();
                        
            if ~isempty(this.refTrajectory)
                % the preview of the reference must cover the whole
                % optimization horizon
                if timestep + this.horizonLength > size(this.refTrajectory, 2)
                    % too short, so append the last value
                    zRef = [this.refTrajectory(:, timestep:end), repmat(this.refTrajectory(:, end), 1, ...
                        timestep + this.horizonLength - size(this.refTrajectory, 2))];
                else
                    zRef = this.refTrajectory(:, timestep:timestep + this.horizonLength);
                end
                solverVars = [{this.scenarioProbabilities}, {this.bufferedControlInputSequences}, {stateMean}, {zRef}];
            else
                solverVars = {this.scenarioProbabilities, this.bufferedControlInputSequences, stateMean};                
            end
  
            % return as column vector            
            [inputs, flag] = this.solver(solverVars);           
            if flag == 1
                warning('ScenarioBasedPredictiveController:DoControlSequenceComputation:ProblemInfeasible', ...
                    '** Optimization problem seems to be infeasible. Returning zero input **');
                inputSequence = zeros(this.sequenceLength * this.dimPlantInput, 1);                
            else                
                numEls = min(this.sequenceLength, this.horizonLength);
                inputSequence = reshape(inputs(:, 1:numEls), numEls * this.dimPlantInput, 1);
                if numEls < this.sequenceLength
                    % take care of the corner case: horizon length < sequence length
                    % if so, repeat last entry
                    addInputs = repmat(inputs(:, end), 1, this.sequenceLength - numEls);
                    inputSequence = [inputSequence; addInputs(:)];
                end                
            end
            
            % update buffer
            this.bufferedControlInputSequences = circshift(this.bufferedControlInputSequences, 1, 3);            
            this.bufferedControlInputSequences(:,:, 1) = reshape(inputSequence, [], this.sequenceLength);     
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, timestep)
            if ~isempty(this.Z)
                % compute performance output and difference to reference
                % trajectory
                performance = this.Z * state - this.refTrajectory(:, timestep);
            else
                performance = state;
            end
            
            stageCosts = Utility.computeStageCosts(performance, input, this.Q, this.R);
        end
        
        %% doCostsComputation
        function costs = doCostsComputation(this, stateTrajectory, appliedInputs)
            numInputs = size(appliedInputs, 2);
            assert(size(stateTrajectory, 2) == numInputs + 1, ...
                'ScenarioBasedPredictiveController:DoCostsComputation:InvalidStateTrajectory', ...
                '** <stateTrajectory> is expected to have %d columns ', numInputs + 1);
            
            if ~isempty(this.Z)
                % compute performance output and difference to reference
                % trajectory
                performance = this.Z * stateTrajectory - this.refTrajectory(:, 1:numInputs + 1);
            else
                performance = stateTrajectory;
            end            
            costs = Utility.computeLQGCosts(numInputs, performance, appliedInputs, this.Q, this.R);
        end
        
        %% doGetDeviationFromRefForState
        function deviation = doGetDeviationFromRefForState(this, state, timestep)
            if ~isempty(this.Z)
                % we track a reference trajectory
                assert(Checks.isScalarIn(timestep, 1, size(this.refTrajectory, 2)) && mod(timestep, 1) == 0, ...
                    'ScenarioBasedPredictiveController:GetDeviationFromRefForState:InvalidTimestep', ...
                    '** Input parameter <timestep> must be in {1, ... %d} **', ...
                    size(this.refTrajectory, 2));

                deviation = this.Z * state - this.refTrajectory(:, timestep);
            else
                % we track the origin
                deviation = state;
            end
        end
    end
    
    methods (Access = private)       
        %% createConstrainedOptimizationProblem
        function solver = createConstrainedOptimizationProblem(this)
            Rsqrt = chol(this.R); % upper Cholesky factor of R, R = Rsqrt' * Rsqrt                    
            inputs = sdpvar(this.dimPlantInput, this.horizonLength, 'full');
            bufferedInputSequences = sdpvar(this.dimPlantInput, this.sequenceLength, this.sequenceLength, 'full');            
            initialState = sdpvar(this.dimPlantState, 1);
                        
            allIdx = this.getAllIndexCombinations();            
            
            numScenarios = size(allIdx, 1);        
            
            scenarioStates = sdpvar(this.dimPlantState, this.horizonLength + 1, numScenarios, 'full'); 
            scenarioProbs = sdpvar(1, numScenarios); 
                        
            if ~isempty(this.Z)                
                % reference tracking
                zRef = sdpvar(this.dimRef, this.horizonLength + 1, 'full');
                % construct for driving differences to the origin
                performance = sdpvar(this.dimRef, this.horizonLength + 1, numScenarios, 'full');
                for j=1:numScenarios
                    for k=1:this.horizonLength + 1
                        performance(:, k, j) = this.Z * scenarioStates(:, k, j) - zRef(:, k);
                    end
                end                
                
                solverVars = {scenarioProbs, bufferedInputSequences, initialState, zRef};                
                
                [P, ~, ~, info] = idare(this.A, this.B, this.Z'*this.Q*this.Z, this.R);
                % compute terminal weighting
                if info.Report > 1
                    % either report == 2: The solution is not finite
                    % or report == 3: No solution found since the symplectic
                    % matrix has eigenvalues on the unit circle
                    P = this.Q;
                    warning('ScenarioBasedPredictiveController:InvalidPlant', ...
                        '** (A,B) seems not stabilizable, cannot use stabilizing solution of associated DARE as terminal weighting Q_N **');
                end  
                P = this.Z * P *this.Z'; % terminal weighting 
            else                 
                % construct for driving states to the origin
                performance = scenarioStates;
                solverVars = {scenarioProbs, bufferedInputSequences, initialState};
                
                [P, ~, ~, info] = idare(this.A, this.B, this.Q, this.R);
                % compute terminal weighting
                if info.Report > 1
                    % either report == 2: The solution is not finite
                    % or report == 3: No solution found since the symplectic
                    % matrix has eigenvalues on the unit circle
                    P = this.Q;
                    warning('ScenarioBasedPredictiveController:InvalidPlant', ...
                        '** (A,B) seems not stabilizable, cannot use stabilizing solution of associated DARE as terminal weighting Q_N **');
                end
            end
            
            defaultInput = zeros(this.dimPlantInput, 1);
            constraints = [];
            costsPerScenario = sdpvar(1, numScenarios);
            % terminal cost, use stabilizing solution of DARE
            for j=1:numScenarios                
                costsPerScenario(j) = performance(:, end, j)' * P * performance(:, end, j);                
                constraints = [constraints, scenarioStates(:, 1, j) == initialState]; %#ok
            end
            numPossibleInputs = [this.sequenceLength + 1:-1:2];
            
            for stage=1:this.horizonLength
                if stage <= this.sequenceLength
                    possibleInputs = cell(1, numPossibleInputs(stage));
                    possibleInputs{1} = inputs(:, stage);                    
                    i = 2;                    
                    for p=1:this.sequenceLength - stage                        
                        possibleInputs{i} = bufferedInputSequences(:, p + stage, p);
                        i= i + 1;
                    end                    
                    possibleInputs{end} = defaultInput;                    
                    for j=1:numScenarios
                        stageIdx = allIdx(j, stage);                                             
                        
                        constraints = [constraints, scenarioStates(:, stage + 1, j) == this.A * scenarioStates(:, stage, j) + this.B * possibleInputs{stageIdx}]; %#ok                        
                        costsPerScenario(j) = costsPerScenario(j) + sum((Rsqrt * possibleInputs{stageIdx}) .^2) + performance(:, stage, j)' * this.Q * performance(:, stage, j);                        
                    end
                else
                    % only one input possible 
                    for j=1:numScenarios     
                        % constraint due to nominal system equation
                        constraints = [constraints, scenarioStates(:, stage + 1, j) == this.A * scenarioStates(:, stage, j) + this.B * inputs(:, stage)]; %#ok                        
                        costsPerScenario(j) = costsPerScenario(j) + sum((Rsqrt * inputs(:, stage)) .^2) + performance(:, stage, j)' * this.Q * performance(:, stage, j);                        
                    end
                end
            end           
    
            options = sdpsettings;
            options.solver = 'quadprog';            
            solver = optimizer(constraints, dot(costsPerScenario, scenarioProbs), options, solverVars, inputs); 
        end
         
        
        %% updateInputProbs
        function updateInputProbs(this, currInputProbs)
            % we can compute the weighting factors alpha_k in advance for
            % the whole horizon
            % for x_k+1 we have numModes possible inputs (incl. default
            % input) -> expected u_k based on alpha_k|k
            % for x_k+2 we have numModes-1 possible inputs (incl. default
            % input) -> expected u_k+1|k based on alpha_k+1|k
            % for x_k+3 we have numModes-2 possible inputs (incl. default
            % input) -> expected u_k+2
            % and so forth
            % the number of factors thus decreases over the horizon   
            
            % update the alphas, based on received inputs probs/mode probs from filter
            this.alphas = cell(1, this.sequenceLength);
            this.alphas{1}  = currInputProbs(:);
            
            % open loop prediction
            for k=1:this.sequenceLength - 1
                this.alphas{k+1} = (this.modeTransitionMatrix')^k*this.alphas{1};
                this.alphas{k+1} = this.alphas{k+1}(k+1:end) / sum(this.alphas{k+1}(k+1:end));
            end
        end
        
        %% updateScenarioProbs
        function updateScenarioProbs(this)
            % compute the scenario probs
            numScenarios = factorial(this.sequenceLength + 1);
            this.scenarioProbabilities = ones(1, numScenarios);
                        
            allIdx = this.getAllIndexCombinations();
            for j=1:numScenarios
                currIdx = allIdx(j, :);
                for stage=1:this.sequenceLength
                    this.scenarioProbabilities(j) = this.scenarioProbabilities(j) * this.alphas{stage}(currIdx(stage));
                end
            end
        end           
        
        %% getAllIndexCombinations
        function allIdx = getAllIndexCombinations(this)             
            sets = cellfun(@(w) 1:numel(w), this.alphas, 'UniformOutput', false);
            v = cell(this.sequenceLength, 1);
            [v{:}] = ndgrid(sets{:});
            allIdx = reshape(cat(this.sequenceLength + 1,v{:}), [], this.sequenceLength);           
        end        
    end      
end

