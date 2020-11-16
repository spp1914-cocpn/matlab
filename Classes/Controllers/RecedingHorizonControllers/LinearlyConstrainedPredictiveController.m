classdef LinearlyConstrainedPredictiveController < SequenceBasedTrackingController
    % Implementation of a predictive controller (based on a quadratic cost function) for linear networked
    % control systems which can deal with linear state and input
    % constraints of the form a'*x_k <= c, b'*u_k <= d.
    % Also, tracking a reference trajectory z_ref according to a
    % performance output z_k = Z*x_k is supported.
    %
    % Literature: 
    %   Achim Hekler, Jörg Fischer, and Uwe D. Hanebeck,
    %   Control over Unreliable Networks Based on Control Input Densities,
    %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
    %   Singapore, Singapore, July 2012.
    %
    %   Alberto Bemporad, Manfred Morari, Vivek Dua, and Efstratios N. Pistikopoulos
    %   The explicit linear quadratic regulator for constrained systems,
    %   Automatica,
    %   Volume 38, Issue 1, pp. 481-485, 2002
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        A;
        B;
        Q;
        R;        
  
        F; % F matrix and its powers F^2, F^3, ..., F^N, with N = sequenceLength -2
        G;
        H;
        alphas;
    end
    
    properties (Access = private)
        etaState;
    end
    
    properties (SetAccess = private, GetAccess = protected)
        stateConstraintWeightings;
        inputConstraintWeightings;
        stateConstraints;
        inputConstraints;
        solver;
    end
    
    properties (SetAccess = private, GetAccess = public)
        horizonLength;
    end   
    
    methods (Access = public)
        %% LinearlyConstrainedPredictiveController
        function this = LinearlyConstrainedPredictiveController(A, B, Q, R, sequenceLength, caDelayProbs, ...
                stateConstraintWeightings, stateConstraints, inputConstraintWeightings, inputConstraints, Z, refTrajectory)
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
            %   >> caDelayProbs (Nonnegative vector with elements summing to 1)
            %      The vector describing the delay distribution of the
            %      CA-network.
            %
            %   >> stateConstraintWeightings (Matrix, dimPlantState-by-s)
            %      The weightings a of the individual
            %      state constraints a'*x_k <= c, column-wise arranged.
            %      If no state constraints are imposed, pass the empty matrix
            %      here.
            %
            %   >> stateConstraints (Vector with s entries)
            %      The individual state constraints c in the form a'*x_k <= c, 
            %      where the i-th entry corresponds to the i-th state constraint.
            %      If no state constraints are imposed, pass the empty matrix
            %      here.
            %
            %   >> inputConstraintWeightings (Matrix, dimPlantInput-by-r)
            %      The weightings b of the individual
            %      input constraints b'*u_k <= d, column-wise arranged.
            %      If no input constraints are imposed, pass the empty matrix
            %      here.
            %
            %   >> inputConstraints (Vector with r entries)
            %      The individual input constraints d in the form b'*u_k <= d, 
            %      where the i-th entry corresponds to the i-th input constraint.
            %      If no input constraints are imposed, pass the empty matrix
            %      here.
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
            %   << this (LinearlyConstrainedPredictiveController)
            %      A new LinearlyConstrainedPredictiveController instance.
            %
            
            % do not require changing setpoints or a reference to be
            % tracked
            if nargin == 10
                Z = [];
                refTrajectory = [];
            elseif nargin ~= 12 
                error('LinearlyConstrainedPredictiveController:InvalidNumberOfArguments', ...
                    '** Constructor must be called with either 10 or 12 arguments **');
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
            
            Validator.validateDiscreteProbabilityDistribution(caDelayProbs);
            probs = Utility.truncateDiscreteProbabilityDistribution(caDelayProbs, sequenceLength + 1);
            
            dimEta = this.dimPlantInput * (sequenceLength * (sequenceLength - 1) / 2);
            this.F = zeros(dimEta, dimEta, sequenceLength - 2);
            % we need F, G, H
            [this.F(:, :, 1), this.G, this.H, ~] = Utility.createActuatorMatrices(sequenceLength, this.dimPlantInput);
          
            % compute the required powers of F
            for i=2:sequenceLength-2
                this.F(:, :, i)  = this.F(:, :, i-1) * this.F(:, :, 1);
            end
            
            this.etaState = zeros(dimEta, 1);
            
            % we can compute the weighting factors alpha_k in advance for
            % the whole horizon
            % for x_k+1 we have numModes possible inputs (incl. default
            % input) -> expected u_k based on alpha_k|k
            % for x_k+2 we have numModes-1 possible inputs (incl. default
            % input) -> expected u_k+1|k based on alpha_k+1|k
            % for x_k+3 we have numModes-2 possible inputs (incl. default
            % input) -> expected u_k+2
            % and so forth            
            P = Utility.calculateDelayTransitionMatrix(probs);
            sums = cumsum(probs);
            this.alphas = cell(1, sequenceLength);            
            this.alphas{1} = zeros(sequenceLength+1, 1);
            
            this.alphas{1}(1) = probs(1);
            for j=2:sequenceLength+1
                this.alphas{1}(j) = prod(1-sums(1:j-1)) * sums(j);
            end
            %this.alphas{1}(sequenceLength) = 1 - sum(this.alphas{1}(1:sequenceLength-1));
            for k=1:sequenceLength-1
                this.alphas{k+1} = (P')^k*this.alphas{1};
                this.alphas{k+1} = this.alphas{k+1}(k+1:end) / sum(this.alphas{k+1}(k+1:end));
            end

            % constraints
            if ~isempty(stateConstraints)
                this.validateStateConstraints(stateConstraintWeightings, stateConstraints);
            end
            if ~isempty(inputConstraints)
                this.validateInputConstraints(inputConstraintWeightings, inputConstraints);
            end
            this.stateConstraintWeightings = stateConstraintWeightings;
            this.inputConstraintWeightings = inputConstraintWeightings;
            this.stateConstraints = stateConstraints;
            this.inputConstraints = inputConstraints;          
            
            this.solver = this.createConstrainedOptimizationProblem();            
        end
        
        %% reset
        function reset(this)
        end
        
        %% getStateConstraints
        function [stateConstraints, constraintWeightings] = getStateConstraints(this)
            % Get the considered state constraints.
            %
            % Returns:
            %   << stateConstraints (Vector with s entries)
            %      The individual state constraints c in the form a'*x_k <= c, 
            %      where the i-th entry corresponds to the i-th state constraint.
            %
            %   << constraintWeightings (Matrix, dimPlantState-by-s)
            %      The weightings a of the individual
            %      state constraints a'*x_k <= c, column-wise arranged.
            %
            stateConstraints = this.stateConstraints;
            constraintWeightings = this.stateConstraintWeightings;
        end
        
        %% getInputConstraints
        function [inputConstraints, constraintWeightings] = getInputConstraints(this)
            % Get the considered input constraints.
            %
            % Returns:
            %   << inputConstraints (Vector with s entries)
            %      The individual input constraints d in the form b'*u_k <= d, 
            %      where the i-th entry corresponds to the i-th input constraint.
            %
            %   << constraintWeightings (Matrix, dimPlantState-by-s)
            %      The weightings b of the individual
            %      input constraints b'*u_k <= d, column-wise arranged.
            %
            inputConstraints = this.inputConstraints;
            constraintWeightings = this.inputConstraintWeightings;
        end
        
        %% changeHorizonLength
        function changeHorizonLength(this, newHorizonLength)
            % Change the length of the optimization horizon to be considered in
            % the future.
            %
            % Parameters:
            %  >> newHorizonLength (Positive integer)
            %     The new length of the optimization horizon to be considered by the controller.
            %
            this.horizonLength = newHorizonLength;
            % the optimization problem to solve changes
            this.solver = this.createConstrainedOptimizationProblem();
        end
        
        %% changeStateConstraints
        function changeStateConstraints(this, newStateConstraintWeightings, newStateConstraints)
            % Change the state constraints to be considered in the future.
            % To remove existing state constraints, simply pass the empty
            % matrix as second paramter. In that case, the value of the
            % first parameter is ignored.
            %
            % Parameters:
            %   >> newStateConstraintWeightings (Matrix, dimPlantState-by-s)
            %      The weightings a of the individual
            %      state constraints a'*x_k <= c, column-wise arranged.
            %
            %   >> newStateConstraints (Vector with s entries)
            %      The individual state constraints c in the form a'*x_k <= c, 
            %      where the i-th entry corresponds to the i-th state constraint.
            %
            if ~isempty(newStateConstraints)
                this.validateStateConstraints(newStateConstraintWeightings, newStateConstraints);
            end 
            this.stateConstraintWeightings = newStateConstraintWeightings;
            this.stateConstraints = newStateConstraints;
            % the optimization problem to solve changes
            this.solver = this.createConstrainedOptimizationProblem();
        end
        
        %% changeInputConstraints
        function changeInputConstraints(this, newInputConstraintWeightings, newInputConstraints)
            % Change the input constraints to be considered in the future.
            % To remove existing input constraints, simply pass the empty
            % matrix as second paramter. In that case, the value of the
            % first parameter is ignored.
            %
            % Parameters:
            %   >> newInputConstraintWeightings (Matrix, dimPlantInput-by-r)
            %      The weightings b of the individual
            %      input constraints b'*u_k <= d, column-wise arranged.
            %
            %   >> newInputConstraints (Vector with r entries)
            %      The individual input constraints d in the form b'*u_k <= d, 
            %      where the i-th entry corresponds to the i-th input constraint.
            %
                                    
            if ~isempty(newInputConstraints)
                this.validateInputConstraints(newInputConstraintWeightings, newInputConstraints);
            end
            this.inputConstraintWeightings = newInputConstraintWeightings;
            this.inputConstraints = newInputConstraints;
            % the optimization problem to solve changes
            this.solver = this.createConstrainedOptimizationProblem();
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
                    'LinearlyConstrainedPredictiveController:SetSequenceLength:InvalidSequenceLength', ...
                    '** <seqLength> must not exceed the current optimization horizon, i.e., {1, ... %d} **', ...
                    this.horizonLength);
            end
        end
        
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, ~, timestep)
            % time step is only required for picking the correct reference
            % trajectory (if present)
            assert(Checks.isPosScalar(timestep) && mod(timestep, 1) == 0, ...
                'LinearlyConstrainedPredictiveController:DoControlSequenceComputation:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be a positive integer **');
            
            [stateMean, ~] = state.getMeanAndCov();
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
                solverVars = {this.etaState, stateMean, zRef};
            else
                solverVars = {this.etaState, stateMean};
            end
            % return as column vector
            [inputs, flag] = this.solver(solverVars);            
            if flag == 1
                warning('LinearlyConstrainedPredictiveController:DoControlSequenceComputation:ProblemInfeasible', ...
                    '** Optimization problem seems to be infeasible. Returning zero input **');
                inputSequence = zeros(this.sequenceLength * this.dimPlantInput, 1);
                
                %eta_k+1=F*eta_k
                this.etaState = this.F(:, :, 1) * this.etaState;
            else 
                inputSequence = reshape(inputs(:, 1:this.sequenceLength), this.sequenceLength * this.dimPlantInput, 1);
                
                %eta_k+1=F*eta_k+G*U_k
                this.etaState = this.F(:, :, 1) * this.etaState + this.G * inputSequence;
            end            
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
                'LinearlyConstrainedPredictiveController:DoCostsComputation:InvalidStateTrajectory', ...
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
                    'LinearlyConstrainedPredictiveController:GetDeviationFromRefForState:InvalidTimestep', ...
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
            dimEta = this.dimPlantInput * (this.sequenceLength * (this.sequenceLength - 1) / 2);
            Rsqrt = chol(this.R); % upper Cholesky factor of R, R = Rsqrt' * Rsqrt
            % ensure that matrices are not symmetric
            states = sdpvar(this.dimPlantState, this.horizonLength + 1, 'full');            
            inputs = sdpvar(this.dimPlantInput, this.horizonLength, 'full');
            eta = sdpvar(dimEta, 1); % changes at runtime, so not fixed
            expectedInputs = sdpvar(this.dimPlantInput, this.horizonLength, 'full');
            
            if ~isempty(this.Z)
                % reference tracking
                zRef = sdpvar(this.dimRef, this.horizonLength + 1);
                % construct for driving differences to the origin
                performance = this.Z * states - zRef;
                
                solverVars = {eta, states(:, 1), zRef};
                
                [K, P, ~] = dlqr(this.A, this.B, this.Z'*this.Q*this.Z, this.R);              
                P = this.Z * P *this.Z'; % terminal weighting 
            else
                % construct for driving states to the origin
                performance = states;
                solverVars = {eta, states(:, 1)};
                
                [K, P, ~] = dlqr(this.A, this.B, this.Q, this.R);
            end           
            
            % we first construct the constraints
            constraints = [];            
            % terminal cost, use stabilizing solution of DARE           
            costs = performance(:, end)' * P * performance(:, end);
            for stage=1:this.horizonLength
                if stage < this.sequenceLength
                    if stage == 1
                        oldInputs = eta;
                    else
                        oldInputs = this.F(:, :, stage-1) * eta;
                    end
                    expectedInputs(:, stage) = inputs(:, stage) * this.alphas{stage}(1);
                    i = 2;
                    for j = stage+1:this.sequenceLength
                        expectedInputs(:, stage) = expectedInputs(:, stage) + this.H(:, :, j) * oldInputs * this.alphas{stage}(i);
                        i = i + 1;
                    end                    
                elseif stage== this.sequenceLength
                    % either the input to be computed or the zero default input
                    expectedInputs(:, this.sequenceLength) = inputs(:, this.sequenceLength) * this.alphas{this.sequenceLength}(1);
                else
                    % switch to unconstrained LQR
                    inputs(:, stage) = -K * states(:, stage);                 
                    expectedInputs(:, stage) = inputs(:, stage); % use computed input directly for prediction 
                end                                
                % constraints due to nominal system equation
                constraints = [constraints, states(:, stage + 1) == this.A * states(:, stage) + this.B * expectedInputs(:, stage)]; 
                costs = costs + sum((Rsqrt * inputs(:, stage)) .^2) + performance(:, stage)' * this.Q * performance(:, stage);
            end
            % add the constraints
            for j=1:numel(this.stateConstraints)
                constraints = [constraints, this.stateConstraintWeightings(:, j)' * states(:, 1:this.sequenceLength+1) <= this.stateConstraints(j)];
            end
            for j=1:numel(this.inputConstraints)
                constraints = [constraints, this.inputConstraintWeightings(:, j)' * inputs(:, 1:this.sequenceLength) <= this.inputConstraints(j)]; 
            end      
            options = sdpsettings;
            options.solver = 'quadprog';
            solver = optimizer(constraints, costs, options, solverVars, inputs);            
        end
                        
        %% validateStateConstraints
        function validateStateConstraints(this, stateWeightings, stateConstraints)
            assert(Checks.isFixedRowMat(stateWeightings, this.dimPlantState), ...
                'LinearlyConstrainedPredictiveController:ValidateStateConstraints:InvalidStateWeightings', ...
                '** Input parameter <stateWeightings> must be a matrix with the weightings (%d-dimensional) column-wise arranged **', ...
                this.dimPlantState);
 
            assert(Checks.isVec(stateConstraints, size(stateWeightings, 2)), ...
                'LinearlyConstrainedPredictiveController:ValidateStateConstraints:InvalidStateConstraints', ...
                '** Input parameter <stateConstraints> must be a vector with %d elements **', ...
                size(stateWeightings, 2));
        end
        
        %% validateInputConstraints
        function validateInputConstraints(this, inputWeightings, inputConstraints)
            assert(Checks.isFixedRowMat(inputWeightings, this.dimPlantInput), ...
                'LinearlyConstrainedPredictiveController:ValidateInputConstraints:InvalidInputWeightings', ...
                '** Input parameter <inputWeightings> must be a matrix with the weightings (%d-dimensional) column-wise arranged **', ...
                this.dimPlantState);

            assert(Checks.isVec(inputConstraints, size(inputWeightings, 2)), ...
                'LinearlyConstrainedPredictiveController:ValidateInputConstraints:InvalidInputConstraints', ...
                '** Input parameter <inputConstraints> must be a vector with %d elements **', ...
                size(inputWeightings, 2));
        end
    end      
end

