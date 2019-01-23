classdef LinearlyConstrainedPredictiveController < SequenceBasedTrackingController
    % Implementation of a predictive controller (based on a quadratic cost function) for linear networked
    % control systems which can deal with linear state and input
    % constraints of the form a'*x_k <= c, b'*u_k <= d.
    % Also, tracking a reference trajectory z_ref according to a
    % performance output z_k = Z*x_k is supported.
  
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        function this = LinearlyConstrainedPredictiveController(A, B, Q, R, sequenceLength, ...
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
            %   >> stateConstraintWeightings (Matrix, dimPlantState-by-s)
            %      The weightings a of the individual
            %      state constraints a'*x_k <= c, column-wise arranged.
            %
            %   >> stateConstraints (Vector with s entries)
            %      The individual state constraints c in the form a'*x_k <= c, 
            %      where the i-th entry corresponds to the i-th state constraint.
            %
            %   >> inputConstraintWeightings (Matrix, dimPlantInput-by-r)
            %      The weightings b of the individual
            %      input constraints b'*u_k <= d, column-wise arranged.
            %
            %   >> inputConstraints (Vector with r entries)
            %      The individual input constraints d in the form b'*u_k <= d, 
            %      where the i-th entry corresponds to the i-th input constraint.
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
            if nargin == 9
                Z = [];
                refTrajectory = [];
            elseif nargin ~= 11 
                error('LinearlyConstrainedPredictiveController:InvalidNumberOfArguments', ...
                    '** Constructor must be called with either 9 or 11 arguments **');
            end
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
                        
            this = this@SequenceBasedTrackingController(dimX, dimU, sequenceLength, Z, refTrajectory, []);
                       
                                  
            Validator.validateCostMatrices(Q, R, this.dimRef, dimU);
            this.Q = Q;
            this.R = R;
            this.A = A;
            this.B = B;
            this.horizonLength = this.sequenceLength;
                                           
            % constraints
            this.validateStateConstraints(stateConstraintWeightings, stateConstraints);
            this.validateInputConstraints(inputConstraintWeightings, inputConstraints);
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
        
        %% changeSequenceLength
        function changeSequenceLength(this, newSequenceLength)
            % Change the length of the control sequences to be created in
            % the future.
            %
            % Parameters:
            %  >> newSequenceLength (Positive integer)
            %     The new length of the input sequences (i.e., the number of
            %     control inputs) to be computed by the controller, which
            %     must not exceed the current length of the optimization
            %     horizon.
            %
            this.sequenceLength = newSequenceLength;           
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
            this.validateStateConstraints(newStateConstraintWeightings, newStateConstraints);
            
            this.stateConstraintWeightings = newStateConstraintWeightings;
            this.stateConstraints = newStateConstraints;
            % the optimization problem to solve changes
            this.solver = this.createConstrainedOptimizationProblem();
        end
        
        %% changeInputConstraints
        function changeInputConstraints(this, newInputConstraintWeightings, newInputConstraints)
            % Change the input constraints to be considered in the future.
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
            this.validateInputConstraints(newInputConstraintWeightings, newInputConstraints);
            
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
            assert(Checks.isPosScalar(timestep)  && mod(timestep, 1) == 0, ...
                'LinearlyConstrainedPredictiveController:DoControlSequenceComputation:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be a positive integer **');
            
            [stateMean, ~] = state.getMeanAndCov();
            if ~isempty(this.refTrajectory)
                zRef = this.refTrajectory(:, timestep:timestep + this.sequenceLength);
                solverVars = {stateMean, zRef};
            else
                solverVars = {stateMean};
            end
            % return as column vector
            [inputs, flag] = this.solver{solverVars};
            if flag == 1
                warning('LinearlyConstrainedPredictiveController:DoControlSequenceComputation:ProblemInfeasible', ...
                    '** Optimization problem seems to be infeasible. Returning zero input **');
                inputSequence = zeros(this.sequenceLength * this.dimPlantInput, 1);
            else
                inputSequence = reshape(inputs(:, 1:this.sequenceLength), this.sequenceLength * this.dimPlantInput, 1);
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
            Rsqrt = chol(this.R); % upper Cholesky factor of R, R = Rsqrt' * Rsqrt
            % ensure that matrices are not symmetric
            states = sdpvar(this.dimPlantState, this.horizonLength + 1, 'full');
            inputs = sdpvar(this.dimPlantInput, this.horizonLength, 'full');
            
            if ~isempty(this.Z)
                % reference tracking
                zRef = sdpvar(this.dimRef, this.horizonLength + 1);
                % construct for driving differences to the origin
                performance = this.Z * states - zRef;
                
                solverVars = {states(:, 1), zRef};
            else
                % construct for driving states to the origin
                performance = states;
                solverVars = {states(:, 1)};
            end
            
            % we first construct the constraints
            constraints = [];
            for j=1:numel(this.stateConstraints)
                constraints = [constraints; this.stateConstraintWeightings(:, j)' * states <= this.stateConstraints(j)];
            end
            for j=1:numel(this.inputConstraints)
                constraints = [constraints; this.inputConstraintWeightings(:, j)' * inputs <= this.inputConstraints(j)]; 
            end
                    
            costs = performance(:, end)' * this.Q * performance(:, end); % terminal costs
            for k=1:this.horizonLength
                % additional constraints due to nominal system equation
                constraints = [constraints; states(:, k + 1) == this.A * states(:, k) + this.B * inputs(:, k)]; 
                costs = costs + sum((Rsqrt * inputs(:, k)) .^2) + performance(:, k)' * this.Q * performance(:, k);
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

