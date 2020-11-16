classdef NominalPredictiveController< SequenceBasedTrackingController
    % Implementation of a linear networked predictive controller based on a
    % nominal linear-quadratic regulator.
    %
    % Literature: 
    %   Guo-Ping Liu,
    %   Predictive Controller Design of Networked Systems With Communication Delays and Data Loss,
    %   IEEE Transactions on Circuits and Systems—II: Express Briefs,
    %   Vol. 57, No. 6, pp. 481-485, 2010.
    %
    %   Huibert Kwakernaak, and Raphael Sivan, 
    %   Linear Optimal Control Systems,
    %   Wiley-Interscience, New York, 1972.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        
    properties (SetAccess = private, GetAccess = protected)
        % (empty or column vector of dimension <dimU>, default zero)
        feedforward = [];
        recomputeGain = false;
        recomputeFeedforward = false;
        A;
        B;
    end
       
    properties (SetAccess = private, GetAccess = public)
        % controller gain matrix (negative of the LQR gain matrix)
        L = [];
        % system state cost matrix;
        % (positive semidefinite matrix of dimension <dimX> x <dimX>)
        Q = [];
        % input cost matrix;
        % (positive definite matrix of dimension <dimU> x <dimU>)
        R = [];

        % steady-state (setpoint to reach)
        % (empty or column vector of dimension <dimX>)
        setpoint = [];
    end
    
    methods (Access = public)
        %% NominalPredictiveController
        function this = NominalPredictiveController(A, B, Q, R, sequenceLength, setpoint)
            % Class constructor.
            %
            % Parameters:
            %  >> A (Square Matrix)
            %     The system matrix of the plant.
            %
            %  >> B (Matrix)
            %     The input matrix of the plant.
            %
            %  >> Q (Positive semi-definite matrix)
            %     The state weighting matrix in the controller's underlying cost function.
            %
            %  >> R (Positive definite matrix)
            %     The input weighting matrix in the controller's underlying cost function.
            %
            %  >> sequenceLength (Positive integer)
            %     The length of the input sequence (i.e., the number of
            %     control inputs) to be computed by the controller.
            %
            %  >> setpoint (Vector, optional)
            %     The set point in the state space that shall be tracked asymptotically instead of the origin.
            %     Note that it is not verified whether the resulting set point tracking problem is solvable. 
            %     For instance, in case dim(x) > dim(u), the problem will generally have no exact solution.
            %     If no set point is provided, the controller will attempt
            %     to drive the plant to the origin instead.
            %
            % Returns:
            %   << obj (NominalPredictiveController)
            %      A new NominalPredictiveController instance.
            %
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedTrackingController(dimX, dimU, sequenceLength, true, [], [], []);
            this.A = A;
            this.B = B;
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
                        
            this.computeAndSetGainMatrix();
            
            if nargin > 5 && ~isempty(setpoint)
                assert(Checks.isVec(setpoint, this.dimPlantState) && all(isfinite(setpoint)), ...
                    [class(this) ':InvalidSetpoint'], ...
                    '** <setpoint> is expected to be a real-valued, %d-dimensional vector', this.dimPlantState);
                this.setpoint = setpoint(:);
                this.computeAndSetFeedforward();
            else
                % no set point provided, so stead state is origin                
                this.feedforward = zeros(dimU, 1);
            end
        end
        
        %% reset
        function reset(~)
            % do nothing
        end
        
        %% changeSequenceLength
        function changeSequenceLength(this, newSequenceLength)
            % Change the length of the control sequences to be created in
            % the future.
            %
            % Parameters:
            %  >> newSequenceLength (Positive integer)
            %     The new length of the input sequences (i.e., the number of
            %     control inputs) to be computed by the controller.
            %
            this.sequenceLength = newSequenceLength;
        end
        
        %% changeCostMatrices
        function changeCostMatrices(this, newQ, newR)
            % Change the state and input weighting matrices in the controller's underlying cost function.
            %
            % Parameters:
            %  >> newQ (Positive semi-definite matrix)
            %     The new state weighting matrix in the controller's underlying cost function.
            %
            %  >> newR (Positive definite matrix)
            %     The new input weighting matrix in the controller's underlying cost function.
            %
            Validator.validateCostMatrices(newQ, newR, this.dimPlantState, this.dimPlantInput);
            this.Q = newQ;
            this.R = newR;
            % requires that gain matrix is recomputed
            this.recomputeGain = true;            
            % and, likewise, the feedforward, if nonzero set point is
            % present
            if ~isempty(this.setpoint)
                this.recomputeFeedforward = true;
            end
        end
        
        %% changeModelParameters
        function changeModelParameters(this, newA, newB)
            Validator.validateSystemMatrix(newA, this.dimPlantState);
            Validator.validateInputMatrix(newB, this.dimPlantState, this.dimPlantInput);                         
            
            this.A = newA;
            this.B = newB;
            
            % requires that gain matrix is recomputed
            this.recomputeGain = true;            
            % and, likewise, the feedforward, if nonzero set point is
            % present
            if ~isempty(this.setpoint)
                this.recomputeFeedforward = true;
            end
        end
        
        %% changeSetPoint
        function changeSetPoint(this, newSetpoint)
            % Change the set point in the state space that shall be tracked asymptotically instead of the origin.
            % Note that it is not verified whether the resulting set point
            % tracking problem is solvable. 
            % For instance, in case dim(x) > dim(u), the problem will
            % generally have no exact solution.
            %
            % Parameters:
            %  >> newSetpoint (Vector)
            %     The new setpoint to be tracked.
            %
            assert(Checks.isVec(newSetpoint, this.dimPlantState) && all(isfinite(newSetpoint)), ...
                'NominalPredictiveController:ChangeSetPoint:InvalidSetpoint', ...
                '** <newSetpoint> is expected to be a real-valued, %d-dimensional vector', this.dimPlantState);

            this.setpoint = newSetpoint(:);
            this.recomputeFeedforward = true;
        end
    end
    
    methods (Access = private)
        %% computeAndSetGainMatrix
        function computeAndSetGainMatrix(this)
            % compute the gain matrix 
            this.L = -dlqr(this.A, this.B, this.Q, this.R);
        end
        
        %% computeAndSetFeedforward
        function computeAndSetFeedforward(this)
            I = speye(this.dimPlantState);
            % we need a stationary gain
            F = this.B \ ((I - this.A - this.B * this.L) / I); % this is the DC gain of the closed loop system, with C = I
            this.feedforward = F * this.setpoint;
        end
    end
    
    methods (Access = protected)
         %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, plantState, varargin)
            if this.recomputeGain
                this.computeAndSetGainMatrix();                
                this.recomputeGain = false;
            end
            if this.recomputeFeedforward
                this.computeAndSetFeedforward();
                this.recomputeFeedforward = false;
            end
            
            [state, ~] = plantState.getMeanAndCov();
            
            inputSequence = zeros(this.sequenceLength * this.dimPlantInput, 1);
            ff = ~isempty(this.setpoint) * this.feedforward;

            for j = 1:this.sequenceLength
                inputSequence((j - 1) * this.dimPlantInput + 1: j * this.dimPlantInput) = this.L * state + ff;
                % predict the system state (i.e., compute estimate since noise is
                % zero-mean)
                % use the computed input u_k for prediction
                state = this.doStatePrediction(state, inputSequence((j - 1) * this.dimPlantInput + 1: j * this.dimPlantInput), j);
%                 state = this.A * state + ...
%                     this.B * inputSequence((j - 1) * this.dimPlantInput + 1: j * this.dimPlantInput);
            end
        end
        
        %% doStatePrediction
        function newState = doStatePrediction(this, state, computedInput, stage)
            % predict the system state (i.e., compute estimate since noise is
            % zero-mean)
            newState = this.A * state + this.B * computedInput;
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, ~)
            if isempty(this.setpoint)
                performance = state;
            else
                performance = state - this.setpoint;
            end
            stageCosts = Utility.computeStageCosts(performance, input, this.Q, this.R);
        end
        
        %% doCostsComputation
        function averageLQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonLength = size(appliedInputs, 2);
            assert(size(stateTrajectory, 2) == horizonLength + 1, ...
                'NominalPredictiveController:DoCostsComputation', ...
                '** <stateTrajectory> is expected to have %d columns ', horizonLength + 1);

            if isempty(this.setpoint)
                performance = stateTrajectory;
            else
                performance = stateTrajectory - this.setpoint;
            end
            averageLQGCosts = Utility.computeLQGCosts(horizonLength, performance, appliedInputs, this.Q, this.R) / horizonLength;
        end
        
        %% doGetDeviationFromRefForState
        function deviation = doGetDeviationFromRefForState(this, trueState, ~)
            if isempty(this.setpoint)
                % we track the origin
                deviation = trueState;
            else
                % we track a nonzero setpoint
                deviation = trueState - this.setpoint;
            end
        end
    end
end

