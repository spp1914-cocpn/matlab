classdef PolePlacementPredictiveController < SequenceBasedTrackingController & ModelParamsChangeable
    % Implementation of a linear networked predictive controller 
    % based on pole placement design.
    %
    % Literature: 
    %   J. Kautsky, N. K. Nichols, and P. van Dooren
    %   Robust Pole Assignment in Linear State Feedback,
    %   International Journal of Control,
    %   Vol. 41, No. 6, pp. 1129-1155, 1985.
    
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
    
    properties (SetAccess = private)
        A;
        B;
        polesCont;
        samplingInterval(1,1) double {mustBePositive} = 0.1;
    end
    
    properties (SetAccess = private, GetAccess = protected)
        % (empty or column vector of dimension <dimU>, default zero)
        feedforward = [];
        recomputeGain = false;
        recomputeFeedforward = false;
    end
    
    properties (SetAccess = private, GetAccess = public)
        % controller gain matrix
        L = []; % negative feedback, u_k = -L*x_k
        poles;

        % steady-state (setpoint to reach)
        % (empty or column vector of dimension <dimX>)
        setpoint = [];
    end
    
    methods (Access = public)
        %% PolePlacementPredictiveController
        function this = PolePlacementPredictiveController(A, B, polesCont, sequenceLength, defaultSamplingInterval, setpoint)            
            % Class constructor.
            %
            % Parameters:
            %  >> A (Square Matrix)
            %     The system matrix of the plant.
            %
            %  >> B (Matrix)
            %     The input matrix of the plant.
            %
            %  >> polesCont (Vector)
            %     The vector containing the desired locations, given in terms of n self-conjugate (complex) numbers of the
            %     closed-loop poles in the s-domain, where n is the state dimension (i.e., the number of rows in A).
            %
            %  >> defaultSamplingInterval (Positive scalar)
            %     The sampling interval to be used for the computation of
            %     the equivalent pole locations in the z-domain.
            %     The actual controller gain is then computed w.r.t. these
            %     poles.
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
            %   << obj (PolePlacementPredictiveController)
            %      A new PolePlacementPredictiveController instance.
            %
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedTrackingController(dimX, dimU, sequenceLength, true, [], [], []);
            this.A = A;
            this.B = B;
            
            % validate the given poles: check if complex poles occur as
            % conjugate pairs
            assert(Checks.isVec(polesCont, dimX), ...
                'PolePlacementPredictiveController:InvalidPoles', ...
                '** <polesCont> is expected to be a vector containing %d self-conjugate elements **', dimX);
            try
                this.polesCont = cplxpair(polesCont(:));
            catch ex
                if strcmp(ex.identifier, 'MATLAB:cplxpair:ComplexValuesPaired')
                    error('PolePlacementPredictiveController:InvalidPoles', ...
                        '** Complex-valued entries in <polesCont> must be self-conjugate pairs **');
                else
                    ex.rethrow();
                end
            end
            % also check the multiplicity of the poles, must not exceed
            % rank(B) (code adapted from Matlab's place funtion)            
            assert(all(diff([0;find(diff(sort(this.polesCont))~=0);dimX]) <= rank(B)), ...
                'PolePlacementPredictiveController:InvalidPoles', ...
                '** <polesCont> must not contain poles with multiplicity > %d=rank(B) **', rank(B));                
            
            if any(real(this.polesCont) >= 0)
                warning('PolePlacementPredictiveController:Unstable', ...
                '** <polesCont> contains at least one unstable pole **'); 
            end
            
            this.samplingInterval = defaultSamplingInterval;
            % compute equivalent poles in z-domain (matched z-transform) and the resulting
            % feedback gain
            this.poles = exp(this.polesCont .* this.samplingInterval);
            this.computeAndSetGainMatrix();
            
            if nargin > 5 && ~isempty(setpoint)
                assert(Checks.isVec(setpoint, this.dimPlantState) && all(isfinite(setpoint)), ...
                    'PolePlacementPredictiveController:InvalidSetpoint', ...
                    '** <setpoint> is expected to be a real-valued, %d-dimensional vector', dimX);
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
        
        %% changeModelParameters
        function changeModelParameters(this, newA, newB, ~)
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
                'PolePlacementPredictiveController:ChangeSetPoint:InvalidSetpoint', ...
                '** <newSetpoint> is expected to be a real-valued, %d-dimensional vector', this.dimPlantState);

            this.setpoint = newSetpoint(:);
            this.recomputeFeedforward = true;
        end
        
        %% changeSamplingInterval
        function changeSamplingInterval(this, newSamplingInterval)
            this.samplingInterval = newSamplingInterval;
            % update  poles in z-domain
            this.poles = exp(this.polesCont .* this.samplingInterval);
            
            % requires that gain matrix is recomputed
            this.recomputeGain = true;            
            % and, likewise, the feedforward, if nonzero set point is
            % present
            if ~isempty(this.setpoint)
                this.recomputeFeedforward = true;
            end
        end
        
    end
    
    methods (Access = private)
        %% computeAndSetGainMatrix
        function computeAndSetGainMatrix(this)
            % compute the gain matrix, based on pole placement   
            this.L = place(this.A, this.B, this.poles);
        end
        
        %% computeAndSetFeedforward
        function computeAndSetFeedforward(this)
            I = speye(this.dimPlantState);
            % we need a stationary gain
            F = this.B \ ((I - this.A + this.B * this.L) / I); % this is the inverse of the DC gain of the closed loop system, with C = I
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
                inputSequence((j - 1) * this.dimPlantInput + 1: j * this.dimPlantInput) = -this.L * state + ff;
                % predict the system state (i.e., compute estimate since noise is
                % zero-mean)
                state = this.A * state + ...
                    this.B * inputSequence((j - 1) * this.dimPlantInput + 1: j * this.dimPlantInput);
            end
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, ~)
%             if isempty(this.setpoint)
%                 performance = state;
%             else
%                 performance = state - this.setpoint;
%             end
%             stageCosts = Utility.computeStageCosts(performance, input, this.Q, this.R);
            stageCosts = 0;
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
            %averageLQGCosts = Utility.computeLQGCosts(horizonLength, performance, appliedInputs, this.Q, this.R) / horizonLength;
            averageLQGCosts = 0;
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

