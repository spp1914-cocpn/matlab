classdef NominalPredictiveController< SequenceBasedController
    % Implementation of a linear networked predictive controller based on a
    % nominal linear-quadratic regulator.
    %
    % Literature: 
    %   Guo-Ping Liu,
    %   Predictive Controller Design of Networked Systems With Communication Delays and Data Loss,
    %   IEEE Transactions on Circuits and Systems—II: Express Briefs,
    %   Vol. 57, No. 6, pp. 481-485, 2010.
    %
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017  Florian Rosenthal <florian.rosenthal@kit.edu>
    %
    %                        Institute for Anthropomatics and Robotics
    %                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
    %                        Karlsruhe Institute of Technology (KIT), Germany
    %
    %                        http://isas.uka.de
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
    
    properties (SetAccess = immutable)
        A;
        B;
        % system state cost matrix;
        % (positive semidefinite matrix of dimension <dimX> x <dimX>)
        Q = [];
        % input cost matrix;
        % (positive definite matrix of dimension <dimU> x <dimU>)
        R = [];
        % controller gain matrix (negative)
        L = [];
    end
    
    methods (Access = public)
        %% NominalPredictiveController
        function this = NominalPredictiveController(A, B, Q, R, sequenceLength)
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
            % Returns:
            %   << obj (NominalPredictiveController)
            %      A new NominalPredictiveController instance.
            %
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength);
            this.A = A;
            this.B = B;
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
            
            % compute the gain matrix 
            [L, ~, ~] = dlqr(A, B, Q, R);
            this.L = -L;
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
    end
    
    methods (Access = protected)
         %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, plantState, varargin)
            [state, ~] = plantState.getMeanAndCovariance();

            inputSequence = zeros(this.sequenceLength * this.dimPlantInput, 1);
            for j = 1:this.sequenceLength
                inputSequence((j - 1) * this.dimPlantInput + 1: j * this.dimPlantInput) = this.L * state;
                % predict the system state (i.e., compute estimate since noise is
                % zero-mean)
                state = this.A * state + ...
                    this.B * inputSequence((j - 1) * this.dimPlantInput + 1: j * this.dimPlantInput);
             end
        end
        
        %% doCostsComputation
        function averageLQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonLength = size(appliedInputs, 2);
            if size(stateTrajectory, 2) ~= horizonLength + 1
                error('NominalPredictiveController:DoCostsComputation', ...
                    '** <stateTrajectory> is expected to have %d columns ', horizonLength + 1);
                    
            end
            averageLQGCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, appliedInputs, this.Q, this.R) / horizonLength;
        end
    end
end

