classdef ExpectedInputPredictiveController< NominalPredictiveController
    % Implementation of a linear networked predictive controller based on a
    % nominal linear-quadratic regulator that uses expected inputs for the
    % open-loop prediction of the state trajectory.
    %
    % Literature: 
    %   Achim Hekler, Jörg Fischer, and Uwe D. Hanebeck,
    %   Control over Unreliable Networks Based on Control Input Densities,
    %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
    %   Singapore, Singapore, July 2012.
    %
    %   Achim Hekler, Jörg Fischer, and Uwe D. Hanebeck,
    %   Sequence-Based Control for Networked Control Systems Based on Virtual Control Inputs,
    %   Proceedings of the IEEE 51st Conference on Decision and Control (CDC 2012),
    %   Maui, Hawaii, USA, December 2012.
    
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
    
    properties (GetAccess = private, SetAccess = immutable)       
        F; % F matrix and its powers F^2, F^3, ..., F^N, with N = sequenceLength -2
        % F^2, F^3, ..., F^N are needed to extract inputs from past
        % sequences, stored in eta_k, that are applicable at stage i of the
        % predction horizon, i.e., at time k+i
        G;
        H;
    end
    
    properties (Access = private)
        caDelayProbs;
        etaState; % eta_k stores inputs from past sequences that are applicable at time k or later, evolves according to eta_k=F*eta_k-1+G*U_k-1
        alphas;
    end
    
    methods (Access = public)
        %% ExpectedInputPredictiveController
        function this = ExpectedInputPredictiveController(A, B, Q, R, sequenceLength, caDelayProbs, setpoint)
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
            %   >> caDelayProbs (Nonnegative vector with elements suming to 1)
            %      The vector describing the delay distribution of the
            %      CA-network.
            %
            %  >> setpoint (Vector, optional)
            %     The set point in the state space that shall be tracked asymptotically instead of the origin.
            %     Note that it is not verified whether the resulting set point tracking problem is solvable. 
            %     For instance, in case dim(x) > dim(u), the problem will generally have no exact solution.
            %     If no set point is provided, the controller will attempt
            %     to drive the plant to the origin instead.
            %
            % Returns:
            %   << obj (ExpectedInputPredictiveController)
            %      A new ExpectedInputPredictiveController instance.
            %
            if nargin < 7
                setpoint = [];
            end
            this = this@NominalPredictiveController(A, B, Q, R, sequenceLength, setpoint);
            Validator.validateDiscreteProbabilityDistribution(caDelayProbs);
            this.caDelayProbs = Utility.truncateDiscreteProbabilityDistribution(caDelayProbs, sequenceLength + 1);
            
            dimEta = this.dimPlantInput * (sequenceLength * (sequenceLength - 1) / 2);
            this.F = zeros(dimEta, dimEta, sequenceLength - 2);
            % we need F, G, H
            [this.F(:, :, 1), this.G, this.H, ~] = Utility.createActuatorMatrices(sequenceLength, this.dimPlantInput);
            
            % compute the required powers of F
            % we use them to compute the expected input u_k+i|k for all
            % stages i of the prediction horizon to obtain estimate
            % x_k+i+1|k
            for i=2:sequenceLength-2
                this.F(:, :, i)  = this.F(:, :, i-1) * this.F(:, :, 1);
            end
            % initially, no old inputs are present
            this.etaState = zeros(dimEta, 1);
            this.computeAndSetInputWeights();
        end
        
        %% reset
        function reset(~)
            % do nothing
        end
        
        %% changeCaDelayProbs
        function changeCaDelayProbs(this, newCaDelayProbs)
            % Change the distribution of the control packet delays to be assumed by the controller.
            %
            % Parameters:
            %  >> newCaDelayProbs (Nonnegative vector)
            %     Vector specifiying the new delay distribution.
            %
            this.caDelayProbs = Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.sequenceLength + 1);
            this.computeAndSetInputWeights();
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
            error('ExpectedInputPredictiveController:ChangeSequenceLength:NotSupported', ...
                '** This operation is not supported **');
        end        
        
    end
    
    methods (Access = private)
        %% computeAndSetInputWeights
        function computeAndSetInputWeights(this)
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
            
            % use the transition probabilities for convenience, not needed
            % a straightforward approach directly computes the alphas based
            % on the delay probabilities
            P = Utility.calculateDelayTransitionMatrix(this.caDelayProbs);
            sums = cumsum(this.caDelayProbs);
            q = cumprod(1 - sums);
            this.alphas = cell(1, this.sequenceLength);            
            this.alphas{1} = zeros(this.sequenceLength + 1, 1);
            
            this.alphas{1}(1) = this.caDelayProbs(1); % this is alpha_0,0=p0
            for j=2:this.sequenceLength + 1
                this.alphas{1}(j) = q(j-1) * sums(j); % alpha_j,0, e.g. alpha_1,0=(1-p0)*(p0+p1) and alpha_2,0=(1-p0)(1-p0-p1)*(p0+p1+p2)
            end
            
            for k=1:this.sequenceLength - 1
                this.alphas{k+1} = (P')^k*this.alphas{1};
                this.alphas{k+1} = this.alphas{k+1}(k+1:end) / sum(this.alphas{k+1}(k+1:end));
            end
        end
    end
    
    methods (Access = protected)        
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, plantState, varargin)
            % plant state is an estimate of the true state at time k: x_k|k
            inputSequence = doControlSequenceComputation@NominalPredictiveController(this, plantState, varargin);
            %eta_k+1=F*eta_k+G*U_k
            this.etaState = this.F(:, :, 1) * this.etaState + this.G * inputSequence;
        end
        
        %% doStatePrediction
        function newState = doStatePrediction(this, state, computedInput, stage)
            % check the stage of the horizon
            % stage starts at 1!
            % stage j -> predict stage for stage j+1 using computed input
            % u_k+j-1
            % compute the expected input first based on the input u_k+j-1|k currently computed
            % open-loop prediction of inputs
            % assumes that default input is zero
            % expected input is simply the mean of a Dirac Mixture (virtual control input)             
            
            if stage == this.sequenceLength % maximum stage
                % either the compute input or the zero default input
                expectedInput = computedInput * this.alphas{this.sequenceLength}(1);
            else
                if stage == 1
                    eta = this.etaState;
                else
                    eta = this.F(:, :, stage-1) * this.etaState;
                end
                expectedInput = computedInput * this.alphas{stage}(1) ...
                        + sum(mtimesx(this.H(:, :, stage+1:end-1), eta) .* reshape(this.alphas{stage}(2:end-1), 1, 1, []), 3);
            end
            % predict the new state for stage j+1 (x_k+j+1|k), based on estimate of current state (x_k+j|k) and
            % expected input u_k+i|k
            newState = this.A * state + this.B * expectedInput;
        end        
    end
end

