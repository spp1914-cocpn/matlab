classdef IMMBasedRecedingHorizonController < SequenceBasedController
    % Implementation of a receding horizon linear sequence-based controller for
    % NCS with networks connecting the controller and actuator, 
    % and sensor and controller, respectively, where application layer ACKs
    % are sent out by the actuator upon reception of applicable control
    % inputs.
    %
    % Literature: 
    %   Florian Rosenthal and Uwe D. Hanebeck,
    %   Sequence-Based Stochastic Receding Horizon Control Using IMM Filtering and Value Function Approximation (submitted),
    %   Proceedings of the IEEE 58th Conference on Decision and Control (CDC 2019),
    %   Nice, France, December 2019.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (SetAccess = immutable, GetAccess = private)
        Q;
        R;
        horizonLength;
                
        numModes;
        transitionMatrix;
        F;
        G;
        
        augA;
        augB;
        R_tilde; % store only the one for the first mode
        Q_tilde; % zero for the first and the last mode
        P; % cost matrices, one for each mode
        
        immf;
        dimEta;
        measurementModel;
        mjls;
    end
    
    properties (Access = private)
        etaState; % the input related part of the state, evolves according to %eta_k+1=F*eta_k + G*U_k
        bufferedInputSequences; % corresponding to the modes
        initialized = true; % flag to indicate whether intial state was set (no update required prior to computation of sequence)
    end
    
    methods (Access = public)
        %% IMMBasedRecedingHorizonController
        function this = IMMBasedRecedingHorizonController(A, B, C, Q, R, caDelayProb, ...
                sequenceLength, maxMeasDelay, W, V, horizonLength, x0, x0Cov)
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength);
            
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
            
            Validator.validateDiscreteProbabilityDistribution(caDelayProb);            
            Validator.validateMeasurementMatrix(C, dimX);
            dimMeas = size(C, 1);
            assert(Checks.isNonNegativeScalar(maxMeasDelay) && mod(maxMeasDelay, 1) == 0, ...
                'IMMBasedRecedingHorizonController:InvalidMaxMeasDelay', ...
                ['** Input parameter <maxMeasDelay> (maximum measurement',...
                'delay (M)) must be a nonnegative integer **']);             
            
            % check the noise covs
            [W_sqrt] = Validator.validateSysNoiseCovarianceMatrix(W, dimX);
            Validator.validateMeasNoiseCovarianceMatrix(V, dimMeas);       
            
            Validator.validateHorizonLength(horizonLength);
            this.horizonLength = horizonLength;
            % check the states
            assert(Checks.isVec(x0, dimX) && all(isfinite(x0)), ...
                'IMMBasedRecedingHorizonController:InvalidX0', ...
                '** Input parameter <x0> (initial estimate of plant state) must be a %d-dimensional vector **', dimX);
            assert(Checks.isCov(x0Cov, dimX), ...
                'IMMBasedRecedingHorizonController:InvalidX0Cov', ...
                '** Input parameter <x0Cov> (initial covariance) must be positive definite and %d-by%d-dimensional **', dimX, dimX);
            
            this.numModes = sequenceLength + 1;
            this.transitionMatrix = Utility.calculateDelayTransitionMatrix(...
                Utility.truncateDiscreteProbabilityDistribution(caDelayProb, this.numModes));
            
            this.dimEta = dimU * (sequenceLength * (sequenceLength - 1) / 2);
            blkdiag(W, zeros(this.dimEta));
            [dimAugState, this.F, this.G, this.augA, this.augB, this.Q_tilde, augR] = ...
                Utility.performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R);
            % augmented state consists of x_k and eta_k
            % augR (R_tilde in the paper) is zero for all modes but the first
            this.R_tilde = augR(:, :, 1);                

            modeFilters = arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, 'UniformOutput', false);
            this.immf = DelayedModeIMMF(modeFilters, this.transitionMatrix, maxMeasDelay);
            this.measurementModel = LinearMeasurementModel(C);
            this.measurementModel.setNoise(Gaussian(zeros(dimMeas, 1), V));
            
            this.mjls = JumpLinearSystemModel(this.numModes, ...
                arrayfun(@(~) LinearPlant(A, B, W), 1:this.numModes, 'UniformOutput', false));

            this.bufferedInputSequences = zeros(dimU, sequenceLength, sequenceLength);

            this.setControllerPlantStateInternal(x0, x0Cov);
            this.P = this.computeCostate(dimAugState);
        end
        
        %% setControllerPlantState
        function setControllerPlantState(this, state)
            % Set the estimate of the plant state.
            %
            % This function is mainly used to set an initial estimate as
            % the controller's dynamics is used for updating the estimate of
            % the plant state.
            %
            % Parameters:
            %   >> state (Subclass of Distribution)
            %      The new system state.
            
            assert(Checks.isClass(state, 'Distribution') && state.getDim() == this.dimPlantState, ...
                'IMMBasedRecedingHorizonController:SetControllerPlantState:InvalidState', ...
                '** <state> must be a %d-dimensional Distribution **', this.dimPlantState);
            
            [x0, x0Cov] = state.getMeanAndCov();
            
            this.setControllerPlantStateInternal(x0, x0Cov);  
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
            [plantState, ~] = this.immf.getStateMeanAndCov();
        end
        
        %% computeControlSequence
        function inputSequence = computeControlSequence(this, measurements, measurementDelays, modeMeas, modeDelays)
            if nargin == 1 || isempty(measurements)
                measurements = [];
                measurementDelays = [];
            end
            if this.initialized
                % at the first time step, no update of the filter required
                this.initialized = false;                
            else
                % update the filter state first, i.e., obtain the posterior for
                % the current time step
                this.mjls.setSystemInput([cell2mat(arrayfun(@(mode) this.bufferedInputSequences(:, mode, mode), ...
                        1:this.sequenceLength, 'UniformOutput', false)) zeros(this.dimPlantInput, 1)]);
                this.immf.step(this.mjls, this.measurementModel, ...
                        measurements, measurementDelays, modeMeas, modeDelays);
            end
            % now we can compute the input sequence
            inputSequence = this.doControlSequenceComputation();
            % predict eta_{k+1} and buffer the freshly created sequence
            this.etaState = this.F * this.etaState + this.G * inputSequence;
            
            this.bufferedInputSequences = circshift(this.bufferedInputSequences, 1, 3);            
            this.bufferedInputSequences(:,:, 1) = inputSequence;  
        end
        
         %% getLastComputationMeasurementData
        function [numUsedMeas, numDiscardedMeas] = getLastComputationMeasurementData(this)
            % Get information on the last performed control sequence computation (due to a call of computeControlSequence).
            %
            % Returns:
            %   << numUsedMeas (Nonnegative integer)
            %      The number of measurements used by the last measurement update.
            %
            %   << numDiscardedMeas (Nonnegative integer)
            %      The number of measurements discarded during the last
            %      measurement update due to their delays.
            
            [numUsedMeas, numDiscardedMeas] = this.immf.getLastUpdateMeasurementData();            
        end
        
        %% reset
        function reset(this)
        
        end
    end
    
    methods (Access = protected)
        %% doControlSequenceComputation
         function inputSequence = doControlSequenceComputation(this)
            % extract the components of the mixture (posterior state)
            [means, ~, probs] = this.immf.getState().getComponents();
                        
            weights = this.transitionMatrix .* probs(:);
            stateMeans = vertcat(reshape(means, this.dimPlantState, 1, this.numModes), ...
                repmat(this.etaState, 1, 1, this.numModes));
            statePart = mtimesx(this.augA, stateMeans);
            % augR is only nonzero for first mode
            part1 = probs(1) * this.R_tilde;
            part2 = zeros(this.dimPlantInput * this.sequenceLength, 1);
            
            for r=1:this.numModes
                BP = mtimesx(this.augB, 'T', this.P(:, :, r));
                part1 = part1 + sum(reshape(weights(:, r), 1, 1, this.numModes) .* mtimesx(BP, this.augB), 3);
                part2 = part2 + sum(reshape(weights(:, r), 1, 1, this.numModes) .* mtimesx(BP, statePart), 3);
            end
            inputSequence = -pinv(part1) * part2;
         end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, ~)                        
            stageCosts = Utility.computeStageCosts(state, input, this.Q, this.R);
        end
        
        %% doCostsComputation
         function lQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonlength = size(appliedInputs, 2);
            assert(size(stateTrajectory, 2) == horizonlength + 1, ...
                'IMMBasedRecedingHorizonController:DoCostsComputation', ...
                '** <stateTrajectory> is expected to have %d columns ', horizonlength + 1);
            
            lQGCosts = Utility.computeLQGCosts(horizonlength, stateTrajectory, appliedInputs, this.Q, this.R);
        end
    end
    
    methods (Access = private)        
        %% computeCostate
        function P_1 = computeCostate(this, dimAugState)
            dimInputSeq = this.dimPlantInput * this.sequenceLength;
            terminalP = repmat(blkdiag(this.Q, zeros(this.dimEta)), 1, 1, this.numModes); % the terminal term, Y_K
            % backward iteration over the horizon
            P_k = terminalP;
            for k=this.horizonLength-1:-1:1
                % compute Y_k+1 using P_k+1
                Y = reshape(cell2mat(arrayfun(@(i) sum(reshape(this.transitionMatrix(i,:), 1, 1, this.numModes) .* P_k,3), ...
                    1:this.numModes, 'UniformOutput', false)), dimAugState, dimAugState, this.numModes);                
               
                %inner term
                BY = mtimesx(this.augB, 'T', Y);
                M = mtimesx(BY, this.augB);
                % this is the expression for M in the paper
                M(:, :, 1) = M(:, :, 1) + this.R_tilde; % augR is only nonzero for the first mode
                M = reshape(cell2mat(arrayfun(@(i) pinv(M(:, :, i)), 1:this.numModes, 'UniformOutput', false)), ...
                    dimInputSeq, dimInputSeq, this.numModes); 
                M = Y - mtimesx(mtimesx(BY, 'T', M), BY);
                P_k = this.Q_tilde + mtimesx(mtimesx(this.augA, 'T', M), this.augA); %P_k
                P_k = (P_k + permute(P_k, [2 1 3])) / 2; % ensure symmetry
            end
            P_1 = P_k; % as desired, P_{k+1} for all modes
        end
        
        %% setControllerPlantStateInternal
        function setControllerPlantStateInternal(this, mean, cov)
            % no checks at all
            % we assume the buffer is initially empty
            this.etaState = zeros(this.dimEta, 1);
            % consequently, the last mode of the MJLS is active
            this.immf.setState(GaussianMixture(repmat(mean(:), 1, this.numModes), ...
                repmat(cov, 1, 1, this.numModes), [zeros(1, this.numModes - 1) 1]));
        end
    end
end

