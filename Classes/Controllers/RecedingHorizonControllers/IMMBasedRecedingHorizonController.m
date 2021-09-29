classdef IMMBasedRecedingHorizonController < SequenceBasedController & ModelParamsChangeable & CaDelayProbsChangeable
    % Implementation of a receding horizon linear sequence-based controller for
    % linear NCS with networks connecting the controller and actuator, 
    % and sensor and controller, respectively, where application layer ACKs
    % are sent out by the actuator upon reception of applicable control
    % inputs.
    %
    % Literature: 
    %   Florian Rosenthal and Uwe D. Hanebeck,
    %   Sequence-Based Stochastic Receding Horizon Control Using IMM Filtering and Value Function Approximation,
    %   Proceedings of the IEEE 58th Conference on Decision and Control (CDC 2019),
    %   Nice, France, December 2019.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        horizonLength double {Validator.validateHorizonLength(horizonLength)} = 1;
                
        numModes;             
    end
    
    properties (SetAccess = immutable, GetAccess = ?IMMBasedRecedingHorizonControllerTest)
        immf;
        dimEta;
        measurementModel;
        mjls;
        % for fixed sequence length, F and G are independent of (A,B)
        F;
        G;        
    end
    
    properties (SetAccess = private, GetAccess = ?IMMBasedRecedingHorizonControllerTest)
        % the input related part of the state, evolves according to %eta_k+1=F*eta_k + G*U_k
        etaState; % holds eta_{k-1}
        % the last computed input sequence
        inputSequence;
    end
    
    properties (GetAccess = private, SetAccess = ?IMMBasedRecedingHorizonControllerTest)
        initialized = true; % flag to indicate whether intial state was set (no update required prior to computation of sequence)
    end
    
    properties (Access = private)
        Q;
        terminalQ; % use unique stabilizing solution of DARE 
        R;
        augA;
        augB;
        R_tilde; % store only the one for the first mode
        Q_tilde;
                
        transitionMatrix;
        P; % cost matrices, one for each mode
        recomputeCostate = false;
    end
    
    properties (SetAccess = immutable, GetAccess = public)        
        useMexImplementation(1,1) logical = true; 
        % by default, we use the C++ (mex) implementation for computation
        % of the costate/cost matrices
        % this is usually faster, but can produce slightly different results
    end
    
    methods (Access = public)
        %% IMMBasedRecedingHorizonController
        function this = IMMBasedRecedingHorizonController(A, B, C, Q, R, modeTransitionMatrix, ...
                sequenceLength, maxMeasDelay, W, measNoise, horizonLength, x0, x0Cov, useMexImplementation)
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
            %   >> modeTransitionMatrix (Stochastic matrix, i.e. a square matrix with nonnegative entries whose rows sum to 1)
            %      The transition matrix of the mode theta_k of the augmented dynamics.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %   >> maxMeasDelay (Nonnegative integer)
            %      The maximum delay a measurement may experience before it
            %      is discarded by the IMM filter.
            %
            %   >> W (Square matrix)
            %      The covariance matrix of the plant noise.
            %
            %   >> measNoise (Distribution)
            %      A distribution denoting the time-invariant measurement noise process v_k.
            %
            %   >> horizonLength (Positive integer)
            %      The horizon length considered by the controller for optimization.
            %
            %   >> x0 (Vector)
            %      The initial estimate of the plant state.
            %
            %   >> x0Cov (Positive definite matrix)
            %      The covariance matrix denoting the uncertainty of the IMM filter's initial estimate.
            %
            %   >> useMexImplementation (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for the computation of the controller
            %      costate matrices which is typically considerably faster than the Matlab
            %      implementation. 
            %      If left out, the default value true is used.
            %
            % Returns:
            %   << this (IMMBasedRecedingHorizonController)
            %      A new IMMBasedRecedingHorizonController instance.
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength, false);
            
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
            
            this.numModes = sequenceLength + 1;
            Validator.validateTransitionMatrix(modeTransitionMatrix, this.numModes);
            this.transitionMatrix = modeTransitionMatrix;
            
            Validator.validateMeasurementMatrix(C, dimX);
            dimMeas = size(C, 1);
            assert(Checks.isNonNegativeScalar(maxMeasDelay) && mod(maxMeasDelay, 1) == 0, ...
                'IMMBasedRecedingHorizonController:InvalidMaxMeasDelay', ...
                ['** Input parameter <maxMeasDelay> (maximum measurement',...
                'delay (M)) must be a nonnegative integer **']);             
            
            % check the noise
            Validator.validateSysNoiseCovarianceMatrix(W, dimX);
            assert(Checks.isClass(measNoise, 'Distribution') && measNoise.getDim() == dimMeas, ...
                    'IMMBasedRecedingHorizonController:InvalidMeasNoise', ...
                    '** Input parameter <measNoise> (measurement noise) must be a %d-dimensional Distribution**', dimMeas);
            
            this.horizonLength = horizonLength;
            % check the states
            assert(Checks.isVec(x0, dimX) && all(isfinite(x0)), ...
                'IMMBasedRecedingHorizonController:InvalidX0', ...
                '** Input parameter <x0> (initial estimate of plant state) must be a %d-dimensional vector **', dimX);
            assert(Checks.isCov(x0Cov, dimX), ...
                'IMMBasedRecedingHorizonController:InvalidX0Cov', ...
                '** Input parameter <x0Cov> (initial covariance) must be positive definite and %d-by-%d-dimensional **', dimX, dimX);
            
            if nargin > 13                
                this.useMexImplementation = useMexImplementation;
            end     
            
            this.dimEta = dimU * (sequenceLength * (sequenceLength - 1) / 2);
            
            [dimAugState, this.F, this.G, this.augA, this.augB, this.Q_tilde, augR] = ...
                Utility.performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R);
            % augmented state consists of x_k and eta_k
            % augR (R_tilde in the paper) is zero for all modes but the first
            this.R_tilde = augR(:, :, 1);                

            modeFilters = arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, 'UniformOutput', false);
            this.immf = DelayedModeIMMF(modeFilters, this.transitionMatrix, maxMeasDelay);
            this.measurementModel = LinearMeasurementModel(C);
            % moment matching to model measurement noise as Gaussian
            [v_mean, V] = measNoise.getMeanAndCov();
            this.measurementModel.setNoise(Gaussian(v_mean, V));
            
            this.mjls = JumpLinearSystemModel(this.numModes, ...
                arrayfun(@(~) LinearPlant(A, B, W), 1:this.numModes, 'UniformOutput', false));
            
            this.inputSequence = zeros(dimU * this.sequenceLength, 1);
            
            this.setControllerPlantStateInternal(x0, x0Cov);
            
            % for stability: use stabilizing solution of DARE as terminal
            % weighting -> only exists in case (A,B) stabilizable + the
            % associated symplectic matrix has no eigenvalues on the unit
            % circle     
            [this.terminalQ,~,~,info] = idare(A, B, Q, R);
            if info.Report > 1
                % either report == 2: The solution is not finite
                % or report == 3: No solution found since the symplectic
                % matrix has eigenvalues on the unit circle
                this.terminalQ = this.Q;
                warning('IMMBasedRecedingHorizonController:InvalidPlant', ...
                    '** (A,B) seems not stabilizable, cannot use stabilizing solution of associated DARE as terminal weighting Q_N **');
            end            
            % computation of costate has to be done only once if model
            % parameters or cost function do not change
            this.P = this.computeCostate();
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
                Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.numModes));
            if ~isequal(newMat, this.transitionMatrix)
                this.transitionMatrix = newMat;
                % also update the matrix of the filter
                this.immf.setModeTransitionMatrix(this.transitionMatrix); 
                % and we must remember recompute the costate, depends on
                % transition properties
                this.recomputeCostate = true;
            end
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

            % for stability: use stabilizing solution of DARE as terminal
            % weighting -> only exists in case (A,B) stabilizable + the
            % associated symplectic matrix has no eigenvalues on the unit
            % circle     
            [this.terminalQ,~,~,info] = idare(this.augA(1:this.dimPlantState, 1:this.dimPlantState, 1), ...
                this.augB(1:this.dimPlantState, 1:this.dimPlantInput, 1), newQ, newR);
            if info.Report > 1
                % either report == 2: The solution is not finite
                % or report == 3: No solution found since the symplectic
                % matrix has eigenvalues on the unit circle
                this.terminalQ = newQ;                
            end  
            
            % requires that cost function augmentation is carried out again
            [this.Q_tilde, augR] = Utility.createAugmentedCostModel(this.sequenceLength, this.Q, this.R);
            % augR (R_tilde in the paper) is zero for all modes but the first
            this.R_tilde = augR(:, :, 1);
            % remember to recompute the costate, depends on cost matrices
            this.recomputeCostate = true;
        end
        
        %% changeModelParameters
        function changeModelParameters(this, newA, newB, newW)
            Validator.validateSystemMatrix(newA, this.dimPlantState);
            Validator.validateInputMatrix(newB, this.dimPlantState, this.dimPlantInput);
            Validator.validateSysNoiseCovarianceMatrix(newW, this.dimPlantState);

            [~, ~, H, J] = Utility.createActuatorMatrices(this.sequenceLength, this.dimPlantInput);
            % in the the augmented model, F, G, H, J do not change
            % so only adapt the "first block row" in augA and augB            
            this.augA(1:this.dimPlantState, 1:this.dimPlantState, 1) = newA;
            for i=2:this.sequenceLength
                % update B*H for all but first and last mode
                % H is empty for first and last mode
                this.augA(1:this.dimPlantState, :, i) = [newA, newB * H(:, :, i)];                
            end
            this.augA(1:this.dimPlantState, 1:this.dimPlantState, this.sequenceLength + 1) = newA;            
            % augB is only different for first mode (multiply B by J)
            % B*J is zero for all modes but first
            %
            this.augB(1:this.dimPlantState, :, 1) = newB * J(:, :, 1);             
            
            % for stability: use stabilizing solution of DARE as terminal
            % weighting -> only exists in case (A,B) stabilizable + the
            % associated symplectic matrix has no eigenvalues on the unit
            % circle
            [this.terminalQ,~,~,info] = idare(newA, newB, this.Q, this.R);  
            if info.Report > 1
                % either report == 2: The solution is not finite
                % or report == 3: No solution found since the symplectic
                % matrix has eigenvalues on the unit circle
                this.terminalQ = this.Q;
                warning('IMMBasedRecedingHorizonController:InvalidPlant', ...
                    '** (A,B) seems not stabilizable, cannot use stabilizing solution of associated DARE as terminal weighting Q_N **');
            end            

            for j=1:this.sequenceLength + 1
                % change the affected parameters of the model
                this.mjls.setSystemMatrixForMode(newA, j);
                this.mjls.setSystemInputMatrixForMode(newB, j);
                this.mjls.setSystemNoiseCovarianceMatrixForMode(newW, j);
            end
            
            % remember to recompute the costate, depends on plant
            % model/augmented model
            this.recomputeCostate = true;
        end
        
         %% setEtaState
        function setEtaState(this, newEta, newInputSeq)
            assert(Checks.isVec(newEta, this.dimEta) && all(isfinite(newEta)), ...
                'IMMBasedRecedingHorizonController:SetEtaState:InvalidEta', ...
                '** <newEta> must be a %d-dimensional vector **', ...
                this.dimEta);
            assert(Checks.isVec(newInputSeq, this.dimPlantInput * this.sequenceLength) && all(isfinite(newInputSeq)), ...
                'IMMBasedRecedingHorizonController:SetEtaState:InvalidInputSeq', ...
                '** <newInputSeq> must be a %d-dimensional vector **', ...
                this.dimPlantInput * this.sequenceLength);
            
            this.etaState = newEta(:);
            this.inputSequence = newInputSeq(:);             
        end
        
        %% setControllerPlantState
        function setControllerPlantState(this, state)
            % Set the estimate of the plant state.
            %
            % This function is mainly used to set an initial estimate as
            % the associated IMM filter is used for updating the estimate of
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
        function inputSeq = computeControlSequence(this, measurements, measurementDelays, modeMeas, modeDelays)
            if nargin == 1
                measurements = [];
                measurementDelays = [];
                modeMeas = [];
                modeDelays = [];
            elseif nargin == 3
                modeMeas = [];
                modeDelays = [];
            elseif nargin ~= 5
                error('IMMBasedRecedingHorizonController:ComputeControlSequence:InvalidNumArgs', ...
                    '** Number of arguments must be 1, 3 or 5 **');             
            end
            
            dimMeas = size(this.measurementModel.measMatrix, 1);            
            assert(isempty(measurements) || Checks.isFixedRowMat(measurements, dimMeas), ...
                    'IMMBasedRecedingHorizonController:ComputeControlSequence:InvalidMeas', ...
                    '** Individual measurements must be %d-dimensional **', dimMeas);
            % validity of given measurement delays is checked by filter
            % likewise, validity of mode observations and corresponding
            % delays
                
            if this.initialized
                % at the first time step, no update of the filter required
                this.initialized = false;                
                newEta = this.etaState;
            else
                % obtain the possible inputs stored in eta_{k-1} and
                % U_{k-1}
                distributedInputs = zeros(this.dimPlantInput, this.numModes);
                % default input zeros for the last mode, so last column
                % remains zero
                % first column is first entry from most recent sequence U_{k-1}
                distributedInputs(:, 1) = this.inputSequence(1:this.dimPlantInput);
                idx = 1;
                for i=2:this.numModes - 1 
                    distributedInputs(:, i) = this.etaState(idx:idx+this.dimPlantInput-1);
                    idx =  idx + (this.sequenceLength - (i - 1)) * this.dimPlantInput;                     
                end                
                % update the filter state first, i.e., obtain the posterior x_k^e for
                % the current time step k
                this.mjls.setSystemInput(distributedInputs);
                this.immf.step(this.mjls, this.measurementModel, measurements, measurementDelays, modeMeas, modeDelays);
                % compute eta_k
                newEta = this.F * this.etaState + this.G * this.inputSequence;
            end
            % now we can compute the input sequence
            inputSeq = this.doControlSequenceComputation(newEta);            
            
            this.etaState = newEta;
            this.inputSequence = inputSeq;            
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
         function inputSequence = doControlSequenceComputation(this, newEta)
            if this.recomputeCostate
                % we need to recompute the costate, at least one parameter
                % (A, B, Q, R, W, transition matrix) has changed
                this.P = this.computeCostate();
                this.recomputeCostate = false;
            end
            % extract the components of the mixture (posterior state)
            [means, ~, probs] = this.immf.getState().getComponents();
            if this.useMexImplementation                
                inputSequence = ...
                    mex_IMMBasedRecedingHorizonController(this.augA, this.augB, this.transitionMatrix, this.R_tilde, ...
                        means, probs, newEta, this.P);
            else
                weights = this.transitionMatrix .* probs(:);
                stateMeans = vertcat(reshape(means, this.dimPlantState, 1, this.numModes), ...
                    repmat(newEta, 1, 1, this.numModes));
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
        function P_1 = computeCostate(this)            
            if this.useMexImplementation            
                P_1 = mex_IMMBasedRecedingHorizonControllerCostate(this.augA, this.augB, this.transitionMatrix, ...
                    this.Q_tilde, this.R_tilde, this.terminalQ, this.horizonLength);
            else
                dimInputSeq = this.dimPlantInput * this.sequenceLength;
                dimAugState = size(this.augA, 2);                
                % backward iteration over the horizon
                P_k = repmat(blkdiag(this.terminalQ, zeros(this.dimEta)), 1, 1, this.numModes); % the terminal term, Y_K
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

