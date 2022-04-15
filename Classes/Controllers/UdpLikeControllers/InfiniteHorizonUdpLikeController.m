classdef InfiniteHorizonUdpLikeController < SequenceBasedController
    % Implementation of an infinite horizon linear sequence-based LQG controller for
    % NCS with a UDP-like networks connecting the controller and actuator, 
    % and sensor and controller, respectively.
    %
    % This implementation is based on the original one by Maxim Dolgov.
    %
    % Literature: 
    %   Maxim Dolgov, Jörg Fischer, and Uwe D. Hanebeck,
    %   Infinite-Horizon Sequence-based Networked Control without Acknowledgments,
    %   Proceedings of the 2015 American Control Conference (ACC 2015),
    %   Chicago, Illinois, USA, July 2015.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Access = protected)
        stationaryModeDistribution;
        
        sysState; % augmented system state
        augmentedMeasurement;
        % availability of current measurement (zero delay) is indicated by
        % first entry of bit vector
        availableMeasurements; % bit vector where 1 indicates that measurement is available
        L; % feedback gain (-L compared to the paper)
        K; % observer gain
    end
    
    properties (GetAccess = protected, SetAccess = immutable)
        % system state cost matrix;
        % (positive semidefinite matrix of dimension <dimX> x <dimX>)
        Q;
        % input cost matrix;
        % (positive definite matrix of dimension <dimU> x <dimU>)
        R; 
        dimState; % dimension of the augmented state
        expAugA;
        expAugB;
        measNoiseCovSqrt; % lower Cholesky factor of the covariance matrix
        measNoiseMean;
        dimMeas;
        measBufferLength;
        dimAugmentedMeas;
        augC;
    end
    
    properties (Access = private)
        lastNumUsedMeas = 0;
        lastNumDiscardedMeas = 0;
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        useMexImplementation(1,1) logical = true; 
        % by default, we use the C++ (mex) implementation for computation of controller gains
        % this is usually faster, but can produce slightly different results
    end
    
    properties (GetAccess = public, Dependent)
        maxMeasurementDelay;
    end
     
    methods
        function maxMeasDelay = get.maxMeasurementDelay(this)
            maxMeasDelay = this.measBufferLength - 1;
        end
    end
    
    methods (Access = public)
        %% InfiniteHorizonUdpLikeController
        function this = InfiniteHorizonUdpLikeController(A, B, C, Q, R, caModeTransitionMatrix, scDelayProb, ...
                sequenceLength, maxMeasDelay, W, V, v_mean, useMexImplementation)
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
            %      The state weighting matrix in the controller's underlying cost function.
            %
            %   >> R (Positive definite matrix)
            %      The input weighting matrix in the controller's underlying cost function.
            %
            %   >> caModeTransitionMatrix (Stochastic matrix, i.e. a square matrix with nonnegative entries whose rows sum to 1)
            %      The transition matrix of the mode theta_k of the augmented dynamics.
            %
            %   >> scDelayProb (Nonnegative vector)
            %      The vector describing the delay distribution of the
            %      SC-network.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %   >> maxMeasDelay (Nonnegative integer)
            %      The maximum delay a measurement may experience before it
            %      is discarded by the controller.
            %
            %   >> W (Square Matrix)
            %      The covariance matrix of the plant noise.
            %
            %   >> V (Square Matrix)
            %      The covariance matrix of the measurement noise.
            %
            %   >> v_mean (Vector)
            %      The mean of the measurement noise.
            %      I the empty matrix is passed, v_mean = 0 is assumed.
            %
            %   >> useMexImplementation (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for the computation of the controller
            %      gains which is usually considerably faster than the Matlab
            %      implementation. 
            %      If left out, the default value true is used.
            %
            % Returns:
            %   << this (InfiniteHorizonUdpLikeController)
            %      A new InfiniteHorizonUdpLikeController instance.
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength, false);
            
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
            
            Validator.validateMeasurementMatrix(C, dimX);
            this.dimMeas = size(C, 1);
            assert(Checks.isNonNegativeScalar(maxMeasDelay) && mod(maxMeasDelay, 1) == 0, ...
                'InfiniteHorizonUdpLikeController:InvalidMaxMeasDelay', ...
                ['** Input parameter <maxMeasDelay> (maximum measurement',...
                'delay (M-1)) must be a nonnegative integer **']);            
                        
            this.measBufferLength = maxMeasDelay + 1; % maxMeasDelay + 1 is buffer length
            
            this.dimAugmentedMeas = this.dimMeas * this.measBufferLength;
            % check the noise covs
            Validator.validateSysNoiseCovarianceMatrix(W, dimX);
            this.measNoiseCovSqrt = Validator.validateMeasNoiseCovarianceMatrix(V, this.dimMeas);
            if  isempty(v_mean)
                this.measNoiseMean = zeros(this.dimMeas, 1);
            else
                assert(Checks.isVec(v_mean, this.dimMeas), ...
                    'InfiniteHorizonUdpLikeController:InvalidMeasNoiseMean', ...
                    '** Input parameter <v_mean> (mean of measurement noise) must be a %d-dimensional vector **',...
                     this.dimMeas);    
                this.measNoiseMean = v_mean(:);
            end         
            
            if nargin > 12
                this.useMexImplementation = useMexImplementation;
            end
            
            Validator.validateTransitionMatrix(caModeTransitionMatrix, sequenceLength + 1);            
            Validator.validateDiscreteProbabilityDistribution(scDelayProb);   
            
            % initially, there is no measurement available (modelled by
            % noise)
            noiseSamples = ...
                Utils.drawGaussianRndSamples(this.measNoiseMean, this.measNoiseCovSqrt, this.measBufferLength); 
            this.augmentedMeasurement = noiseSamples(:);
            this.availableMeasurements = zeros(1, this.measBufferLength);            

            % we need the stationary distribution of the Markov chain
            this.stationaryModeDistribution = Utility.computeStationaryDistribution(caModeTransitionMatrix);
            
            % augmented state consists of x_k, psi_k and eta_k
            dimEta = dimU * (sequenceLength * (sequenceLength - 1) / 2);
            this.dimState = dimX * this.measBufferLength + dimEta;
            this.sysState = zeros(this.dimState, 1);
            [F, G, H, ~] = Utility.createAugmentedPlantModel(sequenceLength, A, B);
            
            % compute the expected augmented sys matrices
            [this.expAugA, this.expAugB, augA, augB] = ...
                this.computeExpectedAugmentedSysMatrices(A, B, F, G, H, dimEta);
            [expAugQ, expAugR] = this.computeExpectedAugmentedCostMatrices(H);
            
            [this.augC, probC] = this.computeAugmentedMeasMatrices(C, ...
                Utility.truncateDiscreteProbabilityDistribution(scDelayProb, this.measBufferLength + 1));
            [this.L, this.K] = this.computeControllerGains(W, V, augA, augB, probC, expAugQ, expAugR);
        end             
        
        %% reset
        function reset(this)
            this.sysState = zeros(this.dimState, 1);
            % initially, there is no measurement available (modelled by
            % noise)
            noiseSamples = ...
                Utils.drawGaussianRndSamples(this.measNoiseMean, this.measNoiseCovSqrt, this.measBufferLength); 
            this.augmentedMeasurement = noiseSamples(:);
            this.availableMeasurements = zeros(1, this.measBufferLength);
            this.lastNumUsedMeas = 0;
            this.lastNumDiscardedMeas = 0;
        end
        
        %% setControllerPlantState
        function setControllerPlantState(this, state)
            % Set the estimate of the plant state.
            %
            % This function is mainly used to set an initial estimate, which is zero by default, as
            % the controller's dynamics is used for updating the estimate of
            % the plant state.
            %
            % Parameters:
            %   >> state (Subclass of Distribution)
            %      The new system state.
            
            assert(Checks.isClass(state, 'Distribution') && state.getDim() == this.dimPlantState, ...
                'InfiniteHorizonUdpLikeController:SetControllerPlantState:InvalidState', ...
                '** <state> must be a %d-dimensional Distribution **', this.dimPlantState);
            
            [x0, ~] = state.getMeanAndCov();
            this.sysState = [x0; zeros(this.dimState - this.dimPlantState, 1)];
        end       
        
        %% getControllerPlantState
        function plantState = getControllerPlantState(this)
            % Get the plant state as currently perceived by the controller, i.e., its internal estimate of the plant state.
            %
            % Returns:
            %   << plantState (Column vector)
            %      The controller's current estimate of the plant state.
            %
            plantState = this.sysState(1:this.dimPlantState);
        end
        
        %% getControllerGains
        function [K, L] = getControllerGains(this)
            % Get the parameters/gains of the linear, mode-independent control law.
            %
            % Returns:
            %   << K (Matrix)
            %      The observer gain K.
            %
            %   << L (Matrix)
            %      The controller gain L.
            %      In contrast to the paper, the control sequence U_k is
            %      given as U_k = L*state_k, so the gain returned is the
            %      negative of the gain computed in the paper.
            %
            K = this.K;
            L = this.L;
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
            
            numUsedMeas = this.lastNumUsedMeas;
            numDiscardedMeas = this.lastNumDiscardedMeas;
        end
        
        %% computeControlSequence
        function inputSequence = computeControlSequence(this, measurements, measurementDelays, ~, ~)
            if nargin == 1 || isempty(measurements)
                applicableMeas = [];
                applicableDelays = [];
            else
                this.checkMeasurementsAndDelays(measurements, measurementDelays);
                [applicableMeas, applicableDelays] = this.getApplicableMeasurements(measurements, measurementDelays);
            end
            
            % absence of measurement is modelled by noise
            noiseSamples = ...
                Utils.drawGaussianRndSamples(this.measNoiseMean, this.measNoiseCovSqrt, this.measBufferLength); 
            this.augmentedMeasurement = noiseSamples(:);
            this.availableMeasurements = zeros(1, this.measBufferLength);
            % now incorporate the received measurements
            for i=1:numel(applicableDelays)
                % store the measurement
                delay = applicableDelays(i);
                this.availableMeasurements(delay + 1) = 1;
                this.augmentedMeasurement(delay * this.dimMeas + 1:(delay + 1) * this.dimMeas) ...
                    = applicableMeas(:, i);
            end
            % now use the measurements to pick the corresponding
            % measurement matrix according to the availability
            idx = sum(2 .^ [0:this.measBufferLength - 1] .* this.availableMeasurements) + 1;
            innovation = this.augmentedMeasurement - this.augC(:, :, idx) * this.sysState;
            
            inputSequence = this.doControlSequenceComputation();
            % update controller state
            this.sysState = this.expAugA * this.sysState + this.K * innovation + this.expAugB * inputSequence;
           
            if nargin == 1
                this.lastNumUsedMeas = 0;
                this.lastNumDiscardedMeas = 0;
            else
                this.lastNumUsedMeas = size(applicableMeas, 2);
                this.lastNumDiscardedMeas = size(measurements, 2) - this.lastNumUsedMeas;
            end
        end
    end
    
    methods (Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this)
            inputSequence = this.L * this.sysState;
        end
         
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, ~)
                        
            stageCosts = Utility.computeStageCosts(state, input, this.Q, this.R);
        end
         
         %% doCostsComputation
         function averageLQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonLength = size(appliedInputs, 2);
            assert(size(stateTrajectory, 2) == horizonLength + 1, ...
                'InfiniteHorizonUdpLikeController:DoCostsComputation', ...
                '** <stateTrajectory> is expected to have %d columns ', horizonLength + 1);
            
            averageLQGCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, appliedInputs, this.Q, this.R) / horizonLength;
        end
    end
    
    methods (Access = private)
        %% computeControllerGains
        function [L,K] = computeControllerGains(this, W, V, augA, augB, probC, expAugQ, expAugR)
            augW = blkdiag(W, zeros(this.dimState - this.dimPlantState));
            augV = Utils.blockDiag(V, this.measBufferLength); % sparse
            
            if this.useMexImplementation                
                 [L, K] = mex_InfiniteHorizonUdpLikeController(augA, augB, probC, expAugQ, expAugR, augW, augV, ...
                    this.expAugA, this.expAugB, this.stationaryModeDistribution, this.augC);                 
                return
            end          
    
            numCombinations = 2 ^ this.measBufferLength;
            numModes = this.sequenceLength + 1;
            
            L = zeros(this.dimPlantInput * this.sequenceLength, this.dimState);
            K = zeros(this.dimState, this.dimAugmentedMeas);
            
            L_old = inf(this.dimPlantInput * this.sequenceLength, this.dimState);
            K_old = inf(this.dimState, this.dimAugmentedMeas);
            
            OverbarPsi = zeros(this.dimState);
            UnderbarPsi = zeros(this.dimState);
            OverbarLambda = zeros(this.dimState);
            UnderbarLambda = zeros(this.dimState);
            
            OverbarPsi_old = inf(this.dimState);
            UnderbarPsi_old = inf(this.dimState);
            
            maxIterNum = 10000;
            convergenceDiff = 1e-10;            
            k = 1;
            while (k <= maxIterNum) && (sum(sum(abs(OverbarPsi_old - OverbarPsi))) > convergenceDiff ...
                    || sum(sum(abs(UnderbarPsi_old - UnderbarPsi))) > convergenceDiff)
                % compute the next iterates
                K_old = K;
                L_old = L;
                OverbarPsi_old = OverbarPsi;
                UnderbarPsi_old = UnderbarPsi;
                OverbarLambda_old = OverbarLambda;
                UnderbarLambda_old = UnderbarLambda;

                BexpL = this.expAugB * L_old;
                part = this.expAugA - BexpL;
                part2 = part * UnderbarPsi_old * part'; % should be symmetric
                part2 = (part2 + part2') / 2;
                part3 = part' * UnderbarLambda_old * part; % should be symmetric
                part3 = (part3 + part3') / 2;
                KVK = K_old * augV * K_old'; % KVK' % should be symmetric
                KVK = (KVK + KVK') / 2;
                LRL = L_old' * expAugR * L_old; % L'RL % should be symmetric
                LRL = (LRL + LRL') / 2;
                
                UnderbarPsi = part2 + KVK;
                OverbarPsi = -part2 + KVK + augW;
                OverbarLambda = LRL + expAugQ - part3;
                UnderbarLambda = LRL;
                
                Lambda_sum = OverbarLambda_old + UnderbarLambda_old;

                % two separate terms required for L
                Bexp_Lambda = -this.expAugB' * UnderbarLambda_old;                              
                % compute the feedback gain L
                B_lambda_Sum = mtimesx(augB, 'T', Lambda_sum);
                L1 = Bexp_Lambda * this.expAugB + expAugR + sum(reshape(this.stationaryModeDistribution, 1, 1, numModes) ...
                    .* mtimesx(B_lambda_Sum, augB), 3);  % should be symmetric
                L2 = Bexp_Lambda * this.expAugA + sum(reshape(this.stationaryModeDistribution, 1, 1, numModes) ...
                    .* mtimesx(B_lambda_Sum, augA) , 3);
                L = pinv((L1 + L1') / 2) * L2;
                %L = lsqminnorm((L1+L1')/2, L2);
                              
                % two separate terms required for K
                PsiC = mtimesx(OverbarPsi_old, this.augC, 'T');
                CPsiC = mtimesx(this.augC, PsiC); % should be symmetric                
                % expectation with respect to C
                K2 = this.expAugA * sum(reshape(probC, 1, 1, numCombinations) .* PsiC, 3);
                K1 = augV + sum(reshape(probC, 1, 1, numCombinations) .* CPsiC, 3); % should be symmetric                 
                K = K2 * pinv((K1 + K1') / 2);                
                %K = lsqminnorm((K1 + K1') / 2, K2')';
 
                K_C = mtimesx(K_old, this.augC); % KC for all i
                
                % expectation with respect to C
                KCOverbarPsiKC = K_old * sum(reshape(probC, 1, 1, numCombinations) .* CPsiC, 3) * K_old'; % should be symmetric
                UnderbarPsi = UnderbarPsi + (KCOverbarPsiKC + KCOverbarPsiKC') / 2; % complete
                
                Aexp_KC = this.expAugA - K_C;
                BexpL_KC = BexpL - K_C;                
                Aexp_KCUnderbarLambdaAexp_KC = mtimesx(Aexp_KC, 'T', mtimesx(UnderbarLambda_old, Aexp_KC)); % should be symmetric
                BexpL_KCUnderbarLambdaBexpL_KC = mtimesx(BexpL_KC, 'T', mtimesx(UnderbarLambda_old, BexpL_KC)); % should be symmetric
                                
                partSum = sum(reshape(probC, 1, 1, numCombinations) .* (Aexp_KCUnderbarLambdaAexp_KC - BexpL_KCUnderbarLambdaBexpL_KC), 3);
                UnderbarLambda = UnderbarLambda + (partSum + partSum') / 2;
                
                BL = mtimesx(augB, L_old); % BL for all modes i
                ABL = augA - BL;
                ABLPsiABL = mtimesx(ABL, mtimesx(UnderbarPsi_old, ABL, 'T')); % should be symmetric
                ABLLambdaABL = mtimesx(ABL, 'T', mtimesx(Lambda_sum, ABL)); % should be symmetric
                BLambdaB = mtimesx(augB, 'T', mtimesx(OverbarLambda_old, augB)); %B'*Lambda*B for all modes i
                
                OverbarPsi = OverbarPsi + sum(reshape(this.stationaryModeDistribution, 1, 1, numModes) ...
                    .* (ABLPsiABL + permute(ABLPsiABL, [2 1 3])) / 2, 3);
                OverbarLambda = OverbarLambda + sum(reshape(this.stationaryModeDistribution, 1, 1, numModes) ...
                    .* (ABLLambdaABL + permute(ABLLambdaABL, [2 1 3])) / 2, 3); % complete
                LBLambdaBL = L_old' * sum(reshape(this.stationaryModeDistribution, 1, 1, numModes) .* BLambdaB, 3) * L_old;
                UnderbarLambda = UnderbarLambda + (LBLambdaBL + LBLambdaBL') / 2;                
                
                for j=1:numModes
                    A_KC = augA(:, :, j) - K_C;
                    BL_KC = BL(:, :, j) - K_C;
                    AKCPsiAKC = mtimesx(A_KC, mtimesx(OverbarPsi_old, A_KC, 'T')); % should be symmetric
                    BLKCLambdaBLKC = mtimesx(BL_KC, 'T', mtimesx(UnderbarLambda_old, BL_KC)); % should be symmetric
                    
                    innerExpectationOverbarPsi = sum(reshape(probC, 1, 1, numCombinations) ...
                        .*  AKCPsiAKC, 3); % should be symmetric
                    innerExpectationUnderbarLambda = sum(reshape(probC, 1, 1, numCombinations) ...
                        .*  BLKCLambdaBLKC, 3); % should be symmetric

                    OverbarPsi = OverbarPsi + ...
                        this.stationaryModeDistribution(j) * (innerExpectationOverbarPsi + innerExpectationOverbarPsi') / 2; % complete
                    UnderbarLambda = UnderbarLambda + ...
                        this.stationaryModeDistribution(j) * (innerExpectationUnderbarLambda + innerExpectationUnderbarLambda') / 2; % complete
                end
                % end compute next iterates
                k = k + 1;
            end
            L = -L; 
          end
        
        %% computeAugmentedMeasMatrices
        function [augC, probC] = computeAugmentedMeasMatrices(this, C, scDelayProbs)
            numCombinations = 2 ^ this.measBufferLength;
            probSums = cumsum(scDelayProbs);
            compProbSums = 1 -probSums;
            
            augC = zeros(this.dimAugmentedMeas, this.dimState, numCombinations);
            probC = ones(1, numCombinations);
                        
            % first slice (augC(:, :, 1)) is zero matrix, i.e., no
            % measurements available
            % last slice (augC(:, :, end)) indicates case that all
            % measurements are available
            for i=0:numCombinations -1
                element = bitget(i, 1:1:this.measBufferLength); % binary vector
                % construct the augmented meas matrix corresponding to the
                % case and the corresponding probability of occurrence
                for j=1:this.measBufferLength
                    if element(j) == 1
                        rowIdx = (j-1) * this.dimMeas + 1:j * this.dimMeas;
                        colIdx = (j-1) * this.dimPlantState + 1: j * this.dimPlantState;
                        augC(rowIdx, colIdx, i+1) = C;
                        probC(i+1) = probC(i+1) * probSums(j);
                    else
                        probC(i+1) = probC(i+1) * compProbSums(j);
                    end
                end
            end
        end
        
        %% computeExpectedAugmentedSysMatrices
        function [expAugA, expAugB, augA, augB] = computeExpectedAugmentedSysMatrices(this, A, B, F, G, H, dimEta)
            numModes = this.sequenceLength + 1;
            E = diag(ones(this.dimPlantState * (this.measBufferLength - 2), 1), -this.dimPlantState);
            D = eye(this.dimPlantState); % first part of the D matrix in the paper
            
            augA = repmat(blkdiag(A, E, F), 1, 1, numModes);
            augB = repmat([zeros(this.dimState - size(G, 1), size(G, 2)); G], 1, 1, numModes);
            % augB is only different for first mode (multiply B by J)
            BJ = B * [speye(this.dimPlantInput) zeros(this.dimPlantInput, this.dimPlantInput * (this.sequenceLength -1))];
            augB(1:this.dimPlantState, :, 1) ...
                = BJ;

            % compute the expected augmented sys matrices
            % only the upper right corner of expAugA is stochastic,
            % mode-dependent, but H is zero for first and last mode
            BH = mtimesx(B, H(:, :, 2:this.sequenceLength)); % compute B * H(:, :, i) in one shot for all i (expect first and last mode)
            augA(1:this.dimPlantState, end - dimEta + 1:end, 2:this.sequenceLength) = BH;
            augA(this.dimPlantState + 1:2 * this.dimPlantState, 1:this.dimPlantState, :) = repmat(D, 1,1, numModes);
            expBH = sum(reshape(this.stationaryModeDistribution(2:this.sequenceLength), 1, 1, numModes - 2) .* BH, 3);

            expAugA = blkdiag(A, E, F);
            expAugA(this.dimPlantState + 1:2 * this.dimPlantState, 1:this.dimPlantState) = D;
            expAugA(1:this.dimPlantState, end - dimEta + 1:end) = expBH;
      
            % as with augB, only the first mode contributes to the expected
            % matrix, as B*J is always zero for all other modes
            expAugB = [ ...
                        this.stationaryModeDistribution(1) * BJ; ...
                        zeros((this.measBufferLength - 1) * this.dimPlantState, size(G, 2)); ...
                        G
                      ];
        end
        
        %% computeExpectedAugmentedCostMatrices
        function [expAugQ, expAugR] = computeExpectedAugmentedCostMatrices(this, H)
            dimZeros = (this.measBufferLength - 1) * this.dimPlantState;
  
            % H is zero for first and last mode
            HRH = mtimesx(H(:, :, 2:this.sequenceLength), 'T', ...
                mtimesx(this.R, H(:, :, 2:this.sequenceLength))); % should be symmetric    
            symmHRH =  (HRH + permute(HRH, [2 1 3])) / 2; % ensure the symmetry of each H'RH
            
            expHRH = sum(reshape(this.stationaryModeDistribution(2:this.sequenceLength), 1, 1, []) .* symmHRH, 3);
            expAugQ = blkdiag(this.Q, zeros(dimZeros), expHRH);
            
            % J matrix is only nonzero for first mode, so only this
            % contributes to expAugR
            J = [speye(this.dimPlantInput) zeros(this.dimPlantInput, this.dimPlantInput * (this.sequenceLength -1))];
            expAugR = this.stationaryModeDistribution(1) * J' * this.R * J;
       end
        
        %% checkMeasurementsAndDelays
        function delays = checkMeasurementsAndDelays(this, measurements, measDelays)
            assert(Checks.isFixedRowMat(measurements, this.dimMeas), ...
                'InfiniteHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeas', ...
                '** Individual measurements must be %d-dimensional **', ...
                this.dimMeas);
            
            numMeas = size(measurements, 2);
            assert(Checks.isNonNegativeVec(measDelays, numMeas) ...
                    && all(arrayfun(@(delay) mod(delay, 1) == 0, measDelays)), ...
                'InfiniteHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeasDelay', ...
                '** Each measurement delay (%d in total) must be a nonnegative integer **', numMeas);
            
            delays = measDelays;
        end
        
        %% getApplicableMeasurements
        function [applicableMeas, applicableDelays] = getApplicableMeasurements(this, measurements, measDelays)
            numMeas = numel(measDelays);
            % the controller assumes only a single measurement per time step
            [~, uniqueIdx, ~] = unique(measDelays, 'stable');
            assert(numel(uniqueIdx) == numMeas, ...
                'InfiniteHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsNotUnique', ...
                '** %s\nIgnoring %d of %d measurements. **', ...
                'Controller assumes only one measurement per time step', ...
                numMeas - numel(uniqueIdx), numMeas);
  
            % find the measurements with valid delays
            idx = find(measDelays(uniqueIdx) <= this.measBufferLength - 1);
            applicableMeas = measurements(:, idx);
            applicableDelays = measDelays(idx);
            if numel(idx) ~= numel(uniqueIdx)
                warning('InfiniteHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsDelayTooLarge', ...
                    '** %s\nIgnoring %d of %d measurements. **', ...
                    'Delay too large', ...
                    numel(uniqueIdx) - numel(idx), numel(uniqueIdx));
            end
        end
    end
end

