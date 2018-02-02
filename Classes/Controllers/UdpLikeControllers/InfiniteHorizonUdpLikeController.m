classdef InfiniteHorizonUdpLikeController < SequenceBasedController
    % Implementation of an infinite horizon linear sequence-based LQG controller for
    % NCS with a UDP-like networks connecting the controller and actuator, 
    % and sensor and controller, respectively.
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
    
    properties (SetAccess = protected)
        % packet delay probability density function of the controller-actuator-
        % link; (vector of dimension >= <sequenceLength>)
        caDelayProb;
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
        C; % measurement matrix
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
        dimMeas;
        measBufferLength;
        dimAugmentedMeas;
        augC;
    end
    
     properties (Access = private)
        lastNumUsedMeas = 0;
        lastNumDiscardedMeas = 0;
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
        function this = InfiniteHorizonUdpLikeController(A, B, C, Q, R, caDelayProb, scDelayProb, ...
                sequenceLength, maxMeasDelay, W, V)
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
            %   >> caDelayProb (Nonnegative vector)
            %      The vector describing the delay distribution of the
            %      CA-network.
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
            % Returns:
            %   << this (InfiniteHorizonUdpLikeController)
            %      A new InfiniteHorizonUdpLikeController instance.
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength);
            
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
            
            Validator.validateMeasurementMatrix(C, dimX);
            this.dimMeas = size(C, 1);
            if ~Checks.isNonNegativeScalar(maxMeasDelay) || mod(maxMeasDelay, 1) ~= 0
                 error('InfiniteHorizonUdpLikeController:InvalidMaxMeasDelay', ...
                    ['** Input parameter <maxMeasDelay> (maximum measurement',...
                     'delay (M-1)) must be a nonnegative integer **']);
            end
                        
            this.measBufferLength = maxMeasDelay + 1; % maxMeasDelay + 1 is buffer length
            
            this.dimAugmentedMeas = this.dimMeas * this.measBufferLength;
            % check the noise covs
            Validator.validateSysNoiseCovarianceMatrix(W, dimX);
            this.measNoiseCovSqrt = Validator.validateMeasNoiseCovarianceMatrix(V, this.dimMeas);
            
            % initially, there is no measurement available (modelled by
            % noise)
            noiseSamples = ...
                Utils.drawGaussianRndSamples(zeros(this.dimMeas, 1), this.measNoiseCovSqrt, this.measBufferLength); 
            this.augmentedMeasurement = noiseSamples(:);
            this.availableMeasurements = bitget(0, 1:this.measBufferLength);
     
            Validator.validateDiscreteProbabilityDistribution(caDelayProb);
            Validator.validateDiscreteProbabilityDistribution(scDelayProb);
            
            elementCount = numel(caDelayProb);
            if elementCount <= sequenceLength
                % fill up with zeros (row vector)
                probHelp = [reshape(caDelayProb, 1, elementCount), zeros(1, sequenceLength + 1 - elementCount)];
            else
                % cut up the distribution and normalize
                probHelp = [caDelayProb(1:sequenceLength), 1 - sum(caDelayProb(1:sequenceLength))];
            end
            % store as row vector, i.e., column-wise arranged
            this.caDelayProb = probHelp;
            
            elementCount = numel(scDelayProb);
            if elementCount <= this.measBufferLength
                % fill up with zeros (row vector)
                trueScDelayProbs = [reshape(scDelayProb, 1, elementCount), zeros(1, this.measBufferLength + 1 - elementCount)];
            else
                % cut up the distribution and normalize
                trueScDelayProbs = [scDelayProb(1:this.measBufferLength), 1 - sum(scDelayProb(1:this.measBufferLength))];
            end
            
            transitionMatrix = Utility.calculateDelayTransitionMatrix(this.caDelayProb);
            % we need the stationary distribution of the Markov chain
            this.stationaryModeDistribution = Utility.computeStationaryDistribution(transitionMatrix);
               
            % augmented state consists of x_k, psi_k and eta_k
            dimEta = dimU * (sequenceLength * (sequenceLength - 1) / 2);
            this.dimState = dimX * this.measBufferLength + dimEta;
            this.sysState = zeros(this.dimState, 1);
            [F, G, H, ~] = Utility.createAugmentedPlantModel(sequenceLength, A, B);
            
            % compute the expected augmented sys matrices
            [this.expAugA, this.expAugB, augA, augB] = ...
                this.computeExpectedAugmentedSysMatrices(A, B, F, G, H, dimEta);
            [expAugQ, expAugR] = this.computeExpectedAugmentedCostMatrices(Q, R, H);
            
            [this.augC, probC] = this.computeAugmentedMeasMatrices(C, trueScDelayProbs);
            [this.L, this.K] = this.computeControllerGains(W, V, augA, augB, probC, expAugQ, expAugR);
        end
        
        %% reset
        function reset(this)
            this.sysState = zeros(this.dimState, 1);
            % initially, there is no measurement available (modelled by
            % noise)
            noiseSamples = ...
                Utils.drawGaussianRndSamples(zeros(this.dimMeas, 1), this.measNoiseCovSqrt, this.measBufferLength); 
            this.augmentedMeasurement = noiseSamples(:);
            this.availableMeasurements = bitget(0, 1:this.measBufferLength);
            this.lastNumUsedMeas = 0;
            this.lastNumDiscardedMeas = 0;
        end
        
        %% getControllerPlantState
        function plantState = getControllerPlantState(this)
            % Get the plant state as currently perceived by the controller, i.e., its internal estimate of the plant state..
            %
            % Returns:
            %   << plantState (Column vector)
            %      The controller's current estimate of the plant state.
            %
            plantState = this.sysState(1:this.dimPlantState);
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
        function inputSequence = computeControlSequence(this, measurements, measurementDelays)
            if nargin == 1 || isempty(measurements)
                applicableMeas = [];
                applicableDelays = [];
            else
                this.checkMeasurementsAndDelays(measurements, measurementDelays);
                [applicableMeas, applicableDelays] = this.getApplicableMeasurements(measurements, measurementDelays);
            end
            noiseSample = Utils.drawGaussianRndSamples(zeros(this.dimMeas, 1), this.measNoiseCovSqrt, 1);
            this.augmentedMeasurement = [noiseSample; this.augmentedMeasurement(1:end-this.dimMeas)]; % shift measurements "one downwards"
            this.availableMeasurements = [0 this.availableMeasurements(1:end-1)];

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
         
         %% doCostsComputation
         function averageLQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonLength = size(appliedInputs, 2);
            if size(stateTrajectory, 2) ~= horizonLength + 1
                error('InfiniteHorizonUdpLikeController:DoCostsComputation', ...
                    '** <stateTrajectory> is expected to have %d columns ', horizonLength + 1);
                    
            end
            averageLQGCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, appliedInputs, this.Q, this.R) / horizonLength;
        end
    end
    
    methods (Access = private)
        %% computeControllerGains
        function [L,K] = computeControllerGains(this, W, V, augA, augB, probC, expAugQ, expAugR)
            augW = blkdiag(W, zeros(this.dimState - this.dimPlantState));
            augV = Utils.blockDiag(V, this.measBufferLength);
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
            
            maxIterNum = 10000;
            convergenceDiff = 1e-10;
            
            k = 1;
            while k <= maxIterNum && (k == 2 || norm (K_old - K) > convergenceDiff ...
                    && norm (L_old - L) > convergenceDiff)
                % compute the next iterates
                K_old = K;
                L_old = L;
                OverbarPsi_old = OverbarPsi;
                UnderbarPsi_old = UnderbarPsi;
                OverbarLambda_old = OverbarLambda;
                UnderbarLambda_old = UnderbarLambda;

                BexpL = this.expAugB * L_old;
                part = this.expAugA - BexpL;
                part2 = part * UnderbarPsi_old * part';
                KVK = K_old * augV * K_old'; % KVK'
                LRL = L_old' * expAugR * L_old; % L'RL

                UnderbarPsi = part2 + KVK;
                OverbarPsi = -part2 + KVK + augW;
                OverbarLambda = LRL + expAugQ -part' * UnderbarLambda_old * part;
                UnderbarLambda = LRL;

                Lambda_sum = OverbarLambda_old + UnderbarLambda_old;

                % two seperate terms required for L
                Bexp_Lambda = -this.expAugB' * UnderbarLambda_old;
                L1 = Bexp_Lambda * this.expAugB + expAugR;
                L2 = Bexp_Lambda * this.expAugA;

                % two seperate terms required for K
                K1 = augV;
                K2 = zeros(this.dimState, this.dimAugmentedMeas);

                K_C = zeros(this.dimState, this.dimState, numCombinations);

                for i=1:numCombinations
                    K_C(:, :, i) = K_old * this.augC(:, :, i);
                    UnderbarPsi = UnderbarPsi +  probC(i) * K_C(:, :, i) * OverbarPsi_old * K_C(:, :, i)'; % complete
                    
                    Aexp_KC = this.expAugA - K_C(:, :, i);
                    BexpL_KC = BexpL - K_C(:, :, i);
                    UnderbarLambda = UnderbarLambda + ...
                        probC(i) * (Aexp_KC' * UnderbarLambda_old * Aexp_KC ...
                                    - BexpL_KC' * UnderbarLambda_old * BexpL_KC);
                    
                    OverbarPsi_C = probC(i) * OverbarPsi_old * this.augC(:, :, i)';
                                        
                    K1 = K1 + this.augC(:, :, i) * OverbarPsi_C;
                    K2 = K2 + OverbarPsi_C;
                end

                for j=1:numModes
                    BL = augB(:, :, j) * L_old;
                    A_BL = augA(:, :, j) - BL;

                    OverbarPsi = OverbarPsi + ...
                        this.stationaryModeDistribution(j) * A_BL * UnderbarPsi_old * A_BL';
                    OverbarLambda = OverbarLambda + ...
                        this.stationaryModeDistribution(j) * (A_BL' * Lambda_sum * A_BL); % complete
                    UnderbarLambda = UnderbarLambda + ...
                            this.stationaryModeDistribution(j) * (BL' * OverbarLambda_old * BL);

                    % compute the feedback gain L
                    B_LambdaSum = this.stationaryModeDistribution(j) * augB(:, :, j)' * Lambda_sum;
                    L1 = L1 + B_LambdaSum * augB(:, :, j);
                    L2 = L2 + B_LambdaSum * augA(:, :, j);

                    innerExpectationOverbarPsi = zeros(this.dimState);
                    innerExpectationUnderbarLambda = zeros(this.dimState);
                    for i=1:numCombinations
                        A_KC = augA(:, :, j) - K_C(:, :, i);
                        BL_KC = BL - K_C(:, :, i);
  
                        innerExpectationOverbarPsi = innerExpectationOverbarPsi + ...
                            probC(i) * A_KC * OverbarPsi_old * A_KC';

                        innerExpectationUnderbarLambda = innerExpectationUnderbarLambda + ...
                             probC(i) * BL_KC' * UnderbarLambda_old * BL_KC;
                    end

                    OverbarPsi = OverbarPsi + ...
                        this.stationaryModeDistribution(j) * innerExpectationOverbarPsi; % complete
                    UnderbarLambda = UnderbarLambda + ...
                        this.stationaryModeDistribution(j) * innerExpectationUnderbarLambda; % complete
                end
                L = pinv(L1) * L2;
                K = this.expAugA * K2 * pinv(K1);
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
            % mode-dependent
            expBH = zeros(this.dimPlantState, dimEta);
            for i = 1:numModes
                augA(1:this.dimPlantState, end - dimEta + 1:end, i) = B * H(:,:,i);
                augA(this.dimPlantState + 1:2 * this.dimPlantState, 1:this.dimPlantState, i) = D;
                
                expBH = expBH + this.stationaryModeDistribution(i) * B * H(:,:,i);
            end
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
        function [expAugQ, expAugR] = computeExpectedAugmentedCostMatrices(this, Q, R, H)
            numModes = this.sequenceLength + 1;
            dimZeros = (this.measBufferLength - 1) * this.dimPlantState;
            
            expHRH = zeros(this.dimState - dimZeros - this.dimPlantState);

            for i = 1:numModes
                expHRH = expHRH + this.stationaryModeDistribution(i) * H(:, :, i)' * R * H(:, :, i);
            end
            expAugQ = blkdiag(Q, zeros(dimZeros), expHRH);
            
            % J matrix is only nonzero for first mode, so only this
            % contributes to expAugR
            J = [speye(this.dimPlantInput) zeros(this.dimPlantInput, this.dimPlantInput * (this.sequenceLength -1))];
            expAugR = this.stationaryModeDistribution(1) * J' * R * J;
       end
        
        %% checkMeasurementsAndDelays
        function delays = checkMeasurementsAndDelays(this, measurements, measDelays)
            if ~Checks.isFixedRowMat(measurements, this.dimMeas)
                error('InfiniteHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeas', ...
                    '** Individual measurements must be %d-dimensional **', this.dimMeas);
            end
            numMeas = size(measurements, 2);
            
            if ~Checks.isNonNegativeVec(measDelays, numMeas) ...
                    || any(arrayfun(@(measDelay) mod(measDelay, 1) ~= 0, measDelays))
                error('InfiniteHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeasDelay', ...
                    '** Each measurement delay (%d in total) must be a nonnegative integer **', numMeas);
            end
            delays = measDelays;
        end
        
        %% getApplicableMeasurements
        function [applicableMeas, applicableDelays] = getApplicableMeasurements(this, measurements, measDelays)
            numMeas = numel(measDelays);
            % the controller assumes only a single measurement per time step
            [~, uniqueIdx, ~] = unique(measDelays, 'stable');
            if numel(uniqueIdx) ~= numMeas
                error('InfiniteHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsNotUnique', ...
                    '** %s\nIgnoring %d of %d measurements. **', ...
                    'Controller assumes only one measurement per time step', ...
                    numMeas - numel(uniqueIdx), numMeas);
            end
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

