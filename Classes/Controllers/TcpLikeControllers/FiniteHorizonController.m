classdef FiniteHorizonController < SequenceBasedController
    % Implementation of the optimal finite horizon linear sequence-based LQG controller for
    % NCS with a TCP-like network connecting the controller and the
    % actuator.
    % In order to perform control tasks subject to state and/or input
    % constraints, this implementation supports linear integral state and input constraints 
    % of the form E{a_0' * x_0 + b_0' * u_0 + a_1' * x_1 + b_1' * u_1 + ... + a_K' * x_K} =< c, 
    % i.e., soft expectation-type integral constraints, where a_i, b_i are state
    % and input weighting, c is a constant and K is the horizon length.
    %
    % While this implementation is to some extent more general than the original one by Maxim Dolgov, 
    % which can be found <a href="matlab:web('http://www.cloudrunner.eu/algorithm/160/constrained-slqg-controller/version/1/')" >here</a>, 
    % some parts are directly based thereof.
    %
    % Literature: 
    %   Jörg Fischer, Achim Hekler, Maxim Dolgov, and Uwe D. Hanebeck,
    %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject to Random Transmission Delays and Packet Losses,
    %   Proceedings of the 2013 American Control Conference (ACC 2013),
    %   Washington D. C., USA, June 2013.
    %      
    %   Maxim Dolgov, Jörg Fischer, and Uwe D. Hanebeck,
    %   Sequence-based LQG Control over Stochastic Networks with Linear Integral Constraints,
    %   Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014), 
    %   Los Angeles, California, USA, December 2014.

    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        % system state cost matrix per time step;
        % (positive semidefinite matrix of dimension <dimX> x <dimX>)
        Q = [];
        % input cost matrix;
        % (positive definite matrix of dimension <dimU> x <dimU>)
        R = [];       
        % CA-network-actuator system
        % mapping for propagation of possible control inputs;
        % (matrix of dimension <dimU*sequenceLength*(sequenceLength-1)/2> x 
        %  <dimU*packetLength*(sequenceLength-1)/2>)
        F = [];
        % matrix mapping control sequence to the possible control inputs;
        % (matrix of dimension <dimU*packetLength*(sequenceLength-1)/2> x
        %  <dimU*sequenceLength>)
        G = [];
        % Markov chain
        % transition matrix of the Markov chain;
        % (matrix of dimension <sequenceLength+1> x <sequenceLength+1>)
        transitionMatrix = [];
        L = [];
        % dimension of the augmented system state; 
        % (positive integer <dimX+dimU*sequenceLength*(sequenceLength-1)/2>)
        dimState = -1;
    end
    
    properties (Access = private)
        feedforward = [];
        multiplier;
        S_0;
        P_0;
        constraintBounds;
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        horizonLength;
        constraintsPresent@logical = false;        
        
        useMexImplementation@logical=true; 
        % by default, we use the C++ (mex) implementation for computation of controller gains
        % this is faster, but can produce slightly different results
    end
    
    properties (SetAccess = private, GetAccess = protected)
        % augmented state space
        % system state of the augmented state space;
        % (column vector of dimension <dimState>)
        sysState = [];
    end
    
    methods (Access = public)
        %% FiniteHorizonController
        function this = FiniteHorizonController(A, B, Q, R, delayProb, sequenceLength, horizonLength, useMexImplementation, ...
                stateConstraintWeightings, inputConstraintWeightings, constraintBounds)
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
            %      The state weighting matrix in the controller's underlying cost function.
            %
            %   >> R (Positive definite matrix)
            %      The input weighting matrix in the controller's underlying cost function.
            %
            %   >> delayProb (Nonnegative vector)
            %      The vector describing the delay distribution of the
            %      CA-network.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %   >> horizonLength (Positive integer)
            %      The length of the considered control horizon.
            %
            %   >> useMexImplementation (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for the computation of the controller
            %      gains which is considerably faster than the Matlab
            %      implementation. 
            %      If left out, the default value true is used.
            %
            %   >> stateConstraintWeightings (3D matrix, optional)
            %      The state weightings of the individual
            %      constraints, slice-wise arranged.
            %
            %   >> inputConstraintWeightings (3D matrix, optional)
            %      The input weightings of the individual
            %      constraints, slice-wise arranged.
            %
            %   >> constraintBounds (Vector, optional)
            %      The bounds of the individual
            %      constraints, given as a row or column vector.
            %
            % Returns:
            %   << obj (FiniteHorizonController)
            %      A new FiniteHorizonController instance.
            
            assert(nargin == 7 || nargin == 8 || nargin == 11, ...
                'FiniteHorizonController:InvalidNumberOfArguments', ...
                '** Constructor must be called with 7, 8 or 11 arguments **');
                        
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength, true);
            
            Validator.validateHorizonLength(horizonLength);
            this.horizonLength = horizonLength;
            
            if nargin > 7
                assert(Checks.isFlag(useMexImplementation), ...
                    'FiniteHorizonController:InvalidUseMexFlag', ...
                    '** <useMexImplementation> must be a flag **');
                this.useMexImplementation = useMexImplementation;
            end
                        
            if nargin == 11
                this.validateLinearIntegralConstraints(stateConstraintWeightings, ...
                    inputConstraintWeightings, constraintBounds);
                this.constraintsPresent = true;
                this.constraintBounds = constraintBounds;
            end
            
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
            
            % delayProb
            Validator.validateDiscreteProbabilityDistribution(delayProb);
            this.transitionMatrix = Utility.calculateDelayTransitionMatrix( ...
                Utility.truncateDiscreteProbabilityDistribution(delayProb, sequenceLength + 1)); 
                        
            % augmented initial system state
            [this.dimState, this.F, this.G, augA, augB, augQ, augR] ...
                = Utility.performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R);
            this.sysState = zeros(this.dimState, 1);
         
            terminalK = blkdiag(this.Q, zeros(this.dimState - dimX)); % K_N for all modes (P_N in second paper mentioned above), symmetric

            if this.constraintsPresent
                terminalStateConstraintWeightings = squeeze(stateConstraintWeightings(:, end, :));
                numConstraints = numel(constraintBounds);
                                
                terminalP = zeros(this.dimState, numConstraints); % P_N for all modes (P_tilde in the paper)
                terminalP(1:dimX, :) = terminalStateConstraintWeightings;

                if this.useMexImplementation
                    [this.L, this.feedforward, this.P_0, this.S_0] = ...
                        mex_FiniteHorizonControllerWithConstraints(augA, augB, augQ, augR, this.transitionMatrix, ...
                            terminalK, terminalP, inputConstraintWeightings, stateConstraintWeightings(:, 1:end-1, :), ...
                            this.horizonLength);
                else                    
                    [this.L, this.feedforward, this.P_0, this.S_0] = ...
                        this.computeGainMatricesWithConstraints(augA, augB, augQ, augR, ...
                            inputConstraintWeightings, stateConstraintWeightings(:, 1:end-1, :), terminalK, terminalP);                 
                end               
            elseif this.useMexImplementation                
                this.L = mex_FiniteHorizonController(augA, augB, augQ, augR, this.transitionMatrix, terminalK, this.horizonLength);                
            else                
                this.L = this.computeGainMatrices(augA, augB, augQ, augR, terminalK);                
            end           
        end
        
        %% reset
        function reset(this)
           this.sysState = zeros(this.dimState,1);
        end % function reset
        
        %% setEtaState
        function setEtaState(this, newEta)
            assert(Checks.isVec(newEta, this.dimState - this.dimPlantState) && all(isfinite(newEta)), ...
                'FiniteHorizonController:SetEtaState:InvalidEta', ...
                '** <newEta> must be a %d-dimensional vector **', ...
                this.dimState - this.dimPlantState);
            
            this.sysState(this.dimPlantState + 1:end) = newEta(:);
        end
        
    end
    
    methods (Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, mode, timestep)
            assert(Checks.isScalarIn(timestep, 1, this.horizonLength) && mod(timestep, 1) == 0, ...
                'FiniteHorizonController:DoControlSequenceComputation:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be in {1, ... %d} **', ...
                this.horizonLength);
            assert(Checks.isScalarIn(mode, 1, this.sequenceLength + 1) && mod(mode, 1) == 0, ...
                'FiniteHorizonController:DoControlSequenceComputation:InvalidMode', ...
                '** Input parameter <mode> (previous plant mode/mode estimate) must be in {1, ... %d} **', ...
                this.sequenceLength + 1);
                        
            [stateMean, ~] = state.getMeanAndCov();
            
            this.sysState(1:this.dimPlantState) = stateMean(:);
            inputSequence = this.L(: , :, mode, timestep) * this.sysState;
            if this.constraintsPresent
                if timestep == 1
                    % laziness for the win: solve optimization problem for
                    % the multiplier just now
                    % then, compute true feedforward terms
                    this.multiplier = this.computeMultiplier(mode, this.sysState);
                end
                inputSequence = inputSequence + this.feedforward(:, :, mode, timestep) * this.multiplier;
            end
            
            this.sysState(this.dimPlantState + 1:end) = ...
                this.F * this.sysState(this.dimPlantState + 1:end) + this.G * inputSequence;
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, timestep)
            assert(Checks.isScalarIn(timestep, 1, this.horizonLength) && mod(timestep, 1) == 0, ...
                'FiniteHorizonController:DoStageCostsComputation:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be in {1, ... %d} **', ...
                this.horizonLength);
            
            stageCosts = Utility.computeStageCosts(state, input, this.Q, this.R);
        end
        
        %% doCostsComputation
        function lQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            assert(size(stateTrajectory, 2) == this.horizonLength + 1, ...
                'FiniteHorizonController:DoCostsComputation:InvalidStateTrajectory', ...
                '** <stateTrajectory> is expected to have %d columns ', this.horizonLength + 1);
            assert(size(appliedInputs, 2) == this.horizonLength, ...
                'FiniteHorizonController:DoCostsComputation:InvalidInputTrajectory', ...
                '** <appliedInputs> is expected to have %d columns ', this.horizonLength);
            
            lQGCosts = Utility.computeLQGCosts(this.horizonLength, stateTrajectory, appliedInputs, this.Q, this.R);
        end
    end
    
    methods (Access = private)
        
        %% validateLinearIntegralConstraints
        function validateLinearIntegralConstraints(this, stateWeightings, inputWeightings, constraints)
            assert(Checks.isMat3D(stateWeightings) && size(stateWeightings, 1) == this.dimPlantState ...
                    && size(stateWeightings, 2) == this.horizonLength + 1, ...
                'FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidStateWeightings', ...
                '** Input parameter <stateWeightings> must be 3d matrix where each slice is %d-by-%d **', ...
                this.dimPlantState, this.horizonLength + 1);
            
            numConstraints  = size(stateWeightings, 3);
            assert(Checks.isMat3D(inputWeightings, this.dimPlantInput, this.horizonLength, numConstraints), ...
                'FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidInputWeightings', ...
                '** Input parameter <inputWeightings> must be 3d matrix with %d slices, where each slice is %d-by-%d **', ...
                numConstraints, this.dimPlantInput, this.horizonLength);
            assert(Checks.isVec(constraints, numConstraints), ...
                'FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidConstraints', ...
                '** Input parameter <contstraints> must be vector with %d elements **', ...
                numConstraints);            
        end
        
        %% computeGainMatrices
        function L = computeGainMatrices(this, augA, augB, augQ, augR, terminalK)
            numModes = this.sequenceLength + 1;
            dim = size(augR, 1);
            
            % we only need K_{k+1} to compute K_k and L_k
            K_k = repmat(terminalK, 1, 1, numModes); % K_N is the same for all modes
            L = zeros(dim, this.dimState, numModes, this.horizonLength);
            for k=this.horizonLength:-1:1
                K_prev = K_k; % K_{k+1}
                KA = mtimesx(K_prev, augA);
                AKA = mtimesx(augA, 'T', KA);
                QAKA = augQ + (AKA + permute(AKA, [2 1 3])) / 2; % ensure symmetry
                BKA = mtimesx(augB, 'T', KA);
                BKB = mtimesx(augB, 'T', mtimesx(K_prev, augB));
                RBKB = augR + (BKB + permute(BKB, [2 1 3])) / 2; % ensure symmetry
                for j=1:numModes
                    P1 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* QAKA, 3);
                    P2 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* BKA, 3);
                    P3 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* RBKB, 3);
                    P4 = pinv(P3) * P2;
                    P5 = P2' * P4;
                    
                    K_k(:,:, j) = P1 - (P5 + P5')/2; % K_k for mode j, ensure that second part is symmetric
                    L(:,:, j, k) = -P4;
                end
            end
        end
%%
%%
        %% computeGainMatricesWithConstraints
        function [L, feedforward, P_0, S_0] = computeGainMatricesWithConstraints(this, augA, augB, augQ, augR, ...
                inputWeightings, stateWeightings, terminalK, terminalP)
            
            numConstraints = size(stateWeightings, 3);
            numModes = this.sequenceLength + 1;
            dim = size(augR, 1);
            
            % we only need S_0 in the end for the computation of the multiplier, so just store S_k and S_{k+1}
            % per mode
            S_k = zeros(numConstraints, numConstraints, numModes);
            % for the computation of L_k and the feedforward part, we
            % require K_{k+1}, so it suffices to store K_k_ and K_{k+1}
            K_k = repmat(terminalK, 1, 1, numModes); % P_N in the paper is the same for all modes
            % we only need P_0 in the end for the computation of the multiplier, so just store P_k and P_{k+1}
            % per mode
            P_k = repmat(terminalP, 1, 1, numModes); % P_tilde_N in the paper is the same for all modes
           
            L = zeros(dim, this.dimState, numModes, this.horizonLength);
            % this one requires a lot of memory
            feedforward = zeros(dim, numConstraints, numModes, this.horizonLength);

            % the follwing three lines are required for the mode-dependent
            % Q_tilde (to avoid explicit calculation of H'*b)
            % cf. Utility.createAugmentedLinearIntegralConstraints
            sums = this.dimPlantInput * cumsum(1:this.sequenceLength -1);
            startIdx = arrayfun(@(i) this.dimState - sums(this.sequenceLength - i + 1) + 1, 2:numModes-1);
            endIdx = startIdx + this.dimPlantInput - 1;
            
            for k = this.horizonLength:-1:1
                S_prev = S_k;
                K_prev = K_k;
                P_prev = P_k;
                
                % R_tilde (augmented input constraint weightings) only nonzero for the first mode
                R_tilde_mode1 = zeros(this.dimPlantInput * this.sequenceLength, numConstraints);
                R_tilde_mode1(1:this.dimPlantInput, :) = inputWeightings(:, k, :);
                %Q_tilde 
                Q_tilde = zeros(this.dimState, numConstraints, numModes);
                Q_tilde(1:this.dimPlantState, :, :) = repmat(squeeze(stateWeightings(:, k, :)), 1, 1, numModes);
                
                % Q_tilde (per mode) depends on H matrix which is
                % zero for first and last mode
                % inline code (with index shift) from Utility.createAugmentedLinearIntegralConstraints
                % to compute the mode-dependent Q_tilde
                for i=2:numModes-1
                    Q_tilde(startIdx(i-1):endIdx(i-1), :, i) = inputWeightings(:, k, :);
                end     
                
                BKB = mtimesx(mtimesx(augB, 'T', K_prev), augB); % shall be symmetric as K_prev is
                RBKB = augR + (BKB + permute(BKB, [2 1 3])) / 2;
                AK = mtimesx(augA, 'T', K_prev);
                AKA = mtimesx(AK, augA); % shall be symmetric
                QAKA = augQ + (AKA + permute(AKA, [2 1 3])) / 2;
                AKB = mtimesx(AK, augB);
                
                APQ_tilde = mtimesx(augA, 'T', P_prev) + Q_tilde;
                % R_tilde (augmented input constraint weightings) only nonzero for the first mode
                BPR_tilde = mtimesx(augB, 'T', P_prev);
                BPR_tilde(:, :, 1) = BPR_tilde(:, :, 1) + R_tilde_mode1; 
                
                for j=1:numModes
                    P1 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* QAKA, 3);
                    P2 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* AKB, 3);
                    P3 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* RBKB, 3);
                    
                    P4 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* APQ_tilde, 3);
                    P5 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* BPR_tilde, 3);
                    
                    S1 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* S_prev, 3);

                    P6 = pinv(P3);
                    P7 = P6 * P5;
                    P8 = P2 * P6 * P2';
                    P9 = P5' * P7;
  
                    K_k(:,:, j) = P1 - (P8 + P8') / 2; % K_k for mode j, shall be symmetric
                    P_k(:, :, j) = P4 - P2 * P7; % P_tilde for mode j
                    S_k(:, :, j) = S1 - (P9 + P9') / 2; % shall be symmetric
                    L(:,:, j, k) = -P6  * P2';
                    feedforward(:, :, j, k) = -P7;
                end                
            end
            S_0 = S_k;
            P_0 = P_k; % P_tilde in the paper
        end
             
        %% computeMultiplier
        function lambda_opt = computeMultiplier(this, initialMode, initialState)
            % directly use matlab's quadprog function (or MOSEK's variant) with interior point
            % algorithm: requires objective to be 1/2x'Hx+f'x
            H = -this.S_0(:, :, initialMode);
            f = this.constraintBounds(:) - transpose(initialState' * this.P_0(:, :, initialMode));
            
            lambda_opt = quadprog(H, f, [], [], [], [], zeros(numel(this.constraintBounds), 1), []);
        end
         
    end
end

