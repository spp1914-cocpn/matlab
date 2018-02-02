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
    
    properties (SetAccess = immutable, GetAccess = protected)
        % system state cost matrix per time step;
        % (positive semidefinite matrix of dimension <dimX> x <dimX>)
        Q = [];
        % input cost matrix;
        % (positive definite matrix of dimension <dimU> x <dimU>)
        R = [];
        % packet delay probability density function of the controller-actuator-
        % link; (vector of dimension >= <sequenceLength>)
        delayProb = [];
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
    end
    
    properties (SetAccess = private, GetAccess = protected)
        % augmented state space
        % system state of the augmented state space;
        % (column vector of dimension <dimState>)
        sysState = [];
    end
    
    methods (Access = public)
        %% FiniteHorizonController
        function this = FiniteHorizonController(A, B, Q, R, delayProb, sequenceLength, horizonLength, ...
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
            
            if nargin ~= 7 && nargin ~= 10
                 error('FiniteHorizonController:InvalidNumberOfArguments', ...
                        '** Constructor must be called with either 7 or 10 arguments **');
            end
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength);
            
            Validator.validateHorizonLength(horizonLength);
            this.horizonLength = horizonLength;
            
             if nargin == 10
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
            elementCount = numel(delayProb);
            if elementCount <= sequenceLength
                % fill up with zeros (row vector)
                probHelp = [reshape(delayProb, 1, elementCount), zeros(1, sequenceLength + 1 - elementCount)];
            else
                % cut up the distribution and normalize
                probHelp = [delayProb(1:sequenceLength), 1 - sum(delayProb(1:sequenceLength))];
            end
            % store as row vector, i.e., column-wise arranged
            this.delayProb = probHelp;

            this.transitionMatrix = Utility.calculateDelayTransitionMatrix(this.delayProb); 
                        
            % augmented initial system state
            [this.dimState, this.F, this.G, augA, augB, augQ, augR] ...
                = Utility.performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R);
            this.sysState = zeros(this.dimState, 1);
         
            terminalK = blkdiag(this.Q, zeros(this.dimState - dimX)); % K_N for all modes (P_N in second paper mentioned above)
            
            if this.constraintsPresent
                terminalStateConstraintWeightings = squeeze(stateConstraintWeightings(:, end, :));
                numConstraints = numel(constraintBounds);
                                
                terminalP = zeros(this.dimState, numConstraints); % P_N for all modes (P_tilde in the paper)
                terminalP(1:dimX, :) = terminalStateConstraintWeightings;

                [this.L, this.feedforward, this.P_0, this.S_0] ...
                    = this.computeGainMatricesWithConstraints(augA, augB, augQ, augR, ...
                    inputConstraintWeightings, stateConstraintWeightings(:, 1:end-1, :), terminalK, terminalP);
            else
                this.L = this.computeGainMatrices(augA, augB, augQ, augR, terminalK);
            end
        end
        
        %% reset
        function reset(this)
           this.sysState = zeros(this.dimState,1);
        end % function reset
    end
    
    methods (Access = protected)
         %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, mode, timestep)
            if ~Checks.isScalarIn(timestep, 1, this.horizonLength) || mod(timestep, 1) ~= 0
              error('FiniteHorizonController:DoControlSequenceComputation:InvalidTimestep', ...
                  '** Input parameter <timestep> (current time step) must be in {1, ... %d} **', ...
                  this.horizonLength);
            end
            if ~Checks.isScalarIn(mode, 1, this.sequenceLength + 1) || mod(mode, 1) ~= 0
              error('FiniteHorizonController:DoControlSequenceComputation:InvalidMode', ...
                  '** Input parameter <mode> (previous plant mode/mode estimate) must be in {1, ... %d} **', ...
                  this.sequenceLength + 1);
            end
            
            [stateMean, ~] = state.getMeanAndCovariance();
            
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
        
        %% doCostsComputation
        function lQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            if size(stateTrajectory, 2) ~= this.horizonLength + 1
                error('FiniteHorizonController:DoCostsComputation:InvalidStateTrajectory', ...
                    '** <stateTrajectory> is expected to have %d columns ', this.horizonLength + 1);
            end
            if size(appliedInputs, 2) ~= this.horizonLength
                 error('FiniteHorizonController:DoCostsComputation:InvalidInputTrajectory', ...
                    '** <appliedInputs> is expected to have %d columns ', this.horizonLength);
            end
            lQGCosts = Utility.computeLQGCosts(this.horizonLength, stateTrajectory, appliedInputs, this.Q, this.R);
        end
    end
    
    methods (Access = private)
        
        %% validateLinearIntegralConstraints
        function validateLinearIntegralConstraints(this, stateWeightings, inputWeightings, constraints)
            if ~Checks.isMat3D(stateWeightings) || size(stateWeightings, 1) ~= this.dimPlantState ...
                        ||  size(stateWeightings, 2) ~= this.horizonLength + 1
                error('FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidStateWeightings', ...
                        '** Input parameter <stateWeightings> must be 3d matrix where each slice is %d-by-%d **', ...
                this.dimPlantState, this.horizonLength + 1);
            end
            numConstraints  = size(stateWeightings, 3);
            if ~Checks.isMat3D(inputWeightings, this.dimPlantInput, this.horizonLength, numConstraints)
                error('FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidInputWeightings', ...
                    '** Input parameter <inputWeightings> must be 3d matrix with %d slices, where each slice is %d-by-%d **', ...
                    numConstraints, this.dimPlantInput, this.horizonLength);
            end
            if ~Checks.isVec(constraints, numConstraints)
                 error('FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidConstraints', ...
                    '** Input parameter <contstraints> must be vector with %d elements **', ...
                    numConstraints);
            end
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
                for j=1:numModes
                    P1 = zeros(this.dimState);
                    P2 = zeros(dim, this.dimState);
                    P3 = zeros(dim);

                    for i=1:numModes
                        modeK = K_prev(:, :, i);
                        
                        p_ji = this.transitionMatrix(j, i);
                         
                        P1 = P1 + p_ji * (augQ(:, :, i) + augA(:, :, i)' * modeK * augA(:, :, i));
                        P2 = P2 + p_ji * augB(:, :, i)' * modeK * augA(:, :, i);
                        P3 = P3 + p_ji * (augR(:, :, i) + augB(:, :, i)' * modeK * augB(:, :, i));
                    end
                    K_k(:,:, j) = P1 - P2' * pinv(P3) * P2; % K_k for mode j
                    L(:,:, j, k) = -pinv(P3) * P2;
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
                
                for j=1:numModes
                    P1 = zeros(this.dimState);
                    P2 = zeros(this.dimState, dim);
                    P3 = zeros(dim);
                    P4 = zeros(this.dimState, numConstraints);
                    P5 = zeros(this.dimPlantInput * this.sequenceLength, numConstraints);
                    S1 = zeros(numConstraints);
                    
                    for i=1:numModes
                        modeK = K_prev(:, :, i);
                        modeP = P_prev(:, :, i);
                        
                        A_iK_i = augA(:, :, i)' * modeK;
                        p_ji = this.transitionMatrix(j, i);
                        P1 = P1 + p_ji * (augQ(:, :, i) + A_iK_i * augA(:, :, i));
                        P2 = P2 + p_ji * A_iK_i * augB(:, :, i);
                        P3 = P3 + p_ji * (augR(:, :, i) + augB(:, :, i)' * modeK * augB(:, :, i));
                      
                        Q_tilde_mode = Q_tilde(:, :, i);
                        % Q_tilde (per mode) depends on H matrix which is
                        % zero for first and last mode
                        % inline code (with index shift) from Utility.createAugmentedLinearIntegralConstraints
                        % to compute the mode-dependent Q_tilde
                        if i ~= numModes && i~= 1
                            Q_tilde_mode(startIdx(i-1):endIdx(i-1), :) = inputWeightings(:, k, :);
                        end
                        P4 = P4 + p_ji * (augA(:, :, i)' * modeP + Q_tilde_mode);
                        P5 = P5 + p_ji * augB(:, :, i)' * modeP; 
                        if i == 1
                            P5 = P5 + p_ji * R_tilde_mode1;
                        end
                        S1 = S1 + p_ji * S_prev(:, :, i);
                    end
                    P6 = pinv(P3);
                    P7 = P6 * P5;
                    K_k(:,:, j) = P1 - P2 * P6 * P2'; % K_k for mode j
                    P_k(:, :, j) = P4 - P2 * P7; % P_tilde for mode j
                    S_k(:, :, j) = S1 - P5' * P7;
                    L(:,:, j, k) = -P6  * P2';
                    feedforward(:, :, j, k) = -P7;
                end
                S_0 = S_k;
                P_0 = P_k;
            end
        end
             
        %% computeMultiplier
        function lambda_opt = computeMultiplier(this, initialMode, initialState)
            % requires MOSEK and yalmip
%             lambda = sdpvar(numel(this.constraintBounds), 1);
%             settings = sdpsettings('solver', 'mosek', 'verbose', 0);
%             objective = dot(this.constraintBounds, lambda) - lambda' * this.S_0(:, :, initialMode) * lambda ...
%                 - initialState' * this.P_0(:, :, initialMode) * lambda;
%             optimize(lambda >= 0, objective, settings);
%             lambda_opt = value(lambda)

            % directly use matlab's quadprog function (or MOSEK's variant) with interior point
            % algorithm: requires objective to be 1/2x'Hx+f'x
            H = -2 * this.S_0(:, :, initialMode);
            f = this.constraintBounds(:) - transpose(initialState' * this.P_0(:, :, initialMode));
            
            lambda_opt = quadprog(H, f, [], [], [], [], zeros(numel(this.constraintBounds), 1), []);
        end
         
    end
end

