classdef FiniteHorizonTrackingController < SequenceBasedTrackingController
    % Implementation of the optimal finite horizon linear sequence-based LQG tracking controller for
    % NCS with a TCP-like network connecting the controller and the
    % actuator.
    %
    % This implementation is based on the original one by Jörg Fischer and Maxim Dolgov.
    %
    % Literature: 
    %   Jörg Fischer, Maxim Dolgov, and Uwe D. Hanebeck,
    %   Optimal Sequence-Based Tracking Control over Unreliable Networks,
    %   Proceedings of the 19th IFAC World Congress (IFAC 2014),
    %   Cape Town, South Africa, August 2014.
    %
    %   Jörg Fischer,
    %   Optimal sequence-based control of networked linear systems,
    %   Karlsruhe series on intelligent sensor-actuator-systems, Volume 15,
    %   KIT Scientific Publishing, 2015.
    
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
        feedforward; % the feedforward constants ff_k, such that U_k = L*x - ff_k
        % dimension of the augmented system state; 
        % (positive integer <dimX+dimU*sequenceLength*(sequenceLength-1)/2>)
        dimState = -1;
    end

    properties (SetAccess = private, GetAccess = protected)
        % augmented state space
        % system state of the augmented state space;
        % (column vector of dimension <dimState>)
        sysState = [];
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        horizonLength;        
        
        useMexImplementation@logical=true; 
        % by default, we use the C++ (mex) implementation for computation of controller gains
        % this is faster, but can produce slightly different results
    end
    
    methods (Access = public)
        %% FiniteHorizonTrackingController
        function this = FiniteHorizonTrackingController(A, B, Q, R, Z, delayProb, ...
                sequenceLength, horizonLength, referenceTrajectory, useMexImplementation)
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
            %      The performance/plant output weighting matrix in the controller's underlying cost function.
            %
            %   >> R (Positive definite matrix)
            %      The input weighting matrix in the controller's underlying cost function.
            %
            %   >> Z (Matrix, n-by-dimPlantState)
            %      The time-invariant plant output (performance) matrix, 
            %      i.e., z_k = Z*x_k
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
            %   >> referenceTrajectory (Matrix, n-by-horizonLength+1)
            %      The reference trajectory to track, given as a matrix
            %      with the reference plant outputs column-wise arranged.
            %            
            %   >> useMexImplementation (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for the computation of the controller
            %      gains which is considerably faster than the Matlab
            %      implementation. 
            %      If left out, the default value true is used.
            %
            % Returns:
            %   << this (FiniteHorizonTrackingController)
            %      A new FiniteHorizonTrackingController instance.
            %
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            Validator.validateHorizonLength(horizonLength);
            this = this@SequenceBasedTrackingController(dimX, dimU, sequenceLength, true, Z, referenceTrajectory, horizonLength + 1);
                        
            this.horizonLength = horizonLength;
                       
            % Q, R
            Validator.validateCostMatrices(Q, R, this.dimRef, dimU);
            this.Q = Q;
            this.R = R;
            
            % delayProb
            Validator.validateDiscreteProbabilityDistribution(delayProb);
            this.transitionMatrix = Utility.calculateDelayTransitionMatrix( ...
                Utility.truncateDiscreteProbabilityDistribution(delayProb, sequenceLength + 1)); 
            
            if nargin > 9
                assert(Checks.isFlag(useMexImplementation), ...
                    'FiniteHorizonTrackingController:InvalidUseMexFlag', ...
                    '** <useMexImplementation> must be a flag **');
                this.useMexImplementation = useMexImplementation;
            end
            
            % augment the model
            [this.dimState, this.F, this.G, augA, augB, augQ, augR] ...
                = Utility.performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R, Z);
            this.sysState = zeros(this.dimState, 1);

            [this.L, this.feedforward] = this.computeControlLaws(augA, augB, augQ, augR);
        end
        
        %% reset
        function reset(this)
           this.sysState = zeros(this.dimState,1);
        end % function reset
        
        %% setEtaState
        function setEtaState(this, newEta)
            assert(Checks.isVec(newEta, this.dimState - this.dimPlantState) && all(isfinite(newEta)), ...
                'FiniteHorizonTrackingController:SetEtaState:InvalidEta', ...
                '** <newEta> must be a %d-dimensional vector **', ...
                this.dimState - this.dimPlantState);
            
            this.sysState(this.dimPlantState + 1:end) = newEta(:);
        end
        
    end
     
    methods (Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, mode, timestep)
            assert(Checks.isScalarIn(timestep, 1, this.horizonLength) && mod(timestep, 1) == 0, ...
                'FiniteHorizonTrackingController:DoControlSequenceComputation:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be in {1, ... %d} **', ...
                this.horizonLength);
            assert(Checks.isScalarIn(mode, 1, this.sequenceLength + 1) && mod(mode, 1) == 0, ...
                'FiniteHorizonTrackingController:DoControlSequenceComputation:InvalidMode', ...
                '** Input parameter <mode> (previous plant mode/mode estimate) must be in {1, ... %d} **', ...
                this.sequenceLength + 1);
             
            [stateMean, ~] = state.getMeanAndCov();
  
            this.sysState(1:this.dimPlantState) = stateMean(:);
            inputSequence = this.L(: , :, mode, timestep) * this.sysState - this.feedforward(:, mode, timestep);
            this.sysState(this.dimPlantState + 1:end) = ...
                this.F * this.sysState(this.dimPlantState + 1:end) + this.G * inputSequence;
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, timestep)
            assert(Checks.isScalarIn(timestep, 1, this.horizonLength) && mod(timestep, 1) == 0, ...
                'FiniteHorizonTrackingController:DoStageCostsComputation:InvalidTimestep', ...
                '** Input parameter <timestep> (current time step) must be in {1, ... %d} **', ...
                this.horizonLength);
            
            performance = this.Z * state - this.refTrajectory(:, timestep);
            stageCosts = Utility.computeStageCosts(performance, input, this.Q, this.R);
        end
        
        function lQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            assert(size(stateTrajectory, 2) == this.horizonLength + 1, ...
                'FiniteHorizonTrackingController:DoCostsComputation:InvalidStateTrajectory', ...
                '** <stateTrajectory> is expected to have %d columns ', this.horizonLength + 1);
            assert(size(appliedInputs, 2) == this.horizonLength, ...
                'FiniteHorizonTrackingController:DoCostsComputation:InvalidInputTrajectory', ...
                '** <appliedInputs> is expected to have %d columns ', this.horizonLength);

            % compute performance output and difference to reference
            % trajectory
            diff = bsxfun(@minus, this.Z * stateTrajectory, this.refTrajectory);
            lQGCosts = Utility.computeLQGCosts(this.horizonLength, diff, appliedInputs, this.Q, this.R);
        end
        
        %% doGetDeviationFromRefForState
        function deviation = doGetDeviationFromRefForState(this, state, timestep)
            assert(Checks.isScalarIn(timestep, 1, this.horizonLength) && mod(timestep, 1) == 0, ...
                'FiniteHorizonTrackingController:GetDeviationFromRefForState:InvalidTimestep', ...
                '** Input parameter <timestep> must be in {1, ... %d} **', ...
                this.horizonLength);
            
            deviation = this.Z * state - this.refTrajectory(:, timestep);
        end
    end
     
    methods (Access = private)
        %% computeControlLaws
        function [L, feedforward] = computeControlLaws(this, augA, augB, augQ, augR)
            numModes = this.sequenceLength + 1;
            dimZero = this.dimState - this.dimPlantState;
            dimR = size(augR, 1);
            stateWeighting = augQ(1:this.dimPlantState, 1:this.dimPlantState, end); % this Z' * Q * Z
            terminalK = blkdiag(stateWeighting, zeros(dimZero));
            % directly compute the transpose
            Qref = [(this.Q * this.Z)'; 
                    zeros(dimZero, this.dimRef)];
            refWeightings = Qref * this.refTrajectory;
            
            if this.useMexImplementation
                [L, feedforward] = mex_FiniteHorizonTrackingController(augA, augB, augQ, augR, ...
                    this.transitionMatrix, terminalK, this.horizonLength, Qref, refWeightings);
            else
                % we only need K_{k+1} to compute K_k and L_k and the
                % feedforward
                K = repmat(terminalK, 1, 1, numModes); % symmetric matrix

                % the same for sigma: only sigma_{k+1} required
                sigma = repmat(refWeightings(:, end), 1, 1, numModes);

                feedforward = zeros(dimR, numModes, this.horizonLength);
                L = zeros(dimR, this.dimState, numModes, this.horizonLength);
                for k = this.horizonLength:-1:1
                    K_prev = K;
                    sigma_prev = sigma;
                    AK = mtimesx(augA, 'T', K_prev);
                    AKA = mtimesx(AK, augA); % this one should be symmetric
                    QAKA = augQ + (AKA + permute(AKA, [2 1 3])) / 2; % ensure symmetry
                    AKB = mtimesx(AK, augB);
                    BKB = mtimesx(augB, 'T', mtimesx(K_prev, augB)); % this one should be symmetric
                    RBKB = augR + (BKB + permute(BKB, [2 1 3])) / 2; % ensure symmetry

                    for j=1:numModes                      
                        P1 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* QAKA, 3);
                        P2 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* AKB, 3);
                        P3 = pinv(sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* RBKB, 3));
                        S1 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* mtimesx(augA, 'T', sigma_prev), 3);
                        S2 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* mtimesx(augB, 'T', sigma_prev), 3);

                        P4 = P2 * P3 * P2'; % should be symmetric

                        K(:, :, j) = P1 - (P4 + P4') / 2; % ensure symmetry
                        sigma(:, 1, j) = refWeightings(:, k) + S1 - P2 * P3 * S2;

                        L(:,:, j, k) = -P3 * P2';
                        feedforward(:, j, k) = P3 * S2;
                    end
                end
            end
        end
     end
end

