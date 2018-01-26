classdef FiniteHorizonTrackingController < SequenceBasedController
    % Implementation of the optimal finite horizon linear sequence-based LQG tracking controller for
    % NCS with a TCP-like network connecting the controller and the
    % actuator.
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
        feedforward; % the feedforward constants ff_k, such that U_k = L*x - ff_k
        % dimension of the augmented system state; 
        % (positive integer <dimX+dimU*sequenceLength*(sequenceLength-1)/2>)
        dimState = -1;
        Z;
        dimRef;
    end

    properties (SetAccess = private, GetAccess = protected)
        % augmented state space
        % system state of the augmented state space;
        % (column vector of dimension <dimState>)
        sysState = [];
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        refTrajectory;
        horizonLength;
    end
    
    methods (Access = public)
        %% FiniteHorizonTrackingController
        function this = FiniteHorizonTrackingController(A, B, Q, R, Z, delayProb, ...
                sequenceLength, horizonLength, referenceTrajectory)
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
            % Returns:
            %   << obj (FiniteHorizonTrackingController)
            %      A new FiniteHorizonTrackingController instance.
            %
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength);
            
            Validator.validateHorizonLength(horizonLength);
            this.horizonLength = horizonLength;
            
            if ~Checks.isFixedColMat(Z, dimX) || any(~isfinite(Z(:)))
                error('FiniteHorizonTrackingController:InvalidZMatrix', ...
                    '** Input parameter <Z> (Plant output/performance maxtrix) must be a real-valued matrix with %d cols **', ...
                    dimX); 
            end
            this.dimRef = size(Z, 1);
            this.Z = Z;
            if ~Checks.isMat(referenceTrajectory, this.dimRef, horizonLength + 1) || any(~isfinite(referenceTrajectory(:)))
                error('FiniteHorizonTrackingController:InvalidReferenceTrajectory', ...
                    '** Reference trajectory <referenceTrajectory> must be a real-valued %d-by-%d matrix **',...
                    this.dimRef, horizonLength + 1);
            end
            this.refTrajectory = referenceTrajectory;
            
            % Q, R
            Validator.validateCostMatrices(Q, R, this.dimRef, dimU);
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
        
        %% getDeviationFromRefForState
        function deviation = getDeviationFromRefForState(this, trueState, timestep)
            % Compute deviation of the performance output for a given true
            % state from the reference at a given time step.
            %
            % Parameters:
            %   >> trueState (Vector dimension dimPlantState)
            %      The plant true state at the given time step
            %
            %   >> timestep (Positive integer)
            %      The time step for which to retrieve the deviation from the
            %      reference trajectory.
            %
            % Returns:
            %   << costs (Nonnegative scalar)
            %      The deviation of the performance output of the given
            %      true state and the reference, i.e., Z * x_true -z_ref.
            %
            if ~Checks.isVec(trueState, this.dimPlantState)
                 error('FiniteHorizonTrackingController:GetDeviationFromRefForState:InvalidTrueState', ...
                   ['** Input parameter <trueState>  must be ' ...
                   'a %d-dimensional vector **'], this.dimPlantState);
            end
            if ~Checks.isScalarIn(timestep, 1, this.horizonLength) || mod(timestep, 1) ~= 0
              error('FiniteHorizonTrackingController:GetDeviationFromRefForState:InvalidTimestep', ...
                  '** Input parameter <timestep> must be in {1, ... %d} **', ...
                  this.horizonLength);
            end
            deviation = this.Z * trueState(:) - this.refTrajectory(:, timestep);
        end
    end
     
    methods (Access = protected)
         %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, mode, timestep)
            if ~Checks.isScalarIn(timestep, 1, this.horizonLength) || mod(timestep, 1) ~= 0
              error('FiniteHorizonTrackingController:DoControlSequenceComputation:InvalidTimestep', ...
                  '** Input parameter <timestep> (current time step) must be in {1, ... %d} **', ...
                  this.horizonLength);
            end
            if ~Checks.isScalarIn(mode, 1, this.sequenceLength + 1) || mod(mode, 1) ~= 0
              error('FiniteHorizonTrackingController:DoControlSequenceComputation:InvalidMode', ...
                  '** Input parameter <mode> (previous plant mode/mode estimate) must be in {1, ... %d} **', ...
                  this.sequenceLength + 1);
            end
            
            [stateMean, ~] = state.getMeanAndCovariance();
  
            this.sysState(1:this.dimPlantState) = stateMean(:);
            inputSequence = this.L(: , :, mode, timestep) * this.sysState - this.feedforward(:, mode, timestep);
            this.sysState(this.dimPlantState + 1:end) = ...
                this.F * this.sysState(this.dimPlantState + 1:end) + this.G * inputSequence;
        end
        
        function lQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            if size(stateTrajectory, 2) ~= this.horizonLength + 1
                error('FiniteHorizonTrackingController:DoCostsComputation:InvalidStateTrajectory', ...
                    '** <stateTrajectory> is expected to have %d columns ', this.horizonLength + 1);
            end
            if size(appliedInputs, 2) ~= this.horizonLength
                 error('FiniteHorizonTrackingController:DoCostsComputation:InvalidInputTrajectory', ...
                    '** <appliedInputs> is expected to have %d columns ', this.horizonLength);
            end
            % compute performance output and difference to reference
            % trajectory
            diff = bsxfun(@minus, this.Z * stateTrajectory, this.refTrajectory);
            lQGCosts = Utility.computeLQGCosts(this.horizonLength, diff, appliedInputs, this.Q, this.R);
        end
    end
     
    methods (Access = private)
        %% computeControlLaws
        function [L, feedforward] = computeControlLaws(this, augA, augB, augQ, augR)
            numModes = this.sequenceLength + 1;
            K = zeros(this.dimState, this.dimState, numModes, this.horizonLength + 1);
            dimZero = this.dimState - this.dimPlantState;
            %K(:,:, :, end) = repmat(blkdiag(this.Z' * this.Q * this.Z, zeros(dimZero)), 1, 1, numModes);
            stateWeighting = augQ(1:this.dimPlantState, 1:this.dimPlantState, end);
            K(:,:, :, end) = repmat(blkdiag(stateWeighting, zeros(dimZero)), 1, 1, numModes);
            
            % directly compute the transpose
            Qref = [(this.Q * this.Z)'; 
                    zeros(dimZero, this.dimRef)];
            refWeightings = Qref * this.refTrajectory;
            sigma = zeros(size(refWeightings, 1), numModes, this.horizonLength + 1);
            sigma(:, :, end) = repmat(refWeightings(:, end), 1, numModes);
            dim = size(augR, 1);
          
            % perform a Riccati-like recursion (backwards) for K and
            % another backward recursion for sigma
            for k = this.horizonLength:-1:1
                for j=1:numModes
                    P1 = zeros(this.dimState);
                    P2 = zeros(this.dimState, dim);
                    P3 = zeros(dim);
                    S1 = zeros(this.dimState, 1); 
                    S2 = zeros(dim, 1);
                    
                    for i =1:numModes
                        modeK = K(:, :, i, k + 1);
                        A_iK_i = augA(:, :, i)' * modeK;
                        p_ji = this.transitionMatrix(j, i);
                        P1 = P1 + p_ji * (augQ(:, :, i) + A_iK_i * augA(:, :, i));
                        P2 = P2 + p_ji * A_iK_i * augB(:, :, i);
                        P3 = P3 + p_ji * (augR(:, :, i) + augB(:, :, i)' * modeK * augB(:, :, i));
                        
                        modeSigma = sigma(:, i, k + 1);
                        S1 = S1 + p_ji * augA(:, :, i)' * modeSigma;
                        S2 = S2 + p_ji * augB(:, :, i)' * modeSigma;
                    end
                    K(:,:, j, k) = P1 - P2 * pinv(P3) * P2';
                    sigma(:, j, k) = refWeightings(:, k) + S1 - P2 * pinv(P3) * S2;
                end
            end
            % no compute feedback and feedforward law
            feedforward = zeros(dim, numModes, this.horizonLength);
            L = zeros(dim, this.dimState, numModes, this.horizonLength);
            for k = 1:this.horizonLength
                for j=1:numModes
                    L1 = zeros(dim);
                    L2 = zeros(dim, this.dimState);
                    S2 = zeros(dim, 1);
                    for i =1:numModes
                        B_iK_i = augB(:, :, i)' * K(:, :, i, k + 1);
                        p_ji = this.transitionMatrix(j, i);
                        
                        L1 = L1 + p_ji * (augR(:, :, i) + B_iK_i * augB(:, :, i));
                        L2 = L2 + p_ji * B_iK_i * augA(:, :, i);
                        S2 = S2 + p_ji * augB(:, :, i)' * sigma(:, i, k + 1);
                    end
                    L(:,:, j, k) = -pinv(L1) * L2;
                    % compute the negative feedforward offsets and
                    % substract from the feedback term
                    feedforward(:, j, k) = pinv(L1) * S2;
                end    
            end
        end 
     end
end

