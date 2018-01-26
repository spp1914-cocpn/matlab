classdef EventTriggeredInfiniteHorizonController < SequenceBasedController
    % Implementation of the infinite horizon event-triggered linear sequence-based LQG controller for
    % NCS with a TCP-like network connecting the controller and the actuator.
    %
    % Literature: 
    %   Maxim Dolgov, Jörg Fischer, and Uwe D. Hanebeck,
    %   Event-based LQG Control over Networks Subject to Random Transmission Delays and Packet Losses,
    %   Proceedings of the 4th IFAC Workshop on Distributed Estimation and Control in Networked Systems (NecSys 2013),
    %   Koblenz, Germany, September 2013.
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
    
    properties (SetAccess = immutable, GetAccess = private)
        %% input properties
        % system state cost matrix;
        % (positive semidefinite matrix of dimension <dimX> x <dimX>)
        Q = [];
        % input cost matrix;
        % (positive definite matrix of dimension <dimU> x <dimU>)
        R = [];
        % packet delay probability density function of the controller-actuator-
        % link; (vector of dimension >= <sequenceLength>)
        delayProb = [];
        %% derived properties
        % CA-network-actuator system
        % mapping for propagation of possible control inputs;
        % (matrix of dimension <dimU*sequenceLength*(sequenceLength-1)/2> x 
        %  <dimU*packetLength*(sequenceLength-1)/2>)
        F = [];
        % matrix mapping control sequence to the possible control inputs;
        % (matrix of dimension <dimU*packetLength*(sequenceLength-1)/2> x
        %  <dimU*sequenceLength>)
        G = [];
        % 3D matrix with history-dependent transition matrices of the
        % Markov chain
        % (matrix of dimension <sequenceLength+1> x <sequenceLength+1> x transmissionCombinations)
        transitionMatrices = [];
        % 4D matrix with controller gains; 
        %  (matrix of dimension <dimU*sequenceLength> x <dimState> x <numModes> x transmissionCombinations)
        L = [];
        % 4D matrix with control covariances; 
        %  (matrix of dimension 
        %   <dimState> x <dimState> x <numModes> x transmissionCombinations)
        P;
        % dimension of the augmented system state; 
        % (positive integer <dimX+dimU*sequenceLength*(sequenceLength-1)/2>)
        dimState = -1;
        % time invariant transmission costs
        transmissionCosts;
       
    end
    
    properties (Access = private)
        % augmented state space
        % system state of the augmented state space;
        % (column vector of dimension <dimState>)
        sysState = []; 
        % bit vector to track the history (1 = sent, 0 = not sent)
        transmissionHistory;
    end
    
    methods (Access = public)
        %% EventTriggeredInfiniteHorizonController
        function this = EventTriggeredInfiniteHorizonController(A, B, Q, R, delayProb, packetLength, transmissionCosts)
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, packetLength);
            
            % Check for stabilizablity of A by B
            if dimX - rank(ctrb(A, B)) > 0
              error('** Plant (A, B) has uncontrollable eigenvalues **');
            end

            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
            
            if ~Checks.isNonNegativeScalar(transmissionCosts)
                 error('EventBasedInfiniteHorizonController:InvalidTransmissionCosts', ...
                    ['** Input parameter <transmissionCosts> ',...
                     'must be a nonegative scalar **']);
            end
            this.transmissionCosts = transmissionCosts;
            
             % delayProb
            Validator.validateDiscreteProbabilityDistribution(delayProb);
            elementCount = numel(delayProb);
            if elementCount <= packetLength
                % fill up with zeros (row vector)
                probHelp = [reshape(delayProb, 1, elementCount), zeros(1, packetLength + 1 - elementCount)];
            else
                % cut up the distribution and normalize
                probHelp = [delayProb(1:packetLength), 1 - sum(delayProb(1:packetLength))];
            end
            % store as row vector, i.e., column-wise arranged
            this.delayProb = probHelp;
            
            this.transitionMatrices = this.computeAllTransitionMatrices();
            this.transmissionHistory = 2 ^ (packetLength + 1); % all sent until now
            
            % augmented initial system state
            [this.dimState, this.F, this.G, augA, augB, augQ, augR] ...
                = Utility.performModelAugmentation(packetLength, dimX, dimU, A, B, Q, R);
            this.sysState = zeros(this.dimState, 1);
            
             % make controller covariance matrices and check result
            this.P = this.makeSteadyStateControlCovarianceMatrices(augA, augB, augQ, augR);
            this.L = this.makeSteadyStateControlGainMatrices(augA, augB, augR);
        end
        
        %% reset
        function reset(this)
           this.sysState = zeros(this.dimState,1);
           this.transmissionHistory = 1; % none sent until now
        end
    end
    
    methods (Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, mode, ~)
            if ~Checks.isScalarIn(mode, 1, this.sequenceLength + 1) || mod(mode, 1) ~= 0
              error('EventBasedInfiniteHorizonController:DoControlSequenceComputation', ...
                  '** Input parameter <mode> (previous plant mode/mode estimate) must be in {1, ... %d} **', ...
                  this.sequenceLength + 1);
            end
            [stateMean, stateCovariance] = state.getMeanAndCovariance();
            this.sysState(1:this.dimPlantState) = stateMean(:);
                      
            shifted = dec2bin(bitshift(this.transmissionHistory - 1, -1), this.sequenceLength);
            futureHistorySend = bin2dec(['1', shifted]) + 1;
            futureHistoryNotSend = bin2dec(['0', shifted]) + 1;
            
            augmentedStateSecondMoment = this.sysState * this.sysState'  ...
                + blkdiag(stateCovariance, zeros(this.dimState - this.dimPlantState));
            sumSend = zeros(this.dimState);
            sumNotSend = zeros(this.dimState);
           
            for i = 1:this.sequenceLength + 1
                sumSend = sumSend + this.transitionMatrices(mode, i, futureHistorySend) * this.P(:, :, i, futureHistorySend);
                sumNotSend = sumNotSend + this.transitionMatrices(mode, i, futureHistoryNotSend) * this.P(:, :, i, futureHistoryNotSend);
            end
            
            if this.transmissionCosts + trace((sumSend - sumNotSend) * augmentedStateSecondMoment) > 0
                % do not send
                this.transmissionHistory = futureHistoryNotSend;
                inputSequence = [];
                this.sysState(this.dimPlantState + 1:end) = ...
                    this.F * this.sysState(this.dimPlantState + 1:end);
            else
                % send
                this.transmissionHistory = futureHistorySend;
                inputSequence = this.L(: , :, mode, futureHistorySend) * this.sysState;
                this.sysState(this.dimPlantState + 1:end) = ...
                    this.F * this.sysState(this.dimPlantState + 1:end) + this.G * inputSequence;
            end
        end
        
        function averageLQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonLength = size(appliedInputs, 2);
            if size(stateTrajectory, 2) ~= horizonLength + 1
                error('InfiniteHorizonController:DoCostsComputation', ...
                    '** <stateTrajectory> is expected to have %d columns ', horizonLength + 1);
            end
            averageLQGCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, appliedInputs, this.Q, this.R) / horizonLength;
        end
    end
    
    methods (Access = private)
        
        function T = computeAllTransitionMatrices(this)
            numModes = this.sequenceLength + 1;
            numTransmissionCombinations = 2 ^ numModes;
            
            T = zeros(numModes, numModes, numTransmissionCombinations);
            for h = 1:numTransmissionCombinations % yields all possible combinations
                history = dec2bin(h-1, numModes);
                
                for i = 1:this.sequenceLength % all but last row
                    for j=1:i
                        T(i,j,h) = bin2dec(history(i+1)) * bin2dec(history(j)) * this.delayProb(j) / (1 - sum(this.delayProb(1:j-1)));
                        for t = 1:j-1
                            T(i,j,h) = T(i,j,h) * (1 - bin2dec(history(t)) ...
                                + bin2dec(history(t)) * (1-sum(this.delayProb(1:t))) / (1-sum(this.delayProb(1:t-1))));
                        end 
                    end
                    T(i,i+1,h) = bin2dec(history(i+1));
                    for t = 1:i
                        T(i,i+1,h) = T(i,i+1,h) * (1 - bin2dec(history(t)) ...
                            + bin2dec(history(t))*(1-sum(this.delayProb(1:t))) / (1-sum(this.delayProb(1:t-1))));
                    end 
                end
                
                for j = 1:this.sequenceLength % last row
                    T(numModes,j,h) = bin2dec(history(j)) * this.delayProb(j)/(1-sum(this.delayProb(1:j-1)));
                    for t = 1:j-1
                        T(numModes,j,h) = T(numModes,j,h) * (1 - bin2dec(history(t)) ...
                            + bin2dec(history(t)) * (1-sum(this.delayProb(1:t))) / (1-sum(this.delayProb(1:t-1))));
                    end 
                end
                
                T(numModes, numModes, h) = 1;
                for t = 1:this.sequenceLength
                    T(numModes, numModes, h) = T(numModes, numModes, h) * (1 - bin2dec(history(t)) ...
                        + bin2dec(history(t)) * (1-sum(this.delayProb(1:t))) / (1-sum(this.delayProb(1:t-1))));
                end
            end
        end
        
        %% makeSteadyStateControlGainMatrices
        function L = makeSteadyStateControlGainMatrices(this, augA, augB, augR)
            numModes = this.sequenceLength + 1;
            numTransmissionCombinations = 2 ^ numModes;
            L = zeros(this.dimPlantInput * this.sequenceLength, this.dimState, numModes, numTransmissionCombinations);
            sizeL = size(L);
            for h = 1:numTransmissionCombinations
                for j = 1:numModes
                    L1 = zeros(this.dimPlantInput * this.sequenceLength);
                    L2 = zeros(sizeL(1:2));
                    for i = 1:numModes
                        B_iP_i = augB(:, :, i)' * this.P(:, :, i, h);
                        L1 = L1 + this.transitionMatrices(j, i, h)* (augR(:, :, i) + B_iP_i * augB(:, :, i));
                        L2 = L2 + this.transitionMatrices(j, i, h) * (B_iP_i * augA(:,:,i));
                    end
                    L(:,:,j, h) = -pinv(L1) * L2;
                end
            end
        end
        
        %% makeSteadyStateControlCovarianceMatrices
        function [P, status] = makeSteadyStateControlCovarianceMatrices(this, augA, augB, augQ, augR)
            numModes = this.sequenceLength + 1;
            numTransmissionCombinations = 2 ^ numModes;
            P = zeros(this.dimState, this.dimState, numModes, numTransmissionCombinations);
  
            previousP = P;

            convergeDiffmin = InfiniteHorizonController.maxIterationDiff;
    
            % initialize P
            for i = 1:numModes
                for h = 1:numTransmissionCombinations
                    P(1:this.dimPlantState, 1:this.dimPlantState, i, h) = this.Q;
                end
            end 
    
            counter = 0; % counts the number of performed iterations
    
            while(1)
                counter = counter + 1;
                % variable to monitor difference between two iterations of currP:
                % we use sum sum of the difference as currP is known that PreviousP >
                % currP in the positive definite sense. Therefore the quadratic form
                % x' currP x with x = [1 1 .... 1 1]'. currP is sum(sum(P)).
                % Therefore, this is an indicator for the change in positive-definite
                % sense. 
                convergeDiff = sum(sum(P(:,:,:,numTransmissionCombinations) - previousP(:,:,:,numTransmissionCombinations)));
                % time region  
                if counter > numModes + 1  
                    % ------- Terminate Conditions ------------------
                    % Terminate Condition 1: stable / converged:
                    % this happens if currP does not change anymore between two
                    % iterations. Instead of zero, we use a lower bound
                    if max(abs(convergeDiff)) < InfiniteHorizonController.minIterationDiff
                        status = 1;
                        break;
                    % Terminate Condition 2: unstable / diverged:
                    % this happens if the closed-loop system is unstable ( controller
                    % not able to stabilize system) or we had numerical problems that
                    % prevented to stop at the converged currP. The cases are hard to
                    % distinguish so we return the best solution found and the caution
                    % flag set.
                    elseif max(abs(convergeDiff)) > InfiniteHorizonController.maxIterationDiff 
                        error('NOT converged or convergence not detected');
                        counter %#ok
                        convergeCounterMin 
                        P = Pmin;
                        status = 0;
                        break;
                    % Terminate Condition 3: Numerical Problems
                    % This should not happen but if it does we had numerical problems.
                    % The best solution found will be returned.
                    elseif min(min(convergeDiff)) < 0
                        counter %#ok<NOPRT>
                        convergeCounterMin %#ok<NOPRT>
                        P = Pmin;
                        status = -1;
                        break;
                    end
        
                    % save best solution obtained yet. In case we diverge due to
                    % numerical problems, we can still use this best solution
                    if sum(abs(convergeDiff)) < convergeDiffmin 
                        convergeCounterMin = counter;
                        convergeDiffmin = convergeDiff;
                        Pmin = P;
                    end
                end

                previousP = P;
                P = zeros(this.dimState, this.dimState, numModes, numTransmissionCombinations);
      
                for h = 1 : numTransmissionCombinations
                    for j = 1 : numModes
                        P1 = 0;
                        P2 = 0;
                        P3 = 0;
                        for i = 1 : numModes
                            hFuture = bin2dec(['1', dec2bin(bitshift(h-1,-1), this.sequenceLength-1)]) + 1;
                            T = this.transitionMatrices(:,:,hFuture);
                            A_iP_i = augA(:,:,i)' * previousP(:,:,i,hFuture);
                            P1 = P1 + T(j,i) * (augQ(:,:,i) + A_iP_i * augA(:,:,i));
                            P2 = P2 + T(j,i) * A_iP_i * augB(:,:,i);
                            P3 = P3 + T(j,i) * (augR(:,:,i) + augB(:,:,i)'* previousP(:,:,i,hFuture) * augB(:,:,i));
                        end
                        P(:,:,j,h) = P1 - P2 * pinv(P3) * P2';
                    end % j: previus Markov mode
                end % h: sending history
            end % while    
        end % function makeSteadyStateControlCovarianceMatrices
    end
end

