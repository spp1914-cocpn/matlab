classdef InfiniteHorizonController < SequenceBasedController
    % Implementation of the optimal infinite horizon linear sequence-based LQG controller for
    % NCS with a TCP-like network connecting the controller and the
    % actuator.
    %
    % Literature: 
    %   Jörg Fischer, Achim Hekler, Maxim Dolgov, and Uwe D. Hanebeck,
    %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject to Random Transmission Delays and Packet Losses,
    %   Proceedings of the 2013 American Control Conference (ACC 2013),
    %   Washington D. C., USA, June 2013.
    %
    %   Jörg Fischer, Maxim Dolgov, and Uwe D. Hanebeck,
    %   On Stability of Sequence-Based LQG Control,
    %   Proceedings of the 52st IEEE Conference on Decision and Control (CDC 2013),
    %   Florence, Italy, December 2013.
    %
    % While this implementation is to some extent more general than the original one of Fischer et al., 
    % which can be found <a href="matlab:
    % web('http://www.cloudrunner.eu/algorithm/142/optimal-sequence-based-lqg-control-over-tcp-like-networks/version/1/')"
    % >here</a>, some parts are directly based thereof.
    
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
        % Markov chain
        % transition matrix of the Markov chain;
        % (matrix of dimension <sequenceLength+1> x <sequenceLength+1>)
        transitionMatrix = [];
        % control gain matrix: U = -L*sysState;
        % (matrix of dimension <dimU*sequenceLength> x <dimState>)
        L = [];
        % dimension of the augmented system state; 
        % (positive integer <dimX+dimU*sequenceLength*(sequenceLength-1)/2>)
        dimState = -1;
        
    end % properties

    properties (SetAccess = immutable, GetAccess = public)
        % status of controller initialization (integer with following values)
        %   1 : successful
        %   0 : maybe successful
        %       (controller gain diverged or its convergence was not
        %       detected; controller works with best gain obtained)
        %  -1 : maybe successful
        %       (numerical problems occured; controller works with best
        %       gain obtained so far)
        %  -2 : not usscessful; 
        %       (system is not stabilizable over the network; no gain
        %       computed)
        status = 0;
    end
    
    properties (SetAccess = private, GetAccess = protected)
        % augmented state space
        % system state of the augmented state space;
        % (column vector of dimension <dimState>)
        sysState = [];
    end
    
    properties (Constant, Access = public)
        % Bounds to determine convergence / divergence of controller gain
        % interation
        
        % lower bound on absolute difference (indicator for convergence)
        minIterationDiff = 1e-6;
        % higher bound on absolute difference (indicator for divergence)
        maxIterationDiff = 1e+11;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %% InfiniteHorizonController
        function this = InfiniteHorizonController(A, B, Q, R, delayProb, sequenceLength)
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
            % Returns:
            %   << obj (InfiniteHorizonController)
            %      A new InfiniteHorizonController instance.
                     
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength);
           
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

            % Check for stabilizablity of A by B
            if dimX - rank(ctrb(A, B)) > 0
              error('InfiniteHorizonController:InvalidPlant', ...
                '** Plant (A, B) has uncontrollable eigenvalues **');
            end
            
            this.transitionMatrix = Utility.calculateDelayTransitionMatrix(this.delayProb);
            
            [this.dimState, this.F, this.G, augA, augB, augQ, augR] ...
                = Utility.performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R);
            this.sysState = zeros(this.dimState, 1);
            
            % make controller covariance matrices and check result
            [P, this.status] = this.computeSteadyStateControlCovarianceMatrices(A, augA, augB, augQ, augR);
            if this.status == -2
              error(['** Stabizilizing infinite time horizon controller', ...
                'does not exist. **']);
            end
            this.L = this.computeSteadyStateControlGainMatrices(augA, augB, augR, P);
        end % function InfiniteHorizonController    
  
        %% reset
        function reset(s)
           s.sysState = zeros(s.dimState,1);
        end % function reset
  
    end % methods public


    methods (Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, mode, ~)
            if ~Checks.isScalarIn(mode, 1, this.sequenceLength + 1) || mod(mode, 1) ~= 0
              error('InfiniteHorizonController:DoControlSequenceComputation:InvalidMode', ...
                  '** Input parameter <mode> (previous plant mode/mode estimate) must be in {1, ... %d} **', ...
                  this.sequenceLength + 1);
            end
            [stateMean, ~] = state.getMeanAndCovariance();
            
            this.sysState(1:this.dimPlantState) = stateMean(:);
            inputSequence = this.L(: , :, mode) * this.sysState;
            this.sysState(this.dimPlantState + 1:end) = ...
                this.F * this.sysState(this.dimPlantState + 1:end) + this.G * inputSequence;
        end
        
        %% doCostsComputation
        function averageLQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonLength = size(appliedInputs, 2);
            if size(stateTrajectory, 2) ~= horizonLength + 1
                error('InfiniteHorizonController:DoCostsComputation:InvalidStateTrajectory', ...
                    '** <stateTrajectory> is expected to have %d columns ', horizonLength + 1);
                    
            end
            averageLQGCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, appliedInputs, this.Q, this.R) / horizonLength;
        end
        
        %% computeSteadyStateControlCovarianceMatrices
        function [Pout, status] = computeSteadyStateControlCovarianceMatrices(s, A, augA, augB, augQ, augR)
            %========================================================================
            % computeSteadyStateControlCovarianceMatrices: computes the steady state
            % solution of the controller by setting the horizon length to infinity
            % (actually done by iterating through the horizon until the difference
            % between two consecutive iteration does not exceed <maxConsecutiveDiff>
            % defined after this comment).
            %
            % Synopsis:
            % computeSteadyStateControlCovarianceMatrices(A, augA, augB, augQ, augR)
            %
            % Output Parameter:
            % - Pout: 3D matrix with control covariances; 
            %      (positive semidefinite matrix of dimension 
            %       <dimState> x <dimState> x <sequenceLength+1>)
            % - status: gives status of the controller gain convergence
            %           1 : successful
            %           0 : maybe successful
            %               (controller gain diverged or its convergence was not
            %               detected; best gain obtained returned that was obtained)
            %          -1 : maybe successful
            %               (numerical problems occured; best gain returned that was
            %               obtained)
            %          -2 : not usscessful; 
            %               (system is not stabilizable over the network; no gain
            %               computed)
            %
            % Pout(:,:,markovMode), markovMode: 1...packetLength+1 (1 = no delay)
            %========================================================================  
            currP = zeros(s.dimState, s.dimState, s.sequenceLength + 1);
            previousP = currP;
            PreviousConvergeDiff = 0;
            %fallingFlag = 0;
            convergeDiffmin = s.maxIterationDiff;
            % counter that counts the number of performed iterations
            counter = 0; 

            % Check if convergence of controller gain is in general possible
            if (s.transitionMatrix(end, end) * (max(abs(eig(A))))^2) > 1
                % no convergence possible
                status = -2;
                Pout = -1;
                return
            end
            % ------------ Interation to calculate currP:
            % for stable closed-loop systems currP should increase every iteration
            % until it is converged --> until difference between iterations is
            % zero. Because of numerical problems in the inversion currP does not
            % totally converge difference only gets small. It can even happen that
            % currP diverges after if the iteration is done too often with too low
            % precision. Therefore, different terminiting criteria for the
            % iteration are implented
            while(1)
                counter = counter + 1;
                % variable to monitor difference between two iterations of currP:
                % we use sum sum of the difference as currP is known that PreviousP >
                % currP in the positive definit sense. Therefore the quadratic form
                % x' currP x with x = [1 1 .... 1 1]'. currP is sum(sum(P)).
                % Therefore, this is an indicator for the change in positive-definite
                % sense. 
                convergeDiff = sum(sum(currP - previousP));
                % variable to monitor change of difference
                convergeAcc = convergeDiff - PreviousConvergeDiff;
                % debuggung stuff
                s.debugging.numIterationsToSteadyStateP = counter;
                s.debugging.convergeDiff(counter + 1, :) = convergeDiff;
                s.debugging.convergeAcc(counter + 1, :) = convergeAcc;
                % evaluate terminate conditions only if we are out of the initial
                % time region  
                if( counter > s.sequenceLength + 2 ) 
                    % ------- Terminate Conditions ------------------
                    % Terminate Condition 1: stable / converged:
                    % this happens if currP does not change anymore between two
                    % iterations. Instead of zero, we use a lower bound
                    if max(abs(convergeDiff)) < s.minIterationDiff
                        %disp('converged');
                        Pout = currP;
                        status = 1;
                        break;
                    % Terminate Condition 2: unstable / diverged:
                    % this happens if the closed-loop system is unstable ( controller
                    % not able to stabilize system) or we had numerical problems that
                    % prevented to stop at the converged currP. The cases are hard to
                    % distinguish so we return the best solution found and the caution
                    % flag set.
                    elseif( max(abs(convergeDiff)) > s.maxIterationDiff )
                        %disp('NOT converged or convergence not detected');
                        counter %#ok<NOPRT>
                        convergeCounterMin %#ok<NOPRT>
                        Pout = Pmin;
                        status = 0;
                        break;
                    % Terminate Condition 3: Numerical Problems
                    % This should not happen but if it does we had numerical problems.
                    % The best solution found will be returned.
                    elseif( min(convergeDiff) < 0 )
                        %disp('Numerical Problems');
                        counter %#ok<NOPRT>
                        convergeCounterMin %#ok<NOPRT>
                        Pout = Pmin;
                        status = -1;
                        break;
                     end
        
                    % save best solution obtained yet. In case we diverge due to
                    % numerical problems, we can still use this best solution
                    if( convergeDiff < convergeDiffmin )
                      convergeCounterMin = counter;
                      convergeDiffmin = convergeDiff;
                      Pmin = currP;
                    end
                end
                PreviousConvergeDiff = convergeDiff;
                previousP = currP;
                currP = zeros(s.dimState, s.dimState, s.sequenceLength+1);
     
                for j = 1 : s.sequenceLength + 1
                    % temporary components of
                    P1 = 0;
                    P2 = 0;
                    P3 = 0;
                    for i = 1 : s.sequenceLength + 1
                        P1 = P1 + s.transitionMatrix(j,i) * ( augQ(:,:,i) + ...
                            augA(:,:,i)' * previousP(:,:,i) * augA(:,:,i) );
                        P2 = P2 + s.transitionMatrix(j,i) * ...
                            (augA(:,:,i)' * previousP(:,:,i) * augB(:,:,i));
                        P3 = P3 + s.transitionMatrix(j,i) * ( augR(:,:,i) ...
                            + augB(:,:,i)' * previousP(:,:,i) * augB(:,:,i) );
                    end % i: current Markov mode
                    currP(:,:,j) = P1 - P2 * pinv(P3) * P2';
                end % j: previus Markov mode
            end % while
        end % function makeSteadyStateControlCovarianceMatrices
  
  
        %% computeSteadyStateControlGainMatrices
        function L = computeSteadyStateControlGainMatrices(this, augA, augB, augR, P)
            %========================================================================
            % computeSteadyStateControlGainMatrices: computes the gain matrices which
            % are used to calculate the control sequence by multiplying with the
            % augmented system state: U = L*x
            %
            % Synopsis:
            % computeSteadyStateControlGainMatrices(augA, augB, augR, P)
            %
            % Output Parameter:
            % - L: 3D matrix of control gains; (matrix of dimension 
            %      <dimU*sequenceLength> x <dimState> x <sequenceLength+1>)
            %
            % L(:,:,markovMode), markovMode: 1...sequenceLength+1 (1 = no delay)
            %========================================================================  
            numModes = this.sequenceLength + 1;
            dim = size(augR, 1);
            
            L = zeros(dim, this.dimState, numModes);
            for j = 1:numModes
                L_1 = zeros(dim);
                L_2 = zeros(dim, this.dimState);
                for i = 1:numModes
                    L_1 = L_1 + this.transitionMatrix(j,i)*...
                        (augR(:,:,i) + augB(:,:,i)'* P(:,:,i)* augB(:,:,i));
                    L_2 = L_2 + this.transitionMatrix(j,i)*...
                        (augB(:,:,i)'* P(:,:,i) * augA(:,:,i));
                end % i: current Markov mode
                L(:,:,j) = - pinv(L_1) * L_2;
            end % i: previous Markov mode
        end % function makeSteadyStateControlGainMatrices
    end
end