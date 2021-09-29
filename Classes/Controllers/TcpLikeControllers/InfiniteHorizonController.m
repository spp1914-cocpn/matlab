classdef InfiniteHorizonController < SequenceBasedController & CaDelayProbsChangeable
    % Implementation of the optimal infinite horizon linear sequence-based LQG controller for
    % NCS with a TCP-like network connecting the controller and the actuator.
    %
    % While this implementation is to some extent more general than the original one of Fischer et al., 
    % which can be found <a href="matlab:web('http://www.cloudrunner.eu/algorithm/142/optimal-sequence-based-lqg-control-over-tcp-like-networks/version/1/')" >here</a>, 
    % some parts are directly based thereof.
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
        
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        %% input properties
        % system state cost matrix;
        % (positive semidefinite matrix of dimension <dimX> x <dimX>)
        Q = [];
        % input cost matrix;
        % (positive definite matrix of dimension <dimU> x <dimU>)
        R = [];
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
        % dimension of the augmented system state; 
        % (positive integer <dimX+dimU*sequenceLength*(sequenceLength-1)/2>)
        dimState = -1;
        
        augA;
        augB;
        augQ;
        augR;
    end

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
        
        useMexImplementation(1,1) logical = true; 
        % by default, we use the C++ (mex) implementation for computation of controller gains
        % this is faster, but can produce slightly different results
    end
    
    properties (SetAccess = private, GetAccess = protected)
        % augmented state space
        % system state of the augmented state space;
        % (column vector of dimension <dimState>)
        sysState = [];
        
        % Markov chain
        % transition matrix of the Markov chain;
        % (matrix of dimension <sequenceLength+1> x <sequenceLength+1>)
        transitionMatrix = [];
        % control gain matrix: U = -L*sysState;
        % (matrix of dimension <dimU*sequenceLength> x <dimState>)
        L = [];
    end
    
    properties (Constant, Access = public)
        % Bounds to determine convergence / divergence of controller gain
        % interation
        
        % lower bound on absolute difference (indicator for convergence)
        minIterationDiff = 1e-6;
        % higher bound on absolute difference (indicator for divergence)
        maxIterationDiff = 1e+11;
    end

    methods (Access = public)
        %% InfiniteHorizonController
        function this = InfiniteHorizonController(A, B, Q, R, modeTransitionMatrix, sequenceLength, useMexImplementation)
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
            %   >> modeTransitionMatrix (Stochastic matrix, i.e. a square matrix with nonnegative entries whose rows sum to 1)
            %      The transition matrix of the mode theta_k of the augmented dynamics.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %   >> useMexImplementation (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for the computation of the controller
            %      gains which is considerably faster than the Matlab
            %      implementation. 
            %      If left out, the default value true is used.
            %
            % Returns:
            %   << obj (InfiniteHorizonController)
            %      A new InfiniteHorizonController instance.
                     
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength, true);        
                       
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
                                     
            % mode transition matrix
            Validator.validateTransitionMatrix(modeTransitionMatrix, sequenceLength + 1);
            this.transitionMatrix = modeTransitionMatrix;
            
             % Check for stabilizablity of A by B
            assert(dimX -rank(ctrb(A, B)) == 0, ...
                'InfiniteHorizonController:InvalidPlant', ...
                '** Plant (A, B) has uncontrollable eigenvalues **');
            
            [this.dimState, this.F, this.G, this.augA, this.augB, this.augQ, this.augR] ...
                = Utility.performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R);
            this.sysState = zeros(this.dimState, 1);
            
            if nargin == 7
                this.useMexImplementation = useMexImplementation;
            end
            
            if this.useMexImplementation
                [this.L, this.status] = mex_InfiniteHorizonController(this.augA, this.augB, this.augQ, this.augR, this.transitionMatrix, A);
                assert(this.status ~= -2, ...
                    'InfiniteHorizonController:ConvergenceImpossible', ...
                    '** Stabizilizing infinite time horizon controller does not exist **');
            else
                
                [P, this.status] = this.computeSteadyStateControlCovarianceMatrices(A);
                assert(this.status ~= -2, ...
                    'InfiniteHorizonController:ConvergenceImpossible', ...
                    '** Stabizilizing infinite time horizon controller does not exist **');
                this.L = this.computeSteadyStateControlGainMatrices(P);
            end
        end
  
        %% reset
        function reset(this)
           this.sysState = zeros(this.dimState,1);
        end 
        
        %% setEtaState
        function setEtaState(this, newEta)
            assert(Checks.isVec(newEta, this.dimState - this.dimPlantState) && all(isfinite(newEta)), ...
                'InfiniteHorizonController:SetEtaState:InvalidEta', ...
                '** <newEta> must be a %d-dimensional vector **', ...
                this.dimState - this.dimPlantState);
            
            this.sysState(this.dimPlantState + 1:end) = newEta(:);
        end
        
         %% changeCaDelayProbs
        function changeCaDelayProbs(this, newCaDelayProbs)
            newMat =  Utility.calculateDelayTransitionMatrix(...
                Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.sequenceLength + 1));
            if ~isequal(newMat, this.transitionMatrix)
                this.transitionMatrix = newMat;
                % and we must recompute the controller gain
                A = this.augA(1:this.dimPlantState, 1:this.dimPlantState, 1); 
                if this.useMexImplementation
                    [this.L, ~] = mex_InfiniteHorizonController(this.augA, this.augB, this.augQ, this.augR, this.transitionMatrix, A);
                else
                    [P, ~] = this.computeSteadyStateControlCovarianceMatrices(A);
                    this.L = this.computeSteadyStateControlGainMatrices(P);           
                end     
            end
        end
        
        %% changeTransitionMatrix
        function changeTransitionMatrix(this, newTransitionMatrix)
            if ~isequal(newTransitionMatrix, this.transitionMatrix)
                this.transitionMatrix = newTransitionMatrix;
                % and we must recompute the controller gain
                A = this.augA(1:this.dimPlantState, 1:this.dimPlantState, 1);
                if this.useMexImplementation
                    [this.L, ~] = mex_InfiniteHorizonController(this.augA, this.augB, this.augQ, this.augR, this.transitionMatrix, A);
                else
                    [P, ~] = this.computeSteadyStateControlCovarianceMatrices(A);
                    this.L = this.computeSteadyStateControlGainMatrices(P);           
                end
            end
        end
    end 


    methods (Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, state, mode, ~)
            assert(Checks.isScalarIn(mode, 1, this.sequenceLength + 1) && mod(mode, 1) == 0, ...
                'InfiniteHorizonController:DoControlSequenceComputation:InvalidMode', ...
                '** Input parameter <mode> (previous plant mode/mode estimate) must be in {1, ... %d} **', ...
                this.sequenceLength + 1);
            
            [stateMean, ~] = state.getMeanAndCov();
            
            this.sysState(1:this.dimPlantState) = stateMean(:);
            inputSequence = this.L(: , :, mode) * this.sysState;
            this.sysState(this.dimPlantState + 1:end) = ...
                this.F * this.sysState(this.dimPlantState + 1:end) + this.G * inputSequence;
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, ~)
                        
            stageCosts = Utility.computeStageCosts(state, input, this.Q, this.R);
        end
        
        %% doCostsComputation
        function averageLQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonLength = size(appliedInputs, 2);
            assert(size(stateTrajectory, 2) == horizonLength + 1, ...
                'InfiniteHorizonController:DoCostsComputation:InvalidStateTrajectory', ...
                '** <stateTrajectory> is expected to have %d columns ', horizonLength + 1);
                       
            averageLQGCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, appliedInputs, this.Q, this.R) / horizonLength;
        end
        
        %% computeSteadyStateControlCovarianceMatrices
        function [Pout, status] = computeSteadyStateControlCovarianceMatrices(this, A)
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
            numModes = this.sequenceLength + 1;
            
            currP = zeros(this.dimState, this.dimState, numModes);
            previousP = currP;
            PreviousConvergeDiff = 0;
            %fallingFlag = 0;
            convergeDiffmin = this.maxIterationDiff;
            % counter that counts the number of performed iterations
            counter = 0; 

            % Check if convergence of controller gain is in general possible
            % Thereom 4d) of CDC paper mentioned above
            if (this.transitionMatrix(end, end) * (max(abs(eig(A))))^2) > 1
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
                this.debugging.numIterationsToSteadyStateP = counter;
                this.debugging.convergeDiff(counter + 1, :) = convergeDiff;
                this.debugging.convergeAcc(counter + 1, :) = convergeAcc;
                % evaluate terminate conditions only if we are out of the initial
                % time region  
                if( counter > this.sequenceLength + 2 ) 
                    % ------- Terminate Conditions ------------------
                    % Terminate Condition 1: stable / converged:
                    % this happens if currP does not change anymore between two
                    % iterations. Instead of zero, we use a lower bound
                    if max(abs(convergeDiff)) < this.minIterationDiff
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
                    elseif( max(abs(convergeDiff)) > this.maxIterationDiff )
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
                   if(convergeDiff < convergeDiffmin)                       
                      convergeCounterMin = counter;
                      convergeDiffmin = convergeDiff;
                      Pmin = currP;
                    end
                end
                PreviousConvergeDiff = convergeDiff;
                previousP = currP;
                currP = zeros(this.dimState, this.dimState, numModes);
                
                AP = mtimesx(this.augA, 'T', previousP);
                BP = mtimesx(this.augB, 'T', previousP);
                APA = mtimesx(AP, this.augA);
                QAPA = this.augQ + (APA + permute(APA, [2 1 3])) / 2; % ensure symmetry
                BPB = mtimesx(BP, this.augB);
                APB = mtimesx(AP, this.augB);
                RBPB = this.augR + (BPB + permute(BPB, [2 1 3])) / 2; % ensure symmetry
                
                for j = 1:numModes
                    P1 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* QAPA, 3);
                    P2 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* APB, 3);
                    P3 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* RBPB, 3);

                    P4 = P2 * pinv(P3) * P2';
                    currP(:,:,j) = P1 - (P4 + P4') / 2;   % ensure that symmetry is maintained                
                end
            end            
        end
  
  
        %% computeSteadyStateControlGainMatrices
        function L = computeSteadyStateControlGainMatrices(this, P)
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
            dim = size(this.augR, 1);
            
            L = zeros(dim, this.dimState, numModes);
            BP = mtimesx(this.augB, 'T', P);
            BPB = mtimesx(BP, this.augB);
            % ensure the symmetry 
            RBPB = this.augR + (BPB + permute(BPB, [2 1 3])) / 2;
            BPA = mtimesx(BP, this.augA);
                        
            for j = 1:numModes
                L_1 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* RBPB, 3);
                L_2 = sum(reshape(this.transitionMatrix(j, :), 1, 1, numModes) .* BPA, 3);
                L(:,:,j) = - pinv(L_1) * L_2;
            end
        end
    end
end
