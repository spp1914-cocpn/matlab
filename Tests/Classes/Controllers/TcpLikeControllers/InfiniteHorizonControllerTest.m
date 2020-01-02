classdef InfiniteHorizonControllerTest< BaseTcpLikeControllerTest
    % Test cases for InfiniteHorizonController.
    
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
    
    properties (Access = private)
        A_ctrb;
        B_ctrb;
    end
    
    methods (Access = protected)
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            % (A, B) must be a controllable pair 
            this.A_ctrb = [0 1 0; 
                      1 0 1;
                      0 0 1];
            this.B_ctrb = [1 0;
                      0 1;
                      1 1];
        end
        
        %% initControllerUnderTest
        function controller = initControllerUnderTest(this)
             controller = InfiniteHorizonController(this.A_ctrb, this.B_ctrb, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength);
        end
        
        %% computeExpectedCosts
        function expectedCosts = computeExpectedCosts(this, states, inputs)
           expectedCosts =  states(:, 1)' * this.Q * states(:, 1) ...
                + states(:, 2)' * this.Q * states(:, 2) + inputs' * this.R * inputs;
            % average costs
            expectedCosts = expectedCosts / size(inputs, 2);
        end
        
        %% computeExpectedStageCosts
        function expectedStageCosts = computeExpectedStageCosts(this, state, input, ~)
             expectedStageCosts =  state' * this.Q * state + input' * this.R * input;
        end
    end
    
    methods (Access = private)
        %% computeExpectedInputs
        function inputs = computeExpectedInputs(this)
            % in this setup, the current input must arrive without delay,
            % or plant evolves open-loop
            % hence, in this setup, G and F are not existent
            % likewise, H
            % due to the structure of the transition matrix,
            % the mode-dependent inputs should not
            % differ
                        
            
            % first mode
            R_1 = this.R;
            Q_1 = this.Q;
            A_1 = this.A_ctrb;
            B_1 = this.B_ctrb;
                       
            % second mode
            R_2 = zeros(this.dimU);
            %Q_2 = this.Q;
            A_2 = this.A_ctrb;
            B_2 = zeros(this.dimX, this.dimU);
            
            % resulting stationary P matrix should be the same for both modes
            % due to the structure of the transition matrix
            P = zeros(this.dimX, this.dimX);
            previousP = P;
            convergeDiffmin = InfiniteHorizonController.maxIterationDiff;
            counter = 0; 
                    
            while(1)
                counter = counter + 1;
                % variable to monitor difference between two iterations of currP:
                % we use sum sum of the difference as currP is known that PreviousP >
                % currP in the positive definit sense. Therefore the quadratic form
                % x' currP x with x = [1 1 .... 1 1]'. currP is sum(sum(P)).
                % Therefore, this is an indicator for the change in positive-definite
                % sense. 
                convergeDiff = sum(sum(P - previousP));
                % variable to monitor change of difference
                              
                % evaluate terminate conditions only if we are out of the initial
                % time region 
                if counter > this.numModes + 1 
                    if max(abs(convergeDiff)) < InfiniteHorizonController.minIterationDiff
                        % convergence
                        break;
                    elseif max(abs(convergeDiff)) > InfiniteHorizonController.maxIterationDiff ...
                            || min(convergeDiff) < 0 
                        % numerical issues
                        P = Pmin;
                        break;
                    end
                    
                    % save best solution obtained yet. In case we diverge due to
                    % numerical problems, we can still use this best solution
                    if convergeDiff < convergeDiffmin 
                        convergeDiffmin = convergeDiff;
                        Pmin = P;
                    end
                end
                  
                previousP = P;
                % same for both modes
                % ensure that matrices are symmetric
                APA = A_1' * previousP * A_1;
                BPB = B_1' * previousP * B_1;
                B2PB2 = B_2' * previousP * B_2;
                P1 = Q_1 + (APA + APA') / 2;
                P2 = this.transitionMatrix(1,1) * A_1' * previousP * B_1 ...
                    + this.transitionMatrix(1,2) * A_2' * previousP * B_2;
                P3 = this.transitionMatrix(1,1) * (R_1 + (BPB + BPB') / 2) ...
                    + this.transitionMatrix(1,2) * (R_2 + (B2PB2 + B2PB2') / 2);
                P4 = P2 * pinv(P3) * P2';
                P = P1 - (P4 + P4') / 2; % should also be symmetric                        
            end
            % input for both modes
            % ensure that matrices are symmetric
            BPB = B_1' * P * B_1;
            B2PB2 = B_2' * P * B_2;
            firstPart = pinv(this.transitionMatrix(1,1) * (R_1 + (BPB + BPB') / 2) ... 
                + this.transitionMatrix(1,2) * (R_2 + (B2PB2 + B2PB2') / 2));
            secondPart = this.transitionMatrix(1,1) * B_1' * P * A_1 + ...
                this.transitionMatrix(1,2) * B_2' * P * A_2;
            
            inputs = -firstPart * secondPart * this.state;
        end
    end
    
    methods (Test)
        %% testInfiniteHorizonControllerInvalidSystemMatrix
        function testInfiniteHorizonControllerInvalidSystemMatrix(this)
             expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
             
             invalidSysMatrix = eye(this.dimX, this.dimX + 1); % not square
             this.verifyError(@() InfiniteHorizonController(invalidSysMatrix, this.B, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength), expectedErrId);
             
             invalidSysMatrix = eye(this.dimX, this.dimX); % square but not finite
             invalidSysMatrix(1, end) = inf;
             this.verifyError(@() InfiniteHorizonController(invalidSysMatrix, this.B, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength), expectedErrId);
        end
        
        %% testInfiniteHorizonControllerInvalidInputMatrix
        function testInfiniteHorizonControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidInputMatrix = eye(this.dimX +1, this.dimU); % invalid dims
            this.verifyError(@() InfiniteHorizonController(this.A, invalidInputMatrix, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength), expectedErrId);
             
            invalidInputMatrix = eye(this.dimX, this.dimU); % correct dims, but not finite
            invalidInputMatrix(1, end) = nan;
            this.verifyError(@() InfiniteHorizonController(this.A, invalidInputMatrix, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength), expectedErrId);
        end
        
                
        %% testInfiniteHorizonControllerInvalidCostMatrices
        function testInfiniteHorizonControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = eye(this.dimX); % Q is not symmetric
            invalidQ(1, end) = 1;
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
            
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, this.Q, invalidR, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, this.Q, invalidR, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, this.Q, invalidR, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
        end
        
        %% testInfiniteHorizonControllerInvalidDelayProbs
        function testInfiniteHorizonControllerInvalidDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength), expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength), expectedErrId);
            
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength), expectedErrId);
        end
        
        %% testInfiniteHorizonControllerInvalidPlant
        function testInfiniteHorizonControllerInvalidPlant(this)
            % the plant (A,B) is uncontrollable/has uncontrollable modes
            expectedErrId = 'InfiniteHorizonController:InvalidPlant';
            
            this.verifyError(@() InfiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
        end
        
        %% testInfiniteHorizonControllerInvalidFlag
        function testInfiniteHorizonControllerInvalidFlag(this)
            expectedErrId = 'InfiniteHorizonController:InvalidUseMexFlag';
            invalidUseMexFlag = 'invalid'; % not a flag
            
            this.verifyError(@() InfiniteHorizonController(this.A_ctrb, this.B_ctrb, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, invalidUseMexFlag), expectedErrId);
        end
        
        %% testInifiniteHorizonControllerConvergenceImpossible
        function testInifiniteHorizonControllerConvergenceImpossible(this)
            % sufficient condition for non-existance of stabilizing controller
            % is given (thereom 4d) of CDC paper
            expectedErrId = 'InfiniteHorizonController:ConvergenceImpossible';
            
            A_invalid = 10 * this.A_ctrb; % the spectral radius of this matrix is 10
            this.verifyError(@() InfiniteHorizonController(A_invalid, this.B_ctrb, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
        end
%%
%%
        %% testInfiniteHorizonController
        function testInfiniteHorizonController(this)
            % successfully initialize controller, should not crash
            controller = InfiniteHorizonController(this.A_ctrb, this.B_ctrb, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength);
            
            this.verifyNotEmpty(controller);
            
            %by default, we use the mex implementation to obtain the
            %controller gains
            this.verifyTrue(controller.useMexImplementation);
            this.verifyTrue(controller.requiresExternalStateEstimate); % controller needs a filter
            
            % status must be integer within [-2, 1]
            this.verifyGreaterThanOrEqual(controller.status, -2);
            this.verifyLessThanOrEqual(controller.status, 1);
            % status should be integer value and 1 (convergence)
            this.verifyTrue(mod(controller.status, 1) == 0);
            this.verifyEqual(controller.status, 1);
        end
%%
%%
        %% testDoControlSequenceComputationInvalidMode
        function testDoControlSequenceComputationInvalidMode(this)
            % call function without a mode given
            expectedErrId = 'MATLAB:minrhs';
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution), ...
                expectedErrId);
            
            expectedErrId = 'InfiniteHorizonController:DoControlSequenceComputation:InvalidMode';
                       
            invalidMode = this; % not a scalar
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode), ...
                expectedErrId);
            
            invalidMode = 0; % scalar, but non-positive
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode), ...
                expectedErrId);
            
            invalidMode = this.numModes + 1; % scalar, but out of bounds
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode), ...
                expectedErrId);
            
            invalidMode = 0.5; % scalar, but not an integer
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode), ...
                expectedErrId);
        end
        
        %% testDoControlSequenceComputationZeroState
        function testDoControlSequenceComputationZeroState(this)
            % perform a sanity check: given state is origin, so computed
            % control sequence (length 1) should be also the zero vector
            % due to the underlying linear control law, independent of the
            % previous mode
            this.assertTrue(this.controllerUnderTest.useMexImplementation);
            expectedSequence = zeros(this.dimU, 1);

            % first mode
            actualSequence = this.controllerUnderTest.computeControlSequence(this.zeroStateDistribution, 1);
            
            this.verifyEqual(actualSequence, expectedSequence);
            
            %second mode
            actualSequence = this.controllerUnderTest.computeControlSequence(this.zeroStateDistribution, 2);
            
            this.verifyEqual(actualSequence, expectedSequence);
        end
        
        %% testDoControlSequenceComputationZeroStateNoMex
        function testDoControlSequenceComputationZeroStateNoMex(this)
            % perform a sanity check: given state is origin, so computed
            % control sequence (length 1) should be also the zero vector
            % due to the underlying linear control law, independent of the
            % previous mode
            useMexImplementation = false;
            controller = InfiniteHorizonController(this.A_ctrb, this.B_ctrb, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, useMexImplementation);
            this.assertFalse(controller.useMexImplementation);
            
            expectedSequence = zeros(this.dimU, 1);

            % first mode
            actualSequence = controller.computeControlSequence(this.zeroStateDistribution, 1);
            
            this.verifyEqual(actualSequence, expectedSequence);
            
            %second mode
            actualSequence = controller.computeControlSequence(this.zeroStateDistribution, 2);
            
            this.verifyEqual(actualSequence, expectedSequence);
        end
        
        %% testDoControlSequenceComputation
        function testDoControlSequenceComputation(this)
            this.assertTrue(this.controllerUnderTest.useMexImplementation);
            expectedInputs = this.computeExpectedInputs();
            
            % check both modes
            % first mode: previous input arrived at plant
            actualInput = this.controllerUnderTest.computeControlSequence(this.stateDistribution, 1);
            this.verifyEqual(actualInput, expectedInputs, 'AbsTol', 1e-12);
            
            % second mode: previous input did not arrive at plant
            actualInput = this.controllerUnderTest.computeControlSequence(this.stateDistribution, 2);
            this.verifyEqual(actualInput, expectedInputs, 'AbsTol', 1e-12);
        end
        
         %% testDoControlSequenceComputationNoMex
        function testDoControlSequenceComputationNoMex(this)
            useMexImplementation = false;
            controller = InfiniteHorizonController(this.A_ctrb, this.B_ctrb, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, useMexImplementation);
            
            this.assertFalse(controller.useMexImplementation);
            expectedInputs = this.computeExpectedInputs();
            
            % check both modes
            % first mode: previous input arrived at plant
            actualInput = controller.computeControlSequence(this.stateDistribution, 1);
            this.verifyEqual(actualInput, expectedInputs, 'AbsTol', 1e-12);
            
            % second mode: previous input did not arrive at plant
            actualInput = controller.computeControlSequence(this.stateDistribution, 2);
            this.verifyEqual(actualInput, expectedInputs, 'AbsTol', 1e-12);
        end
%%
%%
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            horizonLength = 1;
            appliedInputs = ones(this.dimU, horizonLength);
            expectedErrId = 'InfiniteHorizonController:DoCostsComputation:InvalidStateTrajectory';
            
            invalidStateTrajectory = ones(this.dimX, horizonLength + 10); % state trajectory too long
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
            
            invalidStateTrajectory = ones(this.dimX, horizonLength); % state trajectory too short
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
        end
    end
end

