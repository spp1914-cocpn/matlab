classdef NominalPredictiveControllerTest < matlab.unittest.TestCase
    % Test cases for NominalPredictiveController.
    
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
        
    properties (Constant)
        absTol = 1e-5;
    end
    
    properties (Access = protected)
        A;
        B;
        Q;
        R;
        sequenceLength;
        dimX;
        dimU;
        initialPlantState;
        stateTrajectory;
        inputTrajectory;
        horizonLength;
        
        expectedL;
        
        % for reaching a nonzero setpoint
        feedforward;
        setpoint;
        stateTrajectoryNonzeroSetpoint;
        inputTrajectoryNonzeroSetpoint;
        
        controllerUnderTest;
    end
    
    methods (Sealed, TestMethodSetup)
        %% initProperties
        function initProperties(this)
            % use (noise-free) stirred tank example (Example 6.15, p. 500-501) from
            %
            % Huibert Kwakernaak, and Raphael Sivan, 
            % Linear Optimal Control Systems,
            % Wiley-Interscience, New York, 1972.
            %
            
            this.dimX = 2;
            this.dimU = 2;
            
            V = diag([0.01, 1]); % in the book: z = V*x
            this.A = diag([0.9512, 0.9048]);
           
            this.B = [4.877 4.877; -1.1895 3.569];
            this.Q = V * diag([50 0.02]) * V; % R_3 in the book
            this.R = diag([1/3, 3]); % R_2 in the book
                       
            this.expectedL = -[ 0.07125 -0.07029;
                                0.01357 0.04548]; % F in the book
    
            this.initialPlantState = [0.1; 0];
            
            this.sequenceLength = 10;            
            this.horizonLength = this.sequenceLength;
            
            this.initAdditionalProperties();            
            closedLoopSystem = this.initClosedLoopSystem();          
            
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories(closedLoopSystem);
            
            this.setpoint = [42; 42];
            % we require a feedforward term now, computed according to
            % Thereom 6.34 (p. 506), which is applicable as dimX = dimU, from
            %
            % Huibert Kwakernaak, and Raphael Sivan, 
            % Linear Optimal Control Systems,
            % Wiley-Interscience, New York, 1972.
            %
            
            % obtain the feedforward using dc gain  of closed loop system
            % could also be computed directly using transfer matrix of closed loop system: 
            % evaluate H(z) = I * inv(z*I - A -B*L) * B for z=1 (this is the dc gain), then the
            % feedforward is H(z=1)^-1 * setpoint
            this.feedforward = inv(dcgain(closedLoopSystem)) * this.setpoint; %#ok  
            [this.stateTrajectoryNonzeroSetpoint, this.inputTrajectoryNonzeroSetpoint] ...
                = this.computeTrajectoriesNonzeroSetpoint(closedLoopSystem);

            this.controllerUnderTest = this.initControllerUnderTest();
        end
    end
            
    methods (Access = protected)        
        %% callBaseCtor
        function controller = callBaseCtor(this, varargin)
            controller = NominalPredictiveController(varargin{:});
        end

        
        %% initClosedLoopSystem
        function closedLoopSystem = initClosedLoopSystem(this)
            closedLoopSystem = feedback(ss(this.A, this.B, eye(this.dimX), [], []), ss(this.expectedL), +1); % discrete-time model
        end
        
        %% initControllerUnderTest
        function controller = initControllerUnderTest(this)
            controller= NominalPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength); 
        end
        
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            % by default, nothing to do
        end
        
        %% computeTrajectories
        function [stateTrajectory, inputTrajectory] = computeTrajectories(this, closedLoopSystem)
            [~, ~, stateTrajectory] = initial(closedLoopSystem, this.initialPlantState, this.horizonLength);
            stateTrajectory = stateTrajectory';
            
            % compute the correspondig, expected inputs
            inputTrajectory = this.expectedL * stateTrajectory(:, 1:this.horizonLength);
        end
        
        %% computeTrajectoriesNonzeroSetpoint
        function [stateTrajectory, inputTrajectory] = computeTrajectoriesNonzeroSetpoint(this, closedLoopSystem)            
                      
            [~,~, stateTrajectory] = lsim(closedLoopSystem, repmat(this.feedforward', this.horizonLength + 1, 1), ...
                0:this.horizonLength, this.initialPlantState);
            
            stateTrajectory = stateTrajectory';
            % compute the correspondig, expected inputs
            inputTrajectory = this.expectedL * stateTrajectory(:, 1:this.horizonLength) + this.feedforward;            
        end
        
        %% computeCosts
        function costs = computeCosts(this)
            costs =  this.stateTrajectory(:, end)' * this.Q * this.stateTrajectory(:, end);
            for j = 1:this.horizonLength
                costs = costs + this.stateTrajectory(:, j)' * this.Q * this.stateTrajectory(:, j) ...
                    + this.inputTrajectory(:, j)' * this.R * this.inputTrajectory(:, j);
            end
            costs = costs / this.horizonLength;
        end
        
        %% computeCostsNonzeroSetpoint
        function costs = computeCostsNonzeroSetpoint(this)
            diffs = this.stateTrajectoryNonzeroSetpoint - this.setpoint;
            
            costs =  diffs(:, end)' * this.Q * diffs(:, end);
            for j = 1:this.horizonLength
                costs = costs + diffs(:, j)' * this.Q * diffs(:, j) ...
                    + this.inputTrajectoryNonzeroSetpoint(:, j)' * this.R * this.inputTrajectoryNonzeroSetpoint(:, j);
            end
            costs = costs / this.horizonLength;
        end        

    end
    
    methods (Test)        
        %% testCtorInvalidSysMatrix
        function testCtorInvalidSysMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
            
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() this.callBaseCtor(invalidA, this.B, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidA = inf(this.dimX, this.dimX); % square but inf
            this.verifyError(@() this.callBaseCtor(invalidA, this.B, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
        end
        
        %% testCtorInvalidInputMatrix
        function testCtorInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidB = eye(this.dimX + 1, this.dimU); % invalid dims
            this.verifyError(@() this.callBaseCtor(this.A, invalidB, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidB = eye(this.dimX, this.dimU); % correct dims, but nan
            invalidB(1, 1) = nan;
            this.verifyError(@() this.callBaseCtor(this.A, invalidB, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
        end
        
        %% testCtorInvalidCostMatrices
        function testCtorInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() this.callBaseCtor(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() this.callBaseCtor(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() this.callBaseCtor(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() this.callBaseCtor(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() this.callBaseCtor(this.A, this.B, this.Q, invalidR, this.sequenceLength), ...
                expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() this.callBaseCtor(this.A, this.B, this.Q, invalidR, this.sequenceLength), ...
                expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() this.callBaseCtor(this.A, this.B, this.Q, invalidR, this.sequenceLength), ...
                expectedErrId);           
            
        end
        
        %% testCtorInvalidSetpoint
        function testCtorInvalidSetpoint(this)
            expectedErrId = [class(this.controllerUnderTest) ':InvalidSetpoint'];
            
            % not a vector
            invalidSetpoint = this;
            this.verifyError(@() this.callBaseCtor(this.A, this.B, this.Q, this.R, this.sequenceLength, invalidSetpoint), expectedErrId);
            
            % wrong dimension
            invalidSetpoint = [this.setpoint; 42];
            this.verifyError(@() this.callBaseCtor(this.A, this.B, this.Q, this.R, this.sequenceLength, invalidSetpoint), expectedErrId);
            
            % correct dimension, but inf
            invalidSetpoint = this.setpoint;
            invalidSetpoint(2) = inf;
            this.verifyError(@() this.callBaseCtor(this.A, this.B, this.Q, this.R, this.sequenceLength, invalidSetpoint), expectedErrId);
        end
        
        %% testCtor
        function testCtor(this)
            controller = this.callBaseCtor(this.A, this.B, this.Q, this.R, this.sequenceLength);
            
            this.verifyEqual(controller.L, this.expectedL, 'AbsTol', NominalPredictiveControllerTest.absTol);
            this.verifyEmpty(controller.setpoint);
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback
            
            % now pass an empty setpoint
            controller = this.callBaseCtor(this.A, this.B, this.Q, this.R, this.sequenceLength, []);
            
            this.verifyEqual(controller.L, this.expectedL, 'AbsTol', NominalPredictiveControllerTest.absTol);
            this.verifyEmpty(controller.setpoint);
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback
            
            % now a set point
            controller = this.callBaseCtor(this.A, this.B, this.Q, this.R, this.sequenceLength, this.setpoint);
            
            this.verifyEqual(controller.L, this.expectedL, 'AbsTol', NominalPredictiveControllerTest.absTol);
            this.verifyEqual(controller.setpoint, this.setpoint);
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback            
        end
        
        %% testChangeSequenceLength
        function testChangeSequenceLength(this)
            this.assertEqual(this.controllerUnderTest.sequenceLength, this.sequenceLength);
            
            newSeqLength = this.sequenceLength + 2;
            this.controllerUnderTest.changeSequenceLength(newSeqLength);
            
            this.verifyEqual(this.controllerUnderTest.sequenceLength, newSeqLength);
        end
        
        %% testChangeCostMatricesInvalidCostMatrices
        function testChangeCostMatricesInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() this.controllerUnderTest.changeCostMatrices(invalidQ, this.R), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeCostMatrices(invalidQ, this.R), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() this.controllerUnderTest.changeCostMatrices(invalidQ, this.R), ...
                expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() this.controllerUnderTest.changeCostMatrices(invalidQ, this.R), ...
                expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() this.controllerUnderTest.changeCostMatrices(this.Q, invalidR), ...
                expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() this.controllerUnderTest.changeCostMatrices(this.Q, invalidR), ...
                expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() this.controllerUnderTest.changeCostMatrices(this.Q, invalidR), ...
                expectedErrId);
        end
        
        %% testChangeCostMatrices
        function testChangeCostMatrices(this)
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            
            % change Q and R by
            newQ = this.Q + eye(this.dimX);
            newR = this.R + eye(this.dimU);
                        
            this.controllerUnderTest.changeCostMatrices(newQ, newR);
            % this call triggers the recomputation of the gain
            this.controllerUnderTest.computeControlSequence(zeroState);
            
            % compute the new expected gain by solving the DARE
            [~, K, ~] = idare(this.A, this.B, newQ, newR);
            expectedNewL = -K;            
            
            this.verifyEqual(this.controllerUnderTest.L, expectedNewL, 'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testChangeModelParametersInvalidSystemMatrix
        function testChangeModelParametersInvalidSystemMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidDimensions';
            
            invalidA = this.A(1,1); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
                        
            invalidA = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
            
            invalidA = this.A(:, 1); % not square
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
            
            invalidA = this.A; % not finite
            invalidA(end) = inf;
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
        end
        
        %% testChangeModelParametersInvalidInputMatrix
        function testChangeModelParametersInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrixDims';
            
            invalidB = this.B(1,1); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
                        
            invalidB = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
                        
            invalidB = this.B; % not finite
            invalidB(end) = inf;
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
        end
        
        %% testChangeModelParameters
        function testChangeModelParameters(this)
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            
            % change Q and R by
            newA = this.A - 0.1 * eye(this.dimX);
            newB = this.B + 2 * eye(this.dimU);
                        
            this.controllerUnderTest.changeModelParameters(newA, newB);
            % this call triggers the recomputation of the gain
            this.controllerUnderTest.computeControlSequence(zeroState);
            
            % compute the new expected gain by solving the DARE
            [~, K, ~] = idare(newA, newB, this.Q, this.R);
            expectedNewL = -K;            
            
            this.verifyEqual(this.controllerUnderTest.L, expectedNewL, 'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testChangeSetpointInvalidSetpoint
        function testChangeSetpointInvalidSetpoint(this)
            expectedErrId = 'NominalPredictiveController:ChangeSetPoint:InvalidSetpoint';
            
            % not a vector
            invalidSetpoint = this.Q;
            this.verifyError(@() this.controllerUnderTest.changeSetPoint(invalidSetpoint), ...
                expectedErrId);
            
            % wrong dimension
            invalidSetpoint = [this.setpoint; 42];
            this.verifyError(@() this.controllerUnderTest.changeSetPoint(invalidSetpoint), ...
                expectedErrId);
            
            % correct dimension, but inf
            invalidSetpoint = this.setpoint;
            invalidSetpoint(2) = inf;
            this.verifyError(@() this.controllerUnderTest.changeSetPoint(invalidSetpoint), ...
                expectedErrId);
        end
        
        %% testChangeSetpoint
        function testChangeSetpoint(this)
            this.assertEmpty(this.controllerUnderTest.setpoint);
            
            this.controllerUnderTest.changeSetPoint(this.setpoint);
            this.verifyEqual(this.controllerUnderTest.setpoint, this.setpoint);
        end
        
         %% testDoControlSequenceComputationZeroState
        function testDoControlSequenceComputationZeroState(this)
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));

            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            this.verifyEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
            
            % now change the sequence length and compute again
            newSeqLength = this.sequenceLength + 2;
            this.controllerUnderTest.changeSequenceLength(newSeqLength);
            
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            this.verifyEqual(actualSequence, zeros(this.dimU * newSeqLength, 1));
        end
        
         %% testDoControlSequenceComputationZeroStateWithSetpoint
        function testDoControlSequenceComputationZeroStateWithSetpoint(this)
            % assert that the set point is changed
            this.controllerUnderTest.changeSetPoint(this.setpoint);
            this.assertEqual(this.controllerUnderTest.setpoint, this.setpoint);
            
            % perform a sanity check: given state is origin, so first
            % element of the sequence should be the feedforward
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));

            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence(1:this.dimU), this.feedforward, ...
                'AbsTol', 2*1e-4);
            
            % now change the sequence length and compute again
            newSeqLength = this.sequenceLength + 2;
            this.controllerUnderTest.changeSequenceLength(newSeqLength);
            
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * newSeqLength 1];
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence(1:this.dimU), this.feedforward, ...
                'AbsTol', 2*1e-4);
        end
        
        %% testDoControlSequenceComputation
        function testDoControlSequenceComputation(this)
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            actualSequence = this.controllerUnderTest.computeControlSequence(state);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            expectedSequence = this.inputTrajectory(:);
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence, expectedSequence, 'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationWithSetpoint
        function testDoControlSequenceComputationWithSetpoint(this)
            % assert that the set point is changed
            this.controllerUnderTest.changeSetPoint(this.setpoint);
            this.assertEqual(this.controllerUnderTest.setpoint, this.setpoint);
            
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            actualSequence = this.controllerUnderTest.computeControlSequence(state);
             % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            expectedSequence = this.inputTrajectoryNonzeroSetpoint(:);
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence, expectedSequence, 'AbsTol', 2*1e-4);
        end
       
        %% testDoStageCostsComputation
        function testDoStageCostsComputation(this)
            state = this.stateTrajectory(:, end-1);
            input = this.inputTrajectory(:, end);
            timestep = size(this.inputTrajectory, 2);
            
            expectedStageCosts = state' * this.Q * state + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testDoStageCostsComputationWithSetpoint
        function testDoStageCostsComputationWithSetpoint(this)
            % assert that the set point is changed
            this.controllerUnderTest.changeSetPoint(this.setpoint);
            this.assertEqual(this.controllerUnderTest.setpoint, this.setpoint);
            
            state = this.stateTrajectory(:, end-1);
            input = this.inputTrajectory(:, end);
            timestep = size(this.inputTrajectory, 2);
            
            expectedStageCosts = (state - this.setpoint)' * this.Q * (state - this.setpoint) ...
                + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            expectedErrId = 'NominalPredictiveController:DoCostsComputation';
            
            % state trajectory too short
            appliedInputs = this.inputTrajectory;
            invalidStateTrajectory = this.stateTrajectory(:, 1:end-2);
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
            
            % state trajectory too long
            appliedInputs = this.inputTrajectory(:, 1:end-2);
            invalidStateTrajectory = this.stateTrajectory;
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
        end
        
        %% testDoCostsComputation
        function testDoCostsComputation(this)
            expectedCosts = this.computeCosts();
            actualCosts = this.controllerUnderTest.computeCosts(this.stateTrajectory, this.inputTrajectory);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testDoCostsComputationWithSetpoint
        function testDoCostsComputationWithSetpoint(this)
            % assert that the set point is changed
            this.controllerUnderTest.changeSetPoint(this.setpoint);
            this.assertEqual(this.controllerUnderTest.setpoint, this.setpoint);
            
            expectedCosts = this.computeCostsNonzeroSetpoint();
            actualCosts = this.controllerUnderTest.computeCosts(this.stateTrajectoryNonzeroSetpoint, this.inputTrajectoryNonzeroSetpoint);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testDoGetDeviationFromRefForState
        function testDoGetDeviationFromRefForState(this)
            state = this.setpoint;
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(state, 1);
            
            % we want to reach the origin, so there should be a deviation
            this.verifyEqual(actualDeviation, this.setpoint);
        end
        
        %% testDoGetDeviationFromRefForStateWithSetpoint
        function testDoGetDeviationFromRefForStateWithSetpoint(this)
            % assert that the set point is changed
            this.controllerUnderTest.changeSetPoint(this.setpoint);
            this.assertEqual(this.controllerUnderTest.setpoint, this.setpoint);
            
            state = this.setpoint;
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(state, 1);
            
            this.verifyEqual(actualDeviation, zeros(this.dimX, 1));
        end
    end
    
end

