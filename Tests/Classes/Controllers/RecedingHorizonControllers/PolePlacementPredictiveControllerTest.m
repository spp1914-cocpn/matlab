classdef PolePlacementPredictiveControllerTest < matlab.unittest.TestCase
    % Test cases for PolePlacementPredictiveController.
    
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
        absTol = 1e-8;
    end
    
    properties (Access = private)
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
        
        samplingInterval;
        polesDisc;
        polesCont;
        expectedL;
        
        % for reaching a nonzero setpoint
        feedforward;
        setpoint;
        stateTrajectoryNonzeroSetpoint;
        inputTrajectoryNonzeroSetpoint;
        
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
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
                 
            this.polesDisc = [0.1-0.1j; 0.1+0.1j]; % arbitrarily chosen, stable poles
            this.expectedL = place(this.A, this.B, this.polesDisc);

            this.initialPlantState = [0.1; 0];            
            this.samplingInterval = 1; % 1 second 
            this.sequenceLength = 10;            
            this.horizonLength = this.sequenceLength;
 
            closedLoopSystem = feedback(ss(this.A, this.B, eye(this.dimX), [], []), ss(this.expectedL), -1); % discrete-time model
            
            % the corresponding continuous-time pole location
            this.polesCont = log(this.polesDisc) / this.samplingInterval; 

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

            this.controllerUnderTest = PolePlacementPredictiveController(this.A, this.B, this.polesCont, this.sequenceLength, this.samplingInterval);
        end        
    end
    
    methods (Access = private)
        %% computeTrajectories
        function [stateTrajectory, inputTrajectory] = computeTrajectories(this, closedLoopSystem)
            [~, ~, stateTrajectory] = initial(closedLoopSystem, this.initialPlantState, this.horizonLength);
            stateTrajectory = stateTrajectory';
            
            % compute the corresponding, expected inputs
            inputTrajectory = -this.expectedL * stateTrajectory(:, 1:this.horizonLength);
        end
        
        %% computeTrajectoriesNonzeroSetpoint
        function [stateTrajectory, inputTrajectory] = computeTrajectoriesNonzeroSetpoint(this, closedLoopSystem)            
                      
            [~,~, stateTrajectory] = lsim(closedLoopSystem, repmat(this.feedforward', this.horizonLength + 1, 1), ...
                0:this.horizonLength, this.initialPlantState);
            
            stateTrajectory = stateTrajectory';
            % compute the correspondig, expected inputs
            inputTrajectory = -this.expectedL * stateTrajectory(:, 1:this.horizonLength) + this.feedforward;            
        end
    end
    
    methods (Test)
         %% testCtorInvalidSysMatrix
        function testCtorInvalidSysMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
            
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() PolePlacementPredictiveController(invalidA, this.B, this.polesCont, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
            
            invalidA = inf(this.dimX, this.dimX); % square but inf
            this.verifyError(@() PolePlacementPredictiveController(invalidA, this.B, this.polesCont, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
        end
        
        %% testCtorInvalidInputMatrix
        function testCtorInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidB = eye(this.dimX + 1, this.dimU); % invalid dims
            this.verifyError(@() PolePlacementPredictiveController(this.A, invalidB, this.polesCont, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
            
            invalidB = eye(this.dimX, this.dimU); % correct dims, but nan
            invalidB(1, 1) = nan;
            this.verifyError(@() PolePlacementPredictiveController(this.A, invalidB, this.polesCont, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
        end
        
        %% testCtorInvalidPolesCont
        function testCtorInvalidPolesCont(this)
            expectedErrId = 'PolePlacementPredictiveController:InvalidPoles';
            
            invalidPoles = eye(this.dimX + 1, this.dimU); % matrix
            this.verifyError(@() PolePlacementPredictiveController(this.A, this.B, invalidPoles, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
            
            invalidPoles = -1; % incorrect dims
            this.verifyError(@() PolePlacementPredictiveController(this.A, this.B, invalidPoles, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
            
            invalidPoles = [-1+1j; -1]; % not a complex pair
            this.verifyError(@() PolePlacementPredictiveController(this.A, this.B, invalidPoles, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
            
            newB = this.B;
            newB(:, 2) = newB(:, 1); % now the rank of B is 1
            invalidPoles = [-1; -1]; % doubles pole is not allowed
            this.verifyError(@() PolePlacementPredictiveController(this.A, newB, invalidPoles, this.sequenceLength, this.samplingInterval), ...
                expectedErrId);
        end       
       
         %% testCtorUnstablePole
        function testCtorUnstablePole(this)
            expectedWarnId = 'PolePlacementPredictiveController:Unstable';
            
            poles = [-1 eps]; % one pole is unstable, warning expected
            this.verifyWarning(@() PolePlacementPredictiveController(this.A, this.B, poles, this.sequenceLength, this.samplingInterval), ...
                expectedWarnId);          
        end 
                
        %% testCtorInvalidSetpoint
        function testCtorInvalidSetpoint(this)
            expectedErrId = [class(this.controllerUnderTest) ':InvalidSetpoint'];
            
            % not a vector
            invalidSetpoint = this;
            this.verifyError(@() PolePlacementPredictiveController(this.A, this.B, this.polesCont, this.sequenceLength, this.samplingInterval, invalidSetpoint), expectedErrId);
            
            % wrong dimension
            invalidSetpoint = [this.setpoint; 42];
            this.verifyError(@() PolePlacementPredictiveController(this.A, this.B, this.polesCont, this.sequenceLength, this.samplingInterval, invalidSetpoint), expectedErrId);
            
            % correct dimension, but inf
            invalidSetpoint = this.setpoint;
            invalidSetpoint(2) = inf;
            this.verifyError(@() PolePlacementPredictiveController(this.A, this.B, this.polesCont, this.sequenceLength, this.samplingInterval, invalidSetpoint), expectedErrId);
        end
        
        %% testCtor
        function testCtor(this)
            controller = PolePlacementPredictiveController(this.A, this.B, this.polesCont, this.sequenceLength, this.samplingInterval);
            
            this.verifyEqual(controller.poles, this.polesDisc, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
            this.verifyEqual(controller.L, this.expectedL, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
            this.verifyEmpty(controller.setpoint);
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback
            
            % now pass an empty setpoint
            controller = PolePlacementPredictiveController(this.A, this.B, this.polesCont, this.sequenceLength, this.samplingInterval, []);
            
            this.verifyEqual(controller.poles, this.polesDisc, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
            this.verifyEqual(controller.L, this.expectedL, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
            this.verifyEmpty(controller.setpoint);
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback
            
            % now a set point
            controller = PolePlacementPredictiveController(this.A, this.B, this.polesCont, this.sequenceLength, this.samplingInterval, this.setpoint);
            
            this.verifyEqual(controller.poles, this.polesDisc,  'AbsTol', PolePlacementPredictiveControllerTest.absTol);
            this.verifyEqual(controller.L, this.expectedL, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
            this.verifyEqual(controller.setpoint, this.setpoint);
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback            
        end
%%
%%
        %% testChangeSequenceLength
        function testChangeSequenceLength(this)
            this.assertEqual(this.controllerUnderTest.sequenceLength, this.sequenceLength);
            
            newSeqLength = this.sequenceLength + 2;
            this.controllerUnderTest.changeSequenceLength(newSeqLength);
            
            this.verifyEqual(this.controllerUnderTest.sequenceLength, newSeqLength);
        end
%%
%%
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
            
            % change A and B
            newA = this.A - 0.1 * eye(this.dimX);
            newB = this.B + 2 * eye(this.dimU);
                        
            this.controllerUnderTest.changeModelParameters(newA, newB);
            % this call triggers the recomputation of the gain
            this.controllerUnderTest.computeControlSequence(zeroState);
            
            % compute the new expected gain using the same pole locations
            expectedNewL = place(newA, newB, this.polesDisc);            
            
            this.verifyEqual(this.controllerUnderTest.L, expectedNewL, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
        end
%%
%%
        %% testChangeSetpointInvalidSetpoint
        function testChangeSetpointInvalidSetpoint(this)
            expectedErrId = 'PolePlacementPredictiveController:ChangeSetPoint:InvalidSetpoint';
            
            % not a vector
            invalidSetpoint = this.A;
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
%%
%%
        %% testChangeSamplingInterval
        function testChangeSamplingInterval(this)
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            
            % use new sampling interval, impacts the pole locations
            newSamplingInterval = this.samplingInterval / 2;
                        
            this.controllerUnderTest.changeSamplingInterval(newSamplingInterval);
            % this call triggers the recomputation of the gain
            this.controllerUnderTest.computeControlSequence(zeroState);
            
            % compute the new expected gain using the new pole locations
            newPoles = exp(this.polesCont .* newSamplingInterval);
            expectedNewL = place(this.A, this.B, newPoles);            
            
            this.verifyEqual(this.controllerUnderTest.L, expectedNewL, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
        end
%%
%%
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
            % assert that the set point is set
            controller = PolePlacementPredictiveController(this.A, this.B, this.polesCont, ...
                this.sequenceLength, this.samplingInterval, this.setpoint);
            this.assertEqual(controller.setpoint, this.setpoint);
            
            % perform a sanity check: given state is origin, so first
            % element of the sequence should be the feedforward
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));

            actualSequence = controller.computeControlSequence(zeroState);
            
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence(1:this.dimU), this.feedforward, ...
               'AbsTol', PolePlacementPredictiveControllerTest.absTol);
            
            % now change the sequence length and compute again
            newSeqLength = this.sequenceLength + 2;
            controller.changeSequenceLength(newSeqLength);
            
            actualSequence = controller.computeControlSequence(zeroState);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * newSeqLength 1];
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence(1:this.dimU), this.feedforward, ...
                'AbsTol', PolePlacementPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputation
        function testDoControlSequenceComputation(this)
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            actualSequence = this.controllerUnderTest.computeControlSequence(state);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            expectedSequence = this.inputTrajectory(:);
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence, expectedSequence, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
        end
        
         %% testDoControlSequenceComputationWithSetpoint
        function testDoControlSequenceComputationWithSetpoint(this)
            controller = PolePlacementPredictiveController(this.A, this.B, this.polesCont, ...
                this.sequenceLength, this.samplingInterval, this.setpoint);
            this.assertEqual(controller.setpoint, this.setpoint);
            
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            actualSequence = controller.computeControlSequence(state);
             % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            expectedSequence = this.inputTrajectoryNonzeroSetpoint(:);
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence, expectedSequence, 'AbsTol', PolePlacementPredictiveControllerTest.absTol);
        end
%%
%%
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
            controller = PolePlacementPredictiveController(this.A, this.B, this.polesCont, ...
                this.sequenceLength, this.samplingInterval, this.setpoint);
            this.assertEqual(controller.setpoint, this.setpoint);
            
            state = this.setpoint;
            actualDeviation = controller.getDeviationFromRefForState(state, 1);
            
            this.verifyEqual(actualDeviation, zeros(this.dimX, 1));
        end
    end
end

