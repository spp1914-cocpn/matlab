classdef LinearlyConstrainedPredictiveControllerTest < matlab.unittest.TestCase
    % Test cases for LinearlyConstrainedPredictiveController.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties
        A;
        B;
        Q;
        R;
        
        Z;
        refTrajectory;
        Qref;
        stateTrajectoryRef;
        inputTrajectoryRef;
        
        sequenceLength;
        dimX;
        dimU;
        initialPlantState;
        stateTrajectory;
        inputTrajectory;
        horizonLength;
        
        stateConstraintWeightings;
        stateConstraints;
        inputConstraintWeightings;
        inputConstraints;
        
        caDelayProbs;
        
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
        %% initproperties
        function initProperties(this)
            % use (noise-free) stirred tank example (Example 6.15, p. 500-501) from
            %
            % Huibert Kwakernaak, and Raphael Sivan, 
            % Linear Optimal Control Systems,
            % Wiley-Interscience, New York, 1972.
            %
            
            this.sequenceLength = 10;
            this.dimX = 2;
            this.dimU = 2;
            this.horizonLength = this.sequenceLength;
            this.initialPlantState = [0.1; 0];
            
            V = diag([0.01, 1]); % in the book: z = V*x
            this.A = diag([0.9512, 0.9048]);
           
            this.B = [4.877 4.877; -1.1895 3.569];
            this.Q = V * diag([50 0.02]) * V; % R_3 in the book
            this.R = diag([1/3, 3]); % R_2 in the book   
            
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            % this constraints should actually not matter
            this.stateConstraintWeightings = ones(this.dimX, 1);
            this.stateConstraints = 10000;
            this.inputConstraintWeightings = ones(this.dimU, 1);
            this.inputConstraints = 10000;           
           
            % fixed delay is 1 time step
            this.caDelayProbs = [0 1 0 0 0 0];
            
            % only the first state variable is to be driven to the origin
            this.Z = [1 0];
            this.Qref = this.Q(1,1);
            this.refTrajectory = zeros(1, 2 * this.horizonLength);
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
        end
    end
    
    methods (Access = private)        
        %% setupControllerUnderTest
        function setupControllerUnderTest(this, doTrackRef)
            if doTrackRef
                this.controllerUnderTest = LinearlyConstrainedPredictiveController(this.A, this.B, this.Qref, this.R, this.sequenceLength, ...
                    this.caDelayProbs, this.stateConstraintWeightings, this.stateConstraints, this.inputConstraintWeightings, this.inputConstraints, ...
                    this.Z, this.refTrajectory);
            else
                this.controllerUnderTest = LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength, ...
                    this.caDelayProbs, this.stateConstraintWeightings, this.stateConstraints, this.inputConstraintWeightings, this.inputConstraints);
            end
        end
        
        %% computeTrajectories
        function [stateTrajectory, inputTrajectory] = computeTrajectories(this)
            % without effective constraints, problem reduces to LQR                        
            L = dlqr(this.A, this.B, this.Q, this.R);
            
            stateTrajectory = zeros(this.dimX, this.horizonLength + 1);
            inputTrajectory = zeros(this.dimU, this.horizonLength);
            % we have defined the delay probabilities such that every
            % packet is delayed 1 time step
            % so the expected inputs for the prediction of the state are: 
            % at stage 1 to compute x_k+1|k: zero ->"from previous packet"
            % at stage i to compute x_k+i+1|k: u_k+i|k -> "from the current packet"
            stateTrajectory(:, 1) = this.initialPlantState;
            stateTrajectory(:, 2) = this.A * stateTrajectory(:, 1);
            % u_k|k remains zero, since delay is 1 so this entry will never
            % become active
            for j = 2:this.horizonLength
                % compute new input
                inputTrajectory(:, j) = -L * stateTrajectory(:, j);
                % and predict using "expected" input
                stateTrajectory(:, j + 1) = this.A * stateTrajectory(:, j) + this.B * inputTrajectory(:, j);
            end            
        end
        
        %% computeTrajectoriesRef
        function [stateTrajectoryRef, inputTrajectoryRef] = computeTrajectoriesRef(this)
            % without effective constraints, problem reduces to LQR with reference weighting            
            discSys = ss(this.A, this.B, this.Z, [], -1);
            L = lqry(discSys, this.Qref, this.R);
            
            stateTrajectoryRef = zeros(this.dimX, this.horizonLength + 1);
            inputTrajectoryRef = zeros(this.dimU, this.horizonLength);
            % we have defined the delay probabilities such that every
            % packet is delayed 1 time step
            % so the expected inputs for the prediction of the state are: 
            % at stage 1 to compute x_k+1|k: zero ->"from previous packet"
            % at stage i to compute x_k+i+1|k: u_k+i|k -> "from the current packet"
            stateTrajectoryRef(:, 1) = this.initialPlantState;
            stateTrajectoryRef(:, 2) = this.A * stateTrajectoryRef(:, 1);
            % u_k|k remains zero, since delay is 1 so this entry will never
            % become active
            for j = 2:this.horizonLength
                % compute new input
                inputTrajectoryRef(:, j) = -L * stateTrajectoryRef(:, j);
                % and predict using "expected" input
                stateTrajectoryRef(:, j + 1) = this.A * stateTrajectoryRef(:, j) + this.B * inputTrajectoryRef(:, j);
            end      
        end
        
        %% computeCosts
        function costs = computeCosts(this)
            % compute the costs in a straightforward manner            
            costs =  this.stateTrajectory(:, end)' * this.Q * this.stateTrajectory(:, end);
            for j = 1:this.horizonLength
                costs = costs + this.stateTrajectory(:, j)' * this.Q * this.stateTrajectory(:, j) ...
                    + this.inputTrajectory(:, j)' * this.R * this.inputTrajectory(:, j);
            end
        end
        
        %% computeCostsRef
        function costsRef = computeCostsRef(this)
            % compute the costs in a straightforward manner
            performance = this.Z * this.stateTrajectoryRef(:, 1:end) - this.refTrajectory(:, 1:this.horizonLength + 1);
            costsRef =  performance(:, end)' * this.Qref * performance(:, end);
            for j = 1:this.horizonLength
                costsRef = costsRef + performance(:, j)' * this.Qref * performance(:, j) ...
                    + this.inputTrajectoryRef(:, j)' * this.R * this.inputTrajectoryRef(:, j);
            end
        end
    end
    
    methods (Test)
        %% testLinearlyConstrainedPredictiveControllerInvalidSysMatrix
        function testLinearlyConstrainedPredictiveControllerInvalidSysMatrix(this)
            delayProbs = [0.25 0.25 0.25 0.25];
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
                        
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() LinearlyConstrainedPredictiveController(invalidA, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            invalidA = inf(this.dimX, this.dimX); % square but inf
            this.verifyError(@() LinearlyConstrainedPredictiveController(invalidA, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
        end
        
        %% testLinearlyConstrainedPredictiveControllerInvalidInputMatrix
        function testLinearlyConstrainedPredictiveControllerInvalidInputMatrix(this)
            delayProbs = [0.25 0.25 0.25 0.25];
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidB = eye(this.dimX + 1, this.dimU); % invalid dims
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, invalidB, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            invalidB = eye(this.dimX, this.dimU); % correct dims, but nan
            invalidB(1, 1) = nan;
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, invalidB, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
        end
        
        %% testLinearlyConstrainedPredictiveControllerInvalidCostMatrices
        function testLinearlyConstrainedPredictiveControllerInvalidCostMatrices(this)
            delayProbs = [0.25 0.25 0.25 0.25];
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, invalidR, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, invalidR, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, invalidR, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
        end
        
        
        %% testLinearlyConstrainedPredictiveControllerInvalidStateConstraints
        function testLinearlyConstrainedPredictiveControllerInvalidStateConstraints(this) %#ok
            delayProbs = [0.25 0.25 0.25 0.25];
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateStateConstraints:InvalidStateWeightings';
            
            invalidStateConstraintWeightings = ones(this.dimX +1 , 2); % incorrect dims
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, invalidStateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateStateConstraints:InvalidStateConstraints';
            
            invalidStateConstraints = eye(this.dimX); % not a vector
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, invalidStateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            invalidStateConstraints = ones(this.dimX, 1); % incorrect dims
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, invalidStateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
        end
        
        %% testLinearlyConstrainedPredictiveControllerInvalidInputConstraints
        function testLinearlyConstrainedPredictiveControllerInvalidInputConstraints(this) %#ok
            delayProbs = [0.25 0.25 0.25 0.25];
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateInputConstraints:InvalidInputWeightings';
            
            invalidInputConstraintWeightings = ones(this.dimU +1 , 2); % incorrect dims
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                invalidInputConstraintWeightings, this.inputConstraints), ...
                expectedErrId);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateInputConstraints:InvalidInputConstraints';
            
            invalidInputConstraints = eye(this.dimU); % not a vector
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, invalidInputConstraints), ...
                expectedErrId);           
             
            invalidInputConstraints = ones(this.dimU, 1); % incorrect dims
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, invalidInputConstraints), ...
                expectedErrId);
        end
        
        %% testLinearlyConstrainedPredictiveControllerInvalidNumArgs
        function testLinearlyConstrainedPredictiveControllerInvalidNumArgs(this)
            delayProbs = [0.25 0.25 0.25 0.25];
            expectedErrId = 'LinearlyConstrainedPredictiveController:InvalidNumberOfArguments';
            
            % use 11 arguments, we expect either 10 or 12
            this.verifyError(@() LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, delayProbs, this.stateConstraintWeightings, this.stateConstraints, ...
                this.inputConstraintWeightings, this.inputConstraints, []), ...
                expectedErrId);
        end
%%
%%
        %% testLinearlyConstrainedPredictiveController
        function testLinearlyConstrainedPredictiveController(this)
            delayProbs = [0.25 0.25 0.25 0.25];
            % should not crash
            controller = LinearlyConstrainedPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength, delayProbs, ...
                this.stateConstraintWeightings, this.stateConstraints, this.inputConstraintWeightings, this.inputConstraints);
              
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback
            this.verifyEqual(controller.horizonLength, this.sequenceLength); % horizon is by default equal to sequence length
        end
        
        %% testGetStateConstraints
        function testGetStateConstraints(this)
            this.setupControllerUnderTest(false);
            
            [actualConstraints, actualWeightings] = this.controllerUnderTest.getStateConstraints();
            this.verifyEqual(actualConstraints, this.stateConstraints);
            this.verifyEqual(actualWeightings, this.stateConstraintWeightings);
        end
        
        %% testGetInputConstraints
        function testGetInputConstraints(this)
            this.setupControllerUnderTest(false);
            
            [actualConstraints, actualWeightings] = this.controllerUnderTest.getInputConstraints();
            this.verifyEqual(actualConstraints, this.inputConstraints);
            this.verifyEqual(actualWeightings, this.inputConstraintWeightings);
        end
                
        %% testChangeStateConstraintsInvalidConstraints
        function testChangeStateConstraintsInvalidConstraints(this)
            this.setupControllerUnderTest(false);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateStateConstraints:InvalidStateWeightings';
            
            newStateConstraints = this.stateConstraints + 10;
            invalidStateConstraintWeightings = ones(this.dimX +1 , 2); % incorrect dims
            this.verifyError(@() this.controllerUnderTest.changeStateConstraints(invalidStateConstraintWeightings, newStateConstraints), ...
                expectedErrId);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateStateConstraints:InvalidStateConstraints';
            
            newStateConstraintWeightings = this.stateConstraintWeightings + 42;
            invalidStateConstraints = eye(this.dimX); % not a vector
            this.verifyError(@() this.controllerUnderTest.changeStateConstraints(newStateConstraintWeightings, invalidStateConstraints), ...
                expectedErrId);
            
            invalidStateConstraints = ones(this.dimX, 1); % incorrect dims
            this.verifyError(@() this.controllerUnderTest.changeStateConstraints(newStateConstraintWeightings, invalidStateConstraints), ...
                expectedErrId);
        end
        
        %% testChangeStateConstraintsNoConstraints
        function testChangeStateConstraintsNoConstraints(this)
            this.setupControllerUnderTest(false);
            
            newStateConstraints = [];
            newStateConstraintWeightings = this.stateConstraintWeightings + 42; % ignored 
            this.controllerUnderTest.changeStateConstraints(newStateConstraintWeightings, newStateConstraints);
            
            [actualConstraints, ~] = this.controllerUnderTest.getStateConstraints();
            this.verifyEmpty(actualConstraints);            
        end
        
        %% testChangeStateConstraints
        function testChangeStateConstraints(this)
            this.setupControllerUnderTest(false);
            
            newStateConstraints = this.stateConstraints + 10;
            newStateConstraintWeightings = this.stateConstraintWeightings + 42;
            this.controllerUnderTest.changeStateConstraints(newStateConstraintWeightings, newStateConstraints);
            
            [actualConstraints, actualWeightings] = this.controllerUnderTest.getStateConstraints();
            this.verifyEqual(actualConstraints, newStateConstraints);
            this.verifyEqual(actualWeightings, newStateConstraintWeightings);
        end
%%
%%
         %% testChangeInputConstraintsInvalidConstraints
        function testChangeInputConstraintsInvalidConstraints(this)
            this.setupControllerUnderTest(false);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateInputConstraints:InvalidInputWeightings';
            
            invalidInputConstraintWeightings = ones(this.dimU +1 , 2); % incorrect dims
            newInputConstraints = this.inputConstraints + 10;
            this.verifyError(@() this.controllerUnderTest.changeInputConstraints(invalidInputConstraintWeightings, newInputConstraints), ...
                expectedErrId);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:ValidateInputConstraints:InvalidInputConstraints';
            
            newInputConstraintWeightings = this.inputConstraintWeightings + 42;
            invalidInputConstraints = eye(this.dimU); % not a vector
            this.verifyError(@() this.controllerUnderTest.changeInputConstraints(newInputConstraintWeightings, invalidInputConstraints), ...
                expectedErrId);
            
            invalidInputConstraints = ones(this.dimU, 1); % incorrect dims
            this.verifyError(@() this.controllerUnderTest.changeInputConstraints(newInputConstraintWeightings, invalidInputConstraints), ...
                expectedErrId);
        end
%%
%%
        %% testDoControlSequenceComputationInvalidTimestep
        function testDoControlSequenceComputationInvalidTimestep(this)
            this.setupControllerUnderTest(false);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:DoControlSequenceComputation:InvalidTimestep';
            
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            mode = 1; % not needed but has to be passed
            
            invalidTimestep = -1; % not positive
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(zeroState, mode, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = 1.5; % not an integer            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(zeroState, mode, invalidTimestep), ...
                expectedErrId);
        end

        %% testDoControlSequenceComputationZeroState
        function testDoControlSequenceComputationZeroState(this)
            this.setupControllerUnderTest(false);
            
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            expectedSequence = zeros(this.dimU * this.horizonLength, 1);

            actualSequence = this.verifyWarningFree(@() this.controllerUnderTest.computeControlSequence(zeroState, mode, timestep));
            
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', LinearlyConstrainedPredictiveControllerTest.absTol);      
        end
        
         %% testDoControlSequenceComputationZeroStateRef
        function testDoControlSequenceComputationZeroStateRef(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            timestep = 1;
            
            % assert that initial difference to ref is zero
            state = [this.refTrajectory(:, timestep); 42];
            this.assertTrue(this.Z * state == this.refTrajectory(:, timestep));
                        
            mode = 1; % not needed but has to be passed
            expectedSequence = zeros(this.dimU * this.horizonLength, 1);

            actualSequence = this.verifyWarningFree(@() this.controllerUnderTest.computeControlSequence(Gaussian(state, eye(this.dimX)), mode, timestep));
            
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', LinearlyConstrainedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationInfeasibleState
        function testDoControlSequenceComputationInfeasibleState(this)
            this.setupControllerUnderTest(false);
            
            % perform a sanity check: the initial state violates the
            % constraint imposed on the state
            state = [this.stateConstraints; this.stateConstraints];
            this.assertFalse(dot(this.stateConstraintWeightings, state) <= this.stateConstraints);
            
            infeasibleState = Gaussian(state, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            expectedWarnId = 'LinearlyConstrainedPredictiveController:DoControlSequenceComputation:ProblemInfeasible';
            expectedSequence = zeros(this.dimU * this.horizonLength, 1);
            % infeasible problem: issue warning and return zero inputs
            actualSequence = this.verifyWarning(@() this.controllerUnderTest.computeControlSequence(infeasibleState, mode, timestep), ...
                expectedWarnId);
            
            this.verifyEqual(actualSequence, expectedSequence);          
        end
        
        %% testDoControlSequenceComputationInfeasibleStateRef
        function testDoControlSequenceComputationInfeasibleStateRef(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            
            % perform a sanity check: the initial state violates the
            % constraint imposed on the state
            state = [this.stateConstraints; this.stateConstraints];
            this.assertFalse(dot(this.stateConstraintWeightings, state) <= this.stateConstraints);
            
            infeasibleState = Gaussian(state, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            expectedWarnId = 'LinearlyConstrainedPredictiveController:DoControlSequenceComputation:ProblemInfeasible';
            expectedSequence = zeros(this.dimU * this.horizonLength, 1);
            % infeasible problem: issue warning and return zero inputs
            actualSequence = this.verifyWarning(@() this.controllerUnderTest.computeControlSequence(infeasibleState, mode, timestep), ...
                expectedWarnId);
            
            this.verifyEqual(actualSequence, expectedSequence);          
        end
        
        %% testDoControlSequenceComputation
        function testDoControlSequenceComputation(this)
            this.setupControllerUnderTest(false);
            
            % assert that the initial state is feasible
            this.assertTrue(dot(this.stateConstraintWeightings, this.initialPlantState) <= this.stateConstraints);
            
            expectedSequence = this.inputTrajectory(:);
            
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', 2e-5);
        end
        
        %% testDoControlSequenceComputationRef
        function testDoControlSequenceComputationRef(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            
            % assert that the initial state is feasible
            this.assertTrue(dot(this.stateConstraintWeightings, this.initialPlantState) <= this.stateConstraints);
            
            expectedSequence = this.inputTrajectoryRef(:);
            
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', 1.5e-3);
        end
%%
%%
        %% testDoStageCostsComputation
        function testDoStageCostsComputation(this)
            this.setupControllerUnderTest(false);
            
            state = this.stateTrajectory(:, end-1);
            input = this.inputTrajectory(:, end);
            timestep = size(this.inputTrajectory, 2);
            
            expectedStageCosts = state' * this.Q * state + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', LinearlyConstrainedPredictiveControllerTest.absTol);
        end
        
        %% testDoStageCostsComputation
        function testDoStageCostsComputationRef(this)
            % include ref tracking
            this.setupControllerUnderTest(true);
            
            state = this.stateTrajectoryRef(:, end-1);
            input = this.inputTrajectoryRef(:, end);
            timestep = size(this.inputTrajectoryRef, 2);
            
            expectedStageCosts = (this.Z * state - this.refTrajectory(:, timestep))' * this.Qref * (this.Z * state - this.refTrajectory(:, timestep)) ...
                + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', LinearlyConstrainedPredictiveControllerTest.absTol);
        end
%%
%%
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            this.setupControllerUnderTest(false);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:DoCostsComputation:InvalidStateTrajectory';
            
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
            this.setupControllerUnderTest(false);
            
            expectedCosts = this.computeCosts();
            actualCosts = this.controllerUnderTest.computeCosts(this.stateTrajectory, this.inputTrajectory);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', LinearlyConstrainedPredictiveControllerTest.absTol);
        end
        
        %% testDoCostsComputationRef
        function testDoCostsComputationRef(this)
            this.setupControllerUnderTest(true);
            
            expectedCosts = this.computeCostsRef();
            actualCosts = this.controllerUnderTest.computeCosts(this.stateTrajectoryRef, this.inputTrajectoryRef);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', LinearlyConstrainedPredictiveControllerTest.absTol);
        end
%%
%%
        %% testDoGetDeviationFromRefForState
        function testDoGetDeviationFromRefForState(this)
            this.setupControllerUnderTest(false);
            
            state = this.initialPlantState;
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(state, 1);
            
            % we want to reach the origin, so there should be a deviation
            this.verifyEqual(actualDeviation, state);
        end
        
         %% testDoGetDeviationFromRefForStateRef
        function testDoGetDeviationFromRefForStateRef(this)
            this.setupControllerUnderTest(true);
            
                       
            state = this.initialPlantState;
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(state, 1);
            % we want to reach the origin with the first variable, so there should be a deviation
            this.verifyEqual(actualDeviation, state(1));
        end
        
        %% testDoGetDeviationFromRefForStateInvalidTimestep
        function testDoGetDeviationFromRefForStateInvalidTimestep(this)
            this.setupControllerUnderTest(true);
            
            expectedErrId = 'LinearlyConstrainedPredictiveController:GetDeviationFromRefForState:InvalidTimestep';
            
            state = this.initialPlantState;
            invalidTimestep = -1; % not positive
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(state, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = 1.5; % not an integer
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(state, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = size(this.refTrajectory, 2) + 1; % exceed the length of the reference trajectory
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(state, invalidTimestep), ...
                expectedErrId);
        end
    end
    
end

