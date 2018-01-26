classdef FiniteHorizonControllerTest < BaseFiniteHorizonControllerTest
    % Test cases for FiniteHorizonController.
    
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
        numConstraints;
        stateConstraintWeightings;
        inputConstraintWeightings;
        constraintBounds;
    end
    
    methods (Access = private)
        %% computeExpectedInputsWithoutConstraints
        function inputs = computeExpectedInputsWithoutConstraints(this)
            % in this setup, the current input must arrive without delay,
            % or plant evolves open-loop
            % hence, in this setup, G and F are not existent
            % likewise, H
            % due to the structure of the transition matrix,
            % and the one step horizon, the mode-dependent inputs should not
            % differ
                        
            
            % first mode
            R_1 = this.R;
            Q_1 = this.Q;
            A_1 = this.A;
            B_1 = this.B;
                       
            % second mode
            R_2 = zeros(this.dimU);
            Q_2 = this.Q;
            A_2 = this.A;
            B_2 = zeros(this.dimX, this.dimU);
            
            % input for both modes
            firstPart = pinv(this.transitionMatrix(1,1) * (R_1 + B_1' * Q_1 * B_1) ... 
                + this.transitionMatrix(1,2) * (R_2 + B_2' * Q_2 * B_2));
            secondPart = this.transitionMatrix(1,1) * B_1' * Q_1 * A_1 + ...
                this.transitionMatrix(1,2) * B_2' * Q_2 * A_2;
            
            inputs = -firstPart * secondPart * this.state;
        end
        
        %% computeExpectedInputsWithConstraints
        function input = computeExpectedInputsWithConstraints(this)
                       
            % now we must compute the input including the feedforward part due to the
            % constraints
            % due to the structure of the transition matrix,
            % and the one step horizon, the mode-dependent feedforward inputs should not
            % differ
            
            % first mode
            R_1 = this.R;
            Q_1 = this.Q;
            A_1 = this.A;
            B_1 = this.B;
            R_tilde_1 = this.inputConstraintWeightings;
            Q_tilde_1 = this.stateConstraintWeightings;
            
            % second mode
            R_2 = zeros(this.dimU);
            Q_2 = this.Q;
            A_2 = this.A;
            B_2 = zeros(this.dimX, this.dimU);
            R_tilde_2 = zeros(this.dimU, 1, this.numConstraints);
            Q_tilde_2 = this.stateConstraintWeightings;
            
            P_tilde_0 = this.transitionMatrix(1,1) * (A_1' * squeeze(Q_tilde_1(:, 2, :)) + squeeze(Q_tilde_1(:, 1, :))) ...
                + this.transitionMatrix(1,2) * (A_2' * squeeze(Q_tilde_2(:, 2, :)) + squeeze(Q_tilde_2(:, 1, :))) ...
                - (this.transitionMatrix(1,1) * A_1' * Q_1 * B_1 + this.transitionMatrix(1,2) * A_2' * Q_2 * B_2) ...
                * pinv(this.transitionMatrix(1,1) * (R_1 + B_1' * Q_1 * B_1) + this.transitionMatrix(1,2) * (R_2 + B_2' * Q_2 * B_2)) ...
                * (this.transitionMatrix(1,1) * (B_1' * squeeze(Q_tilde_1(:, 2, :)) + squeeze(R_tilde_1(:, 1, :)))...
                + this.transitionMatrix(1,2) * (B_2' * squeeze(Q_tilde_2(:, 2, :)) + squeeze(R_tilde_2(:, 1, :))));
            
            firstPart =  (this.transitionMatrix(1,1) * (B_1' * squeeze(Q_tilde_1(:, 2, :)) + squeeze(R_tilde_1(:, 1, :)))...
                + this.transitionMatrix(1,2) * (B_2' * squeeze(Q_tilde_2(:, 2, :)) + squeeze(R_tilde_2(:, 1, :))))' ...
                * pinv(this.transitionMatrix(1,1) * (R_1 + B_1' * Q_1 * B_1) + this.transitionMatrix(1,2) * (R_2 + B_2' * Q_2 * B_2)) ...
                * (this.transitionMatrix(1,1) * (B_1' * squeeze(Q_tilde_1(:, 2, :)) + squeeze(R_tilde_1(:, 1, :)))...
                + this.transitionMatrix(1,2) * (B_2' * squeeze(Q_tilde_2(:, 2, :)) + squeeze(R_tilde_2(:, 1, :))));
            
            S_0 = - firstPart;
            lambda_0 = ones(this.numConstraints, 1);
            lb = zeros(this.numConstraints, 1);
            ub = inf(this.numConstraints, 1);
            options = optimoptions(@fmincon, 'Algorithm', 'interior-point', ...
                'OptimalityTolerance', 1e-21, 'Display', 'off');
            lambda = fmincon(@(l) -l'* S_0 * l - this.state' * P_tilde_0 * l + reshape(this.constraintBounds, 1, this.numConstraints) * l, ...
                lambda_0, [], [], [], [], lb, ub, [], options);
            
            feedforward = -pinv(this.transitionMatrix(1,1) * (R_1 + B_1' * Q_1 * B_1) + this.transitionMatrix(1,2) * (R_2 + B_2' * Q_2 * B_2)) ...
                * (this.transitionMatrix(1,1) * (B_1' * squeeze(Q_tilde_1(:, 2, :)) + squeeze(R_tilde_1(:, 1, :)))...
                + this.transitionMatrix(1,2) * (B_2' * squeeze(Q_tilde_2(:, 2, :)) + squeeze(R_tilde_2(:, 1, :))));
            input = this.computeExpectedInputsWithoutConstraints() + feedforward * lambda;
        end
             
    end
    
    methods (Access = protected)
        %% initControllerUnderTest
        function controller = initControllerUnderTest(this)
            controller = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength);
        end
        
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            this.initAdditionalProperties@BaseFiniteHorizonControllerTest();
            this.numConstraints = 5;
            this.stateConstraintWeightings = ones(this.dimX, 2, this.numConstraints);
            this.inputConstraintWeightings = ones(this.dimU, 1, this.numConstraints);
            this.constraintBounds = ones(1, this.numConstraints);
            for j = 1:this.numConstraints
                this.stateConstraintWeightings(:, :, j) ...
                    = this.stateConstraintWeightings(:, :, j) * j; 
                this.inputConstraintWeightings(:, :, j) ...
                    = this.inputConstraintWeightings(:, :, j) * j;
                this.constraintBounds(j) =  this.constraintBounds(j) * j;
            end
        end
        
        %% computeExpectedCosts
        function expectedCosts = computeExpectedCosts(this, states, inputs)
             expectedCosts =  states(:, 1)' * this.Q * states(:, 1) ...
                + states(:, 2)' * this.Q * states(:, 2) + inputs' * this.R * inputs;
        end
    end
%%
%%
    methods (Test)
        %% testFiniteHorizonControllerInvalidNumberOfArgs
        function testFiniteHorizonControllerInvalidNumberOfArgs(this)
            expectedErrId = 'FiniteHorizonController:InvalidNumberOfArguments';
            
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength), expectedErrId);
                 
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings), expectedErrId);
            
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, this.inputConstraintWeightings), expectedErrId);
        end
    
        %% testFiniteHorizonControllerInvalidSystemMatrix
        function testFiniteHorizonControllerInvalidSystemMatrix(this)
             expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
             
             invalidSysMatrix = eye(this.dimX, this.dimX + 1); % not square
             this.verifyError(@() FiniteHorizonController(invalidSysMatrix, this.B, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId);
             
             invalidSysMatrix = eye(this.dimX, this.dimX); % square but not finite
             invalidSysMatrix(1, end) = inf;
             this.verifyError(@() FiniteHorizonController(invalidSysMatrix, this.B, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId);
        end
        
        %% testFiniteHorizonControllerInvalidInputMatrix
        function testFiniteHorizonControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidInputMatrix = eye(this.dimX +1, this.dimU); % invalid dims
            this.verifyError(@() FiniteHorizonController(this.A, invalidInputMatrix, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId);
             
            invalidInputMatrix = eye(this.dimX, this.dimU); % correct dims, but not finite
            invalidInputMatrix(1, end) = nan;
            this.verifyError(@() FiniteHorizonController(this.A, invalidInputMatrix, this.Q, this.R, ...
                 this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId);
            
        end
        
        %% testFiniteHorizonControllerInvalidHorizonLength
        function testFiniteHorizonControllerInvalidHorizonLength(this)
            expectedErrId = 'Validator:ValidateHorizonLength:InvalidHorizonLength';
            
            invalidHorizonLength = eye(3); % not a scalar
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength), expectedErrId);
            
            invalidHorizonLength = 0; % not a positive scalar
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength), expectedErrId);
            
            invalidHorizonLength = 1.5; % not an integer
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength), expectedErrId);
            
            invalidHorizonLength = inf; % not finite
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength), expectedErrId);
        end
        
        %% testFiniteHorizonControllerInvalidCostMatrices
        function testFiniteHorizonControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() FiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() FiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() FiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = eye(this.dimX); % Q is not symmetric
            invalidQ(1, end) = 1;
            this.verifyError(@() FiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() FiniteHorizonController(this.A, this.B, invalidQ, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, invalidR, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, invalidR, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, invalidR, ...
                this.delayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
        end
        
        %% testFiniteHorizonControllerInvalidDelayProbs
        function testFiniteHorizonControllerInvalidDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength, this.horizonLength), expectedErrId)
        end
%%
%%
        %% testFiniteHorizonControllerInvalidLinearIntegralConstraints
        function testFiniteHorizonControllerInvalidLinearIntegralConstraints(this)
            expectedErrId = 'FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidStateWeightings';
            
            invalidStateWeightings = ones(this.dimX, this.horizonLength + 1, this.numConstraints, this.numConstraints);
            % not a 3D matrix
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                invalidStateWeightings, this.inputConstraintWeightings, this.constraintBounds), expectedErrId);
            
            invalidStateWeightings = ones(this.dimX, this.horizonLength, this.numConstraints);
            % not 3D matrix, but incorrect dimension
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                invalidStateWeightings, this.inputConstraintWeightings, this.constraintBounds), expectedErrId);
            
            invalidStateWeightings = ones(this.dimX + 2, this.horizonLength + 1, this.numConstraints);
            % not 3D matrix, but incorrect dimension
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                invalidStateWeightings, this.inputConstraintWeightings, this.constraintBounds), expectedErrId);
            
            expectedErrId = 'FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidInputWeightings';
            
            invalidInputWeightings = ones(this.dimX, this.horizonLength, this.numConstraints);
            % wrong dimension
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, invalidInputWeightings, this.constraintBounds), expectedErrId);
            
            invalidInputWeightings = ones(this.dimU, this.horizonLength + 1, this.numConstraints);
            % wrong dimension
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, invalidInputWeightings, this.constraintBounds), expectedErrId);
            
            invalidInputWeightings = ones(this.dimU, this.horizonLength, this.numConstraints + 1);
            % wrong number of slices
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, invalidInputWeightings, this.constraintBounds), expectedErrId);
            
            expectedErrId = 'FiniteHorizonController:ValidateLinearIntegralConstraints:InvalidConstraints';
            
            invalidConstraintBounds = ones(this.dimX, this.horizonLength, this.numConstraints);
            % not a vector
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, this.inputConstraintWeightings, invalidConstraintBounds), expectedErrId);
            
            invalidConstraintBounds = ones(1, this.numConstraints + 2);
            % wrong dimension
            this.verifyError(@() FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, this.inputConstraintWeightings, invalidConstraintBounds), expectedErrId);
        end
%%
%%
        %% testFiniteHorizonController
        function testFiniteHorizonController(this)
            % create a controller without constraints
            controller = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength);
            
            this.verifyEqual(controller.horizonLength, this.horizonLength);
            this.verifyFalse(controller.constraintsPresent);
            
            % now one with constraints given
            controller = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, this.inputConstraintWeightings, this.constraintBounds);
            
            this.verifyEqual(controller.horizonLength, this.horizonLength);
            this.verifyTrue(controller.constraintsPresent);
            
            % now the border case where only a single constraint is present
            stateWeightings = squeeze(this.stateConstraintWeightings(:, :, 1));
            inputWeightings = squeeze(this.inputConstraintWeightings(:, :, 1));
            controller = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                stateWeightings, inputWeightings, this.constraintBounds(1));
            
            this.verifyEqual(controller.horizonLength, this.horizonLength);
            this.verifyTrue(controller.constraintsPresent);
        end
%%        
%%               
        %% testDoControlSequenceComputationZeroStateWithoutConstraints
        function testDoControlSequenceComputationZeroStateWithoutConstraints(this)
            % perform a sanity check: given state is origin, so compute
            % control sequence (length 1) should be also the zero vector
            % due to the underlying linear control law, independent of the
            % previous mode
            expectedSequence = zeros(this.dimU, 1);

            % first mode
            actualSequence = this.controllerUnderTest.computeControlSequence(this.zeroStateDistribution, 1, ...
                    this.horizonLength);
            
            this.verifyEqual(actualSequence, expectedSequence);
            
            %second mode
            actualSequence = this.controllerUnderTest.computeControlSequence(this.zeroStateDistribution, 2, ...
                    this.horizonLength);
            
            this.verifyEqual(actualSequence, expectedSequence);
        end
        
        %% testDoControlSequenceComputationZeroStateWithConstraints
        function testDoControlSequenceComputationZeroStateWithConstraints(this)
            % perform a sanity check: given state is origin,
            % and weightings for constraints are zero
            % so computed control sequence (length 1) should be also the zero vector
            % due to the underlying linear control law, independent of the
            % previous mode
            expectedSequence = zeros(this.dimU, 1);
                       
            controllerWithConstraints = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                zeros(this.dimX, 2, this.numConstraints), ...
                zeros(this.dimU, 1, this.numConstraints), this.constraintBounds);
            
            % first mode
            actualSequence = controllerWithConstraints.computeControlSequence(this.zeroStateDistribution, 1, ...
                    this.horizonLength);
            
            this.verifyEqual(actualSequence, expectedSequence);
            
            controllerWithConstraints = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                zeros(this.dimX, 2, this.numConstraints), ...
                zeros(this.dimU, 1, this.numConstraints), this.constraintBounds);
            
            %second mode
            actualSequence = controllerWithConstraints.computeControlSequence(this.zeroStateDistribution, 2, ...
                    this.horizonLength);
            
            this.verifyEqual(actualSequence, expectedSequence);
        end
%%
%%
        %% testDoControlSequenceComputationWithoutConstraints
        function testDoControlSequenceComputationWithoutConstraints(this)
            expectedInputs = this.computeExpectedInputsWithoutConstraints();
            
            % check both modes
            % first mode: previous input arrived at plant
            actualInput = this.controllerUnderTest.computeControlSequence(this.stateDistribution, 1, ...
                    this.horizonLength);
            this.verifyEqual(actualInput, expectedInputs);
            
            % second mode: previous input did not arrive at plant
            actualInput = this.controllerUnderTest.computeControlSequence(this.stateDistribution, 2, ...
                    this.horizonLength);
            this.verifyEqual(actualInput, expectedInputs);
        end
        
        %% testDoControlSequenceComputationWithConstraints
        function testDoControlSequenceComputationWithConstraints(this)
            expectedInputs = this.computeExpectedInputsWithConstraints();
            
            controllerWithConstraints = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, this.inputConstraintWeightings, this.constraintBounds);
            
            % check both modes
            % first mode: previous input arrived at plant
            actualInput = controllerWithConstraints.computeControlSequence(this.stateDistribution, 1, ...
                    this.horizonLength);
            this.verifyEqual(actualInput, expectedInputs, 'AbsTol', 1e-5);
            
            controllerWithConstraints = FiniteHorizonController(this.A, this.B, this.Q, this.R, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, ...
                this.stateConstraintWeightings, this.inputConstraintWeightings, this.constraintBounds);
            
            % second mode: previous input dot not arrive at plant
            actualInput = controllerWithConstraints.computeControlSequence(this.stateDistribution, 2, ...
                    this.horizonLength);
            this.verifyEqual(actualInput, expectedInputs, 'AbsTol', 1e-5);
        end
    end
end

