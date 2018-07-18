classdef FiniteHorizonTrackingControllerTest < BaseFiniteHorizonControllerTest
    % Test cases for FiniteHorizonTrackingController.
    
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
    
    properties
        % performance matrix is simply the identity, and reference is
        % matrix of 1
        Z;
        dimRef;
        refTrajectory;
    end
    
    methods (Access = private)
        function feedback = computeFeedback(this)
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
            
            feedback = -firstPart * secondPart * this.state;
    
        end
        
        function feedforward = computeFeedforward(this)
            % in this setup, the current input must arrive without delay,
            % or plant evolves open-loop
            % hence, in this setup, G and F are not existent
            % likewise, H
            
            % first mode
            R_1 = this.R;
            B_1 = this.B;
                       
            % second mode
            R_2 = zeros(this.dimU);
            B_2 = zeros(this.dimX, this.dimU);
            
            sigma = this.Z * this.Q * this.refTrajectory(:, end);
            K = this.Z' * this.Q * this.Z;
            
            firstPart = -pinv(this.transitionMatrix(1,1) * (R_1 + B_1' * K * B_1) ... 
                + this.transitionMatrix(1,2) * (R_2 + B_2' * K * B_2));
            secondPart = this.transitionMatrix(1,1) * B_1' * sigma + this.transitionMatrix(1,2) * B_2' * sigma;
            
            % the same for both modes
            feedforward = firstPart * secondPart;
        end
    end
    
    methods (Access = protected)
        %% initControllerUnderTest
        function controller = initControllerUnderTest(this)
            controller = FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory);
        end
        
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            this.initAdditionalProperties@BaseFiniteHorizonControllerTest();
            this.dimRef = this.dimX;
            this.Z = eye(this.dimRef);
            this.refTrajectory = ones(this.dimX, this.horizonLength + 1);
        end
        
        %% computeExpectedCosts
        function expectedCosts = computeExpectedCosts(this, states, inputs)
            expectedCosts = (states(:, 1) - this.refTrajectory(:, 1))' * this.Q * (states(:, 1) - this.refTrajectory(:, 1)) ...
                + (states(:, 2) - this.refTrajectory(:, 2))' * this.Q * (states(:, 2) - this.refTrajectory(:, 2)) ...
                + inputs' * this.R * inputs;
        end
        
        %% computeExpectedStageCosts
        function expectedStageCosts = computeExpectedStageCosts(this, state, input, timestep)
             expectedStageCosts =  (state - this.refTrajectory(:, timestep))' * this.Q ...
                 * (state - this.refTrajectory(:, timestep)) + input' * this.R * input;
        end
    end
    
    methods (Test)
        %% testFiniteHorizonTrackingControllerInvalidSystemMatrix
        function testFiniteHorizonTrackingControllerInvalidSystemMatrix(this)
             expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
             
             invalidSysMatrix = eye(this.dimX, this.dimX + 1); % not square
             this.verifyError(@() FiniteHorizonTrackingController(invalidSysMatrix, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
             
             invalidSysMatrix = eye(this.dimX, this.dimX); % square but not finite
             invalidSysMatrix(1, end) = inf;
             this.verifyError(@() FiniteHorizonTrackingController(invalidSysMatrix, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
        end
        
        %% testFiniteHorizonTrackingControllerInvalidInputMatrix
        function testFiniteHorizonTrackingControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidInputMatrix = eye(this.dimX +1, this.dimU); % invalid dims
            this.verifyError(@() FiniteHorizonTrackingController(this.A, invalidInputMatrix, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            invalidInputMatrix = eye(this.dimX, this.dimU); % correct dims, but not finite
            invalidInputMatrix(1, end) = nan;
            this.verifyError(@() FiniteHorizonTrackingController(this.A, invalidInputMatrix, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
        end
        
        %% testFiniteHorizonTrackingControllerInvalidHorizonLength
        function testFiniteHorizonTrackingControllerInvalidHorizonLength(this)
            expectedErrId = 'Validator:ValidateHorizonLength:InvalidHorizonLength';
            
            invalidHorizonLength = eye(3); % not a scalar
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength, this.refTrajectory), expectedErrId);
            
            invalidHorizonLength = 0; % not a positive scalar
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength, this.refTrajectory), expectedErrId);
            
            invalidHorizonLength = 1.5; % not an integer
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength, this.refTrajectory), expectedErrId);
            
            invalidHorizonLength = inf; % not finite
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, invalidHorizonLength, this.refTrajectory), expectedErrId);
        end
                                        
        %% testFiniteHorizonTrackingControllerInvalidCostMatrices
        function testFiniteHorizonTrackingControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, invalidQ, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, invalidQ, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, invalidQ, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = eye(this.dimX); % Q is not symmetric
            invalidQ(1, end) = 1;
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, invalidQ, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, invalidQ, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, invalidR, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, invalidR, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, invalidR, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
        end
        
        %% testFiniteHorizonTrackingControllerInvalidDelayProbs
        function testFiniteHorizonTrackingControllerInvalidDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                invalidDelayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                invalidDelayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                invalidDelayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory), expectedErrId);
        end
        
        %% testFiniteHorizonTrackingController
        function testFiniteHorizonTrackingController(this)
            controller =  FiniteHorizonTrackingController(this.A, this.B, this.Q, this.R, this.Z, ...
                this.delayProbs, this.sequenceLength, this.horizonLength, this.refTrajectory);
            
            this.verifyEqual(controller.horizonLength, this.horizonLength);
            this.verifyEqual(controller.refTrajectory, this.refTrajectory);
        end
                    
        %% testGetDeviationFromRefForStateInvalidTimestep
        function testGetDeviationFromRefForStateInvalidTimestep(this)
            trueState = zeros(this.dimX, 1);
            expectedErrId = 'FiniteHorizonTrackingController:GetDeviationFromRefForState:InvalidTimestep';
            
            invalidTimestep = this; % not a scalar
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(trueState, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = -1; % scalar, but negative
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(trueState, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = this.horizonLength + 1; % scalar, but out of bounds
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(trueState, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = 0.5; % scalar, but not an integer
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(trueState, invalidTimestep), ...
                expectedErrId);
        end

        %% testGetDeviationFromRefForState
        function testGetDeviationFromRefForState(this)
            timestep = 1;
            trueState = 42 * ones(this.dimX, 1);
            
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(trueState, timestep);
            expectedDeviation = this.Z * trueState - this.refTrajectory(:, timestep);
            
            this.verifyEqual(actualDeviation, expectedDeviation);
        end
        
%%
%%
        %% testDoControlSequenceComputationZeroState
        function testDoControlSequenceComputationZeroState(this)
            % perform a sanity check: given state is origin so computed
            % control sequence (length 1) should only consist of a
            % feedforward term due to the underlying linear feedback control law, independent of the
            % previous mode
            expectedSequence = this.computeFeedforward();

            % first mode
            actualSequence = this.controllerUnderTest.computeControlSequence(this.zeroStateDistribution, 1, ...
                    this.horizonLength);
            
            this.verifyEqual(actualSequence, expectedSequence);
            
            %second mode
            actualSequence = this.controllerUnderTest.computeControlSequence(this.zeroStateDistribution, 2, ...
                    this.horizonLength);
            
            this.verifyEqual(actualSequence, expectedSequence);
        end
        
        %% testDoControlSequenceComputation
        function testDoControlSequenceComputation(this)
            % now the state is not the origin
            expectedSequence = this.computeFeedback() + this.computeFeedforward();
            
            % check both modes
            % first mode: previous input arrived at plant
            actualSequence = this.controllerUnderTest.computeControlSequence(this.stateDistribution, 1, ...
                    this.horizonLength);
            this.verifyEqual(actualSequence, expectedSequence);
            
            % second mode: previous input dot not arrive at plant
            actualSequence = this.controllerUnderTest.computeControlSequence(this.stateDistribution, 2, ...
                    this.horizonLength);
            this.verifyEqual(actualSequence, expectedSequence);
        end
    end
end

