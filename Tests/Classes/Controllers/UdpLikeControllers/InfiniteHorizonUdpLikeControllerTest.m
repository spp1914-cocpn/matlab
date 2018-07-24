classdef InfiniteHorizonUdpLikeControllerTest < matlab.unittest.TestCase
    % Test cases for InfiniteHorizonUdpLikeController.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Constant)
        absTol = 1e-8;
    end
    
     properties (Access = private)
        A;
        B;
        Q;
        R;
        C;
        sequenceLength;
        dimX;
        dimU;
        dimY;
        
        W;
        V;
        v_mean;
        
        caDelayProbs;
        scDelayProbs;
        maxMeasDelay;
        
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.dimX = 3;
            this.dimU = 2;
            this.dimY = 2;
            
            this.sequenceLength = 1; % so we have two modes
            this.caDelayProbs = ones(1, 5) / 5;
            this.scDelayProbs = ones(1, 6) / 6;
                                   
            % from the delay probs and the 2 modes we get the following
            % transition matrix
            %this.transitionMatrix = [1/5 4/5; 1/5 4/5];
                
            this.A = eye(this.dimX);
            this.B = ones(this.dimX, this.dimU);
            this.C = 2 * ones(this.dimY, this.dimX);
            
            this.Q = eye(this.dimX);
            this.R = eye(this.dimU);
              
            this.maxMeasDelay = 1;
            
            this.W = eye(this.dimX) * 0.01;
            this.V = eye(this.dimY) * 0.001;
            this.v_mean = repmat(0.5, 1, this.dimY);
            
            this.controllerUnderTest = InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.v_mean);
        end
    end
    
    methods (Test)
        %% testInfiniteHorizonUdpLikeControllerInvalidSystemMatrix
        function testInfiniteHorizonUdpLikeControllerInvalidSystemMatrix(this)
             expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
             
             invalidSysMatrix = eye(this.dimX, this.dimX + 1); % not square
             this.verifyError(@() InfiniteHorizonUdpLikeController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
             
             invalidSysMatrix = eye(this.dimX, this.dimX); % square but not finite
             invalidSysMatrix(1, end) = inf;
             this.verifyError(@() InfiniteHorizonUdpLikeController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidInputMatrix
        function testInfiniteHorizonUdpLikeControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidInputMatrix = eye(this.dimX +1, this.dimU); % invalid dims
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
             
            invalidInputMatrix = eye(this.dimX, this.dimU); % correct dims, but not finite
            invalidInputMatrix(1, end) = nan;
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
        end
        
        
        %% testInfiniteHorizonUdpLikeControllerInvalidCostMatrices
        function testInfiniteHorizonUdpLikeControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = eye(this.dimX); % Q is not symmetric
            invalidQ(1, end) = 1;
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidMeasurementMatrix
        function testInfiniteHorizonUdpLikeControllerInvalidMeasurementMatrix(this)
            expectedErrId = 'Validator:ValidateMeasurementMatrix:InvalidMeasMatrix';
            
            invalidMeasMatrix = eye(this.dimY, this.dimX + 1); % invalid dims
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
             
            invalidMeasMatrix = eye(this.dimY, this.dimX); % correct dims, but not finite
            invalidMeasMatrix(1, end) = nan;
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidMaxMeasDelay
        function testInfiniteHorizonUdpLikeControllerInvalidMaxMeasDelay(this)
            expectedErrId = 'InfiniteHorizonUdpLikeController:InvalidMaxMeasDelay';
            
            invalidMaxMeasDelay = [1 2]; % not a scalar
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidMaxMeasDelay = -1; % negative scalar
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidMaxMeasDelay = 1.5; % not an integer
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidMaxMeasDelay = inf; % not finite
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V), ...
                expectedErrId);
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidSysNoiseCovariance
        function testInfiniteHorizonUdpLikeControllerInvalidSysNoiseCovariance(this)
            expectedErrId = 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCovDim';
            
            invalidW = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, invalidW, this.V), ...
                expectedErrId);
                       
            invalidW = ones(this.dimU); % W is not pd
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, invalidW, this.V), ...
                expectedErrId); 
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidMeasNoiseCovariance
        function testInfiniteHorizonUdpLikeControllerInvalidMeasNoiseCovariance(this)
            expectedErrId = 'Validator:ValidateMeasNoiseCovarianceMatrix:InvalidCovDim';
            
            invalidV = eye(this.dimY + 1, this.dimY); % not square
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, invalidV), ...
                expectedErrId);
                       
            invalidV = ones(this.dimY); % V is not pd
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, invalidV), ...
                expectedErrId);
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidMeasNoiseMean
        function testInfiniteHorizonUdpLikeControllerInvalidMeasNoiseMean(this)
            expectedErrId = 'InfiniteHorizonUdpLikeController:InvalidMeasNoiseMean';
            
            invalidNoiseMean = this; % not a vector
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidNoiseMean), ...
                expectedErrId);
            
            invalidNoiseMean = ones(this.dimY + 2, 1); % invalid dims
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidNoiseMean), ...
                expectedErrId);       
            
            tmp = {1};
            invalidNoiseMean = repmat(tmp, 1, this.dimY); % not a numeric vector
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidNoiseMean), ...
                expectedErrId); 
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidCaDelayProbs
        function testInfiniteHorizonUdpLikeControllerInvalidCaDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
        end
        
        %% testInfiniteHorizonUdpLikeControllerInvalidScDelayProbs
        function testInfiniteHorizonUdpLikeControllerInvalidScDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V), ...
                expectedErrId);
        end
%%
%%
        %% testInfiniteHorizonUdpLikeController
        function testInfiniteHorizonUdpLikeController(this)
            controller = InfiniteHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V);
            
            expectedControllerPlantState = zeros(this.dimX, 1);
            this.verifyEqual(controller.getControllerPlantState(), expectedControllerPlantState);
            this.verifyEqual(controller.maxMeasurementDelay, this.maxMeasDelay);
            this.verifyEqual(controller.sequenceLength, this.sequenceLength);
        end
%%
%%
        %% testReset
        function testReset(this)
            this.controllerUnderTest.computeControlSequence();
            zeroState = zeros(this.dimX, 1);
            % assert that the state has changed
            this.assertNotEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
            
            this.controllerUnderTest.reset();
            
            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
        end
%%
%%
        %% testGetLastComputationMeasurementDataNoMeasurements
        function testGetLastComputationMeasurementDataNoMeasurements(this)
            % compute a sequence without using measurements
            this.controllerUnderTest.computeControlSequence();
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testGetLastComputationMeasurementDataNoMeasDiscarded
        function testGetLastComputationMeasurementDataNoMeasDiscarded(this)
            measurements = ones(this.dimY, 2); % two meas, column-wise arranged
            delays = [0 1];
            
            this.controllerUnderTest.computeControlSequence(measurements, delays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testGetLastComputationMeasurementData
        function testGetLastComputationMeasurementData(this)
            measurements = ones(this.dimY, 4); % four meas, column-wise arranged
            delays = [0 2 1 3]; % two meas are too old
            
            warning off;
            this.controllerUnderTest.computeControlSequence(measurements, delays);
            warning on;
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 2);
        end
%%
%%
        %% testComputeControlSequenceNoMeasurements
        function testComputeControlSequenceNoMeasurements(this)
            zeroState = zeros(this.dimX, 1);
            % the initial state is the origin, and we have a linear control
            % law
            expectedInputSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence();
            
            this.verifyEqual(actualInputSequence, expectedInputSequence);
            % the controller state should have changed
            this.verifyNotEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
        end
        
        %% testComputeControlSequenceInvalidMeasurements
        function testComputeControlSequenceInvalidMeasurements(this)
            expectedErrId = 'InfiniteHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeas';
            
            invalidMeasurements = ones(this.dimY + 1, 2); % two meas, column-wise arranged, but of wrong dimension
            delays = [0 1];            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(invalidMeasurements, delays), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceInvalidMeasDelays
        function testComputeControlSequenceInvalidMeasDelays(this)
            expectedErrId = 'InfiniteHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeasDelay';
            
            measurements = ones(this.dimY, 3); % three meas, column-wise arranged
            invalidDelays = [0 -1 1]; % negative delay            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 0.5 1]; % fractional delay            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 1]; % not enough entries           
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceInvalidMeasDelaysNotUnique
        function testComputeControlSequenceInvalidMeasDelaysNotUnique(this)
            expectedErrId = 'InfiniteHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsNotUnique';
            
            measurements = ones(this.dimY, 3); % three meas, column-wise arranged
            invalidDelays = [0 1 0]; % only one measurement per delay allowed          
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays), ...
                expectedErrId);            
     
        end
        
        %% testComputeControlSequenceIssueWarningDiscardMeasurement
        function testComputeControlSequenceIssueWarningDiscardMeasurement(this)
            expectedWarnId = 'InfiniteHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsDelayTooLarge';
            
            zeroState = zeros(this.dimX, 1);
            % the initial state is the origin, and we have a linear control
            % law
            expectedInputSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            measurements = ones(this.dimY, 4); % four meas, column-wise arranged
            delays = [0 2 1 3]; % two meas are too old     
            actualInputSequence = this.verifyWarning(@() this.controllerUnderTest.computeControlSequence(measurements, delays), ...
                expectedWarnId);            
     
            this.verifyEqual(actualInputSequence, expectedInputSequence);
            % the controller state should have changed
            this.verifyNotEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
        end
%%
%%
        %% testDoStageCostsComputation
        function testDoStageCostsComputation(this)
            state = ones(this.dimX, 1);
            input = 1.5 * ones(this.dimU, 1);
            timestep = 1;
            
            expectedStageCosts = state' * this.Q * state + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', InfiniteHorizonUdpLikeControllerTest.absTol);
        end
%%
%%
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            horizonLength = 10;
            inputs = ones(this.dimU, horizonLength);
            expectedErrId = 'InfiniteHorizonUdpLikeController:DoCostsComputation';
            
            invalidStates = ones(this.dimX, horizonLength + 2); %  trajectory too long
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStates, inputs), ...
                expectedErrId);
            
            invalidStates = ones(this.dimX, horizonLength); %  trajectory too short
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStates, inputs), ...
                expectedErrId);
        end
        
        %% testDoCostsComputation
        function testDoCostsComputation(this)
            horizonLength = 100;  
            inputs = ones(this.dimU, horizonLength);
            states = 2.25 * ones(this.dimX, horizonLength + 1);
            
            expectedCosts = states(:, end)' * this.Q * states(:, end);
            for j=1:horizonLength
                expectedCosts = expectedCosts + states(:, j)' * this.Q * states(:, j) ...
                    + inputs(:, j)' * this.R * inputs(:, j);
            end
            expectedCosts = expectedCosts / horizonLength;
            
            actualCosts = this.controllerUnderTest.computeCosts(states, inputs);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', NominalPredictiveControllerTest.absTol);
        end
    end
end

