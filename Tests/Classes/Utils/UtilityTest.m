classdef UtilityTest < matlab.unittest.TestCase
    % Test cases for Utility.
    
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
    
    properties (Constant, Access = private)
        absTol = 1e-8;
    end
    
    properties (Access = private)
        dimX;
        dimU;
        dimZ;
        
        A;
        B;
        Z;
        
        Q; 
        R;
        Q_z;
    end
    
    methods (Access = private)
        %% verifyEqualWithAbsTol
        function verifyEqualWithAbsTol(this, actual, expected)
            this.verifyEqual(actual, expected, 'AbsTol', UtilityTest.absTol);
        end
        
        %% verifyAugmentedPlantMatricesOneElementSequence
        function verifyAugmentedPlantMatricesOneElementSequence(this, actualF, actualG, actualAugA, actualAugB)
             % we expected F, G, to be empty
            this.verifyEmpty(actualF);
            this.verifyEmpty(actualG);
                           
            % augA should be simply be two A matrices
            this.verifySize(actualAugA, [this.dimX, this.dimX, 2]);
            this.verifyEqual(actualAugA(:, :, 1), this.A);
            this.verifyEqual(actualAugA(:, :, 2), this.A);
            
            % augB should be simply consist of B and the zero matrix
            this.verifySize(actualAugB, [this.dimX, this.dimU, 2]);
            this.verifyEqual(actualAugB(:, :, 1), this.B);
            this.verifyEqual(actualAugB(:, :, 2), zeros(this.dimX, this.dimU));
        end
        
        %% verifyAugmentedPlantMatrices
        function verifyAugmentedPlantMatrices(this, actualF, actualG, actualAugA, actualAugB, controlSequenceLength, ...
                actualH, actualJ)
            numModes = controlSequenceLength + 1;
            
            dimAugState = this.dimX + this.dimU * controlSequenceLength * (controlSequenceLength - 1) / 2;
            
            expectedF = zeros(18); % dimU * sum(1:controlSequenceLength)
            expectedF(10:15,4:9) = eye(2 * this.dimU);
            expectedF(16:18, 13:15) = eye(this.dimU);
            
            expectedG = zeros(18,12);
            expectedG(1:9, 4:end) = eye(3 * this.dimU);
            
            expectedH = zeros(this.dimU, 18, numModes);
            expectedH(:, 1:this.dimU, 2) = eye(this.dimU);
            expectedH(:, 10:12, 3) = eye(this.dimU);
            expectedH(:, end-2:end, 4) = eye(this.dimU);
            
            expectedJ = zeros(this.dimU, 12, numModes);
            expectedJ(:, 1:this.dimU, 1) = eye(this.dimU);
                        
            expectedAugA = zeros(dimAugState, dimAugState, numModes);
            expectedAugB = zeros(dimAugState, controlSequenceLength * this.dimU, numModes);
            for j=1:numModes
                expectedAugA(:, :, j) = [   this.A, this.B * expectedH(:, :, j);
                                            zeros(dimAugState - this.dimX,this.dimX), expectedF
                                            ];
                expectedAugB(:, :, j) = [this.B * expectedJ(:, :, j); expectedG];
            end
            
            this.verifyEqual(actualF, expectedF);
            this.verifyEqual(actualG, expectedG);
            if nargin == 8
                this.verifyEqual(actualH, expectedH);
                this.verifyEqual(actualJ, expectedJ);
            end
            this.verifyEqual(actualAugA, expectedAugA);
            this.verifyEqual(actualAugB, expectedAugB);
        end
        
        %% verifyAugmentedCostMatricesOneElementSequence
        function verifyAugmentedCostMatricesOneElementSequence(this, actualAugQ, actualAugR, useZ)
            
            this.verifySize(actualAugQ, [this.dimX this.dimX 2]);
            if useZ
                this.verifyEqual(actualAugQ(:, :, 1), this.Z' * this.Q_z * this.Z);
                this.verifyEqual(actualAugQ(:, :, 2), this.Z' * this.Q_z * this.Z);
            else
                this.verifyEqual(actualAugQ(:, :, 1), this.Q);
                this.verifyEqual(actualAugQ(:, :, 2), this.Q);
            end
            this.verifySize(actualAugR, [this.dimU this.dimU 2]);
            this.verifyEqual(actualAugR(:, :, 1), this.R);
            this.verifyEqual(actualAugR(:, :, 2), zeros(this.dimU));
        end
        
        %% verifyAugmentedCostMatrices
        function verifyAugmentedCostMatrices(this, actualAugQ, actualAugR, useZ, controlSequenceLength)
            numModes = controlSequenceLength + 1;
            dimAugState = this.dimX + this.dimU * controlSequenceLength * (controlSequenceLength - 1) / 2;
            
            H = zeros(this.dimU, 18, numModes);
            H(:, 1:this.dimU, 2) = eye(this.dimU);
            H(:, 10:12, 3) = eye(this.dimU);
            H(:, end-2:end, 4) = eye(this.dimU);
            J = zeros(this.dimU, 12, numModes);
            J(:, 1:this.dimU, 1) = eye(this.dimU);
                       
            expectedAugQ = zeros(dimAugState, dimAugState, numModes);
            expectedAugR = zeros(12, 12, numModes);
            for j=1:numModes
                if useZ
                    expectedAugQ(:, :, j) = blkdiag(this.Z' * this.Q_z * this.Z, H(:, :, j)' * this.R * H(:, :, j));
                else
                    expectedAugQ(:, :, j) = blkdiag(this.Q, H(:, :, j)' * this.R * H(:, :, j));
                end
                expectedAugR(:, :, j) = J(:, :, j)' * this.R * J(:, :, j);
            end
            
            this.verifyEqual(actualAugQ, expectedAugQ);
            this.verifyEqual(actualAugR, expectedAugR);
        end
        
        %% verifyEqualAugmentedLinearIntegralConstraints
        function verifyEqualAugmentedLinearIntegralConstraints(this, ...
                augmentedStateWeightings, augmentedInputWeightings, controlSequenceLength, stateWeights, inputWeights)
            
            numModes = controlSequenceLength + 1;
            horizonLength = size(stateWeights, 2);
            numConstraints = size(stateWeights, 3);
            
            dimAugmentedState = this.dimX + this.dimU * controlSequenceLength * (controlSequenceLength - 1) / 2;
            dimSequence = this.dimU * controlSequenceLength;
            H = zeros(this.dimU, 18, numModes);
            H(:, 1:this.dimU, 2) = eye(this.dimU);
            H(:, 10:12, 3) = eye(this.dimU);
            H(:, end-2:end, 4) = eye(this.dimU);
            
            J = zeros(this.dimU, 12, numModes);
            J(:, 1:this.dimU, 1) = eye(this.dimU);
            
            expectedAugmentedStateWeights = zeros(dimAugmentedState, horizonLength, numConstraints, numModes);
            expectedAugmentedInputWeights = zeros(dimSequence, horizonLength, numConstraints, numModes);
            for c = 1:numConstraints
                for j=1:numModes
                    expectedAugmentedStateWeights(:, :, c, j) = [stateWeights(:, :, c); H(:, :, j)' * inputWeights(:, :, c)];
                    expectedAugmentedInputWeights(:, :, c, j) = J(:, :, j)' * inputWeights(:, :, c);
                end
            end
            this.verifyEqual(augmentedStateWeightings, expectedAugmentedStateWeights);
            this.verifyEqual(augmentedInputWeightings, expectedAugmentedInputWeights);
        end
        
        %% computeCostsAndTrajectories
        function [stateTrajectory, inputTrajectory, expectedCosts] = computeCostsAndTrajectories(this)
            
            % use horizon length of 11
            horizonLength = 11;
            stateTrajectory = repmat(gallery('moler', this.dimX), 1, 6);
            inputTrajectory = repmat(gallery('minij', this.dimU), 1, 4);
            inputTrajectory = inputTrajectory(:, 1:end-1);
            
            expectedCosts = stateTrajectory(:, end)' * this.Q * stateTrajectory(:, end);
            for j = 1:horizonLength
                expectedCosts = expectedCosts + stateTrajectory(:, j)' * this.Q * stateTrajectory(:, j);
                expectedCosts = expectedCosts + inputTrajectory(:, j)' * this.R * inputTrajectory(:, j);
            end
        end
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.dimX = 2;
            this.dimU = 3;
            this.dimZ = 4;
            this.A = [1 1; 0 1];
            this.B = [0 1 0; 1 0 1];
            
            this.Q = gallery('moler', this.dimX);
            this.R = gallery('moler', this.dimU);
            this.Z = ones(this.dimZ, this.dimX);
            this.Q_z = gallery('minij', this.dimZ);
        end
    end
    
    methods (Test)
%%
%%      
        %% testComputeStageCosts
        function testComputeStageCosts(this)
            state = [1 3]';
            input = [-2 -3 5]';
            
            actualStageCosts = Utility.computeStageCosts(state, input, this.Q, this.R);
            expectedStateCosts = state' * this.Q * state + input' * this.R * input;
            
            this.verifyEqualWithAbsTol(actualStageCosts, expectedStateCosts);
        end

        %% testComputeLQGCosts
        function testComputeLQGCosts(this)
            [stateTrajectory, inputTrajectory, expectedCosts] = this.computeCostsAndTrajectories(); 

            horizonLength = size(inputTrajectory, 2);
            actualCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, inputTrajectory, this.Q, this.R);
            this.verifyEqualWithAbsTol(actualCosts, expectedCosts);
        end
        
        %% testComputeLGQCostsMultiDim
        function testComputeLGQCostsMultiDim(this)
            % now check if we pass 3d arrays, as if multiple state and
            % input trajectories are to processed
            
            [stateTrajectory, inputTrajectory, expectedCosts] = this.computeCostsAndTrajectories();
            
            states = repmat(stateTrajectory, 1, 1, 3);
            inputs = repmat(inputTrajectory, 1, 1, 3);
            horizonLength = size(inputTrajectory, 2);
            % now the actual costs
            actualCosts = Utility.computeLQGCosts(horizonLength, states, inputs, this.Q, this.R);
            
            this.verifySize(actualCosts, [1 3]);
            this.verifyEqualWithAbsTol(actualCosts, repmat(expectedCosts, 1, 3));
        end
        
        %% testComputeLQGCostsInvalidDims
        function testComputeLQGCostsInvalidDims(this)
            [stateTrajectory, inputTrajectory, ~] = this.computeCostsAndTrajectories();
            
            horizonLength = size(inputTrajectory, 2);
            % first, the state trajectory is too short
            actualCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory(:, 1), inputTrajectory, this.Q, this.R);
            this.verifyEmpty(actualCosts);
            
            % now, horizon is too short
            horizonLength = 2;
            actualCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, inputTrajectory, this.Q, this.R);
            this.verifyEmpty(actualCosts);
            
            % finally, dimensions do not match: states or inputs are 3d array
            horizonLength = size(inputTrajectory, 2);
            actualCosts = Utility.computeLQGCosts(horizonLength, repmat(stateTrajectory, 1, 1, 3), inputTrajectory, this.Q, this.R);
            this.verifyEmpty(actualCosts);

            actualCosts = Utility.computeLQGCosts(horizonLength, stateTrajectory, repmat(inputTrajectory, 1, 1, 3), this.Q, this.R);
            this.verifyEmpty(actualCosts);
        end
%%
%%        
        %% testCreateAugmentedPlantModelOneElementSequence
        function testCreateAugmentedPlantModelOneElementSequence(this)
                        
            controlSequenceLength = 1; % no sequence, just the input for the current time step
            
            [actualF, actualG, actualH, actualJ, actualAugA, actualAugB] = ...
                Utility.createAugmentedPlantModel(controlSequenceLength, this.A, this.B);
            
            this.verifyAugmentedPlantMatricesOneElementSequence(actualF, actualG, actualAugA, actualAugB);
            
            % H should be empty
            this.verifyEmpty(actualH);
            % J should be two matrices of dimension dimU, first one I and
            % second one 0
            this.verifySize(actualJ, [this.dimU, this.dimU, 2]);
            this.verifyEqual(actualJ(:, :, 1), eye(this.dimU));
            this.verifyEqual(actualJ(:, :, 2), zeros(this.dimU));
        end
        
        %% testCreateAugmentedPlantModel
        function testCreateAugmentedPlantModel(this)
            
            controlSequenceLength = 4;
                        
            [actualF, actualG, actualH, actualJ, actualAugA, actualAugB] = ...
                Utility.createAugmentedPlantModel(controlSequenceLength, this.A, this.B);
  
            this.verifyAugmentedPlantMatrices(actualF, actualG, actualAugA, actualAugB, ...
                controlSequenceLength, actualH, actualJ);
            
        end
%%
%%        
        %% testCreateAugmentedCostModelOneElementSequence
        function testCreateAugmentedCostModelOneElementSequence(this)

            controlSequenceLength = 1; % no sequence, just the input for the current time step
            
            % first, do not use the Z matrix
            [actualAugQ, actualAugR] = Utility.createAugmentedCostModel(controlSequenceLength, this.Q, this.R);
            
            this.verifyAugmentedCostMatricesOneElementSequence(actualAugQ, actualAugR, false);
             
            % now utilize the Z matrix
            [actualAugQ, actualAugR] = Utility.createAugmentedCostModel(controlSequenceLength, this.Q_z, this.R, this.Z);
   
            this.verifyAugmentedCostMatricesOneElementSequence(actualAugQ, actualAugR, true);
        end
        
        %% testCreateAugmentedCostModel
        function testCreateAugmentedCostModel(this)
            
            controlSequenceLength = 4;
            
            % first, do not use the Z matrix
            [actualAugQ, actualAugR] = Utility.createAugmentedCostModel(controlSequenceLength, this.Q, this.R);
   
            this.verifyAugmentedCostMatrices(actualAugQ, actualAugR, false, controlSequenceLength);
            
            [actualAugQ, actualAugR] = Utility.createAugmentedCostModel(controlSequenceLength, this.Q_z, this.R, this.Z);
 
            this.verifyAugmentedCostMatrices(actualAugQ, actualAugR, true, controlSequenceLength);
        end
%%
%%        
        %% testPerformModelAugmentationOneElementSequence
        function testPerformModelAugmentationOneElementSequence(this)
            controlSequenceLength = 1; % no sequence, just the input for the current time step
            
            % first, do not use the Z matrix
            [dimAugmentedState, F, G, augA, augB, augQ, augR] = ...
                Utility.performModelAugmentation(controlSequenceLength, this.dimX, this.dimU, ...
                this.A, this.B, this.Q, this.R);
            
            this.verifyEqual(dimAugmentedState, this.dimX);
            this.verifyAugmentedPlantMatricesOneElementSequence(F, G, augA, augB);
            
            this.verifyAugmentedCostMatricesOneElementSequence(augQ, augR, false);
            
            % now also use the Z matrix
             [dimAugmentedState, F, G, augA, augB, augQ, augR] = ...
                Utility.performModelAugmentation(controlSequenceLength, this.dimX, this.dimU, ...
                this.A, this.B, this.Q_z, this.R, this.Z);
            
            this.verifyEqual(dimAugmentedState, this.dimX);
            this.verifyAugmentedPlantMatricesOneElementSequence(F, G, augA, augB);
            this.verifyAugmentedCostMatricesOneElementSequence(augQ, augR, true);
            
        end

        %% testPerformModelAugmentation
        function testPerformModelAugmentation(this)
            controlSequenceLength = 4;
            expectedDimAugmentedState = this.dimX + this.dimU * controlSequenceLength * (controlSequenceLength - 1) / 2;
            
             % first, do not use the Z matrix
            [dimAugmentedState, F, G, augA, augB, augQ, augR] = ...
                Utility.performModelAugmentation(controlSequenceLength, this.dimX, this.dimU, ...
                this.A, this.B, this.Q, this.R);
            
            this.verifyEqual(dimAugmentedState, expectedDimAugmentedState);
            this.verifyAugmentedPlantMatrices(F, G, augA, augB, controlSequenceLength);
            this.verifyAugmentedCostMatrices(augQ, augR, false, controlSequenceLength);
            
            % now also use the Z matrix
            [dimAugmentedState, F, G, augA, augB, augQ, augR] = ...
                Utility.performModelAugmentation(controlSequenceLength, this.dimX, this.dimU, ...
                this.A, this.B, this.Q_z, this.R, this.Z);
            
            this.verifyEqual(dimAugmentedState, expectedDimAugmentedState);
            this.verifyAugmentedPlantMatrices(F, G, augA, augB, controlSequenceLength);
            this.verifyAugmentedCostMatrices(augQ, augR, true, controlSequenceLength);
        end
%%
%%        
        %% testCreateAugmentedLinearIntegralConstraintsOneElementSequence
        function testCreateAugmentedLinearIntegralConstraintsOneElementSequence(this)
            controlSequenceLength = 1; % no sequence, just the input for the current time step
            horizonLength = 100;
            numConstraints = 5;
            stateWeights = ones(this.dimX, horizonLength, numConstraints);
            inputWeights = 2 * ones(this.dimU, horizonLength, numConstraints);
            
            [augmentedStateWeightings, augmentedInputWeightings] ...
                = Utility.createAugmentedLinearIntegralConstraints(controlSequenceLength, stateWeights, inputWeights);
            
            this.verifySize(augmentedStateWeightings, [this.dimX, horizonLength, numConstraints, 2]);
            this.verifyEqual(augmentedStateWeightings(:, :, :, 1), stateWeights);
            this.verifyEqual(augmentedStateWeightings(:, :, :, 2), stateWeights);
            
            this.verifySize(augmentedInputWeightings, [this.dimU, horizonLength, numConstraints, 2]);
            this.verifyEqual(augmentedInputWeightings(:, :, :, 1), inputWeights);
            this.verifyEqual(augmentedInputWeightings(:, :, :, 2), zeros(this.dimU, horizonLength, numConstraints));
        end

        %% testCreateAugmentedLinearIntegralConstraints
        function testCreateAugmentedLinearIntegralConstraints(this)
            controlSequenceLength = 4;
            horizonLength = 100;
            numConstraints = 5;
            
            stateWeights = ones(this.dimX, horizonLength, numConstraints);
            inputWeights = 2 * ones(this.dimU, horizonLength, numConstraints);
            
            [augmentedStateWeightings, augmentedInputWeightings] ...
                = Utility.createAugmentedLinearIntegralConstraints(controlSequenceLength, stateWeights, inputWeights);
            
            this.verifyEqualAugmentedLinearIntegralConstraints(...
                augmentedStateWeightings, augmentedInputWeightings, controlSequenceLength, stateWeights, inputWeights);
  
        end
        
        %% testCreateAugmentedLinearIntegralConstraintsSingleConstraint
        function testCreateAugmentedLinearIntegralConstraintsSingleConstraint(this)
            % same as above, but now test the border case that only one
            % constraint is present
            controlSequenceLength = 4;
            horizonLength = 100;
            numConstraints = 1;
            
            stateWeights = ones(this.dimX, horizonLength, numConstraints);
            inputWeights = 2 * ones(this.dimU, horizonLength, numConstraints);
            
            [augmentedStateWeightings, augmentedInputWeightings] ...
                = Utility.createAugmentedLinearIntegralConstraints(controlSequenceLength, stateWeights, inputWeights);
            
            this.verifyEqualAugmentedLinearIntegralConstraints(...
                augmentedStateWeightings, augmentedInputWeightings, controlSequenceLength, stateWeights, inputWeights);
        end
%%        
%%      
        %% testCalculateDelayTransitionMatrix
        function testCalculateDelayTransitionMatrix(this)
            delayProbs = [0.1 0.3 0.4 0.1 0.1];
            
            % compute transition matrix straightforwardly as described in
            %
            % Jörg Fischer, Achim Hekler, Maxim Dolgov, and Uwe D. Hanebeck,
            % Optimal Sequence-Based LQG Control over TCP-like Networks Subject to Random Transmission Delays and Packet Losses,
            % Proceedings of the 2013 American Control Conference (ACC 2013),
            % Washington D. C., USA, June 2013
            %
            expectedTransitionMatrix = zeros(5);
            
            expectedTransitionMatrix(1,1) = delayProbs(1);
            expectedTransitionMatrix(1,2) = 1 - delayProbs(1);
            expectedTransitionMatrix(2,1) = delayProbs(1);
            expectedTransitionMatrix(2,2) = delayProbs(2);
            expectedTransitionMatrix(2,3) = 1 - sum(delayProbs(1:2));
            expectedTransitionMatrix(3,1) = delayProbs(1);
            expectedTransitionMatrix(3,2) = delayProbs(2);
            expectedTransitionMatrix(3,3) = delayProbs(3);
            expectedTransitionMatrix(3,4) = 1 - sum(delayProbs(1:3));
            expectedTransitionMatrix(4,1) = delayProbs(1);
            expectedTransitionMatrix(4,2) = delayProbs(2);
            expectedTransitionMatrix(4,3) = delayProbs(3);
            expectedTransitionMatrix(4,4) = delayProbs(4);
            expectedTransitionMatrix(4,5) = 1 - sum(delayProbs(1:4));
            expectedTransitionMatrix(5,1) = delayProbs(1);
            expectedTransitionMatrix(5,2) = delayProbs(2);
            expectedTransitionMatrix(5,3) = delayProbs(3);
            expectedTransitionMatrix(5,4) = delayProbs(4);
            expectedTransitionMatrix(5,5) = 1 - sum(delayProbs(1:4));
            
            actualTransitionMatrix = Utility.calculateDelayTransitionMatrix(delayProbs);
            this.verifyEqualWithAbsTol(actualTransitionMatrix, expectedTransitionMatrix);
        end
        
        %% testComputeStationaryDistribution
        function testComputeStationaryDistributionInvalidTransitionMatrix(this)
            expectedErrId = 'Utility:ComputeStationaryDistribution:InvalidTransitionMatrix';
            
            invalidTransitionMatrix = [ 1/2 1/3 1/3;
                                        1/2 1/6 1/3;
                                        1/4 1/12 2/3
                                    ];
                                
            this.verifyError(@() Utility.computeStationaryDistribution(invalidTransitionMatrix), expectedErrId);
        end
        
        %% testComputeStationaryDistribution
        function testComputeStationaryDistribution(this)
            transitionMatrix = [1/3 1/3 1/3;
                                1/2 1/6 1/3;
                                1/4 1/12 2/3
                                ];
            
            expectedStationaryDist = [9/28 5/28 1/2]';
            % compute the stationary distribution using eigenvalue
            % decomposition
            stationaryDist = Utility.computeStationaryDistribution(transitionMatrix);
            this.verifyEqualWithAbsTol(stationaryDist(:), expectedStationaryDist);
            
            % compute the stationary distribution by iterating
            stationaryDist = Utility.computeStationaryDistribution(transitionMatrix, false);
            this.verifyEqualWithAbsTol(stationaryDist(:), expectedStationaryDist);
        end
        
        %% testCalculateRMSE
        function testCalculateRMSE(this)
            % perform a sanity check:
            % exact match between true states and estimates at every time
            % step, assume two runs
            trueStates = repmat(gallery('prolate', 20, 0.1), 1, 1, 2);
            estimates = trueStates;
            
            actualRMSE = Utility.calculateRMSE(trueStates, estimates);
            this.verifySize(actualRMSE, [1 20]);
            this.verifyEqual(actualRMSE, zeros(1, 20));
            
            % now we have difference between estimates and true states
            % 4d state, 200 time steps, 10 runs
            trueStates = ones(this.dimX + 2, 200, 10);
            estimates = trueStates .* 2;
            % differenence vector is [1 1 1 1] -> norm 2 per run and time
            % step -> expected rmse is 10 * 2 / 10 = 2 per time step
            expectedRMSE = ones(1, 200) * 2;
            actualRMSE = Utility.calculateRMSE(trueStates, estimates);
            this.verifySize(actualRMSE, [1 200]);
            this.verifyEqualWithAbsTol(actualRMSE, expectedRMSE);
        end
        
        %% testCalculateRMSEInvalidDims
        function testCalculateRMSEInvalidDims(this)
            % dimensions do not match
           trueStates = ones(this.dimX, 200, 10);
           estimates = ones(this.dimX, 200);
           
           actualRMSE = Utility.calculateRMSE(trueStates, estimates);
           this.verifyEmpty(actualRMSE);
        end
    end
    
end

