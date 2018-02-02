classdef NominalPredictiveControllerTest < matlab.unittest.TestCase
    % Test cases for NominalPredictiveController.
    
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
        
    properties (Constant)
        absTol = 1e-5;
    end
    
    properties
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
        
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            % use (noise-free) stirred tank example (Example 6.15, p. 500-501) from
            %
            % Huibert Kwakernaak, and Raphael Sivan, 
            % Linear Optimal Control Systems,
            % Wiley-Interscience, New York, 1972.
            %
            
            V = diag([0.01, 1]); % in the book: z = V*x
            this.A = diag([0.9512, 0.9048]);
           
            this.B = [4.877 4.877; -1.1895 3.569];
            this.Q = V * diag([50 0.02]) * V; % R_3 in the book
            this.R = diag([1/3, 3]); % R_2 in the book
                       
            this.expectedL = -[ 0.07125 -0.07029;
                                0.01357 0.04548]; % F in the book
    
            this.sequenceLength = 10;
            this.dimX = 2;
            this.dimU = 2;
            this.horizonLength = this.sequenceLength;
            this.initialPlantState = [0.1; 0];
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            this.controllerUnderTest = NominalPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength);
        end
    end
    
    methods (Access = private)
        %% computeTrajectories
        function [stateTrajectory, inputTrajectory] = computeTrajectories(this)
            stateTrajectory = zeros(this.dimX, this.horizonLength + 1);
            inputTrajectory = zeros(this.dimU, this.horizonLength);
            
            stateTrajectory(:, 1) = this.initialPlantState;
            for j = 1:this.horizonLength
                % compute new input
                inputTrajectory(:, j) = this.expectedL * stateTrajectory(:, j);
                stateTrajectory(:, j + 1) = this.A * stateTrajectory(:, j) + this.B * inputTrajectory(:, j);
            end
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
    end
    
    methods (Test)
        
        %% testNominalPredictiveControllerInvalidSysMatrix
        function testNominalPredictiveControllerInvalidSysMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
            
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() NominalPredictiveController(invalidA, this.B, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidA = inf(this.dimX, this.dimX); % square but inf
            this.verifyError(@() NominalPredictiveController(invalidA, this.B, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
        end
        
        %% testNominalPredictiveControllerInvalidInputMatrix
        function testNominalPredictiveControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidB = eye(this.dimX + 1, this.dimU); % invalid dims
            this.verifyError(@() NominalPredictiveController(this.A, invalidB, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidB = eye(this.dimX, this.dimU); % correct dims, but nan
            invalidB(1, 1) = nan;
            this.verifyError(@() NominalPredictiveController(this.A, invalidB, this.Q, this.R, this.sequenceLength), ...
                expectedErrId);
        end
        
        %% testNominalPredictiveControllerInvalidCostMatrices
        function testNominalPredictiveControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() NominalPredictiveController(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() NominalPredictiveController(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() NominalPredictiveController(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() NominalPredictiveController(this.A, this.B, invalidQ, this.R, this.sequenceLength), ...
                expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() NominalPredictiveController(this.A, this.B, this.Q, invalidR, this.sequenceLength), ...
                expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() NominalPredictiveController(this.A, this.B, this.Q, invalidR, this.sequenceLength), ...
                expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() NominalPredictiveController(this.A, this.B, this.Q, invalidR, this.sequenceLength), ...
                expectedErrId);
        end
        
        %% testNominalPredictiveController
        function testNominalPredictiveController(this)
            controller = NominalPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength);
            
            this.verifyEqual(controller.L, this.expectedL, 'AbsTol', NominalPredictiveControllerTest.absTol);
        end
        
        %% testChangeSequenceLength
        function testChangeSequenceLength(this)
            this.assertEqual(this.controllerUnderTest.sequenceLength, this.sequenceLength);
            
            newSeqLength = this.sequenceLength + 2;
            this.controllerUnderTest.changeSequenceLength(newSeqLength);
            
            this.verifyEqual(this.controllerUnderTest.sequenceLength, newSeqLength);
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
    end
    
end
