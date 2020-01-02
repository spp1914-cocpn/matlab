classdef IMMBasedRecedingHorizonControllerTest < matlab.unittest.TestCase
    % Test cases for IMMBasedRecedingHorizonController.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        C;
        sequenceLength;
        dimX;
        dimU;
        dimY;
        
        modeTransitionMatrix;
                        
        W;
        V;
        measNoise; 
        x0;
        x0Cov;
        
        caDelayProbs;
        maxMeasDelay;
        
        horizonLength;
        
        augA;
        augB;
        augC;
        augQ;
        augR;
        augX0;
        augX0Cov;         
        
        F;
        G;
        
        controllerUnderTest;
  
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.dimX = 3;
            this.dimU = 2;
            this.dimY = 2;
            
            this.horizonLength = 3;
            
            this.sequenceLength = 2; % so we have 3 modes
            this.caDelayProbs = ones(1, 5) / 5;
                        
            % from the delay probs and the 3 modes we get the following
            % transition matrix
            this.modeTransitionMatrix = [1/5 4/5 0; 1/5 1/5 3/5; 1/5 1/5 3/5];
            
            this.A = eye(this.dimX);
            this.B = ones(this.dimX, this.dimU);
            this.C = 2 * ones(this.dimY, this.dimX);
            
            this.Q = eye(this.dimX);
            this.R = eye(this.dimU) * 2;
              
            this.maxMeasDelay = 1;
            
            this.W = eye(this.dimX) * 0.01;
            this.V = eye(this.dimY) * 0.001;
            
            this.measNoise = Gaussian(zeros(this.dimY, 1), this.V);
            
            this.x0 = ones(3,1);
            this.x0Cov = eye(this.dimX) * 0.2;
            
            this.initAdditionalProperties();
            
            this.controllerUnderTest = IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov);
        end
    end
    
    methods (Access = private)
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            [dimAugmentedState, this.F, this.G, this.augA, this.augB, this.augQ, this.augR] = ...
                Utility.performModelAugmentation(this.sequenceLength, this.dimX, this.dimU, this.A, this.B, this.Q, this.R);
            
            dimEta = dimAugmentedState - this.dimX;
            
            this.augX0 = [this.x0(:); zeros(dimEta, 1)];
            this.augX0Cov = blkdiag(this.x0Cov, zeros(dimEta));
        end
        
        %% computeCostate
        function P = computeCostate(this, transitionMatrix)
            dimAugState = size(this.augA, 1);
            numModes = this.sequenceLength + 1;
            
            P = zeros(dimAugState, dimAugState, numModes, this.horizonLength);
            P(1:this.dimX, 1:this.dimX, :, end) = repmat(this.Q, 1, 1, numModes); %P_k+K
                        
            for n=this.horizonLength-1:-1:1
                for r=1:numModes
                    Y = zeros(dimAugState);
                    for i=1:numModes
                        Y = Y + transitionMatrix(r,i) * P(:, :, i, n+1);
                    end
                    M = this.augR(:, :, r) + transpose(this.augB(:, :, r)) * Y * this.augB(:, :, r);
                    M = (M + M') / 2;
                    P_tmp = this.augQ(:, :, r) + transpose(this.augA(:, :, r)) * Y * this.augA(:, :, r) ...
                        - transpose(this.augA(:, :, r)) * Y * this.augB(:, :, r) * pinv(M) * transpose(this.augB(:, :, r)) ...
                        * Y * this.augA(:, :, r);
                    % ensure symmetry
                    P(:, :, r, n) = (P_tmp + P_tmp') / 2;
                end
            end
        end
        
        %% computeExpectedInputForState
        function inputSeq = computeExpectedInputForState(this, modeStates, modeProbs, transitionMatrix)
            numModes = this.sequenceLength + 1;
                                    
            P = this.computeCostate(transitionMatrix);
            
            L = modeProbs(1) * this.augR(:, :, 1);
            S = zeros(this.sequenceLength * this.dimU, 1);
            for r =1:numModes
                for i=1:numModes
                    BPB = transpose(this.augB(:, :, i)) * P(:, :, r, 1) * this.augB(:, :, i);
                    % ensure symmetry
                    L = L + transitionMatrix(i,r) * modeProbs(i) * (BPB + BPB') / 2;
                    S = S + transitionMatrix(i,r) * modeProbs(i) * transpose(this.augB(:, :, i)) * P(:, :, r, 1) * this.augA(:, :, i) * modeStates(:, i);
                end
            end
            inputSeq = -pinv(L) * S;
        end
    end
    
    methods(Test)
        %% testIMMBasedRecedingHorizonControllerInvalidSystemMatrix
        function testIMMBasedRecedingHorizonControllerInvalidSystemMatrix(this)
             expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';             
             invalidSysMatrix = eye(this.dimX, this.dimX + 1); % not square
             this.verifyError(@() IMMBasedRecedingHorizonController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
             
             invalidSysMatrix = eye(this.dimX, this.dimX); % square but not finite
             invalidSysMatrix(1, end) = inf;
             this.verifyError(@() IMMBasedRecedingHorizonController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidInputMatrix
        function testIMMBasedRecedingHorizonControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidInputMatrix = eye(this.dimX +1, this.dimU); % invalid dims
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
             
            invalidInputMatrix = eye(this.dimX, this.dimU); % correct dims, but not finite
            invalidInputMatrix(1, end) = nan;
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidMeasurementMatrix
        function testIMMBasedRecedingHorizonControllerInvalidMeasurementMatrix(this)
            expectedErrId = 'Validator:ValidateMeasurementMatrix:InvalidMeasMatrix';
            
            invalidMeasMatrix = eye(this.dimY, this.dimX + 1); % invalid dims
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
             
            invalidMeasMatrix = eye(this.dimY, this.dimX); % correct dims, but not finite
            invalidMeasMatrix(1, end) = nan;
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidCostMatrices
        function testIMMBasedRecedingHorizonControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = eye(this.dimX); % Q is not symmetric
            invalidQ(1, end) = 1;
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidCaDelayProbs
        function testIMMBasedRecedingHorizonControllerInvalidCaDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidMaxMeasDelay
        function testIMMBasedRecedingHorizonControllerInvalidMaxMeasDelay(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:InvalidMaxMeasDelay';
            
            invalidMaxMeasDelay = [1 2]; % not a scalar
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidMaxMeasDelay = -1; % negative scalar
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidMaxMeasDelay = 1.5; % not an integer
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidMaxMeasDelay = inf; % not finite
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
         %% testIMMBasedRecedingHorizonControllerInvalidSysNoiseCovariance
        function testIMMBasedRecedingHorizonControllerInvalidSysNoiseCovariance(this)
            expectedErrId = 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCovDim';
            
            invalidW = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, invalidW, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
                       
            invalidW = ones(this.dimU); % W is not pd
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, invalidW, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        
        %% testIMMBasedRecedingHorizonControllerInvalidMeasNoise
        function testIMMBasedRecedingHorizonControllerInvalidMeasNoise(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:InvalidMeasNoise';
            
            invalidNoise = eye(this.dimY + 1, this.dimY); % not a Distribution
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, invalidNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
                       
            invalidNoise = Gaussian(); % wrong dimension
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, invalidNoise, this.horizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidHorizonLegth
        function testIMMBasedRecedingHorizonControllerInvalidHorizonLegth(this)
            expectedErrId = 'Validator:ValidateHorizonLength:InvalidHorizonLength';
            
            invalidHorizonLength = [1 2]; % not a scalar
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, invalidHorizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidHorizonLength = -1; % negative scalar
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, invalidHorizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidHorizonLength = 0; 
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, invalidHorizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidHorizonLength = 1.5; % not an integer
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, invalidHorizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
            
            invalidHorizonLength = inf; % not finite
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, invalidHorizonLength, ...
                this.x0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidX0
        function testIMMBasedRecedingHorizonControllerInvalidX0(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:InvalidX0';
            
            invalidX0 = ones(this.dimX +1); % not a vector
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                invalidX0, this.x0Cov), expectedErrId);
            
            invalidX0 = ones(this.dimX +1, 1); % wrong dimensions
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                invalidX0, this.x0Cov), expectedErrId);
            
            invalidX0 = ones(this.dimX, 1); % not finite
            invalidX0(1) = nan;
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                invalidX0, this.x0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidX0Cov
        function testIMMBasedRecedingHorizonControllerInvalidX0Cov(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:InvalidX0Cov';
            
            invalidX0Cov = ones(this.dimX +1, 1); % not a matrix
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, invalidX0Cov), expectedErrId);
            
            invalidX0Cov = ones(this.dimX +1, this.dimX); % wrong dimensions
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, invalidX0Cov), expectedErrId);
            
            invalidX0Cov = ones(this.dimX); % not pd
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, invalidX0Cov), expectedErrId);
        end
        
        %% testIMMBasedRecedingHorizonControllerInvalidFlag
        function testIMMBasedRecedingHorizonControllerInvalidFlag(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:InvalidUseMexFlag';
            
            invalidUseMexFlag = 'invalid'; % not a flag            
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov, invalidUseMexFlag), expectedErrId);
            
            invalidUseMexFlag = [false true]; % not a flag            
            this.verifyError(@() IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov, invalidUseMexFlag), expectedErrId);
        end
%%
%%
        %% testIMMBasedRecedingHorizonController
        function testIMMBasedRecedingHorizonController(this)
            controller = IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov);
                        
            this.verifyEqual(controller.getControllerPlantState(), this.x0); 
            this.verifyFalse(controller.requiresExternalStateEstimate); % does not require a filter or state feedback
            
            % by default, we use the mex implementation to obtain the
            % controller costate matrices
            this.verifyTrue(controller.useMexImplementation);
        end
%%
%%
    %% testChangeCaDelayProbs
    function testChangeCaDelayProbs(this)
        this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
               
        numModes = this.sequenceLength + 1;
        % expected input sequence, buffer is initially empty
        eta = zeros(size(this.F, 1), 1);
        modeStates = repmat([this.x0; eta], 1, numModes);
        modeProbs = zeros(1, numModes);
        modeProbs(end) = 1;
        expectedInputSequence = this.computeExpectedInputForState(modeStates, modeProbs, this.modeTransitionMatrix);

        actualInputSequence = this.controllerUnderTest.computeControlSequence();
        this.assertEqual(actualInputSequence, expectedInputSequence, 'AbsTol', IMMBasedRecedingHorizonControllerTest.absTol);
        this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
        
        % compute the new expected controller state
        newModeStates = zeros(size(this.augA, 1), numModes);
        for i=1:numModes
            newModeStates(:, i) = this.augA(:, :, i) * modeStates(:, i) + this.augB(:, :, i) * expectedInputSequence;
        end
        
        % now change the ca delay probs
        newCaDelayProbs = ones(1, 10) / 10;
        newTransitionMatrix = Utility.calculateDelayTransitionMatrix(...
                Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, numModes));
        this.controllerUnderTest.changeCaDelayProbs(newCaDelayProbs);            
                    
        expectedInputSequence = this.computeExpectedInputForState(newModeStates, ...
            this.modeTransitionMatrix' * modeProbs(:), newTransitionMatrix);
        actualInputSequence = this.controllerUnderTest.computeControlSequence();
        this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', IMMBasedRecedingHorizonControllerTest.absTol); 
    end
%%
%%
        %% testSetEtaStateInvalidEta
        function testSetEtaStateInvalidEta(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:SetEtaState:InvalidEta';
            
            etaDim = this.dimU * this.sequenceLength * (this.sequenceLength - 1) / 2;
            inputSeq = ones(this.dimU * this.sequenceLength, 1);
            
            invalidEta = this; % not a vector
            this.verifyError(@() this.controllerUnderTest.setEtaState(invalidEta, inputSeq), ...
                expectedErrId);
            
            invalidEta = ones(etaDim, etaDim); % not a vector
            this.verifyError(@() this.controllerUnderTest.setEtaState(invalidEta, inputSeq), ...
                expectedErrId);            
            
            invalidEta = ones(etaDim +1, 1); % invalid dimension
            this.verifyError(@() this.controllerUnderTest.setEtaState(invalidEta, inputSeq), ...
                expectedErrId);            
            
            invalidEta = ones(etaDim, 1);
            invalidEta(1) = -inf; % not finite
            this.verifyError(@() this.controllerUnderTest.setEtaState(invalidEta, inputSeq), ...
                expectedErrId);
        end
        
        %% testSetEtaStateInvalidInputSeq
        function testSetEtaStateInvalidInputSeq(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:SetEtaState:InvalidInputSeq';
            
            etaDim = this.dimU * this.sequenceLength * (this.sequenceLength - 1) / 2;
            eta = ones(etaDim, 1);
            
            invalidInputSeq = this; % not a vector
            this.verifyError(@() this.controllerUnderTest.setEtaState(eta, invalidInputSeq), ...
                expectedErrId);
            
            invalidInputSeq = ones(this.dimU * this.sequenceLength); % not a vector
            this.verifyError(@() this.controllerUnderTest.setEtaState(eta, invalidInputSeq), ...
                expectedErrId);            
            
            invalidInputSeq = ones(this.dimU * this.sequenceLength +1, 1); % invalid dimension
            this.verifyError(@() this.controllerUnderTest.setEtaState(eta, invalidInputSeq), ...
                expectedErrId);            
            
            invalidInputSeq = ones(this.dimU * this.sequenceLength, 1);
            invalidInputSeq(1) = -inf; % not finite
            this.verifyError(@() this.controllerUnderTest.setEtaState(eta, invalidInputSeq), ...
                expectedErrId);
        end
%%
%%
        %% testSetControllerPlantStateInvalidState
        function testSetControllerPlantStateInvalidState(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:SetControllerPlantState:InvalidState';
            
            invalidState = this; % not a Distribution
            this.verifyError(@() this.controllerUnderTest.setControllerPlantState(invalidState), expectedErrId);
            
            invalidState = Gaussian(0, 1); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.setControllerPlantState(invalidState), expectedErrId);
        end
        
        %% testSetControllerPlantState
        function testSetControllerPlantState(this)
            stateMean = 0.2 * ones(this.dimX, 1);                        
            initialState = Gaussian(stateMean, 42 * eye(this.dimX));
            
            this.controllerUnderTest.setControllerPlantState(initialState);
            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), stateMean);
        end    
%%
%%
        %% testGetLastComputationMeasurementDataNoMeasurements
        function testGetLastComputationMeasurementDataNoMeasurements(this)
            % compute a sequence without using measurements or mode
            % observations
            this.controllerUnderTest.computeControlSequence();
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
            
            % compute a sequence without using measurements but with mode
            % observations
            modeObservation = 1;
            modeDelay = 1;
            this.controllerUnderTest.computeControlSequence([], [], modeObservation, modeDelay);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testGetLastComputationMeasurementDataNoMeasDiscarded
        function testGetLastComputationMeasurementDataNoMeasDiscarded(this)
            this.controllerUnderTest.computeControlSequence(); % ensures that controller/filter is initialized
            
            measurements = ones(this.dimY, 2); % two meas, column-wise arranged
            delays = [0 1];
                                    
            this.controllerUnderTest.computeControlSequence(measurements, delays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 0);
            
            % now the same but with a single mode observation
            % should not result in a change
            modeObservation = 1;
            modeDelay = 1;
            this.controllerUnderTest.computeControlSequence(measurements, delays, modeObservation, modeDelay);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testGetLastComputationMeasurementData
        function testGetLastComputationMeasurementData(this)                        
            this.controllerUnderTest.computeControlSequence(); % ensures that controller/filter is initialized
                        
            measurements = ones(this.dimY, 4); % four meas, column-wise arranged
            delays = [0 2 1 3]; % two meas are too old
            
            this.controllerUnderTest.computeControlSequence(measurements, delays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 2);
            
            % now the same but with a single mode observation
            % should not result in a change
            modeObservation = 1;
            modeDelay = 1;
            this.controllerUnderTest.computeControlSequence(measurements, delays, modeObservation, modeDelay);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 2);
        end
%%
%%
        %% testComputeControlSequenceInvalidNumArgs
        function testComputeControlSequenceInvalidNumArgs(this)
            expectedErrId = 'IMMBasedRecedingHorizonController:ComputeControlSequence:InvalidNumArgs';
            
            measurement = ones(this.dimY, 1); % only 2 arguments given
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurement), ...
                expectedErrId);
            
            delay = 1;
            mode = 1; % 4 arguments given
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurement, delay, mode), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceInvalidMeasurements
        function testComputeControlSequenceInvalidMeasurements(this)
            this.controllerUnderTest.computeControlSequence(); % ensures that controller/filter is initialized
            
            expectedErrId = 'IMMBasedRecedingHorizonController:ComputeControlSequence:InvalidMeas';
            
            % two measurements, column-wise arranged, but of wrong dimension
            invalidMeasurements = ones(this.dimY + 1, 2);
            delays = [0 1];            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(invalidMeasurements, delays), ...
                expectedErrId);          
      
            % existence of mode observation shoud not change anything
            modeObservation = 2;
            modeDelay = 1;
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(invalidMeasurements, delays, modeObservation, modeDelay), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceInvalidMeasDelays
        function testComputeControlSequenceInvalidMeasDelays(this)
            this.controllerUnderTest.computeControlSequence(); % ensures that controller/filter is initialized
            
            expectedErrId = 'Filter:InvalidMeasDelay';
            
            modeObservation = 2;
            modeDelay = 1;
            measurements = ones(this.dimY, 3); % three measurements, column-wise arranged
            invalidDelays = [0 -1 1]; % negative delay            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 0.5 1]; % fractional delay            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays, modeObservation, modeDelay), ...
                expectedErrId);
            
            invalidDelays = [0 1]; % not enough entries           
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays, modeObservation, modeDelay), ...
                expectedErrId);
        end    
        
         %% testComputeControlSequenceInvalidModes
        function testComputeControlSequenceInvalidModes(this)            
            this.controllerUnderTest.computeControlSequence(); % ensures that controller/filter is initialized
            
            expectedErrId = 'Filter:InvalidModeObservations';
                        
            % two mode observations, column-wise arranged, but of wrong dimension
            invalidModes = ones(2, 2);
            delays = [0 1];        
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([],[], invalidModes, delays), ...
                expectedErrId);
            
            % two mode observations, but out of scope
            invalidModes = [-1 1];
            delays = [0 1];
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], invalidModes, delays), ...
                expectedErrId);
            
            % two mode observations, but out of scope
            invalidModes = [1 this.sequenceLength + 2];
            delays = [0 1];
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], invalidModes, delays), ...
                expectedErrId);
            
            % two mode observations, but out of scope
            invalidModes = [0 1];
            delays = [0 1];
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], invalidModes, delays), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceInvalidModeDelays
        function testComputeControlSequenceInvalidModeDelays(this)
            this.controllerUnderTest.computeControlSequence(); % ensures that controller/filter is initialized
            
            expectedErrId = 'Filter:InvalidModeDelays';
                        
            modes = [1 2];
            invalidDelays = 1; % too few        
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([],[], modes, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [-1 0]; % negative
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], modes, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 0.5]; % fractional
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], modes, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 inf]; % not finite
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], modes, invalidDelays), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceZeroStateNoMeasurementsNoModes
        function testComputeControlSequenceZeroStateNoMeasurementsNoModes(this)
            zeroState = zeros(this.dimX, 1);
            this.controllerUnderTest.setControllerPlantState(Gaussian(zeroState, this.x0Cov));
            
            this.assertTrue(this.controllerUnderTest.useMexImplementation);
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
            % the initial state is the origin, and we have a linear control
            % law, so perform a sanity check
            expectedInputSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence();
            
            this.verifyEqual(actualInputSequence, expectedInputSequence);
            % the controller state should not have changed, as we have no
            % measurements that could be used to correct it
            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);           
        end
        
        %% testComputeControlSequenceZeroStateNoMeasurementsNoModesNoMex
        function testComputeControlSequenceZeroStateNoMeasurementsNoModesNoMex(this)
             controller = IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov, false);
            this.assertFalse(controller.useMexImplementation);            
            
            zeroState = zeros(this.dimX, 1);
            
            controller.setControllerPlantState(Gaussian(zeroState, this.x0Cov));
                        
            this.assertEqual(controller.getControllerPlantState(), zeroState);
            % the initial state is the origin, and we have a linear control
            % law, so perform a sanity check
            expectedInputSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            actualInputSequence = controller.computeControlSequence();
            
            this.verifyEqual(actualInputSequence, expectedInputSequence);
            % the controller state should not have changed, as we have no
            % measurements that could be used to correct it
            this.verifyEqual(controller.getControllerPlantState(), zeroState);           
        end
        
        %% testComputeControlSequenceNoMeasurementsNodeModes
        function testComputeControlSequenceNoMeasurementsNodeModes(this)
            this.assertTrue(this.controllerUnderTest.useMexImplementation);            
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
            
            numModes = this.sequenceLength + 1;
            % expected input sequence, buffer is initially empty
            eta = zeros(size(this.F, 1), 1);
            modeStates = repmat([this.x0; eta], 1, numModes);
            modeProbs = zeros(1, numModes);
            modeProbs(end) = 1;
            expectedInputSequence = this.computeExpectedInputForState(modeStates, modeProbs, this.modeTransitionMatrix);
           
            actualInputSequence = this.controllerUnderTest.computeControlSequence();
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', IMMBasedRecedingHorizonControllerTest.absTol);
            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
            
            % compute the new expected controller state
            newModeStates = zeros(size(this.augA, 1), numModes);
            for i=1:numModes
                newModeStates(:, i) = this.augA(:, :, i) * modeStates(:, i) + this.augB(:, :, i) * expectedInputSequence;
            end            
            expectedInputSequence = this.computeExpectedInputForState(newModeStates, ...
                this.modeTransitionMatrix' * modeProbs(:), this.modeTransitionMatrix);
            actualInputSequence = this.controllerUnderTest.computeControlSequence();
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', IMMBasedRecedingHorizonControllerTest.absTol);            
        end
        
        %% testComputeControlSequenceNoMeasurementsNodeModesNoMex
        function testComputeControlSequenceNoMeasurementsNodeModesNoMex(this)
            controller = IMMBasedRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.caDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.measNoise, this.horizonLength, ...
                this.x0, this.x0Cov, false);
            
            this.assertFalse(controller.useMexImplementation);            
            this.assertEqual(controller.getControllerPlantState(), this.x0);
            
            numModes = this.sequenceLength + 1;
            % expected input sequence, buffer is initially empty
            eta = zeros(size(this.F, 1), 1);
            modeStates = repmat([this.x0; eta], 1, numModes);
            modeProbs = zeros(1, numModes);
            modeProbs(end) = 1;
            expectedInputSequence = this.computeExpectedInputForState(modeStates, modeProbs, this.modeTransitionMatrix);
           
            actualInputSequence = controller.computeControlSequence();
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', IMMBasedRecedingHorizonControllerTest.absTol);
            this.verifyEqual(controller.getControllerPlantState(), this.x0);
            
            % compute the new expected controller state
            newModeStates = zeros(size(this.augA, 1), numModes);
            for i=1:numModes
                newModeStates(:, i) = this.augA(:, :, i) * modeStates(:, i) + this.augB(:, :, i) * expectedInputSequence;
            end            
            expectedInputSequence = this.computeExpectedInputForState(newModeStates, ...
                this.modeTransitionMatrix' * modeProbs(:), this.modeTransitionMatrix);
            actualInputSequence = controller.computeControlSequence();
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', IMMBasedRecedingHorizonControllerTest.absTol);            
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
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
        end
%%
%%
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            trajectoryLength = 10;
            inputs = ones(this.dimU, trajectoryLength);
            expectedErrId = 'IMMBasedRecedingHorizonController:DoCostsComputation';
            
            invalidStates = ones(this.dimX, trajectoryLength + 2); %  trajectory too long
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStates, inputs), ...
                expectedErrId);
            
            invalidStates = ones(this.dimX, trajectoryLength); %  trajectory too short
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStates, inputs), ...
                expectedErrId);
        end
        
         %% testDoCostsComputation
        function testDoCostsComputation(this)
            trajectoryLength = 100;  
            inputs = ones(this.dimU, trajectoryLength);
            states = 2.25 * ones(this.dimX, trajectoryLength + 1);
            
            % cost function is simple a quadratic
            expectedCosts = states(:, end)' * this.Q * states(:, end);
            for j=1:trajectoryLength
                expectedCosts = expectedCosts + states(:, j)' * this.Q * states(:, j) ...
                    + inputs(:, j)' * this.R * inputs(:, j);
            end
                        
            actualCosts = this.controllerUnderTest.computeCosts(states, inputs);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
        end
    end
end

