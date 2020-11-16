classdef ValidatorTest< matlab.unittest.TestCase
    % Test cases for Validator.
    
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
        dimX = 2;
        dimU = 2;
        dimY = 4;
    end
    
    properties (Access = private)
        infMat;
        nanMat;
        dimXEye;
        dimYXEye;
        negativeDefiniteMatrix;
        
        pdMatrixDimX;
        pdMatrixDimU;
        pdMatrixDimY;
    end
    
    
    methods (TestClassSetup)
        function initProperties(this)
            this.infMat = inf(this.dimX, this.dimX);
            this.nanMat = nan(this.dimY, this.dimX); % nan
            this.dimXEye = eye(this.dimX);
            this.dimYXEye = eye(this.dimY, this.dimX);
            
            this.negativeDefiniteMatrix = [-1 0; 0 0];
            this.pdMatrixDimX = gallery('moler', this.dimX);
            this.pdMatrixDimU = gallery('minij', this.dimU);
            this.pdMatrixDimY = gallery('moler', this.dimY);
        end
    end
    
    methods (Test)
        
        %% testValidateSystemMatrix
        function testValidateSystemMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix'; 
            
            A = this.dimYXEye; % not square
            this.verifyError(@() Validator.validateSystemMatrix(A), expectedErrId);
            A = this.infMat; % square, but inf
            this.verifyError(@() Validator.validateSystemMatrix(A), expectedErrId);
            
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidDimensions';
            A = eye(this.dimX * 2); % square but wrong dimensions
            this.verifyError(@() Validator.validateSystemMatrix(A, this.dimX), expectedErrId);
            
            % finally, a successful run
            A = this.dimXEye;
            Validator.validateSystemMatrix(A, this.dimX);
        end
        
        %% testValidateInputMatrix
        function testValidateInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix'; 
           
            B = eye(1, 2); % wrong number of rows
            this.verifyError(@() Validator.validateInputMatrix(B, this.dimX), expectedErrId);
            B = [1 NaN; 0 0]; % nan
            this.verifyError(@() Validator.validateInputMatrix(B, this.dimX), expectedErrId);
            
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrixDims'; 
            
            B = eye(1, 2); % wrong number of cols
            this.verifyError(@() Validator.validateInputMatrix(B, 1, 3), expectedErrId);
            B = [1 NaN; 0 0]; % nan
            this.verifyError(@() Validator.validateInputMatrix(B, 2, 2), expectedErrId);
           
            % finally, two successful run, B matrix is supposed to be
            % square
            B = this.dimXEye;
            Validator.validateInputMatrix(B, this.dimX);
            Validator.validateInputMatrix(B, this.dimX, this.dimX);
        end
        
        %% testValidateMeasurementMatrix
        function testValidateMeasurementMatrix(this)
            expectedErrId = 'Validator:ValidateMeasurementMatrix:InvalidMeasMatrix'; 
           
            C = eye(this.dimY, this.dimX + 1); % wrong number of cols
            this.verifyError(@() Validator.validateMeasurementMatrix(C, this.dimX), expectedErrId);
            C = this.nanMat; % nan
            this.verifyError(@() Validator.validateMeasurementMatrix(C, this.dimX), expectedErrId);
           
            expectedErrId = 'Validator:ValidateMeasurementMatrix:InvalidMeasMatrixDims'; 
            C = eye(this.dimY + 1, this.dimX); % wrong dimensions
            this.verifyError(@() Validator.validateMeasurementMatrix(C, this.dimX, this.dimY), expectedErrId);
            C = this.nanMat; % nan
            this.verifyError(@() Validator.validateMeasurementMatrix(C, this.dimX, this.dimY), expectedErrId);
            
            % finally, a successful run
            C = this.dimYXEye;
            Validator.validateMeasurementMatrix(C, this.dimX, this.dimY);
        end
        
        %% testValidateMeasNoiseCovarianceMatrix
        function testValidateMeasNoiseCovarianceMatrix(this)
            expectedErrId = 'Validator:ValidateMeasNoiseCovarianceMatrix:InvalidCovDim';
            
            V = -this.pdMatrixDimY; % matrix is not a valid cov matrix
            this.verifyError(@() Validator.validateMeasNoiseCovarianceMatrix(V, this.dimY), expectedErrId);
            
            V = this.pdMatrixDimU; % cov, but wrong dimensions
            this.verifyError(@() Validator.validateMeasNoiseCovarianceMatrix(V, this.dimY), expectedErrId);
            
            % now a successful run
            V = this.pdMatrixDimY;
            expected = chol(this.pdMatrixDimY, 'Lower');
            actual = Validator.validateMeasNoiseCovarianceMatrix(V, this.dimY);
            this.verifyEqual(actual, expected);
            
            expectedErrId = 'Validator:ValidateMeasNoiseCovarianceMatrix:InvalidCov';
            V = -this.pdMatrixDimY; % matrix is not a valid cov matrix
            this.verifyError(@() Validator.validateMeasNoiseCovarianceMatrix(V), expectedErrId);
        end
        
        %% testValidateSysNoiseCovarianceMatrix
        function testValidateSysNoiseCovarianceMatrix(this)
            expectedErrId = 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCovDim';
            
            W = -this.pdMatrixDimX; % matrix is not a valid cov matrix
            this.verifyError(@() Validator.validateSysNoiseCovarianceMatrix(W, this.dimX), expectedErrId);
            
            W = this.pdMatrixDimY; % cov, but wrong dimensions
            this.verifyError(@() Validator.validateSysNoiseCovarianceMatrix(W, this.dimX), expectedErrId);
            
            % now a successful run
            W = this.pdMatrixDimX;
            expected = chol(this.pdMatrixDimX, 'Lower');
            actual = Validator.validateSysNoiseCovarianceMatrix(W, this.dimX);
            this.verifyEqual(actual, expected);
            
            expectedErrId = 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCov';
            W = -this.pdMatrixDimX; % matrix is not a valid cov matrix
            this.verifyError(@() Validator.validateSysNoiseCovarianceMatrix(W), expectedErrId);
        end
        
        %% testValidateCostMatrices
        function testValidateCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            R = eye(this.dimU);
            Q = this.dimYXEye; % matrix is not square
            this.verifyError(@() Validator.validateCostMatrices(Q, R, this.dimX, this.dimU), expectedErrId);
            
            Q = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() Validator.validateCostMatrices(Q, R, this.dimX, this.dimU), expectedErrId);
            Q = [1 0; Inf 0];
            this.verifyError(@() Validator.validateCostMatrices(Q, R, this.dimX, this.dimU), expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            Q = this.negativeDefiniteMatrix; % Q is not psd
            this.verifyError(@() Validator.validateCostMatrices(Q, R, this.dimX, this.dimU), expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            Q = this.dimXEye;
            
            R = this.dimYXEye; % matrix is not square
            this.verifyError(@() Validator.validateCostMatrices(Q, R, this.dimX, this.dimU), expectedErrId);
            R = [1 0; NaN 0];
            this.verifyError(@() Validator.validateCostMatrices(Q, R, this.dimX, this.dimU), expectedErrId);
            R = this.negativeDefiniteMatrix; % R is not pd
            this.verifyError(@() Validator.validateCostMatrices(Q, R, this.dimX, this.dimU), expectedErrId);
            
            % finally, a successful run
            Q = this.pdMatrixDimX;
            R = this.pdMatrixDimU;
            Validator.validateCostMatrices(Q, R, this.dimX, this.dimU)
        end
        
        %% testValidateDiscreteProbabilityDistribution
        function testValidateDiscreteProbabilityDistribution(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            prob = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() Validator.validateDiscreteProbabilityDistribution(prob), expectedErrId);
            
            prob = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() Validator.validateDiscreteProbabilityDistribution(prob), expectedErrId);
           
            prob = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() Validator.validateDiscreteProbabilityDistribution(prob), expectedErrId);
            
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbsNum';
            % sums up to one but wrong number of elements
            uniformProbs = 1 / 10 * ones(10, 1);
            expectedNumel = 11;
            this.verifyError(@() Validator.validateDiscreteProbabilityDistribution(uniformProbs, expectedNumel), ...
                expectedErrId);
            
            % finally, a successful run
            uniformProbs = 1 / 100 * ones(100, 1);
            Validator.validateDiscreteProbabilityDistribution(uniformProbs);
        end
        
        %% testValidateTransitionMatrix
        function testValidateTransitionMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrix';
            
            % invalid transition matrix: not square
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1];
            this.verifyError(@() Validator.validateTransitionMatrix(invalidTransitionMatrix), expectedErrId); 
            
            % invalid transition matrix: square, but row sums not 1
            invalidTransitionMatrix = [0.7 0.3; 0.2 0.5];
            this.verifyError(@() Validator.validateTransitionMatrix(invalidTransitionMatrix), expectedErrId); 
            
            % invalid transition matrix: square, but elements not all in [0,1]
            invalidTransitionMatrix = [0.7 0.3; 1.1 -0.1];
            this.verifyError(@() Validator.validateTransitionMatrix(invalidTransitionMatrix), expectedErrId); 

            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';            
            % invalid transition matrix: square, but wrong dimension
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.3 0.5];
            expectedDim = 4;
            this.verifyError(@() Validator.validateTransitionMatrix(invalidTransitionMatrix, expectedDim), expectedErrId); 
            
            % finally, a successful run
            expectedDim = 3;
            validTransitionMatrix = eye(3);
            Validator.validateTransitionMatrix(validTransitionMatrix, expectedDim)
        end
        
        %% testValidateHorizonLength
        function testValidateHorizonLength(this)
            expectedErrId = 'Validator:ValidateHorizonLength:InvalidHorizonLength';
            
            invalidHorizonLength = eye(3); % not a scalar
            this.verifyError(@() Validator.validateHorizonLength(invalidHorizonLength), expectedErrId);
            
            invalidHorizonLength = 0; % not a positive scalar
            this.verifyError(@() Validator.validateHorizonLength(invalidHorizonLength), expectedErrId);
            
            invalidHorizonLength = 1.5; % not an integer
            this.verifyError(@() Validator.validateHorizonLength(invalidHorizonLength), expectedErrId);
            
            invalidHorizonLength = inf; % not finite
            this.verifyError(@() Validator.validateHorizonLength(invalidHorizonLength), expectedErrId);
            
            % finally, a successful run
            validHorizonLength = 3;
            Validator.validateHorizonLength(validHorizonLength);
        end
        
        %% testValidateSequenceLength
        function testValidateSequenceLength(this)
            expectedErrId = 'Validator:ValidateSequenceLength:InvalidSequenceLength';
            
            invalidSequenceLength = eye(3); % not a scalar
            this.verifyError(@() Validator.validateSequenceLength(invalidSequenceLength), expectedErrId);
            
            invalidSequenceLength = 0; % not a positive scalar
            this.verifyError(@() Validator.validateSequenceLength(invalidSequenceLength), expectedErrId);
            
            invalidSequenceLength = 1.5; % not an integer
            this.verifyError(@() Validator.validateSequenceLength(invalidSequenceLength), expectedErrId);
            
            invalidSequenceLength = inf; % not finite
            this.verifyError(@() Validator.validateSequenceLength(invalidSequenceLength), expectedErrId);
            
            % finally, a successful run
            validSequenceLength = 3;
            Validator.validateSequenceLength(validSequenceLength);
        end
        
        %% testValidateMaxPacketDelay
        function testMaxPacketDelay(this)
            expectedErrId = 'Validator:ValidateMaxPacketDelay:InvalidMaxPacketDelay';
            
            invalidMaxPacketDelay = eye(3); % not a scalar
            this.verifyError(@() Validator.validateMaxPacketDelay(invalidMaxPacketDelay), expectedErrId);
            
            invalidMaxPacketDelay = -1; % not a nonnegative scalar
            this.verifyError(@() Validator.validateMaxPacketDelay(invalidMaxPacketDelay), expectedErrId);
            
            invalidMaxPacketDelay = 1.5; % not an integer
            this.verifyError(@() Validator.validateMaxPacketDelay(invalidMaxPacketDelay), expectedErrId);
            
            invalidMaxPacketDelay = nan; % not finite nor inf
            this.verifyError(@() Validator.validateMaxPacketDelay(invalidMaxPacketDelay), expectedErrId);
            
            % finally, a successful run
            validMaxPacketDelay = 3;
            Validator.validateMaxPacketDelay(validMaxPacketDelay);
            validMaxPacketDelay = inf;
            Validator.validateMaxPacketDelay(validMaxPacketDelay);
        end
    end
   
end

