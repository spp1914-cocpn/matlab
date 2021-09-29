classdef LinearPlantTest < matlab.unittest.TestCase
    % Test cases for LinearPlant.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Access = private)
        A; % sysMatrix
        B;
        dimX;
        dimU;
        W; % process noise cov
        
        emptySysInput;
        sysInput;
        invalidSysInput;
        invalidDimSysInput;
        
        invalidAMatrix;
        invalidDimsAMatrix;
        invalidBMatrix;
        invalidDimsBMatrix;
        
        plantNoise;
        plantUnderTest;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.dimX = 3;
            this.dimU = 2;
            this.A = eye(this.dimX);
            this.B = ones(this.dimX, this.dimU);
            
            this.W = gallery('moler', this.dimX);
            this.plantNoise = Gaussian(zeros(this.dimX, 1), this.W);
            
            this.plantUnderTest = LinearPlant(this.A, this.B, this.W);
            
            this.emptySysInput = [];
            this.sysInput = [4 2]';
            this.invalidSysInput = [1 nan];
            this.invalidDimSysInput = 42;
            
            this.invalidAMatrix = inf(this.dimX, this.dimX);
            this.invalidDimsAMatrix = this.B;
            this.invalidBMatrix = inf(this.dimX, this.dimU);
            this.invalidDimsBMatrix = this.B'; % flipped dimensions
        end
    end
    
    methods (Test)
        
        %% testLinearPlant
        function testLinearPlant(this)
            linearPlant = LinearPlant(this.A, this.B, this.W);
            
            this.verifyEqual(linearPlant.dimState, this.dimX);
            this.verifyEqual(linearPlant.inputMatrix, this.B);
            
            noise = linearPlant.noise;
            this.verifyClass(noise, ?Gaussian);
            this.verifyEqual(noise, this.plantNoise);
        end
         
        %% testLinearPlantSysNoiseMatrix
        function testLinearPlantSysNoiseMatrix(this)
            subW = this.W(1:2, 1:2);
            G = ones(this.dimX, 2);
            actualNoise = Gaussian(zeros(2,1), subW);
            
            linearPlant = LinearPlant(this.A, this.B, subW, G);
            this.verifyEqual(linearPlant.dimState, this.dimX);
            this.verifyEqual(linearPlant.inputMatrix, this.B);
            noise = linearPlant.noise;
            
            this.verifyClass(noise, ?Gaussian);
            this.verifyEqual(noise, actualNoise);
            
        end
        
        %% testSetNoise
        function testSetNoise(this)
            W_new = this.W * 0.75; % less noise now
            expectedNoise = Gaussian(zeros(this.dimX, 1), W_new);
            
            this.plantUnderTest.setNoise(W_new);
            
            newNoise = this.plantUnderTest.noise;
            this.verifyClass(newNoise, ?Gaussian);
            this.verifyEqual(newNoise, expectedNoise);
            
            % now specify covariance in terms of a vector
            W_new = 1:1:this.dimX;
            expectedNoise = Gaussian(zeros(this.dimX, 1), W_new);
            
            this.plantUnderTest.setNoise(W_new);
            
            newNoise = this.plantUnderTest.noise;
            this.verifyClass(newNoise, ?Gaussian);
            this.verifyEqual(newNoise, expectedNoise);
        end
        
        %% testSetNoiseInvalidCov
        function testSetNoiseInvalidCov(this)
            expectedErrId = 'LinearPlant:SetNoise:InvalidCovariance';
                        
            this.assertEmpty(this.plantUnderTest.sysNoiseMatrix); % no G matrix
            
            subW = this.W(1:2, 1:2);
            this.verifyError(@() this.plantUnderTest.setNoise(subW), expectedErrId);
                        
            % now specify covariance in terms of a vector
            subW = [1 1];
            this.verifyError(@() this.plantUnderTest.setNoise(subW), expectedErrId);
        end
        
        %% testSetNoiseSysNoiseMatrixInvalidCov
        function testSetNoiseSysNoiseMatrixInvalidCov(this)
            expectedErrId = 'LinearPlant:SetNoise:InvalidCovariance';
            
            G = ones(this.dimX, 2);
            this.plantUnderTest.setSystemNoiseMatrix(G);
            this.assertEqual(this.plantUnderTest.sysNoiseMatrix, G);
                       
            this.verifyError(@() this.plantUnderTest.setNoise(this.W), expectedErrId);
                        
            % now specify covariance in terms of a vector
            W_new = 1:1:this.dimX;
            this.verifyError(@() this.plantUnderTest.setNoise(W_new), expectedErrId);
        end
        
        %% testSetNoiseSysNoiseMatrix
        function testSetNoiseSysNoiseMatrix(this)
            subW = this.W(1:2, 1:2);
            G = ones(this.dimX, 2);
            this.plantUnderTest.setSystemNoiseMatrix(G);
            this.assertEqual(this.plantUnderTest.sysNoiseMatrix, G);
                        
            expectedNoise = Gaussian(zeros(2, 1), subW);
            
            this.plantUnderTest.setNoise(subW);
            
            newNoise = this.plantUnderTest.noise;
            this.verifyClass(newNoise, ?Gaussian);
            this.verifyEqual(newNoise, expectedNoise);
            
            % now specify covariance in terms of a vector
            W_new = [1 1];
            expectedNoise = Gaussian(zeros(2, 1), W_new);
            
            this.plantUnderTest.setNoise(W_new);
            
            newNoise = this.plantUnderTest.noise;
            this.verifyClass(newNoise, ?Gaussian);
            this.verifyEqual(newNoise, expectedNoise);
        end
        
        %% testSetSystemInputInvalidInput
        function testSetSystemInputInvalidInput(this)
            expectedErrId = 'LinearPlant:InvalidInput';
            
            % input of wrong dimension
            this.verifyError(@() this.plantUnderTest.setSystemInput(this.invalidDimSysInput), expectedErrId);
            % invalid input vector
            this.verifyError(@() this.plantUnderTest.setSystemInput(this.invalidSysInput), expectedErrId);
        end
        
        %% testSetSystemInputEmptyInput
        function testSetSystemInputEmptyInput(this)
            this.plantUnderTest.setSystemInput(this.emptySysInput);
            
            actualInput = this.plantUnderTest.getSystemInput();
            this.verifyEmpty(actualInput);
        end
        
        %% testSetSystemInput
        function testSetSystemInput(this)
            expectedInput = this.sysInput;
            this.plantUnderTest.setSystemInput(expectedInput);
            
            actualInput = this.plantUnderTest.getSystemInput();
            this.verifyEqual(actualInput, expectedInput);
    
        end
        
        %% testGetSystemInput
        function testGetSystemInput(this)
            % by default, no input should be returned
            actualInput = this.plantUnderTest.getSystemInput();
            this.verifyEmpty(actualInput);
            
            expectedInput = this.sysInput;
            this.plantUnderTest.setSystemInput(expectedInput);
            actualInput = this.plantUnderTest.getSystemInput();
            this.verifyEqual(actualInput, expectedInput);
        end
        
        %% testSetSystemInputMatrixInvalidMatrices
        function testSetSystemInputMatrixInvalidMatrices(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrixDims';
             
            % input matrix of wrong dimensions
            this.verifyError(@() this.plantUnderTest.setSystemInputMatrix(this.invalidDimsBMatrix), expectedErrId);
            % invalid input matrix
            this.verifyError(@() this.plantUnderTest.setSystemInputMatrix(this.invalidBMatrix), expectedErrId);
            
            % change of input dimension is not allowed
            this.verifyError(@() this.plantUnderTest.setSystemInputMatrix([this.B, this.B]), expectedErrId);
        end
        
        %% testSetSystemInputMatrixNoInput
        function testSetSystemInputMatrixNoInput(this)
            this.assertEmpty(this.plantUnderTest.getSystemInput());
            expectedB = this.B * 2;
            
            this.plantUnderTest.setSystemInputMatrix(expectedB);
            this.verifyEqual(this.plantUnderTest.inputMatrix, expectedB);
            
            % check the side effect
            plantState = ones(this.dimX, 1);
            noiseSample = zeros(this.dimX, 1); % for simplicity
            expectedState = this.A * plantState;
            actualState = this.plantUnderTest.systemEquation(plantState, noiseSample);
            this.verifyEqual(actualState, expectedState);
        end
        
        
        %% testSetSystemInputMatrixInput
        function testSetSystemInputMatrixInput(this)
            this.assertEmpty(this.plantUnderTest.getSystemInput());
                        
            expectedB = this.B * 2;
            
            this.plantUnderTest.setSystemInputMatrix(expectedB);
            this.plantUnderTest.setSystemInput(this.sysInput);
            
            this.verifyEqual(this.plantUnderTest.inputMatrix, expectedB);
            
            % check the side effect
            plantState = ones(this.dimX, 1);
            noiseSample = zeros(this.dimX, 1); % for simplicity
            expectedState = this.A * plantState + expectedB * this.sysInput;
            actualState = this.plantUnderTest.systemEquation(plantState, noiseSample);
            this.verifyEqual(actualState, expectedState);
                        
            this.assertEqual(this.plantUnderTest.getSystemInput(), this.sysInput);
            
            % now change the input matrix again, this time an input is
            % already given
            expectedB = this.B * 3;
            this.plantUnderTest.setSystemInputMatrix(expectedB);
            % check the side effect
            plantState = ones(this.dimX, 1);
            noiseSample = zeros(this.dimX, 1); % for simplicity
            expectedState = this.A * plantState + expectedB * this.sysInput;
            actualState = this.plantUnderTest.systemEquation(plantState, noiseSample);
            this.verifyEqual(actualState, expectedState);
        end
        
        %% testSetSystemMatrixInvalidMatrices
        function testSetSystemMatrixInvalidMatrices(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidDimensions';
            
            % matrix of wrong dimensions, not square
            this.verifyError(@() this.plantUnderTest.setSystemMatrix(this.invalidDimsAMatrix), expectedErrId);
            % invalid matrix
            this.verifyError(@() this.plantUnderTest.setSystemMatrix(this.invalidAMatrix), expectedErrId);
            
            % change of state dimension is not allowed
            this.verifyError(@() this.plantUnderTest.setSystemMatrix(blkdiag(this.A, 2)), expectedErrId);
        end
        
        %% testSetSystemMatrix
        function testSetSystemMatrix(this)
            expectedA = this.A * 2;
            
            this.plantUnderTest.setSystemMatrix(expectedA);
            this.verifyEqual(this.plantUnderTest.sysMatrix, expectedA);
            this.verifyEqual(this.plantUnderTest.dimState, this.dimX);            
        end
        
        %% testSetStateConstraintsInvalidLowerBound
        function testSetStateConstraintsInvalidLowerBound(this)
            expectedErrId = 'LinearPlant:SetStateConstraints:InvalidLowerBound';
            
            upperBound = inf(this.dimX, 1);
            invalidLowerBound = this; % not a vector
            this.verifyError(@() this.plantUnderTest.setStateConstraints(invalidLowerBound, upperBound), ...
                expectedErrId);
            
            invalidLowerBound = -inf(this.dimX +1, 1); %invalid dim
            this.verifyError(@() this.plantUnderTest.setStateConstraints(invalidLowerBound, upperBound), ...
                expectedErrId);            
            
            invalidLowerBound = -inf(this.dimX +1, 1); %contains nan
            invalidLowerBound(end-1) = nan;
            this.verifyError(@() this.plantUnderTest.setStateConstraints(invalidLowerBound, upperBound), ...
                expectedErrId);
        end
        
        %% testSetStateConstraintsInvalidUpperBound
        function testSetStateConstraintsInvalidUpperBound(this)
            expectedErrId = 'LinearPlant:SetStateConstraints:InvalidUpperBound';
            
            lowerBound = inf(this.dimX, 1);
            invalidUpperBound = this; % not a vector
            this.verifyError(@() this.plantUnderTest.setStateConstraints(lowerBound, invalidUpperBound), ...
                expectedErrId);
            
            invalidUpperBound = inf(this.dimX +1, 1); %invalid dim
            this.verifyError(@() this.plantUnderTest.setStateConstraints(lowerBound, invalidUpperBound), ...
                expectedErrId);            
            
            invalidUpperBound = inf(this.dimX +1, 1); %contains nan
            invalidUpperBound(end-1) = -nan;
            this.verifyError(@() this.plantUnderTest.setStateConstraints(lowerBound, invalidUpperBound), ...
                expectedErrId);
        end
        
        %% setStateConstraints
        function setStateConstraints(this)
            lowerBound = ones(this.dimX, 1);
            lowerBound(end) = -inf;
            upperBound = 10 * ones(1, this.dimX); % use row vector here
            upperBound(end-1) = inf;
            
            this.plantUnderTest.setStateConstraints(lowerBound, upperBound);
            this.verifyEqual(this.plantUnderTest.stateConstraints, [lowerBound upperBound(:)]);
        end
        
        %% testIsValidStateError
        function testIsValidStateError(this)
            expectedErrId = 'LinearPlant:IsValidState:InvalidSystemState';
            
            errState = this; % not a vector
            this.verifyError(@() this.plantUnderTest.isValidState(errState), expectedErrId);
            
            errState = ones(this.dimX-1, 1); % wrong dimension
            this.verifyError(@() this.plantUnderTest.isValidState(errState), expectedErrId);
            
            errState = ones(1, this.dimX); % not a col vector
            this.verifyError(@() this.plantUnderTest.isValidState(errState), expectedErrId);
            
            errState = ones(this.dimX, 1);
            errState(end-1) = -inf; % not finite
            this.verifyError(@() this.plantUnderTest.isValidState(errState), expectedErrId);
        end
        
        %% testIsValidStateNoConstraints
        function testIsValidStateNoConstraints(this)
            this.assertEmpty(this.plantUnderTest.stateConstraints);
            
            validState = ones(this.dimX, 1);
            this.verifyTrue(this.plantUnderTest.isValidState(validState));
        end
        
        %% testIsValidStateOnlyLowerBound
        function testIsValidStateOnlyLowerBound(this)
            lowerBound = ones(this.dimX, 1);
            upperBound = inf(1, this.dimX);
            this.plantUnderTest.setStateConstraints(lowerBound, upperBound);
            
            this.assertNotEmpty(this.plantUnderTest.stateConstraints);
            
            validState = 42 * ones(this.dimX, 1);
            invalidState = -validState;
            this.verifyTrue(this.plantUnderTest.isValidState(validState));
            this.verifyFalse(this.plantUnderTest.isValidState(invalidState));
        end
        
        %% testIsValidStateOnlyUpperBound
        function testIsValidStateOnlyUpperBound(this)
            lowerBound = -inf(this.dimX, 1);
            upperBound = ones(1, this.dimX);
            this.plantUnderTest.setStateConstraints(lowerBound, upperBound);
            
            this.assertNotEmpty(this.plantUnderTest.stateConstraints);
            
            validState = -42 * ones(this.dimX, 1);
            invalidState = -validState;
            this.verifyTrue(this.plantUnderTest.isValidState(validState));
            this.verifyFalse(this.plantUnderTest.isValidState(invalidState));
        end
        
        %% testIsValidState
        function testIsValidState(this)
            lowerBound = ones(this.dimX, 1);
            lowerBound(end) = -inf;
            upperBound = 10 * ones(1, this.dimX); % use row vector here
            upperBound(end-1) = inf;
            this.plantUnderTest.setStateConstraints(lowerBound, upperBound);
            
            this.assertNotEmpty(this.plantUnderTest.stateConstraints);
            
            validState = ones(this.dimX, 1);
            validState2 = ones(this.dimX, 1);
            validState2(end-1) = 1e52;
            invalidState = -validState;
            this.verifyTrue(this.plantUnderTest.isValidState(validState));
            this.verifyTrue(this.plantUnderTest.isValidState(validState2));
            this.verifyFalse(this.plantUnderTest.isValidState(invalidState));
        end
        
        %% testCopy
        function testCopy(this)            
            lowerBound = ones(this.dimX, 1);
            lowerBound(end) = -inf;
            upperBound = 10 * ones(1, this.dimX); % use row vector here
            upperBound(end-1) = inf;
            this.plantUnderTest.setStateConstraints(lowerBound, upperBound);
            this.plantUnderTest.setSystemInput(this.sysInput);
            
            plantCopy = this.plantUnderTest.copy();
            this.verifyEqual(this.plantUnderTest, plantCopy);
            
            plantCopy.setSystemMatrix(this.A / 2);
            this.verifyNotEqual(this.plantUnderTest, plantCopy);
            this.verifyEqual(this.plantUnderTest.sysMatrix, this.A);
            
            plantCopy.setSystemMatrix(this.A);
            this.verifyEqual(this.plantUnderTest, plantCopy); % equal again
            
            % now change the system input matrix
            plantCopy.setSystemInputMatrix(this.B / 2);
            this.verifyNotEqual(this.plantUnderTest, plantCopy);
            this.verifyEqual(this.plantUnderTest.inputMatrix, this.B);
            
            plantCopy.setSystemInputMatrix(this.B);
            this.verifyEqual(this.plantUnderTest, plantCopy); % equal again
                        
            % now change the system noise matrix
            plantCopy.setSystemNoiseMatrix(this.A);
            this.verifyNotEqual(this.plantUnderTest, plantCopy);
            this.verifyEmpty(this.plantUnderTest.sysNoiseMatrix);
            
            plantCopy.setSystemNoiseMatrix([]);
            this.verifyEqual(this.plantUnderTest, plantCopy); % equal again
            
            % now change the input
            plantCopy.setSystemInput([]);
            this.verifyNotEqual(this.plantUnderTest, plantCopy);
            this.verifyEmpty(plantCopy.getSystemInput());
            
            plantCopy.setSystemInput(this.sysInput);
            this.verifyEqual(this.plantUnderTest, plantCopy); % equal again
            
            this.plantUnderTest.setStateConstraints(ones(this.dimX, 1), 10 * ones(1, this.dimX)); % slightly different
            this.verifyNotEqual(this.plantUnderTest, plantCopy);
            this.verifyEqual(plantCopy.stateConstraints, [lowerBound(:), upperBound(:)]);
            
            this.plantUnderTest.setStateConstraints(lowerBound, upperBound);
            this.verifyEqual(this.plantUnderTest, plantCopy); % equal again
            
            % finally, change the noise
            plantCopy.setNoise(this.W * 1.5); % change noise covariance
            [~,newCov] = plantCopy.noise.getMeanAndCov();
            [~, oldCov] = this.plantUnderTest.noise.getMeanAndCov();
            this.verifyEqual(newCov, this.W * 1.5);
            this.verifyNotEqual(oldCov, newCov);            
            this.verifyNotEqual(this.plantUnderTest, plantCopy);
            
            plantCopy.setNoise(this.W);
            this.verifyEqual(this.plantUnderTest, plantCopy); % equal again
        end
    end
    
end

