classdef LinearPlantTest < matlab.unittest.TestCase
    % Test cases for LinearPlant.
    
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
        
        %% testSetSystemInputMatrix(this)
        function testSetSystemInputMatrixInvalidMatrices(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
             
            % input matrix of wrong dimensions
            this.verifyError(@() this.plantUnderTest.setSystemInputMatrix(this.invalidDimsBMatrix), expectedErrId);
            % invalid input matrix
            this.verifyError(@() this.plantUnderTest.setSystemInputMatrix(this.invalidBMatrix), expectedErrId);
        end
        
        %% testSetSystemInputMatrix
        function testSetSystemInputMatrix(this)
            expectedB = this.B * 2;
            
            this.plantUnderTest.setSystemInputMatrix(expectedB);
            this.verifyEqual(this.plantUnderTest.inputMatrix, expectedB);
        end
        
        %% testSetSystemMatrixInvalidMatrices
        function testSetSystemMatrixInvalidMatrices(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
            
            % matrix of wrong dimensions, not square
            this.verifyError(@() this.plantUnderTest.setSystemMatrix(this.invalidDimsAMatrix), expectedErrId);
            % invalid matrix
            this.verifyError(@() this.plantUnderTest.setSystemMatrix(this.invalidAMatrix), expectedErrId);
        end
        
        %% testSetSystemMatrix
        function testSetSystemMatrix(this)
            expectedA = this.A * 2;
            
            this.plantUnderTest.setSystemMatrix(expectedA);
            this.verifyEqual(this.plantUnderTest.sysMatrix, expectedA);
            this.verifyEqual(this.plantUnderTest.dimState, this.dimX);
            
            % change the dimension of the state
            expectedA = eye(this.dimX * 3);
            this.plantUnderTest.setSystemMatrix(expectedA);
            this.verifyEqual(this.plantUnderTest.sysMatrix, expectedA);
            this.verifyEqual(this.plantUnderTest.dimState, this.dimX * 3);
        end
    end
    
end

