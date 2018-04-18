classdef DelayedKFSystemModelTest < matlab.unittest.TestCase
    % Test cases for DelayedKFSystemModel.
    
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
        zeroMeasDelay = 0;
    end
    
    properties (Access = private)
        dimX;
        dimU;
        
        A;
        B;
        W;
        plantNoise;
        numPossibleInputs;
        maxMeasDelay;
        inputProbs;
        uncertainInputs;
        invalidUncertainInputs;
        nanUncertainInputs;
        
        modelUnderTest;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.numPossibleInputs = 5;
            this.maxMeasDelay = 5;
            this.dimX = 2;
            this.dimU = 1;
            
            this.A = [1 1; 0 1];
            this.B = [0; 1];
            this.W = 0.5 * eye(this.dimX);
            this.plantNoise = Gaussian(zeros(this.dimX, 1), this.W);
            this.inputProbs = ones(1, this.numPossibleInputs) / this.numPossibleInputs;
            
            this.uncertainInputs = ones(this.dimU, this.numPossibleInputs) .* [1:this.numPossibleInputs];
            this.invalidUncertainInputs = eye(this.dimX);
            this.nanUncertainInputs = this.uncertainInputs;
            this.nanUncertainInputs(1) = nan;
            
            this.modelUnderTest = DelayedKFSystemModel(this.A, this.B, this.plantNoise, ...
                this.numPossibleInputs, this.maxMeasDelay, this.inputProbs);
          end
    end
    
    methods (Test)
        %% testDelayedKFSystemModelInvalidMaxMeasDelay
        function testDelayedKFSystemModelInvalidMaxMeasDelay(this)
            expectedErrId = 'DelayedKFSystemModel:InvalidMaxMeasDelay';
            
            invalidMaxMeasDelay = -this.maxMeasDelay; % negative value
            this.verifyError(@() DelayedKFSystemModel(this.A, this.B, this.plantNoise, this.numPossibleInputs, ...
                invalidMaxMeasDelay, this.inputProbs), expectedErrId);
            
            invalidMaxMeasDelay = this.maxMeasDelay / 2; % fractional
            this.verifyError(@() DelayedKFSystemModel(this.A, this.B, this.plantNoise, this.numPossibleInputs, ...
                invalidMaxMeasDelay, this.inputProbs), expectedErrId);
        end
        
        %% testDelayedKFSystemModelInvalidNumModes
        function testDelayedKFSystemModelInvalidNumModes(this)
            expectedErrId = 'DelayedKFSystemModel:InvalidNumModes';
            
            invalidNumPossibleInputs = 0; % not a positive value
            this.verifyError(@() DelayedKFSystemModel(this.A, this.B, this.plantNoise, invalidNumPossibleInputs, ...
                this.maxMeasDelay, this.inputProbs), expectedErrId);
            
            invalidNumPossibleInputs = this.maxMeasDelay / 2; % fractional
            this.verifyError(@() DelayedKFSystemModel(this.A, this.B, this.plantNoise, invalidNumPossibleInputs, ...
                this.maxMeasDelay, this.inputProbs), expectedErrId);
            
        end
        
        %% testDelayedKFSystemModelInvalidNoise
        function testDelayedKFSystemModelInvalidNoise(this)
            expectedErrId = 'DelayedKFSystemModel:InvalidSystemNoise';
            
            invalidNoise = this.W; % directly use the noise covariance
            this.verifyError(@() DelayedKFSystemModel(this.A, this.B, invalidNoise, this.numPossibleInputs, ...
                this.maxMeasDelay, this.inputProbs), expectedErrId);
        end
        
        %% testDelayedKFSystemModel
        function testDelayedKFSystemModel(this)
            expectedAugSysMatrix = [this.A, zeros(this.dimX, this.dimX * this.maxMeasDelay); ...
                Utils.blockDiag(eye(this.dimX), this.maxMeasDelay), zeros(this.dimX * this.maxMeasDelay, this.dimX)];
            expectedAugNoiseMatrix = [eye(this.dimX); zeros(this.dimX * this.maxMeasDelay, this.dimX)];
            
            model = DelayedKFSystemModel(this.A, this.B, this.plantNoise, ...
                this.numPossibleInputs , this.maxMeasDelay, this.inputProbs);
            
            % check the constructed augmented matrices
            augSysMatrix = model.sysMatrix;
            augNoiseMatrix = model.sysNoiseMatrix;
            
            this.verifyEqual(augSysMatrix, expectedAugSysMatrix);
            this.verifyEqual(augNoiseMatrix, expectedAugNoiseMatrix);
        end
        
        %% testDelayedKFSystemModel
        function testDelayedKFSystemModelZeroMeasDelay(this)
            % in case of zero measurement delay, this model reduces to a
            % simple, linear model
            expectedAugSysMatrix = this.A;
            expectedAugNoiseMatrix = eye(this.dimX);
            
            zeroMeasDelayModel = DelayedKFSystemModel(this.A, this.B, this.plantNoise, ...
                this.numPossibleInputs, DelayedKFSystemModelTest.zeroMeasDelay, this.inputProbs);
            
            % check the constructed augmented matrices: remove sparsity
            augSysMatrix = full(zeroMeasDelayModel.sysMatrix);
            augNoiseMatrix = full(zeroMeasDelayModel.sysNoiseMatrix);
            
            this.verifyEqual(augSysMatrix, expectedAugSysMatrix);
            this.verifyEqual(augNoiseMatrix, expectedAugNoiseMatrix);
        end
        
        %% testReset
        function testReset(this)
            % set uncertain inputs, then reset the model
            expectedPlantNoise = this.plantNoise;
            [expectedMean, expectedCov] = this.plantNoise.getMeanAndCov();
            
            this.uncertainInputs = ones(this.dimU, this.numPossibleInputs) .* [1:this.numPossibleInputs];
            this.modelUnderTest.setSystemInput(this.uncertainInputs);
            
            this.verifyNotEqual(this.modelUnderTest.noise, expectedPlantNoise);
            
            this.modelUnderTest.reset();
            
            actualPlantNoise = this.modelUnderTest.noise;
            [actualMean, actualCov] = actualPlantNoise.getMeanAndCov();
            
            this.verifyClass(actualPlantNoise, ?Gaussian);
            this.verifyEqual(actualMean, expectedMean);
            this.verifyEqual(actualCov, expectedCov);
        end
        
        %% testSetSystemMatrixInvalidMatrix
        function testSetSystemMatrixInvalidMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidDimensions';
            
            newSysMatrix = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() this.modelUnderTest.setSystemMatrix(newSysMatrix), expectedErrId);
            newSysMatrix = inf(this.dimX, this.dimX); % square, but inf
            this.verifyError(@() this.modelUnderTest.setSystemMatrix(newSysMatrix), expectedErrId);
        end
        
        %% testSetSystemMatrix
        function testSetSystemMatrix(this)
            expectedNewSysMatrix = this.A' * this.A;
            augmentedSysMatrix = full(this.modelUnderTest.sysMatrix);
            
            this.modelUnderTest.setSystemMatrix(expectedNewSysMatrix);
            
            actualAugmentedNewSysMatrix = full(this.modelUnderTest.sysMatrix);
            actualNewSysMatrix = actualAugmentedNewSysMatrix(1:this.dimX, 1:this.dimX);
            
            this.verifyEqual(actualNewSysMatrix, expectedNewSysMatrix);
            % compare the remaining blocks
            this.verifyEqual(size(actualAugmentedNewSysMatrix), size(augmentedSysMatrix));
            this.verifyEqual(actualAugmentedNewSysMatrix(this.dimX+1:end, :), ...
                augmentedSysMatrix(this.dimX+1:end, :));
            this.verifyEqual(actualAugmentedNewSysMatrix(:, this.dimX+1:end), ...
                augmentedSysMatrix(:, this.dimX+1:end));
        end
        
        %% testSetSystemMatrixEmptyMatrix
        function testSetSystemMatrixEmptyMatrix(this)
            expectedNewSysMatrix = eye(this.dimX);
            augmentedSysMatrix = full(this.modelUnderTest.sysMatrix);
            
            this.modelUnderTest.setSystemMatrix([]);
            
            actualAugmentedNewSysMatrix = full(this.modelUnderTest.sysMatrix);
            actualNewSysMatrix = actualAugmentedNewSysMatrix(1:this.dimX, 1:this.dimX);
            
            this.verifyEqual(actualNewSysMatrix, expectedNewSysMatrix);
            % compare the remaining blocks
            this.verifyEqual(size(actualAugmentedNewSysMatrix), size(augmentedSysMatrix));
            this.verifyEqual(actualAugmentedNewSysMatrix(this.dimX+1:end, :), ...
                augmentedSysMatrix(this.dimX+1:end, :));
            this.verifyEqual(actualAugmentedNewSysMatrix(:, this.dimX+1:end), ...
                augmentedSysMatrix(:, this.dimX+1:end));
        end
        
        %% testSetSystemNoiseMatrixInvalidMatrix
        function testSetSystemNoiseMatrixInvalidMatrix(this)
            expectedErrId = 'DelayedKFSystemModel:SetSystemNoiseMatrix:InvalidDimensions';
            
            newSysNoiseMatrix = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() this.modelUnderTest.setSystemNoiseMatrix(newSysNoiseMatrix), expectedErrId);
            newSysNoiseMatrix = inf(this.dimX, this.dimX); % square, but inf
            this.verifyError(@() this.modelUnderTest.setSystemNoiseMatrix(newSysNoiseMatrix), expectedErrId);
        end
        
        %% testSetSystemNoiseMatrix
        function testSetSystemNoiseMatrix(this)
            expectedNewSysNoiseMatrix = this.A' * this.A;
            augmentedSysNoiseMatrix = full(this.modelUnderTest.sysNoiseMatrix);
            
            this.modelUnderTest.setSystemNoiseMatrix(expectedNewSysNoiseMatrix);
            
            actualAugmentedNewSysNoiseMatrix = full(this.modelUnderTest.sysNoiseMatrix);
            actualNewSysNoiseMatrix = actualAugmentedNewSysNoiseMatrix(1:this.dimX, :);
            
            this.verifyEqual(actualNewSysNoiseMatrix, expectedNewSysNoiseMatrix);
            % compare the remaining blocks
            this.verifyEqual(size(actualAugmentedNewSysNoiseMatrix), size(actualAugmentedNewSysNoiseMatrix));
            this.verifyEqual(actualAugmentedNewSysNoiseMatrix(this.dimX+1:end, :), ...
                augmentedSysNoiseMatrix(this.dimX+1:end, :));
            this.verifyEqual(actualAugmentedNewSysNoiseMatrix(:, this.dimX+1:end), ...
                augmentedSysNoiseMatrix(:, this.dimX+1:end));
        end
        
        %% testSetSystemNoiseMatrixEmptyMatrix
        function testSetSystemNoiseMatrixEmptyMatrix(this)
            expectedNewSysNoiseMatrix = eye(this.dimX);
            augmentedSysNoiseMatrix = full(this.modelUnderTest.sysNoiseMatrix);
            
            this.modelUnderTest.setSystemNoiseMatrix([]);
            
            actualAugmentedNewSysNoiseMatrix = full(this.modelUnderTest.sysNoiseMatrix);
            actualNewSysNoiseMatrix = actualAugmentedNewSysNoiseMatrix(1:this.dimX, :);
            
            this.verifyEqual(actualNewSysNoiseMatrix, expectedNewSysNoiseMatrix);
            % compare the remaining blocks
            this.verifyEqual(size(actualAugmentedNewSysNoiseMatrix), size(actualAugmentedNewSysNoiseMatrix));
            this.verifyEqual(actualAugmentedNewSysNoiseMatrix(this.dimX+1:end, :), ...
                augmentedSysNoiseMatrix(this.dimX+1:end, :));
            this.verifyEqual(actualAugmentedNewSysNoiseMatrix(:, this.dimX+1:end), ...
                augmentedSysNoiseMatrix(:, this.dimX+1:end));
        end
        
        %% testSetSystemInputInvalidInput
        function testSetSystemInputInvalidInput(this)
            expectedErrId = 'DelayedKFSystemModel:InvalidInput';
            
            % inputs as matrix of wrong dimension
            invalidInputs = this.invalidUncertainInputs;
            this.verifyError(@() this.modelUnderTest.setSystemInput(invalidInputs), expectedErrId);
            % inputs as matrix of correct dimension, but with nan
            invalidInputs = this.nanUncertainInputs;
            this.verifyError(@() this.modelUnderTest.setSystemInput(invalidInputs), expectedErrId);
        end
        
        %% testSetSystemInput
        function testSetSystemInput(this)
            % compute the additional uncertainty
            inputMean = mean(this.uncertainInputs);
            inputCov = cov(this.uncertainInputs, 1); % normalize by number of samples
            expectedMean = this.B * inputMean;
            expectedCov = this.W + this.B * inputCov * this.B';
                  
            this.modelUnderTest.setSystemInput(this.uncertainInputs);
            
            actualPlantNoise = this.modelUnderTest.noise;
            [actualMean, actualCov] = actualPlantNoise.getMeanAndCov();
            
            this.verifyClass(actualPlantNoise, ?Gaussian);
            this.verifyEqual(actualMean, expectedMean);
            this.verifyEqual(actualCov, expectedCov);
        end
    end
    
end

