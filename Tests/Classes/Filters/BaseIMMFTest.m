classdef (Abstract) BaseIMMFTest < matlab.unittest.TestCase
    % This class contains test cases and aditional functions to facilitate
    % testing of IMM-like filters.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Constant, Access = protected)
        absTol = 1e-8;
    end
    
    properties (SetAccess = private, GetAccess = protected)
        numModes;
        modeTransitionMatrix;
        dimX;
        dimU;
        dimY;
        
        A1;
        A2;
        B;
        C;
        W;
        V;
        input;
        measurement;
        
        modePlantModels;
        jumpLinearPlantModel;
        
        measModel;
        
        modeFilters;
        
        stateGaussian;
        stateGaussianMean;
        stateGaussianCov;
        
        stateGaussianMixture;
        stateGaussianMixtureMean;
        stateGaussianMixtureCov;
        mixtureWeights;
        mixtureMeans;
        mixtureCovs;
        
        predictedMixtureWeights;
        predictedMixtureMeans;
        predictedMixtureCovs;
        
        stateInvalidGaussianMixture;
        
        filterUnderTest;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.numModes = 2;
            this.modeTransitionMatrix = ones(this.numModes, this.numModes) / this.numModes;
            
            this.modeFilters = arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, ...
                'UniformOutput', false);
                        
            this.dimX = 3;
            this.dimU = 2;
            this.dimY = 2;
            
            this.A1 = eye(this.dimX);
            this.A2 = eye(this.dimX) * 2;
            this.B = ones(this.dimX, this.dimU);
            this.C = ones(this.dimY, this.dimX); % meas matrix
            this.input = ones(this.dimU, 1);
            this.measurement = [1/2; 1];
            this.W = 0.01 * eye(this.dimX); % sys noise cov
            this.V = 0.1 * eye(this.dimY); % meas noise cov
            
            this.modePlantModels = {    LinearPlant(this.A1, this.B, this.W) ,...
                                        LinearPlant(this.A2, this.B, this.W)
                                    };
                                   
            this.jumpLinearPlantModel = JumpLinearSystemModel(this.numModes, this.modePlantModels);
            this.jumpLinearPlantModel.setSystemInput(this.input);
            
            this.measModel = LinearMeasurementModel(this.C);
            this.measModel.setNoise(Gaussian(zeros(this.dimY, 1), this.V));
            
            this.initDistributions();
            this.computePredictionGaussianMixture();
            
            this.filterUnderTest = this.initFilterUnderTest();
        end
    end
    
    methods (Access = protected, Abstract)
        filter = initFilterUnderTest(this);
    end
    
    methods (Access = protected)
         
        %% verifyEqualWithAbsTol
        function verifyEqualWithAbsTol(this, actual, expected)
            this.verifyEqual(actual, expected, 'AbsTol', BaseIMMFTest.absTol);
        end
        
        %% verifyPredictionGaussianMixture
        function verifyPredictionGaussianMixture(this)
            expectedWeights = this.predictedMixtureWeights;
            expectedMean = this.predictedMixtureWeights(1) * this.predictedMixtureMeans(:, 1) ...
                + this.predictedMixtureWeights(2) * this.predictedMixtureMeans(:, 2);
            expectedCov = this.predictedMixtureWeights(1) * this.predictedMixtureCovs(:, :, 1) ...
                + this.predictedMixtureWeights(2) * this.predictedMixtureCovs(:, :, 2);
            expectedCov = expectedCov + this.predictedMixtureWeights(1) ...
                * (expectedMean - this.predictedMixtureMeans(:, 1)) * transpose(expectedMean - this.predictedMixtureMeans(:, 1));
            expectedCov = expectedCov ...
                + (this.predictedMixtureWeights(2) * ...
                + (expectedMean - this.predictedMixtureMeans(:, 2)) * transpose(expectedMean - this.predictedMixtureMeans(:, 2)));
            
            % first check the individual filter states
            actualModeFilterStates = cellfun(@getState, this.modeFilters, 'UniformOutput', false);
            for j=1:this.numModes
                this.verifyClass(actualModeFilterStates{j}, ?Gaussian);
                % if suffices to compare mean and cov
                [actualMean, actualCov] = actualModeFilterStates{j}.getMeanAndCov();
                this.verifyEqualWithAbsTol(actualMean, this.predictedMixtureMeans(:, j));
                this.verifyEqualWithAbsTol(actualCov, this.predictedMixtureCovs(:, :, j));
            end
            
             % now the resulting Gaussian mixture
            predictedState = this.filterUnderTest.getState();
            this.verifyClass(predictedState, ?GaussianMixture);
            
            [~, ~, actualWeights] = predictedState.getComponents();
            [actualMean, actualCov] = this.filterUnderTest.getStateMeanAndCov();
                   
            % first, a sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualCov - this.stateGaussianMixtureCov), 0);
            
            this.verifyEqual(actualWeights, expectedWeights);
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqual(actualCov, actualCov'); % shall be symmetric
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
        end
        
        %% verifyPredictionGaussian
        function verifyPredictionGaussian(this)
            % init the filter with a Gaussian, and use same plant model for
            % all modes; equivalent to simply using a KF
            
            expectedWeights = this.predictedMixtureWeights;
            expectedMean = this.A1 * this.stateGaussianMean + this.B * this.input;
            expectedCov = this.A1 * this.stateGaussianCov * this.A1' + this.W;
            
            % first check the individual filter states; should be all the
            % same
            actualModeFilterStates = cellfun(@getState, this.modeFilters, 'UniformOutput', false);
            for j=1:this.numModes
                this.verifyClass(actualModeFilterStates{j}, ?Gaussian);
                % if suffices to compare mean and cov
                [actualMean, actualCov] = actualModeFilterStates{j}.getMeanAndCov();
                this.verifyEqualWithAbsTol(actualMean, expectedMean);
                this.verifyEqualWithAbsTol(actualCov, expectedCov);
            end
            
            % now the resulting Gaussian mixture
            predictedState = this.filterUnderTest.getState();
            this.verifyClass(predictedState, ?GaussianMixture);
            
            [~, ~, actualWeights] = predictedState.getComponents();
            [actualMean, actualCov] = this.filterUnderTest.getStateMeanAndCov();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualCov - this.stateGaussianCov), 0);
            
            % resulting Gaussian mixture consists of two equally weighted
            % Gaussians with same mean and cov
            this.verifyEqual(actualWeights, expectedWeights);
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqual(actualCov, actualCov'); % shall be symmetric
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
        end
    end
    
    methods (Access = private)
        %% initDistributions
        function initDistributions(this)
            this.stateGaussianMean = zeros(this.dimX, 1);
            this.stateGaussianCov = 0.5 * eye(this.dimX);
            this.stateGaussian = Gaussian(this.stateGaussianMean, this.stateGaussianCov);
            
            this.mixtureMeans = [zeros(this.dimX, 1) ones(this.dimX, 1)];
            this.mixtureCovs(:, :, 2) = 2 * eye(this.dimX);
            this.mixtureCovs(:, :, 1) = 0.5 * eye(this.dimX);
            this.mixtureWeights = [2/3 1/3];
            this.stateGaussianMixture = GaussianMixture(this.mixtureMeans, this.mixtureCovs, this.mixtureWeights);
            this.stateGaussianMixtureMean = 1/3 * ones(this.dimX, 1);
            this.stateGaussianMixtureCov = eye(this.dimX) + 2/9 * ones(this.dimX); % cov of means + mean of covs
                 
            invalidMixtureMeans = [this.mixtureMeans 0.75 * ones(this.dimX, 1)];
            invalidMixtureCovs(:, :, 3) = eye(this.dimX);
            invalidMixtureCovs(:, :, 2) = 2 * eye(this.dimX);
            invalidMixtureCovs(:, :, 1) = 0.5 * eye(this.dimX);
            invalidMixtureWeights = [1/3 1/2 1/6];
            this.stateInvalidGaussianMixture = ...
                GaussianMixture(invalidMixtureMeans, invalidMixtureCovs, invalidMixtureWeights);
            
        end
        
        %% computePredictionGaussianMixture
        function computePredictionGaussianMixture(this)
            % use straightforward implementation of the algorithm given
            % in (Table II) in
            %
            %   X. Rong Li and Vesselin P. Jilkov,
            %   Survey of maneuvering target tracking - Part V: Multiple-Model Methods,
            %   IEEE Transactions on Aerospace and Electronic Systems 41.4 (2005): 1255-1321.
            
            this.predictedMixtureWeights = [1/2 1/2]; 
            
            mixingWeights = [   2/3 2/3;
                                1/3 1/3];
            
            mixingEstimates(:, 2) = this.mixtureMeans(:, 1) * mixingWeights(1,2) ...
                + this.mixtureMeans(:, 2) * mixingWeights(2,2);
            mixingEstimates(:, 1) = this.mixtureMeans(:, 1) * mixingWeights(1,1) ...
                + this.mixtureMeans(:, 2) * mixingWeights(2,1);
            
            
            this.predictedMixtureMeans(:, 2) = this.A2 * mixingEstimates(:, 2) + this.B * this.input;
            this.predictedMixtureMeans(:, 1) = this.A1 * mixingEstimates(:, 1) + this.B * this.input;
            
            mixingCovs(:, :, 2) = mixingWeights(1,2) * (this.mixtureCovs(:, :, 1) ...
                + (mixingEstimates(:, 2) - this.mixtureMeans(:, 1)) * transpose(mixingEstimates(:, 2) - this.mixtureMeans(:, 1))) ...
                + mixingWeights(2,2) * (this.mixtureCovs(:, :, 2) ...
                + (mixingEstimates(:, 2) - this.mixtureMeans(:, 2)) * transpose(mixingEstimates(:, 2) - this.mixtureMeans(:, 2)));
            
             mixingCovs(:, :, 1) = mixingWeights(1,1) * (this.mixtureCovs(:, :, 1) ...
                + (mixingEstimates(:, 1) - this.mixtureMeans(:, 1)) * transpose(mixingEstimates(:, 1) - this.mixtureMeans(:, 1))) ...
                + mixingWeights(2,1) * (this.mixtureCovs(:, :, 2) ...
                + (mixingEstimates(:, 1) - this.mixtureMeans(:, 2)) * transpose(mixingEstimates(:, 1) - this.mixtureMeans(:, 2)));
            
            this.predictedMixtureCovs(:, :, 2) = this.A2 * mixingCovs(:, :, 2) * this.A2' + this.W;
            this.predictedMixtureCovs(:, :, 1) = this.A1 * mixingCovs(:, :, 1) * this.A1' + this.W;
        end
    end
    
    methods (Test)
        %% testSetModeTransitionMatrixInvalidMatrix
        function testSetModeTransitionMatrixInvalidMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';
            
            % transition matrix must be square
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1];
            this.verifyError(@() this.filterUnderTest.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
            
            % transition matrix must be square and have appropriate
            % dimensions
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1; 0.5 0.25 0.25];
            this.verifyError(@() this.filterUnderTest.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
            
            % row sums must be 1
            invalidTransitionMatrix = [0.7 0.3; 0.6 0.5];
            this.verifyError(@() this.filterUnderTest.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
            
            % entries must be in [0,1]
            invalidTransitionMatrix = [1.3 -0.3; 0.6 0.4];
            this.verifyError(@() this.filterUnderTest.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
        end
              
                
        %% testSetModeTransitionMatrix
        function testSetModeTransitionMatrix(this)
            expectedNewTransitionMatrix = this.modeTransitionMatrix';
            
            this.filterUnderTest.setModeTransitionMatrix(expectedNewTransitionMatrix);
            
            this.verifyEqual(this.filterUnderTest.modeTransitionProbs, expectedNewTransitionMatrix);
        end
        
        %% testSetStateInvalidDistribution
        function testSetStateInvalidDistribution(this)
            expectedErrId = 'Filter:InvalidSystemState';
            
            invalidState = eye(this.dimX); % not a Distribution
            this.verifyError(@() this.filterUnderTest.setState(invalidState), expectedErrId);
        end
        
        %% testSetStateInvalidGaussianMixture
        function testSetStateInvalidGaussianMixture(this)
            expectedErrId = 'Filter:InvalidStateGaussianMixture';
            
            invalidState = this.stateInvalidGaussianMixture; % mixture with too many components
            this.verifyError(@() this.filterUnderTest.setState(invalidState), expectedErrId);
        end
        
        %% testSetStateDistribution
        function testSetStateDistribution(this)
            
            state = this.stateGaussian;
            this.filterUnderTest.setState(state);
            
            actualModeFilterStates = cellfun(@getState, this.modeFilters, 'UniformOutput', false);
            [expectedMean, expectedCov] = this.stateGaussian.getMeanAndCov();
            for j=1:this.numModes
                this.verifyClass(actualModeFilterStates{j}, ?Gaussian);
                % if suffices to compare mean and cov
                [actualMean, actualCov] = actualModeFilterStates{j}.getMeanAndCov();
                this.verifyEqual(actualMean, expectedMean);
                this.verifyEqual(actualCov, expectedCov);
            end
        end
        
        %% testSetStateGaussianMixture
        function testSetStateGaussianMixture(this)
           
            state = this.stateGaussianMixture;
            this.filterUnderTest.setState(state);
            
            actualModeFilterStates = cellfun(@getState, this.modeFilters, 'UniformOutput', false);
 
            for j=1:this.numModes
                this.verifyClass(actualModeFilterStates{j}, ?Gaussian);
                % if suffices to compare mean and cov
                [actualMean, actualCov] = actualModeFilterStates{j}.getMeanAndCov();
                this.verifyEqual(actualMean, this.mixtureMeans(:, j));
                this.verifyEqual(actualCov, this.mixtureCovs(:, :, j));
            end
        end
        
        %% testGetStateGaussianMixture
        function testGetStateGaussianMixture(this)
            expectedNewState = this.stateGaussianMixture;
            expectedMean = this.stateGaussianMixtureMean;
            expectedCov = this.stateGaussianMixtureCov;
            expectedWeights = this.mixtureWeights;
            
            this.filterUnderTest.setState(expectedNewState);
            actualNewState = this.filterUnderTest.getState();
            [actualMean, actualCov] = actualNewState.getMeanAndCov();
                        
            this.verifyClass(actualNewState, ?GaussianMixture);
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
            
            [~, ~, actualWeights] = actualNewState.getComponents();
            this.verifyEqual(actualWeights, expectedWeights);
        end
        
        %% testGetStateGaussian
        function testGetStateGaussian(this)
            expectedNewState = this.stateGaussian;
            expectedMean = this.stateGaussianMean;
            expectedCov = this.stateGaussianCov;
            % we expect equally weighted components, given as a row vector
            expectedWeights = [1/2 1/2];
            
            this.filterUnderTest.setState(expectedNewState);
            actualNewState = this.filterUnderTest.getState();
            [actualMean, actualCov] = actualNewState.getMeanAndCov();
                        
            this.verifyClass(actualNewState, ?GaussianMixture);
            this.verifyEqual(actualMean, expectedMean, 'AbsTol', 1e-8);
            this.verifyEqual(actualCov, expectedCov, 'AbsTol', 1e-8);
            
            [~, ~, actualWeights] = actualNewState.getComponents();
            this.verifyEqual(actualWeights, expectedWeights);
        end
        
        %% testGetStateMeanAndCovGaussianMixture
        function testGetStateMeanAndCovGaussianMixture(this)
            newState = this.stateGaussianMixture;
            expectedMean = this.stateGaussianMixtureMean;
            expectedCov = this.stateGaussianMixtureCov;
            expectedCovSqrt = chol(this.stateGaussianMixtureCov)';
            
            % set the new state and retrieve the point estimate
            this.filterUnderTest.setState(newState);
             [actualMean, actualCov, actualCovSqrt] = this.filterUnderTest.getStateMeanAndCov();
            
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
            this.verifyEqualWithAbsTol(actualCovSqrt, expectedCovSqrt);
            
        end
        
        %% testGetStateMeanAndCovGaussian
        function testGetStateMeanAndCovGaussian(this)
            newState = this.stateGaussian;
            expectedMean = this.stateGaussianMean;
            expectedCov = this.stateGaussianCov;
            expectedCovSqrt = chol(this.stateGaussianCov)';
            
            % set the new state and retrieve the point estimate
            this.filterUnderTest.setState(newState);
            [actualMean, actualCov, actualCovSqrt] = this.filterUnderTest.getStateMeanAndCov();
            
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
            this.verifyEqualWithAbsTol(actualCovSqrt, expectedCovSqrt);
        end
        
        %% testGetModeEstimateInitialMode
        function testGetModeEstimateInitialMode(this)
            % setState has not been called, so mode should be empty
            [maxMode, prob] = this.filterUnderTest.getModeEstimate();
            
            this.verifyEmpty(maxMode);
            this.verifyEmpty(prob);
        end
        
        %% testGetModeEstimateGaussianMixture
        function testGetModeEstimateGaussianMixture(this)
            newState = this.stateGaussianMixture;
            expectedMaxMode = 1;
            expectedMaxModeProb = 2 / 3;
            
            this.filterUnderTest.setState(newState);
            [maxMode, prob] = this.filterUnderTest.getModeEstimate();
            
            this.verifyEqual(maxMode, expectedMaxMode);
            this.verifyEqual(prob, expectedMaxModeProb);
        end
        
        %% testGetModeEstimateGaussian
        function testGetModeEstimateGaussian(this)
            newState = this.stateGaussian;
            % mode probability is uniform, first maxima should be returned
            expectedMaxMode = 1;
            expectedMaxModeProb = 1/ 2;
            
            this.filterUnderTest.setState(newState);
            [maxMode, prob] = this.filterUnderTest.getModeEstimate();
            
            this.verifyEqual(maxMode, expectedMaxMode);
            this.verifyEqual(prob, expectedMaxModeProb);
        end
    end
    
end

