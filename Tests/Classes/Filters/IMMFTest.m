classdef IMMFTest < BaseIMMFTest;
    % Test cases for IMMF.
    
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
    
    properties (Constant, Access = private)
        defaultFiltername = 'Interacting Multiple Model Filter';
    end
    
    properties (Access = private)
        filterName;
        
        updatedMixtureWeights;
        updatedMixtureMeans;
        updatedMixtureCovs;
        
    end
    
    methods (Access = private)
        %% computeUpdateGaussianMixture        
        function computeUpdateGaussianMixture(this)
            % Use straightforward implementation of the algorithm given
            % in (Table II) in
            %
            %   X. Rong Li and Vesselin P. Jilkov,
            %   A survey of maneuvering target tracking - Part V: Multiple-Model Methods,
            %   IEEE Transactions on Aerospace and Electronic Systems 41.4 (2005): 1255-1321.
            
            % we use a single measurement model for all modes
            % first the innovations/ residuals
            residuals(:, 2) = this.measurement - this.C * this.mixtureMeans(:, 2);
            residuals(:, 1) = this.measurement - this.C * this.mixtureMeans(:, 1);
            residualCovs(:, :, 2) = this.C * this.mixtureCovs(:, :, 2) * this.C' + this.V;
            residualCovs(:, :, 1) = this.C * this.mixtureCovs(:, :, 1) * this.C' + this.V;
            
            % compute the individual gains for all filters
            gains(:, :, 2) = this.mixtureCovs(:, :, 2) * this.C' * inv(residualCovs(:, :, 2));
            gains(:, :, 1) = this.mixtureCovs(:, :, 1) * this.C' * inv(residualCovs(:, :, 1));
             
            % now perform the update
            this.updatedMixtureMeans(:, 2) = this.mixtureMeans(:, 2) + gains(:, :, 2) * residuals(:, 2);
            this.updatedMixtureMeans(:, 1) = this.mixtureMeans(:, 1) + gains(:, :, 1) * residuals(:, 1);
            this.updatedMixtureCovs(:, :, 2) = this.mixtureCovs(:, :, 2) - gains(:, :, 2) * residualCovs(:, :, 2) * gains(:, :, 2)';
            this.updatedMixtureCovs(:, :, 1) = this.mixtureCovs(:, :, 1) - gains(:, :, 1) * residualCovs(:, :, 1) * gains(:, :, 1)';
            
            % for the weight updates we need the Gaussian-assumed likelihoods
            likelihoods(2) = mvnpdf(residuals(:, 2), zeros(this.dimY, 1), residualCovs(:, :, 2));
            likelihoods(1) = mvnpdf(residuals(:, 1), zeros(this.dimY, 1), residualCovs(:, :, 1));

            normalizationConstant = this.mixtureWeights(1) * likelihoods(1) + this.mixtureWeights(2) * likelihoods(2);
            this.updatedMixtureWeights(2) = this.mixtureWeights(2) * likelihoods(2) / normalizationConstant;
            this.updatedMixtureWeights(1) = this.mixtureWeights(1) * likelihoods(1) / normalizationConstant;
        end
                     
        %% verifyUpdateGaussianMixture
        function verifyUpdateGaussianMixture(this)
            % compute the expected parameters of the posterior mixture
            expectedWeights = this.updatedMixtureWeights;
            expectedMean = this.updatedMixtureWeights(1) * this.updatedMixtureMeans(:, 1) ...
                + this.updatedMixtureWeights(2) * this.updatedMixtureMeans(:, 2);
            expectedCov = this.updatedMixtureWeights(1) * this.updatedMixtureCovs(:, :, 1) ...
                + this.updatedMixtureWeights(2) * this.updatedMixtureCovs(:, :, 2);
            expectedCov = expectedCov + this.updatedMixtureWeights(1) ...
                * (expectedMean - this.updatedMixtureMeans(:, 1)) * transpose(expectedMean - this.updatedMixtureMeans(:, 1));
            expectedCov = expectedCov ...
                + (this.updatedMixtureWeights(2) * ...
                + (expectedMean - this.updatedMixtureMeans(:, 2)) * transpose(expectedMean - this.updatedMixtureMeans(:, 2)));
            
              
            % first check the individual filter states
            actualModeFilterStates = cellfun(@getState, this.modeFilters, 'UniformOutput', false);
            for j=1:this.numModes
                this.verifyClass(actualModeFilterStates{j}, ?Gaussian);
                % if suffices to compare mean and cov
                [actualMean, actualCov] = actualModeFilterStates{j}.getMeanAndCov();
                this.verifyEqualWithAbsTol(actualMean, this.updatedMixtureMeans(:, j));
                this.verifyEqualWithAbsTol(actualCov, this.updatedMixtureCovs(:, :, j));
            end
            
            % now the resulting Gaussian mixture
            updatedState = this.filterUnderTest.getState();
            this.verifyClass(updatedState, ?GaussianMixture);
            
            [~, ~, actualWeights] = updatedState.getComponents();
            [actualMean, actualCov] = this.filterUnderTest.getStateMeanAndCov();
            
            % first, a sanity check: posterior cov must be smaller than prior
            % cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(this.stateGaussianMixtureCov - actualCov), 0);
            
            this.verifyEqualWithAbsTol(actualWeights, expectedWeights);
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqual(actualCov, actualCov'); % shall be symmetric
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
        end
        
        %% verifyUpdateGaussian
        function verifyUpdateGaussian(this)
            % first check the individual filter states
            actualModeFilterStates = cellfun(@getState, this.modeFilters, 'UniformOutput', false);
            
            residualCov = this.C * this.stateGaussianCov * this.C' + this.V;
            gain = this.stateGaussianCov * this.C' * inv(residualCov);
            % mean should be unchanged as we perform update with expected measurement
            expectedMean = this.stateGaussianMean; 
            expectedCov = this.stateGaussianCov - gain * residualCov * gain';
            for j=1:this.numModes
                this.verifyClass(actualModeFilterStates{j}, ?Gaussian);
                % if suffices to compare mean and cov
                [actualMean, actualCov] = actualModeFilterStates{j}.getMeanAndCov();
                
                this.verifyEqualWithAbsTol(actualMean, expectedMean);
                this.verifyEqualWithAbsTol(actualCov, expectedCov);
            end
            
            % now the resulting Gaussian mixture
            updatedState = this.filterUnderTest.getState();
            this.verifyClass(updatedState, ?GaussianMixture);
            
            [~, ~, actualWeights] = updatedState.getComponents();
            [actualMean, actualCov] = this.filterUnderTest.getStateMeanAndCov();
           
            % mode probabilities should be as before the update; 
            expectedWeights = [1/2 1/2];
            this.verifyEqualWithAbsTol(actualWeights, expectedWeights);
            % posterior estimate shall be unchanged
            this.verifyEqualWithAbsTol(actualMean, this.stateGaussianMean);
            this.verifyEqual(actualCov, actualCov'); % shall be symmetric
            % mixture covariance is equal to that of the individual Gaussians
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
        end
    end
    
    methods (Access = protected)
        function filter = initFilterUnderTest(this)
            filter = IMMF(this.modeFilters, this.modeTransitionMatrix, this.filterName);
        end
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.filterName = 'TestIMMFilter';
            
            initProperties@BaseIMMFTest(this);
            this.computeUpdateGaussianMixture();
        end
    end
       
    methods (Test)
        
        %% testIMMFNoName
        function testIMMFNoName(this)
            immf = IMMF(this.modeFilters, this.modeTransitionMatrix);
            
            this.verifyEqual(immf.getName(), IMMFTest.defaultFiltername);
        end
        
        %% testIMMFInvalidModeFilters
        function testIMMFInvalidModeFilters(this)
            expectedErrId = 'MATLAB:class:RequireClass';
            
            % no cell passed
            invalidModeFilters = eye(this.dimX);
            this.verifyError(@() IMMF(invalidModeFilters, this.modeTransitionMatrix, this.filterName), expectedErrId);
            
            expectedErrId = 'Filter:InvalidModeFilters:InvalidFilterType';
            
            % cell with invalid filter passed
            invalidModeFilters = {EKF('KF'); this.filterUnderTest};
            this.verifyError(@() IMMF(invalidModeFilters, this.modeTransitionMatrix, this.filterName), expectedErrId);
        end
        
        %% testIMMFInvalidModeTransitionMatrix
        function testIMMFInvalidModeTransitionMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';
            
            % transition matrix must be square
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1];
            this.verifyError(@() IMMF(this.modeFilters, invalidTransitionMatrix, this.filterName), expectedErrId);
            
            % transition matrix must be square and have appropriate
            % dimensions
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1; 0.5 0.25 0.25];
            this.verifyError(@() IMMF(this.modeFilters, invalidTransitionMatrix, this.filterName), expectedErrId);
            
            % row sums must be 1
            invalidTransitionMatrix = [0.7 0.3; 0.6 0.5];
            this.verifyError(@() IMMF(this.modeFilters, invalidTransitionMatrix, this.filterName), expectedErrId);
            
            % entries must be in [0,1]
            invalidTransitionMatrix = [1.3 -0.3; 0.6 0.4];
            this.verifyError(@() IMMF(this.modeFilters, invalidTransitionMatrix, this.filterName), expectedErrId);
        end
        
        %% testIMMF
        function testIMMF(this)
            immf = IMMF(this.modeFilters, this.modeTransitionMatrix, this.filterName);
            
            this.verifyEqual(immf.getName(), this.filterName);
            this.verifyEqual(immf.numModes, this.numModes);
            this.verifyEqual(immf.modeTransitionProbs, this.modeTransitionMatrix);
        end
        
        %% testPerformUpdateInvalidMeasModel
        function testPerformUpdateInvalidMeasModel(this)
            expectedErrId = 'Filter:UnsupportedMeasurementModel';
                      
            % invalid model: number of meas models in cell array does not
            % match
            invalidMeasModel = {this.measModel, this.measModel, this.measModel};
            this.verifyError(@() this.filterUnderTest.update(invalidMeasModel, this.measurement), expectedErrId);
            
            % invalid model: not a valid model
            % the error now should actually be thrown by the GaussianFilter
            % class
            expectedErrId = 'Filter:InvalidUnobservableStateDimension';
            invalidMeasModel = eye(this.dimY);
            this.verifyError(@() this.filterUnderTest.update(invalidMeasModel, this.measurement), expectedErrId);
        end
        
        %% testPerformUpdateGaussianMixtureOneMeasModel
        function testPerformUpdateGaussianMixtureOneMeasModel(this)
            % set the prior state and update to obtain the posterior
            this.filterUnderTest.setState(this.stateGaussianMixture);
            this.filterUnderTest.update(this.measModel, this.measurement);
            
            this.verifyUpdateGaussianMixture();
        end
        
        %% testPerformUpdateGaussianMixtureMeasModelsPerMode
        function testPerformUpdateGaussianMixtureMeasModelsPerMode(this)
            % set the prior state and update to obtain the posterior
            this.filterUnderTest.setState(this.stateGaussianMixture);
            % pass a cell array of 2 meas models
            measModels = {this.measModel this.measModel};
            this.filterUnderTest.update(measModels, this.measurement);
            
            this.verifyUpdateGaussianMixture();
        end
         
        
        %% testPerformUpdateGaussianOneMeasModel
        function testPerformUpdateGaussianOneMeasModel(this)
            % check whether the prior estimate remains unchanged if
            % expected measurement is passed to update
            
            % set the prior state and update to obtain the posterior
            this.filterUnderTest.setState(this.stateGaussian);
            expectedMeas = this.C * this.stateGaussianMean;
            this.filterUnderTest.update(this.measModel, expectedMeas);
                                    
            this.verifyUpdateGaussian();
        end
                     
        %% testPerformUpdateGaussianMeasModelsPerMode
        function testPerformUpdateGaussianMeasModelsPerMode(this)
            % check whether the prior estimate remains unchanged if
            % expected measurement is passed to update
            
            % set the prior state and update to obtain the posterior
            this.filterUnderTest.setState(this.stateGaussian);
            expectedMeas = this.C * this.stateGaussianMean;
            % pass a cell array of 2 meas models
            measModels = {this.measModel this.measModel};
            this.filterUnderTest.update(measModels, expectedMeas);
                        
            this.verifyUpdateGaussian();
        end
        
%%
%%
        %% testPerformPredictionInvalidSysModel
        function testPerformPredictionInvalidSysModel(this)
            expectedErrId = 'Filter:UnsupportedSystemModel';
                    
            % invalid sys model: neither a JumpLinearSystem nor an
            % arbitrary SystemModel
            invalidSysModel = eye(this.dimU);
            this.verifyError(@() this.filterUnderTest.predict(invalidSysModel), expectedErrId);
            
            % now an invalid JumpLinearSystem: number of modes does not match
            invalidSysModel = JumpLinearSystemModel(1, this.modePlantModels(1));
            this.verifyError(@() this.filterUnderTest.predict(invalidSysModel), expectedErrId);
        end
        
        %% testPerformPredictionGaussianMixture
        function testPerformPredictionGaussianMixture(this)
                       
            this.filterUnderTest.setState(this.stateGaussianMixture);
            this.filterUnderTest.predict(this.jumpLinearPlantModel);
            
            this.verifyPredictionGaussianMixture();
        end
                       
        %% testPerformPredictionGaussian
        function testPerformPredictionGaussian(this)
            % init the filter with a Gaussian, and use same plant model for
            % all modes; equivalent to simply using a KF
            plant = LinearPlant(this.A1, this.B, this.W);
            plant.setSystemInput(this.input);
            
            
            this.filterUnderTest.setState(this.stateGaussian);
            this.filterUnderTest.predict(plant);
            
            this.verifyPredictionGaussian();
        end
    end
    
end

