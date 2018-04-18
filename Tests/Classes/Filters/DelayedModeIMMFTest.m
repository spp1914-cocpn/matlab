classdef DelayedModeIMMFTest < BaseIMMFTest
    % Test cases for DelayedModeIMMF.
    
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
        defaultFiltername = 'IMM (with Delayed Measurements and Modes)';
        zeroMaxMeasDelay = 0;
    end
    
    properties (Access = private)
        filterName;
        zeroDelayFiltername
        maxMeasDelay;
        
        updatedMixtureWeights;
        updatedMixtureMeans;
        updatedMixtureCovs;
    end
    
    methods (Access = private)
        %% computeZeroMeasDelayUpdateGaussianMixture
        function computeZeroMeasDelayUpdateGaussianMixture(this)
            % assumes that computePredictionGaussianMixture from base class has already
            % been called
            
            % Use straightforward implementation of the algorithm given
            % in (Table II) in
            %
            %   X. Rong Li and Vesselin P. Jilkov,
            %   A survey of maneuvering target tracking - Part V: Multiple-Model Methods,
            %   IEEE Transactions on Aerospace and Electronic Systems 41.4 (2005): 1255-1321.
            
            % we use a single measurement model for all modes
            % first the innovations/ residuals
            residuals(:, 2) = this.measurement - this.C * this.predictedMixtureMeans(:, 2);
            residuals(:, 1) = this.measurement - this.C * this.predictedMixtureMeans(:, 1);
            residualCovs(:, :, 2) = this.C * this.predictedMixtureCovs(:, :, 2) * this.C' + this.V;
            residualCovs(:, :, 1) = this.C * this.predictedMixtureCovs(:, :, 1) * this.C' + this.V;
            
            % compute the individual gains for all filters
            gains(:, :, 2) = this.predictedMixtureCovs(:, :, 2) * this.C' * inv(residualCovs(:, :, 2));
            gains(:, :, 1) = this.predictedMixtureCovs(:, :, 1) * this.C' * inv(residualCovs(:, :, 1));
             
            % now perform the update
            this.updatedMixtureMeans(:, 2) = this.predictedMixtureMeans(:, 2) + gains(:, :, 2) * residuals(:, 2);
            this.updatedMixtureMeans(:, 1) = this.predictedMixtureMeans(:, 1) + gains(:, :, 1) * residuals(:, 1);
            this.updatedMixtureCovs(:, :, 2) = this.predictedMixtureCovs(:, :, 2) - gains(:, :, 2) * residualCovs(:, :, 2) * gains(:, :, 2)';
            this.updatedMixtureCovs(:, :, 1) = this.predictedMixtureCovs(:, :, 1) - gains(:, :, 1) * residualCovs(:, :, 1) * gains(:, :, 1)';
            
            % for the weight updates we need the Gaussian-assumed likelihoods
            likelihoods(2) = mvnpdf(residuals(:, 2), zeros(this.dimY, 1), residualCovs(:, :, 2));
            likelihoods(1) = mvnpdf(residuals(:, 1), zeros(this.dimY, 1), residualCovs(:, :, 1));

            normalizationConstant = this.predictedMixtureWeights(1) * likelihoods(1) + this.predictedMixtureWeights(2) * likelihoods(2);
            this.updatedMixtureWeights(2) = this.predictedMixtureWeights(2) * likelihoods(2) / normalizationConstant;
            this.updatedMixtureWeights(1) = this.predictedMixtureWeights(1) * likelihoods(1) / normalizationConstant;
        end
        
        %% verifyZeroMeasDelayUpdateGaussianMixture
        function verifyZeroMeasDelayUpdateGaussianMixture(this)
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
                
                % sanity checks for the individual covs:
                % posterior cov should be smaller than prior cov
                this.verifyGreaterThan(eig(this.predictedMixtureCovs(:, :, j) - actualCov), -sqrt(eps));
            end
            
            % now the resulting Gaussian mixture
            updatedState = this.filterUnderTest.getState();
            this.verifyClass(updatedState, ?GaussianMixture);
            
            [~, ~, actualWeights] = updatedState.getComponents();
            [actualMean, actualCov] = this.filterUnderTest.getStateMeanAndCov();
                      
            this.verifyEqualWithAbsTol(actualWeights, expectedWeights);
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqual(actualCov, actualCov'); % shall be symmetric
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
        end
    end
    
    methods (Access = protected)
        function filter = initFilterUnderTest(this)
            filter = DelayedModeIMMF(this.modeFilters, this.modeTransitionMatrix, this.maxMeasDelay, this.filterName);
        end
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.filterName = 'TestDelayedIMMF';
            this.zeroDelayFiltername = 'TestZeroDelayedKF';
                      
            this.maxMeasDelay = 3;
            
            initProperties@BaseIMMFTest(this);
            this.computeZeroMeasDelayUpdateGaussianMixture();
        end
    end
   
    methods (Test)
        %% testDelayedModeIMMFNoName
        function testDelayedModeIMMFNoName(this)
            filter = DelayedModeIMMF(this.modeFilters, this.modeTransitionMatrix, this.maxMeasDelay);
            
            this.verifyEqual(filter.getName(), DelayedModeIMMFTest.defaultFiltername);
        end
        
        %% testDelayedModeIMMFZeroMaxMeasDelay
        function testDelayedModeIMMFZeroMaxMeasDelay(this)
            filter = DelayedModeIMMF(this.modeFilters, this.modeTransitionMatrix, ...
                DelayedModeIMMFTest.zeroMaxMeasDelay, this.zeroDelayFiltername);
            
            this.verifyEqual(filter.getMaxMeasurementDelay(), DelayedModeIMMFTest.zeroMaxMeasDelay);
            this.verifyEqual(filter.getName(), this.zeroDelayFiltername);
        end
        
        %% testDelayedModeIMMFInvalidModeFilters
        function testDelayedModeIMMFInvalidModeFilters(this)
            expectedErrId = 'MATLAB:class:RequireClass';
            
            % no cell passed
            invalidModeFilters = eye(this.dimX);
            this.verifyError(@() DelayedModeIMMF(invalidModeFilters, this.modeTransitionMatrix, this.maxMeasDelay, this.filterName), ...
                expectedErrId);
            
            expectedErrId = 'Filter:InvalidModeFilters:InvalidFilterType';
            
           % cell with invalid filter passed
            invalidModeFilters = {EKF('KF'); this.filterUnderTest};
            this.verifyError(@() DelayedModeIMMF(invalidModeFilters, this.modeTransitionMatrix, this.maxMeasDelay, this.filterName), ...
                expectedErrId);
        end
        
         %% testDelayedModeIMMFInvalidModeTransitionMatrix
        function testDelayedModeIMMFInvalidModeTransitionMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';
            
            % transition matrix must be square
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1];
            this.verifyError(@() DelayedModeIMMF(this.modeFilters, invalidTransitionMatrix, this.maxMeasDelay, this.filterName),...
                expectedErrId);
            
            % transition matrix must be square and have appropriate
            % dimensions
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1; 0.5 0.25 0.25];
            this.verifyError(@() DelayedModeIMMF(this.modeFilters, invalidTransitionMatrix, this.maxMeasDelay, this.filterName),...
                expectedErrId);
            
            % row sums must be 1
            invalidTransitionMatrix = [0.7 0.3; 0.6 0.5];
            this.verifyError(@() DelayedModeIMMF(this.modeFilters, invalidTransitionMatrix, this.maxMeasDelay, this.filterName),...
                expectedErrId);
            
            % entries must be in [0,1]
            invalidTransitionMatrix = [1.3 -0.3; 0.6 0.4];
            this.verifyError(@() DelayedModeIMMF(this.modeFilters, invalidTransitionMatrix, this.maxMeasDelay, this.filterName),...
                expectedErrId);
        end
        
        %% testDelayedModeIMMF
        function testDelayedModeIMMF(this)
            filter = DelayedModeIMMF(this.modeFilters, this.modeTransitionMatrix, this.maxMeasDelay, this.filterName); 
           
            this.verifyEqual(filter.getMaxMeasurementDelay(), this.maxMeasDelay);
            this.verifyEqual(filter.getName(), this.filterName);
        end
%%
%%
        %% testGetPreviousModeEstimateInvalidFlag
        function testGetPreviousModeEstimateInvalidFlag(this)
            expectedErrId = 'Filter:GetPreviousModeEstimate:InvalidFlag';
            newState = this.stateGaussianMixture;
                        
            this.filterUnderTest.setState(newState);
            invalidFlag = eye(2); % not a logical
            this.verifyError(@() this.filterUnderTest.getPreviousModeEstimate(invalidFlag), expectedErrId);
            invalidFlag = [true false]; % not a scalar
            this.verifyError(@() this.filterUnderTest.getPreviousModeEstimate(invalidFlag), expectedErrId);
        end
        
        %% testGetPreviousModeEstimateGaussianMixture
        function testGetPreviousModeEstimateGaussianMixture(this)
            newState = this.stateGaussianMixture;
            expectedMaxMode = 1;
            expectedMaxModeProb = 2 / 3;
            
            this.filterUnderTest.setState(newState);
            afterStep = true;
            [maxMode, prob] = this.filterUnderTest.getPreviousModeEstimate(afterStep);
            
            this.verifyEqual(maxMode, expectedMaxMode);
            this.verifyEqual(prob, expectedMaxModeProb);
            
            afterStep = false;
            [maxMode, prob] = this.filterUnderTest.getPreviousModeEstimate(afterStep);
            
            this.verifyEqual(maxMode, expectedMaxMode);
            this.verifyEqual(prob, expectedMaxModeProb);
        end
        
        %% testGetPreviousModeEstimateGaussian
        function testGetPreviousModeEstimateGaussian(this)
            newState = this.stateGaussian;
            % mode probability is uniform, first maxima should be returned
            expectedMaxMode = 1;
            expectedMaxModeProb = 1/ 2;

            this.filterUnderTest.setState(newState);
             afterStep = true;
            [maxMode, prob] = this.filterUnderTest.getPreviousModeEstimate(afterStep);
            
            this.verifyEqual(maxMode, expectedMaxMode);
            this.verifyEqual(prob, expectedMaxModeProb);
            
            afterStep = false;
            [maxMode, prob] = this.filterUnderTest.getPreviousModeEstimate(afterStep);
            
            this.verifyEqual(maxMode, expectedMaxMode);
            this.verifyEqual(prob, expectedMaxModeProb);
        end
%%
%%
       %% testSetModeTransitionMatrix
        function testSetModeTransitionMatrix(this)
            % override the test method from the base class
            expectedNewTransitionMatrix = this.modeTransitionMatrix';
            
            this.filterUnderTest.setModeTransitionMatrix(expectedNewTransitionMatrix);
            
            this.verifyEqual(this.filterUnderTest.getModeTransitionMatrix(), expectedNewTransitionMatrix);
        end
%%
%%
        %% testPredict
        function testPredict(this)
            expectedErrId = 'Filter:InvalidPredictionStep';
            
            this.verifyError(@() this.filterUnderTest.predict(this.jumpLinearPlantModel), expectedErrId);
        end
%%
%%
        %% testUpdate
        function testUpdate(this)
            expectedErrId = 'Filter:InvalidUpdateStep';
        
            this.verifyError(@() this.filterUnderTest.update(this.measModel, this.measurement), expectedErrId);
        end
%%
%%
        %% testStepInvalidNumberOfArguments
        function testStepInvalidNumberOfArguments(this)
            expectedErrId = 'Filter:InvalidStep:InvalidNumArgs';
            
            % pass 6 arguments, i.e, an observed mode, but not a
            % corresponding delay
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, 0, this.numModes), ...
                expectedErrId);
            
        end
        
        %% testStepInvalidSystemModel
        function testStepInvalidSystemModel(this)
            expectedErrId = 'Filter:UnsupportedSystemModel';
            
            % try to use a single linear model for all modes
            invalidSysModel = LinearPlant(this.A1, this.B, this.W);
            this.verifyError(@() this.filterUnderTest.step(invalidSysModel, this.measModel, this.measurement, 0), ...
                expectedErrId);
        end
        
        %% testStepInvalidMeasurements
        function testStepInvalidMeasurements(this)
            expectedErrId = 'Filter:InvalidMeasurements';
            
            invalidMeasurements = ones(this.dimY, this.dimY, this.dimY); % not a matrix
            
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, invalidMeasurements, 0), ...
                expectedErrId);
        end
        
        %% testStepInvalidMeasDelays
        function testStepInvalidMeasDelays(this)
            expectedErrId = 'Filter:InvalidMeasDelay';
            
            % delays must not be a negative scalar
            invalidDelay = -1;
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, invalidDelay), ...
                expectedErrId);
            
             % delays must not be fractional
            invalidDelay = 1.5;
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, invalidDelay), ...
                expectedErrId);
            
            % if delays are given as nonnegative vector, must consist of
            % <numMeas> elements
            invalidDelays = [1 2]; % vector is too large
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, invalidDelays), ...
                expectedErrId);
        end
        
        %% function testStepInvalidModes(this)
        function testStepInvalidModes(this)
            expectedErrId = 'Filter:InvalidModeObservations';
            
            measDelay = 0;

            % we observe an invalid mode
            invalidMode = 0;
            modeDelay = 0;
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, ...
                measDelay, invalidMode, modeDelay), expectedErrId);
            
            % we observe a valid and an invalid mode
            invalidModes = [1 this.numModes + 1];
            modeDelays = [0 1];
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, ...
                measDelay, invalidModes, modeDelays), expectedErrId);
        end
        
        %% testStepInvalidModeDelays
        function testStepInvalidModeDelays(this)
            expectedErrId = 'Filter:InvalidModeDelays';
            
            measDelay = 0;
            observedModes = [1 this.numModes];
            
            invalidModeDelays = 1; % a single delay
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, ...
                measDelay, observedModes, invalidModeDelays), expectedErrId);
            
            invalidModeDelays = [1 1.5]; % a valid and an invalid delay
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, ...
                measDelay, observedModes, invalidModeDelays), expectedErrId);
            
            invalidModeDelays = [1 2 3]; % vector of delays is too large
            this.verifyError(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, ...
                measDelay, observedModes, invalidModeDelays), expectedErrId);
        end
        
        %% testStepIssueWarningDiscardModeObservation
        function testStepIssueWarningDiscardModeObservation(this)
            % we just observe a mode, don't have a measurmeent at hand
            % mode observation ispected to be discarded with warning, but
            % result should be a successful prediction
             this.filterUnderTest.setState(this.stateGaussianMixture);
            expectedWarningId = 'Filter:IgnoringModeObservations';
            
            modeDelay = this.maxMeasDelay + 2;
            modeObservation = 1;
            
            this.verifyWarning(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, [], ...
                [], modeObservation, modeDelay), expectedWarningId);
            
            % check if the prediction was done correctly
            this.verifyPredictionGaussianMixture();
             % finally check if the update data was stored corerctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testStepNoMeasurementsGaussianMixture
        function testStepNoMeasurementsGaussianMixture(this)
            this.filterUnderTest.setState(this.stateGaussianMixture);
            % in case no measurements are used, we just expect a
            % prediction
                      
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, [], []);
            
            this.verifyPredictionGaussianMixture();
            
            % finally check if the update data was stored corerctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testStepNoMeasurementsGaussian
        function testStepNoMeasurementsGaussian(this)
            this.filterUnderTest.setState(this.stateGaussian);
            
            % init the filter with a Gaussian, and use same plant model for
            % all modes; equivalent to simply using a KF for prediction
            plant = LinearPlant(this.A1, this.B, this.W);
                        
            modeModels = {  plant ,...
                           plant };
                                   
            jumpModel = JumpLinearSystemModel(this.numModes, modeModels);
            jumpModel.setSystemInput(this.input);
                                   
            this.filterUnderTest.step(jumpModel, this.measModel, [], []);
            
            this.verifyPredictionGaussian();
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testStepDiscardMeasurement
        function testStepDiscardMeasurement(this)
            this.filterUnderTest.setState(this.stateGaussianMixture);
            % in case all measurements are discarded, we just expect a
            % prediction
                        
            measDelay = this.maxMeasDelay + 1; % so the measurement will be discarded
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, this.measurement, measDelay);
            
            this.verifyPredictionGaussianMixture();
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 1);
        end
       
        %% testStepOverwriteExistingModeObservation
        function testStepOverwriteExistingModeObservation(this)
            this.filterUnderTest.setState(this.stateGaussianMixture);
            
            modeObservation = 1;
            modeDelay = 1;
            measDelay = 0;
            % first invoke step, should be warning-free
            this.assertWarningFree(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                this.measurement, measDelay, modeObservation, modeDelay));
            
            expectedWarningId = 'Filter:OverwritingModeObservation';
            % now invoke step again
            this.verifyWarning(@() this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                this.measurement, measDelay, modeObservation, modeDelay + 1), expectedWarningId);
        end
        
        %% testStepZeroMeasDelayGaussianMixtureOneMeasModel
        function testStepZeroMeasDelayGaussianMixtureOneMeasModel(this)
            this.filterUnderTest.setState(this.stateGaussianMixture);
            
            measDelay = DelayedModeIMMFTest.zeroMaxMeasDelay;
            % calling step without mode observations and a non-delayed
            % measurement should reduce to a combinted prediction + update
            % from an ordinary IMMF
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                this.measurement, measDelay)
            
            this.verifyZeroMeasDelayUpdateGaussianMixture();
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 1);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testStepZeroMeasDelayGaussianMixtureMeasModelsPerMode
        function testStepZeroMeasDelayGaussianMixtureMeasModelsPerMode(this)
            this.filterUnderTest.setState(this.stateGaussianMixture);
            
            measDelay = DelayedModeIMMFTest.zeroMaxMeasDelay;
            % calling step without mode observations and a non-delayed
            % measurement should reduce to a combinted prediction + update
            % from an ordinary IMMF
            % pass a cell array of 2 meas models
            measModels = {this.measModel this.measModel};
            this.filterUnderTest.step(this.jumpLinearPlantModel, measModels, ...
                this.measurement, measDelay)
            
            this.verifyZeroMeasDelayUpdateGaussianMixture();
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 1);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testStepOnlyModeObservationGaussianMixture
        function testStepOnlyModeObservationGaussianMixture(this)
            this.filterUnderTest.setState(this.stateGaussianMixture);
            
            modeDelay = 0; % we know the current mode
            modeObservation = 1;
            % calling step now without a measurement should be simply a
            % prediction from the IMMF as the current mode does not affect
            % the prediction from the last to the current step
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                [], [], modeObservation, modeDelay)

            this.verifyPredictionGaussianMixture();
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
    end
    
end

