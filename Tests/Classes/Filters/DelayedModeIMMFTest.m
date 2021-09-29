classdef DelayedModeIMMFTest < BaseIMMFTest
    % Test cases for DelayedModeIMMF.
    
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
            this.verifyGaussianMixture(this.updatedMixtureMeans, this.updatedMixtureCovs, this.updatedMixtureWeights);
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
            this.zeroDelayFiltername = 'TestZeroDelayedIMMF';
                      
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
            if verLessThan('matlab', '9.8')
                % Matlab R2018 or R2019
                expectedErrId = 'MATLAB:UnableToConvert';
            else
                expectedErrId = 'MATLAB:validation:UnableToConvert';
            end
            
            % no cell passed
            invalidModeFilters = eye(this.dimX);
            this.verifyError(@() DelayedModeIMMF(invalidModeFilters, this.modeTransitionMatrix, this.maxMeasDelay, this.filterName), ...
                expectedErrId);
            
            expectedErrId = 'Filter:InvalidModeFilters:InvalidFilterType';
            
           % cell with invalid filter passed
            invalidModeFilters = {EKF('KF'); this.filterUnderTest};
            this.verifyError(@() DelayedModeIMMF(invalidModeFilters, this.modeTransitionMatrix, this.maxMeasDelay, this.filterName), ...
                expectedErrId);
            
            % cell with invalid objects passed
            invalidModeFilters = {1; this};
            this.verifyError(@() DelayedModeIMMF(invalidModeFilters, this.modeTransitionMatrix, this.maxMeasDelay, this.filterName), ...
                expectedErrId);
            
            if verLessThan('matlab', '9.8')
                % Matlab R2018 or R2019
                expectedErrId = 'MATLAB:type:InvalidInputSize';
            else
                expectedErrId = 'MATLAB:validation:IncompatibleSize';
            end
            
            % cell with invalid dimensions passed
            invalidModeFilters = cell(2,2);
            invalidModeFilters{1} = {EKF('KF'); EKF('KF2')};
            invalidModeFilters{end} = {EKF('KF'); EKF('KF2')};
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
        
        %% testSetState
        function testSetState(this)
            this.filterUnderTest.setState(this.stateGaussianMixture);
            % test the side effect
            for j=1:this.maxMeasDelay + 1
                this.verifyEqual(this.filterUnderTest.stateHistory{j}, this.stateGaussianMixture);
            end
            this.verifyEqual(this.filterUnderTest.getState(), this.stateGaussianMixture);
            
            this.filterUnderTest.setState(this.stateGaussian);
            % test the side effect; estimate is stored as Gaussian mixture,
            % so check the moments
            for j=1:this.maxMeasDelay + 1
                [mean, cov] = this.filterUnderTest.stateHistory{j}.getMeanAndCov();
                this.verifyEqual(mean, this.stateGaussianMean);
                this.verifyEqual(cov, this.stateGaussianCov);
            end
            
            [mean, cov] = this.filterUnderTest.getState().getMeanAndCov();
            this.verifyEqual(mean, this.stateGaussianMean);
            this.verifyEqual(cov, this.stateGaussianCov);
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
            mixinClass = ?ModeTransitionMatrixChangeable;
            this.assertTrue(ismember(mixinClass.Name, superclasses(this.filterUnderTest)));
            
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
        
        %% testStepInvalidModes
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
            
            % check the history
            this.verifyEmpty(this.filterUnderTest.measurementHistory{3});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{2});
            this.verifyEqual(this.filterUnderTest.measurementHistory{1}, this.measurement);
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{4});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{3});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{2});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{1});
            [a,b,c] = this.filterUnderTest.stateHistory{4}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights)
            [a,b,c] = this.filterUnderTest.stateHistory{3}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights)           
            [a,b,c] = this.filterUnderTest.stateHistory{2}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights);
            this.verifyEqual(this.filterUnderTest.stateHistory{1}, this.filterUnderTest.getState());
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
            
            % check the history
            this.verifyEmpty(this.filterUnderTest.measurementHistory{3});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{2});
            this.verifyEqual(this.filterUnderTest.measurementHistory{1}, this.measurement);
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{4});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{3});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{2});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{1});
            [a,b,c] = this.filterUnderTest.stateHistory{4}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights)
            [a,b,c] = this.filterUnderTest.stateHistory{3}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights)           
            [a,b,c] = this.filterUnderTest.stateHistory{2}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights);
            this.verifyEqual(this.filterUnderTest.stateHistory{1}, this.filterUnderTest.getState());
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
                [], [], modeObservation, modeDelay);

            this.verifyPredictionGaussianMixture();
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.filterUnderTest.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
            
            % check the history
            this.verifyEmpty(this.filterUnderTest.measurementHistory{3});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{2});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{1});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{4});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{3});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{2});
            this.verifyEqual(this.filterUnderTest.trueModeHistory{1}, modeObservation);
            
            [a,b,c] = this.filterUnderTest.stateHistory{4}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights)
            [a,b,c] = this.filterUnderTest.stateHistory{3}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights)           
            [a,b,c] = this.filterUnderTest.stateHistory{2}.getComponents();
            this.verifyEqual(a, this.mixtureMeans);
            this.verifyEqual(b, this.mixtureCovs);
            this.verifyEqual(c, this.mixtureWeights);
            this.verifyEqual(this.filterUnderTest.stateHistory{1}, this.filterUnderTest.getState());
        end
        
        %% testStepDelayedModeObservationsMeasurementsGaussianMixture
        function testStepDelayedModeObservationsMeasurementsGaussianMixture(this)
            % a setup similar to TCP-like: we have the previous modes
            % available, but not all measurements
                        
            modeDelays = [1 3 2]; % we know the modes from the previous three time steps
            modeObservations = [1 1 2];
            measDelay = 2; 
            measurements = this.measurement;            
            
            this.filterUnderTest.setState(this.stateGaussianMixture.copy());
            this.jumpLinearPlantModel.setSystemInput(this.input);
            % now perform the steps to incorporate delayed mode observations
            % and measurement
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                measurements, measDelay, modeObservations, modeDelays);
                      
            % another measurement appears, and we know the previous mode
            % also the input changes
            newMeasDelay = 2;
            newMeasurement = 2 * this.measurement;
            newMode = 2;
            newModeDelay = 1;
            newInput = 2 * this.input;
            this.jumpLinearPlantModel.setSystemInput(newInput);
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                newMeasurement, newMeasDelay, newMode, newModeDelay);
                      
            % use an IMMF as reference
            % check that the state history is correctly maintained
            immf = IMMF(arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, ...
                'UniformOutput', false), this.modeTransitionMatrix);
            
            immf.setState(this.stateGaussianMixture.copy()); % time k-5
            plant = JumpLinearSystemModel(this.numModes, this.modePlantModels);
            plant.setSystemInput([]);
            
            immf.predict(plant); % time k-4
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0]));
                       
            immf.step(plant, this.measModel, measurements); % time k-3
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1]));
            
            [a,b,c]=this.filterUnderTest.stateHistory{4}.getComponents();
            this.verifyEqual(a, intermediateMeans);
            this.verifyEqual(b, intermediateCovs);
            this.verifyEqual(c, [0 1]);
            
            immf.step(plant, this.measModel, newMeasurement); % time k-2      
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0]));
            
            [a,b,c]=this.filterUnderTest.stateHistory{3}.getComponents();
            this.verifyEqual(a, intermediateMeans);
            this.verifyEqual(b, intermediateCovs);
            this.verifyEqual(c, [1 0]);
            
            plant.setSystemInput(this.input);
            immf.predict(plant); % state at time k-1
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1]));
            
            [a,b,c]=this.filterUnderTest.stateHistory{2}.getComponents();
            this.verifyEqual(a, intermediateMeans);
            this.verifyEqual(b, intermediateCovs);
            this.verifyEqual(c, [0 1]);
            
            plant.setSystemInput(newInput);
            immf.predict(plant); % state at time k
            
            [expectedMixtureMeans, expectedMixtureCovs, expectedWeights] = immf.getState().getComponents();
            
            this.verifyGaussianMixture(expectedMixtureMeans, expectedMixtureCovs, expectedWeights);            
            this.verifyEqual(this.filterUnderTest.stateHistory{1}, immf.getState());
        end
        
        %% testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffSysMatrix
        function testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffSysMatrix(this)
            % a setup similar to TCP-like: we have the previous modes
            % available, but not all measurements
            % test the case that after some filter invocations, the
            % underlying system dynamics changes
                        
            modeDelays = [1 3 2]; % we know the modes from the previous three time steps
            modeObservations = [1 1 2];
            measDelay = 2; 
            measurements = this.measurement;            
            
            this.filterUnderTest.setState(this.stateGaussianMixture.copy());
            this.jumpLinearPlantModel.setSystemInput(this.input);
            % now perform the steps to incorporate delayed mode observations
            % and measurement
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                measurements, measDelay, modeObservations, modeDelays);
            % another measurement appears, and we know the previous mode
            % also the input changes
            % and the system matrices
            newMeasDelay = 2;
            newMeasurement = 2 * this.measurement;
            newMode = 2;
            newModeDelay = 1;
            newInput = 2 * this.input;
            
            % copy the plant model before changing its properties
            plantCopy = this.jumpLinearPlantModel.copy();
            
            this.jumpLinearPlantModel.setSystemInput(newInput);
            this.jumpLinearPlantModel.setSystemMatrixForMode(this.A1 / 2, 1);
            this.jumpLinearPlantModel.setSystemMatrixForMode(this.A2 * 4, 2);
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                newMeasurement, newMeasDelay, newMode, newModeDelay);
            
            % use an IMMF as reference
            immf = IMMF(arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, ...
                'UniformOutput', false), this.modeTransitionMatrix);
            
            immf.setState(this.stateGaussianMixture.copy());
            plantCopy.setSystemInput([]);
            
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            immf.step(plantCopy, this.measModel, measurements); 
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1])); % mode observation
            
            immf.step(plantCopy, this.measModel, newMeasurement);            
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            plantCopy.setSystemInput(this.input);
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1]));
            plantCopy.setSystemMatrixForMode(this.A1 / 2, 1);
            plantCopy.setSystemMatrixForMode(this.A2 * 4, 2);
            plantCopy.setSystemInput(newInput);
            immf.predict(plantCopy);
            
            [expectedMixtureMeans, expectedMixtureCovs, expectedWeights] = immf.getState().getComponents();
            
            this.verifyGaussianMixture(expectedMixtureMeans, expectedMixtureCovs, expectedWeights);
            % check the side effect
            this.verifyEmpty(this.filterUnderTest.sysModelsHistory{1});
            this.verifyEqual(this.filterUnderTest.sysModelsHistory{2}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{3}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{4}, this.jumpLinearPlantModel);
            
            this.verifyEmpty(this.filterUnderTest.measurementHistory{1});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{2});
            this.verifyEqual(this.filterUnderTest.measurementHistory{3}, newMeasurement);            
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{1});
            this.verifyEqual(this.filterUnderTest.trueModeHistory{2}, newMode);
            this.verifyEqual(this.filterUnderTest.trueModeHistory{3}, modeObservations(1));
            this.verifyEqual(this.filterUnderTest.trueModeHistory{4}, modeObservations(3));
        end
        
        %% testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffInputMatrix
        function testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffInputMatrix(this)
            % a setup similar to TCP-like: we have the previous modes
            % available, but not all measurements
            % test the case that after some filter invocations, the
            % underlying system dynamics (input matrix B) changes
                        
            modeDelays = [1 3 2]; % we know the modes from the previous three time steps
            modeObservations = [1 1 2];
            measDelay = 2; 
            measurements = this.measurement;            
            
            this.filterUnderTest.setState(this.stateGaussianMixture.copy());
            this.jumpLinearPlantModel.setSystemInput(this.input);
            % now perform the steps to incorporate delayed mode observations
            % and measurement
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                measurements, measDelay, modeObservations, modeDelays);
            % another measurement appears, and we know the previous mode
            % also the input changes
            % and the system matrices
            newMeasDelay = 2;
            newMeasurement = 2 * this.measurement;
            newMode = 2;
            newModeDelay = 1;
            newInput = 2 * this.input;
            
            % copy the plant model before changing its properties
            plantCopy = this.jumpLinearPlantModel.copy();
                        
            this.jumpLinearPlantModel.setSystemInput(newInput);
            this.jumpLinearPlantModel.setSystemInputMatrixForMode(ones(this.dimX, this.dimU), 1);
            this.jumpLinearPlantModel.setSystemInputMatrixForMode(zeros(this.dimX, this.dimU), 2);
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                newMeasurement, newMeasDelay, newMode, newModeDelay);
            
            % and finally, another prediction step (no measurements given),
            % but with another input
            this.jumpLinearPlantModel.setSystemInput(newInput * 1.5);
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, [], []);
            
            % use an IMMF as reference
            immf = IMMF(arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, ...
                'UniformOutput', false), this.modeTransitionMatrix);
            
            immf.setState(this.stateGaussianMixture.copy());
            plantCopy.setSystemInput([]);
            
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            immf.step(plantCopy, this.measModel, measurements); 
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1])); % mode observation
            
            immf.step(plantCopy, this.measModel, newMeasurement);            
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            plantCopy.setSystemInput(this.input);
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1]));
            plantCopy.setSystemInputMatrixForMode(ones(this.dimX, this.dimU), 1);
            plantCopy.setSystemInputMatrixForMode(zeros(this.dimX, this.dimU), 2);
            plantCopy.setSystemInput(newInput);
            immf.predict(plantCopy);
            
            plantCopy.setSystemInput(newInput * 1.5);
            immf.predict(plantCopy);
            
            [expectedMixtureMeans, expectedMixtureCovs, expectedWeights] = immf.getState().getComponents();
            
            this.verifyGaussianMixture(expectedMixtureMeans, expectedMixtureCovs, expectedWeights);
            % check the side effect
            this.verifyEmpty(this.filterUnderTest.sysModelsHistory{1});
            this.verifyEqual(this.filterUnderTest.sysModelsHistory{2}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{3}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{4}, this.jumpLinearPlantModel);
            
            this.verifyEmpty(this.filterUnderTest.measurementHistory{1});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{2});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{3});            
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{1});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{2});
            this.verifyEqual(this.filterUnderTest.trueModeHistory{3}, newMode);
            this.verifyEqual(this.filterUnderTest.trueModeHistory{4}, modeObservations(1));
        end
        
        %% testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffNoise
        function testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffNoise(this)
            % a setup similar to TCP-like: we have the previous modes
            % available, but not all measurements
            % test the case that after some filter invocations, the
            % system noise changes
                        
            modeDelays = [1 3 2]; % we know the modes from the previous three time steps
            modeObservations = [1 1 2];
            measDelay = 2; 
            measurements = this.measurement;            
            
            this.filterUnderTest.setState(this.stateGaussianMixture.copy());
            this.jumpLinearPlantModel.setSystemInput(this.input);
            % now perform the steps to incorporate delayed mode observations
            % and measurement
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                measurements, measDelay, modeObservations, modeDelays);
            % another measurement appears, and we know the previous mode
            % also the input changes
            % and the process noise changes
            newMeasDelay = 2;
            newMeasurement = 2 * this.measurement;
            newMode = 2;
            newModeDelay = 1;
            newInput = 2 * this.input;
            
            % copy the plant model before changing its properties
            plantCopy = this.jumpLinearPlantModel.copy();
            
            this.jumpLinearPlantModel.setSystemInput(newInput);
            this.jumpLinearPlantModel.setSystemNoiseCovarianceMatrixForMode(eye(this.dimX), 1);
            this.jumpLinearPlantModel.setSystemNoiseCovarianceMatrixForMode(this.W / 2, 2);
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                newMeasurement, newMeasDelay, newMode, newModeDelay);
            
            % and finally, another prediction step (no measurements given),
            % but with another input, and new noise cov
            this.jumpLinearPlantModel.setSystemInput(newInput * 1.5);
            % simply switch the noise covs
            this.jumpLinearPlantModel.setSystemNoiseCovarianceMatrixForMode(eye(this.dimX), 2);
            this.jumpLinearPlantModel.setSystemNoiseCovarianceMatrixForMode(this.W / 2, 1);
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, [], []);
            
            % use an IMMF as reference
            immf = IMMF(arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, ...
                'UniformOutput', false), this.modeTransitionMatrix);
            
            immf.setState(this.stateGaussianMixture.copy());
            plantCopy.setSystemInput([]);
            
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            immf.step(plantCopy, this.measModel, measurements); 
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1])); % mode observation
            
            immf.step(plantCopy, this.measModel, newMeasurement);            
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            plantCopy.setSystemInput(this.input);
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1]));
            plantCopy.setSystemNoiseCovarianceMatrixForMode(eye(this.dimX), 1);
            plantCopy.setSystemNoiseCovarianceMatrixForMode(this.W / 2, 2);
            plantCopy.setSystemInput(newInput);
            immf.predict(plantCopy);
            
            plantCopy.setSystemInput(newInput * 1.5);
            plantCopy.setSystemNoiseCovarianceMatrixForMode(eye(this.dimX), 2);
            plantCopy.setSystemNoiseCovarianceMatrixForMode(this.W / 2, 1);
            immf.predict(plantCopy);
            
            [expectedMixtureMeans, expectedMixtureCovs, expectedWeights] = immf.getState().getComponents();
            
            this.verifyGaussianMixture(expectedMixtureMeans, expectedMixtureCovs, expectedWeights);
            % check the side effect
            this.verifyEmpty(this.filterUnderTest.sysModelsHistory{1});
            this.verifyEqual(this.filterUnderTest.sysModelsHistory{2}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{3}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{4}, this.jumpLinearPlantModel);
            
            this.verifyEmpty(this.filterUnderTest.measurementHistory{1});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{2});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{3});            
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{1});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{2});
            this.verifyEqual(this.filterUnderTest.trueModeHistory{3}, newMode);
            this.verifyEqual(this.filterUnderTest.trueModeHistory{4}, modeObservations(1));
        end
%%
%%
        %% testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffTransitionMatrices
        function testStepDelayedModeObservationsMeasurementsGaussianMixtureDiffTransitionMatrices(this)
            % a setup similar to TCP-like: we have the previous modes
            % available, but not all measurements
            % test the case that after some filter invocations, the
            % mode transition matrix changes
                        
            modeDelays = [1 3 2]; % we know the modes from the previous three time steps
            modeObservations = [1 1 2];
            measDelay = 2; 
            measurements = this.measurement;            
            
            this.filterUnderTest.setState(this.stateGaussianMixture.copy());
            this.jumpLinearPlantModel.setSystemInput(this.input);
            % now perform the steps to incorporate delayed mode observations
            % and measurement
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                measurements, measDelay, modeObservations, modeDelays);
            % another measurement appears, and we know the previous mode
            % also the input changes
            % and the mode transition matrix changes
            newMeasDelay = 2;
            newMeasurement = 2 * this.measurement;
            newMode = 2;
            newModeDelay = 1;
            newInput = 2 * this.input;
            newModeTransitionMatrix = [1/3 2/3; 3/4 1/4];
            
            % copy the plant model
            plantCopy = this.jumpLinearPlantModel.copy();
            
            this.jumpLinearPlantModel.setSystemInput(newInput);
            this.filterUnderTest.setModeTransitionMatrix(newModeTransitionMatrix);
            % this step is not affected by the changed transition matrix
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, ...
                newMeasurement, newMeasDelay, newMode, newModeDelay);
            
            % and finally, another two prediction steps (no measurements given),
            % now the impact of the new transition matrix becomes visible
            this.jumpLinearPlantModel.setSystemInput(newInput * 1.5);
            this.filterUnderTest.setModeTransitionMatrix(fliplr(newModeTransitionMatrix));
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, [], []);
            this.jumpLinearPlantModel.setSystemInput([]);
            this.filterUnderTest.step(this.jumpLinearPlantModel, this.measModel, [], []);
            
            % use an IMMF as reference
            immf = IMMF(arrayfun(@(mode) EKF(sprintf('KF for mode %d', mode)), 1:this.numModes, ...
                'UniformOutput', false), this.modeTransitionMatrix);
            
            immf.setState(this.stateGaussianMixture.copy());
            plantCopy.setSystemInput([]);
            
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            immf.step(plantCopy, this.measModel, measurements); 
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1])); % mode observation
            
            immf.step(plantCopy, this.measModel, newMeasurement);            
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [1 0])); % mode observation
                        
            plantCopy.setSystemInput(this.input);
            immf.predict(plantCopy);
            [intermediateMeans, intermediateCovs, ~] = immf.getState().getComponents();
            immf.setState(GaussianMixture(intermediateMeans, intermediateCovs, [0 1]));
            plantCopy.setSystemInput(newInput);
            immf.predict(plantCopy);
            
            plantCopy.setSystemInput(newInput * 1.5);
            immf.setModeTransitionMatrix(newModeTransitionMatrix);
            immf.predict(plantCopy);
            
            plantCopy.setSystemInput([]);
            immf.setModeTransitionMatrix(fliplr(newModeTransitionMatrix));
            immf.predict(plantCopy);
            
            [expectedMixtureMeans, expectedMixtureCovs, expectedWeights] = immf.getState().getComponents();
            
            this.verifyGaussianMixture(expectedMixtureMeans, expectedMixtureCovs, expectedWeights);
            % check the side effect
            this.verifyEmpty(this.filterUnderTest.sysModelsHistory{1});
            this.verifyEqual(this.filterUnderTest.sysModelsHistory{2}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{3}, this.jumpLinearPlantModel);
            this.verifyNotEqual(this.filterUnderTest.sysModelsHistory{4}, this.jumpLinearPlantModel);
            
            this.verifyEmpty(this.filterUnderTest.measurementHistory{1});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{2});
            this.verifyEmpty(this.filterUnderTest.measurementHistory{3});            
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{1});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{2});
            this.verifyEmpty(this.filterUnderTest.trueModeHistory{3});
            this.verifyEqual(this.filterUnderTest.trueModeHistory{4}, newMode);
            
        end
    end    
end

