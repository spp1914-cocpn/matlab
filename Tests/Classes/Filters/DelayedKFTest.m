classdef DelayedKFTest < matlab.unittest.TestCase
    % Test cases for DelayedKF.
    
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
        defaultFiltername = 'Delayed KF';
        absTol = 1e-12;
        zeroMeasDelay = 0;
    end
    
    properties (Access = private)
        zeroDelayFilter;
        zeroDelayFilterModel;
        delayFilter;
        delayFilterModel;
        
        dimX;
        dimU;
        dimY;
        numPossibleInputs;
        inputProbs;
        maxMeasDelay;
        filtername;
        zeroDelayFiltername;
        A;
        B;
        C; % meas matrix
        W; % sys noise cov
        V; % meas noise cov
        W_red; % subspace noise
        G; % the corresponding sys noise matrix
        
        uncertainInputs;
        inputMean;
        inputCov;
        
        stateGaussianMean;
        stateGaussianCov;
        stateGaussian;
        
        stateGaussianMixtureMean;
        stateGaussianMixtureCov;
        stateGaussianMixture;
        
        numMeas;
        delayedMeasurements;
        applicableMeasurements;
        notApplicableMeasurements;
        zeroDelayMeasurement;
        delays;
        applicableDelays;
        notApplicableDelays;
        measModel;
        measNoise;
        
        modeTransitionProbabilities;
    end
    
    methods (TestMethodSetup)
         %% initProperties
        function initProperties(this)
            this.maxMeasDelay = 4;
            this.numPossibleInputs = 5;
            % all inputs can be "active" with equal probability
            this.inputProbs = ones(1, this.numPossibleInputs) / this.numPossibleInputs;
            this.modeTransitionProbabilities = repmat(this.inputProbs, this.numPossibleInputs, 1); % stationary distribution is uniform
            this.filtername = 'TestDelayedKF';
            this.zeroDelayFiltername = 'TestZeroDelayedKF';
            
            this.dimX = 3;
            this.dimU = 2;
            this.dimY = 2;
            
            this.initStateDistributions();
            this.initMeasurements();
            this.initInputs();
            
            this.A = 2 * eye(this.dimX);
            this.B = ones(this.dimX, this.dimU);
            this.C = ones(this.dimY, this.dimX) * 3;
            this.W = gallery('moler', this.dimX);
            this.W_red = gallery('moler', this.dimX -1); 
            this.G = ones(this.dimX, this.dimX - 1);            
            this.V = gallery('minij', this.dimY); 
                            
            this.measModel = LinearMeasurementModel(this.C);
            this.measNoise = Gaussian(zeros(this.dimY, 1), this.V);
            this.measModel.setNoise(this.measNoise);
                   
            this.delayFilterModel = LinearPlant(this.A, this.B, this.W);
            this.zeroDelayFilterModel = LinearPlant(this.A, this.B, this.W_red, this.G);
                                   
            this.delayFilter = DelayedKF(this.maxMeasDelay, this.modeTransitionProbabilities, this.filtername);
            this.zeroDelayFilter = DelayedKF(DelayedKFTest.zeroMeasDelay, this.modeTransitionProbabilities, this.zeroDelayFiltername);
            
            this.delayFilter.setPossibleInputs(this.uncertainInputs);
            this.zeroDelayFilter.setPossibleInputs(this.uncertainInputs);
        end
    end
    
    methods (Access = private)
        %% initMeasurements
        function initMeasurements(this)
            this.numMeas = 4;
            this.delayedMeasurements = ones(this.dimY, this.numMeas);
            this.delays = [0 1 5 3];
            % extract the expected measurements, i.e, those with delay <= maxMeasDelay
            this.applicableDelays = this.delays([1 2 4]);
            this.applicableMeasurements = this.delayedMeasurements(:, [1 2 4]);
              
            this.notApplicableMeasurements = this.delayedMeasurements(:, 3);
            this.notApplicableDelays = this.delays(3);
            
            this.zeroDelayMeasurement = this.delayedMeasurements(:, 1);
        end
        
        %% initInputs
        function initInputs(this)
             this.uncertainInputs = [1 2 3 4 5;
                                    2 3 4 5 1
                                    ];
            
            this.inputMean = [3; 3];
            diffs = bsxfun(@minus, this.inputMean, this.uncertainInputs);
            this.inputCov = (diffs * diffs') / this.numPossibleInputs;
        end
        
        %% initStateDistributions
        function initStateDistributions(this)
            this.stateGaussianMean = 0.25 * ones(this.dimX, 1);
            this.stateGaussianCov = 0.5 * eye(this.dimX);
            this.stateGaussian = Gaussian(this.stateGaussianMean, this.stateGaussianCov);
        
            mixtureMeans = [zeros(this.dimX, 1) ones(this.dimX, 1)];
            mixtureCovs(:, :, 2) = 2 * eye(this.dimX);
            mixtureCovs(:, :, 1) = 0.5 * eye(this.dimX);
            mixtureWeights = [2/3 1/3];
            this.stateGaussianMixture = GaussianMixture(mixtureMeans, mixtureCovs, mixtureWeights);
            this.stateGaussianMixtureMean = 1/3 * ones(this.dimX, 1);
            this.stateGaussianMixtureCov = eye(this.dimX) + 2/9 * ones(this.dimX); % cov of means + mean of covs
        end
        
        %% computeUpdatedMoments
        function [expectedUpdatedMean, expectedUpdatedCov] = computeUpdatedMoments(this)
            % first, compute the prediction step
            zeroMatrix = zeros(this.dimX, this.dimX * this.maxMeasDelay);
            % augmented system noise matrix
            augG = [eye(this.dimX); zeroMatrix'];
            augA = [this.A, zeroMatrix; ...
                Utils.blockDiag(eye(this.dimX), this.maxMeasDelay), zeroMatrix'];
            augCov = repmat(this.stateGaussianCov, this.maxMeasDelay + 1, this.maxMeasDelay + 1);
            
            predictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean;
            predictedCov = augA * augCov * augA' + augG * (this.W + this.B * this.inputCov * this.B') * augG';
            
            % Use straightforward implementation of the filter algorithm given
            % in section 4.2.1 in
            %
            %   Maryam Moayedi, Yung Kuan Foo and Yeng Chai Soha
            %   Filtering for networked control systems with single/multiple measurement packets
            %   subject to multiple-step measurement delays and multiple packet dropouts
            %   International Journal of Systems Science (2011)
            %
            gammaF = [eye(this.dimX) repmat(zeros(this.dimX), 1, this.maxMeasDelay)];
            augV = Utils.blockDiag(this.V, this.maxMeasDelay + 1);
            augMean = [predictedMean; repmat(this.stateGaussianMean, this.maxMeasDelay, 1)];  
            % process the delayed measurements one after another
            expectedUpdatedMean = predictedMean;
            expectedUpdatedCov = predictedCov;
            for j=1:numel(this.applicableDelays)
                delay = this.applicableDelays(j);
                meas = this.applicableMeasurements(:, j);
                augD = zeros(this.dimY, (this.maxMeasDelay + 1) * this.dimY);
                augD(:, this.dimY * delay + 1:this.dimY * (delay + 1)) = eye(this.dimY);
                augC = zeros(this.dimY, (this.maxMeasDelay + 1) * this.dimX);
                augC(:, this.dimX * delay + 1:this.dimX * (delay + 1)) = this.C;
                
                L = gammaF * expectedUpdatedCov * augC' * inv(augC * expectedUpdatedCov * augC' + augD * augV * augD');
                expectedUpdatedMean = expectedUpdatedMean + L * (meas - augC * augMean);
                F = [L' zeros(this.dimY, this.dimX * this.maxMeasDelay)]';
                expectedUpdatedCov = (eye((this.maxMeasDelay + 1) * this.dimX) - F * augC) * expectedUpdatedCov ...
                    * (eye((this.maxMeasDelay + 1) * this.dimX) - F * augC)' + F * augD * augV * augD' * F';
            end
            expectedUpdatedCov = expectedUpdatedCov(1:this.dimX, 1:this.dimX);
        end
        
        %% verifyEqualWithAbsTol
        function verifyEqualWithAbsTol(this, actual, expected)
            this.verifyEqual(actual, expected, 'AbsTol', DelayedKFTest.absTol);
        end
    end
    
    methods (Test)
        %% testDelayedKFNoName
        function testDelayedKFNoName(this)            
            filter = DelayedKF(this.maxMeasDelay, this.modeTransitionProbabilities);
            [mode, probability] = filter.getPreviousModeEstimate();
            
            this.verifyEqual(filter.getName(), DelayedKFTest.defaultFiltername);
            this.verifyEqual(filter.getStateDim(), 0);
            
            this.verifyEqual(mode, this.numPossibleInputs);
            this.verifyEqual(probability, 1);
        end
        
        %% testDelayedKFZeroMaxMeasDelay
        function testDelayedKFZeroMaxMeasDelay(this)
            filter = DelayedKF(DelayedKFTest.zeroMeasDelay, this.modeTransitionProbabilities, this.zeroDelayFiltername);            
            
            this.verifyEqual(filter.getMaxMeasurementDelay(), DelayedKFTest.zeroMeasDelay);
            this.verifyEqual(filter.getName(), this.zeroDelayFiltername);
            this.verifyEqual(filter.getStateDim(), 0);
        end
        
        %% testDelayedKFInvalidModeTransitionMatrix
        function testDelayedKFInvalidModeTransitionMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrix';
            
            % transition matrix must be square
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1];
            this.verifyError(@() DelayedKF(DelayedKFTest.zeroMeasDelay, invalidTransitionMatrix, this.filtername),...
                expectedErrId);    
            
            % row sums must be 1
            invalidTransitionMatrix = [0.7 0.3; 0.6 0.5];
            this.verifyError(@() DelayedKF(DelayedKFTest.zeroMeasDelay, invalidTransitionMatrix, this.filtername),...
                expectedErrId);
            
            % entries must be in [0,1]
            invalidTransitionMatrix = [1.3 -0.3; 0.6 0.4];
            this.verifyError(@() DelayedKF(DelayedKFTest.zeroMeasDelay, invalidTransitionMatrix, this.filtername),...
                expectedErrId);
        end
        
        %% testDelayedKF
        function testDelayedKF(this)
            filter = DelayedKF(this.maxMeasDelay, this.modeTransitionProbabilities, this.filtername);
            [mode, probability] = filter.getPreviousModeEstimate();
            
            this.verifyEqual(filter.getMaxMeasurementDelay(), this.maxMeasDelay);
            this.verifyEqual(filter.getName(), this.filtername);
            this.verifyEqual(filter.getStateDim(), 0);
            
            this.verifyEqual(mode, this.numPossibleInputs);
            this.verifyEqual(probability, 1);
        end
        
        %% testSetStateInvalidState
        function testSetStateInvalidState(this)
            expectedErrId = 'Filter:InvalidSystemState'; % issued from super class GaussianFilter
            
            invalidState = this; % not a distribution
            this.verifyError(@() this.delayFilter.setState(invalidState), expectedErrId);
        end
        
        %% testSetStateGaussian
        function testSetStateGaussian(this)
            this.delayFilter.setState(this.stateGaussian);
            actualState = this.delayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
    
            [actualMean, actualCov] = actualState.getMeanAndCov();
               
            this.verifyEqualWithAbsTol(actualMean, this.stateGaussianMean);
            this.verifyEqualWithAbsTol(actualCov, this.stateGaussianCov);
            
        end
        
        %% testSetStateGaussianMixture
        function testSetStateGaussianMixture(this)
            % use a Distribution that is not Gaussian
            this.delayFilter.setState(this.stateGaussianMixture);
            actualState = this.delayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
            
            [actualMean, actualCov] = actualState.getMeanAndCov();
                        
            this.verifyEqualWithAbsTol(actualMean, this.stateGaussianMixtureMean);
            this.verifyEqualWithAbsTol(actualCov, this.stateGaussianMixtureCov);
            
        end
%%
%%
        %% testGetPreviousModeEstimate
        function testGetPreviousModeEstimate(this)
            this.delayFilter.setState(this.stateGaussianMixture);
            % first assert, that the previous mode estimate is correct
            % (we assume that initially we are in last mode)
            [mode, probability] = this.delayFilter.getPreviousModeEstimate();            
            this.assertEqual(mode, this.numPossibleInputs);
            this.assertEqual(probability, 1);
            
            % now perform a step, without measurements or mode observations
            this.delayFilter.step(this.delayFilterModel, this.measModel, [], []);
            
            % the mode should have been predicted correctly
            % we end up in the limit distribution
            [newMode, newProbability] = this.delayFilter.getPreviousModeEstimate(); 
            this.verifyEqual(newProbability, this.inputProbs(1));
            this.verifyEqual(newMode, 1);
        end
%%
%%
        %% testGetModeEstimate
        function testGetModeEstimate(this)
            % use a different mode transition matrix in this test case
            transitionMatrix = eye(this.numPossibleInputs);
            filterUnderTest = DelayedKF(this.maxMeasDelay, transitionMatrix);
            
            % first assert, that the previous mode estimate is correct
            % (we assume that initially we are in last mode)
            [mode, probability] = filterUnderTest.getPreviousModeEstimate();            
            this.assertEqual(mode, this.numPossibleInputs);
            this.assertEqual(probability, 1);
            
            % with this transition matrix, we should always remain in that
            % mode
            [currentMode, currentModeProbability] = filterUnderTest.getModeEstimate();            
            this.verifyEqual(currentMode, this.numPossibleInputs);
            this.verifyEqual(currentModeProbability, 1);
        end
%%
%%
        %% testSetModeTransitionMatrixInvalidMatrix
        function testSetModeTransitionMatrixInvalidMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';
            
            mixinClass = ?ModeTransitionMatrixChangeable;
            this.assertTrue(ismember(mixinClass.Name, superclasses(this.delayFilter)));
            
            % transition matrix must be square
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1];
            this.verifyError(@() this.delayFilter.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
            
            % transition matrix must be square and have appropriate
            % dimensions (5x5 expected)
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1; 0.5 0.25 0.25];
            this.verifyError(@() this.delayFilter.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
            
            % row sums must be 1
            invalidTransitionMatrix = [0.7 0.3; 0.6 0.5];
            this.verifyError(@() this.delayFilter.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
            
            % entries must be in [0,1]
            invalidTransitionMatrix = [1.3 -0.3; 0.6 0.4];
            this.verifyError(@() this.delayFilter.setModeTransitionMatrix(invalidTransitionMatrix), expectedErrId);
        end
        
        %% testSetModeTransitionMatrix
        function testSetModeTransitionMatrix(this)
            mixinClass = ?ModeTransitionMatrixChangeable;
            this.assertTrue(ismember(mixinClass.Name, superclasses(this.delayFilter)));
            
            % first assert, that the previous mode estimate is correct
            % (we assume that initially we are in last mode)
            [mode, probability] = this.delayFilter.getPreviousModeEstimate();            
            this.assertEqual(mode, this.numPossibleInputs);
            this.assertEqual(probability, 1);
            
            expectedNewTransitionMatrix = eye(this.numPossibleInputs);
            
            this.delayFilter.setModeTransitionMatrix(expectedNewTransitionMatrix);
            % check a side effect: new transition matrix shouldn't have
            % impact on current true mode which is a prediction based on
            % previous estimate
            % we expect a uniform distribution which is the limit
            % distribution
            [currentMode, currentModeProbability] = this.delayFilter.getModeEstimate();            
            this.verifyEqual(currentMode, 1);
            this.verifyEqual(currentModeProbability, this.inputProbs(1));
        end
%%
%%
        %% testSetPossibleInputsInvalidInputs
        function testSetPossibleInputsInvalidInputs(this)
            expectedErrorId = 'DelayedKF:SetPossibleInputs:InvalidInputs';
            
            invalidInputs = this; % not a matrix
            this.verifyError(@() this.delayFilter.setPossibleInputs(invalidInputs), expectedErrorId);
                        
            invalidInputs = ones(this.dimU, this.numPossibleInputs + 1); % wrong number of cols
            this.verifyError(@() this.delayFilter.setPossibleInputs(invalidInputs), expectedErrorId);
            
            invalidInputs = ones(this.dimU, this.numPossibleInputs); % not finite
            invalidInputs(1) = nan;
            this.verifyError(@() this.delayFilter.setPossibleInputs(invalidInputs), expectedErrorId);
        end
        
        %% testSetPossibleInputs
        function testSetPossibleInputs(this)
            newInputs = ones(this.dimU, this.numPossibleInputs);
            
            this.delayFilter.setPossibleInputs(newInputs);
            this.verifyEqual(this.delayFilter.possibleInputs, newInputs);
            
            this.delayFilter.setPossibleInputs([]); % empty matrix is allowed
            this.verifyEmpty(this.delayFilter.possibleInputs);
        end
%%
%%
        %% testPredict
        function testPredict(this)
            expectedErrId = 'Filter:InvalidPredictionStep';
            
            this.verifyError(@() this.delayFilter.predict(this.delayFilterModel), expectedErrId);
        end
%%
%%
        %% testUpdate
        function testUpdate(this)
            expectedErrId = 'Filter:InvalidUpdateStep';
        
            this.verifyError(@() this.delayFilter.update(this.measModel, this.delayedMeasurements, this.delays), expectedErrId);
        end
%%
%%            
        %% testStepInvalidNumberOfArguments
        function testStepInvalidNumberOfArguments(this)
            expectedErrId = 'Filter:InvalidStep:InvalidNumArgs';
            
            % pass 6 arguments, i.e, an observed mode, but not a
            % corresponding delay
            observedMode = 1;
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, 0, observedMode), ...
                expectedErrId);            
        end
  
        %% testStepInvalidMeasurements
        function testStepInvalidMeasurements(this)
            expectedErrId = 'Filter:InvalidMeasurements';
            
            invalidMeasurements = ones(this.dimY, this.dimY, this.dimY); % not a matrix
            
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, invalidMeasurements, 0), ...
                expectedErrId);
        end

        %% testStepInvalidMeasDelays
        function testStepInvalidMeasDelays(this)
            expectedErrId = 'Filter:InvalidMeasDelay';
            
            % delays must not be a negative scalar
            invalidDelay = [1 -1 3];
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, invalidDelay), ...
                expectedErrId);
            
             % delays must not be fractional
            invalidDelay = [1.5 1 2];
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, invalidDelay), ...
                expectedErrId);
            
            % if delays are given as nonnegative vector, must consist of
            % <numMeas> elements
            invalidDelays = [1 2]; % vector is too short
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, invalidDelays), ...
                expectedErrId);
        end   
        
        %% testStepInvalidModes(this)
        function testStepInvalidModes(this)
            expectedErrId = 'Filter:InvalidModeObservations';

            % we observe an invalid mode
            invalidMode = 0;
            modeDelay = 0;
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, ...
                this.delays, invalidMode, modeDelay), expectedErrId);
            
            % we observe a valid and an invalid mode
            invalidModes = [1 this.numPossibleInputs + 1];
            modeDelays = [0 1];
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, ...
                this.delays, invalidModes, modeDelays), expectedErrId);
        end
        
        %% testStepInvalidModeDelays
        function testStepInvalidModeDelays(this)
            expectedErrId = 'Filter:InvalidModeDelays';
            
            observedModes = [1 this.numPossibleInputs];
            
            invalidModeDelays = 1; % a single delay
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, ...
                this.delays, observedModes, invalidModeDelays), expectedErrId);
            
            invalidModeDelays = [1 1.5]; % a valid and an invalid delay
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, ...
                this.delays, observedModes, invalidModeDelays), expectedErrId);
            
            invalidModeDelays = [1 2 3]; % vector of delays is too large
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, ...
                this.delays, observedModes, invalidModeDelays), expectedErrId);
        end        
        
        %% testStepIssueWarningDiscardModeObservation
        function testStepIssueWarningDiscardModeObservation(this)
            % we just observe a mode, don't have a measurmeent at hand
            % mode observation is expected to be discarded with warning, but
            % result should be a successful prediction            
            expectedWarningId = 'Filter:IgnoringModeObservations';
            this.delayFilter.setState(this.stateGaussian)
            
            modeDelay = 100; % way too old
            modeObservation = 1;
            
            this.verifyWarning(@() this.delayFilter.step(this.delayFilterModel, this.measModel, [], ...
                [], modeObservation, modeDelay), expectedWarningId);
            
            % should reduce to a normal prediction of a Kalman filter
            expectedPredictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean;
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.W ...
                + this.B * this.inputCov * this.B';
            
            % the resulting mode distribution should be the stationary one
            [mode, probability] = this.delayFilter.getPreviousModeEstimate();        
            actualPredictedState = this.delayFilter.getState();
            
            this.verifyClass(actualPredictedState, ?Gaussian);
            
            [actualPredictedMean, actualPredictedCov] = actualPredictedState.getMeanAndCov();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualPredictedCov - this.stateGaussianCov), 0);
            
            this.verifyEqualWithAbsTol(actualPredictedMean, expectedPredictedMean);
            this.verifyEqualWithAbsTol(actualPredictedCov, expectedPredictedCov);
            this.verifyEqualWithAbsTol(actualPredictedCov, actualPredictedCov'); % should be symmetric
            
            this.verifyEqual(probability, this.inputProbs(1));
            this.verifyEqual(mode, 1);
        end        
%%
%%
        %% testPerformUpdateInvalidMeasModel
        function testPerformUpdateInvalidMeasModel(this)
            expectedErrId = 'Filter:UnsupportedMeasurementModel';
            
            this.delayFilter.setState(this.stateGaussian);
            % not a LinearMeasurementModel
            invalidMeasModel = this.delayFilterModel;
            this.verifyError(@() this.delayFilter.step(this.delayFilterModel, invalidMeasModel, this.delayedMeasurements, this.delays), expectedErrId);
        end
        
        %% testPerformUpdateNoApplicableMeas
        function testPerformUpdateNoApplicableMeas(this)                       
            this.delayFilter.setState(this.stateGaussian);
            % try to perform a measurement update with a measurement that is too old
            this.delayFilter.step(this.delayFilterModel, this.measModel, this.notApplicableMeasurements, this.notApplicableDelays);
            
            % should reduce to a normal prediction of a Kalman filter
            expectedPredictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean;
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.W ...
                + this.B * this.inputCov * this.B';
            
            % the resulting mode distribution should be the stationary one
            [mode, probability] = this.delayFilter.getPreviousModeEstimate();        
            actualPredictedState = this.delayFilter.getState();
            
            this.verifyClass(actualPredictedState, ?Gaussian);
            
            [actualPredictedMean, actualPredictedCov] = actualPredictedState.getMeanAndCov();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualPredictedCov - this.stateGaussianCov), 0);
            
            this.verifyEqualWithAbsTol(actualPredictedMean, expectedPredictedMean);
            this.verifyEqualWithAbsTol(actualPredictedCov, expectedPredictedCov);
            this.verifyEqualWithAbsTol(actualPredictedCov, actualPredictedCov'); % should be symmetric
            
            this.verifyEqual(probability, this.inputProbs(1));
            this.verifyEqual(mode, 1);
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.delayFilter.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 1);
        end
        
        %% testPerformUpdateZeroMaxMeasDelay
        function testPerformUpdateZeroMaxMeasDelay(this)
            % we expect that only the first measurement is processed
            this.zeroDelayFilter.setState(this.stateGaussian);
            this.zeroDelayFilter.step(this.zeroDelayFilterModel, this.measModel, this.delayedMeasurements, this.delays);
            
            % prediction
            expectedPredictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean;
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.G * this.W_red * this.G' ...
                + this.B * this.inputCov * this.B';
            
            % update should be a usual Kalman filter update step
            residualCov = this.C * expectedPredictedCov * this.C' + this.V;
            gain = expectedPredictedCov * this.C' * inv(residualCov);
            expectedMean = expectedPredictedMean + gain * (this.zeroDelayMeasurement - this.C * expectedPredictedMean);
            expectedCov = expectedPredictedCov - gain * residualCov * gain';
            
            % the resulting mode distribution should be the stationary one
            [mode, probability] = this.zeroDelayFilter.getPreviousModeEstimate();
            
            actualState = this.zeroDelayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
            
            [actualMean, actualCov] = actualState.getMeanAndCov();
            
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
            this.verifyEqualWithAbsTol(actualCov, actualCov'); % should be symmetric
            
            this.verifyEqual(probability, this.inputProbs(1));
            this.verifyEqual(mode, 1);
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.zeroDelayFilter.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 1);
            this.verifyEqual(actualNumDiscardedMeas, this.numMeas - 1);
        end
        
        %% testPerformUpdate
        function testPerformUpdate(this)
            % two measurements should be processed, one discarded
            this.delayFilter.setState(this.stateGaussian);
            this.delayFilter.step(this.delayFilterModel, this.measModel, this.delayedMeasurements, this.delays);
            
            [expectedUpdatedMean, expectedUpdatedCov] = this.computeUpdatedMoments();
            
            actualState = this.delayFilter.getState();
            % the resulting mode distribution should be the stationary one
            [mode, probability] = this.delayFilter.getPreviousModeEstimate();
            
            this.verifyClass(actualState, ?Gaussian);
            
            [actualMean, actualCov] = actualState.getMeanAndCov();
                      
            this.verifyEqualWithAbsTol(actualCov, actualCov'); % should be symmetric
            this.verifyEqualWithAbsTol(actualMean, expectedUpdatedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedUpdatedCov);            
            
            this.verifyEqual(probability, this.inputProbs(1));
            this.verifyEqual(mode, 1);
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.delayFilter.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, numel(this.applicableDelays));
            this.verifyEqual(actualNumDiscardedMeas, numel(this.notApplicableDelays));
        end        
%%
%%        
        %% testPerformPredictionInvalidSysModel
        function testPerformPredictionInvalidSysModel(this)
            expectedErrId = 'Filter:UnsupportedSystemModel';
            
            this.delayFilter.setState(this.stateGaussian);
            % not a LinearPlant
            invalidSysModel = InvertedPendulum(1, 1, 1, 1, 1);
            this.verifyError(@() this.delayFilter.step(invalidSysModel, this.measModel, this.delayedMeasurements, this.delays), expectedErrId);
        end
        
        %% testPerformPrediction
        function testPerformPrediction(this)
            % should reduce to a normal prediction of a Kalman filter
            expectedPredictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean;
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.W ...
                + this.B * this.inputCov * this.B';
            
            % set the state and perform a prediction
            this.delayFilter.setState(this.stateGaussian);
            this.delayFilter.step(this.delayFilterModel, this.measModel, [], []);
            actualPredictedState = this.delayFilter.getState();
            % the resulting mode distribution should be the stationary one
            [mode, probability] = this.delayFilter.getPreviousModeEstimate();
            
            this.verifyClass(actualPredictedState, ?Gaussian);
            
            [actualPredictedMean, actualPredictedCov] = actualPredictedState.getMeanAndCov();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualPredictedCov - this.stateGaussianCov), 0);
            
            this.verifyEqualWithAbsTol(actualPredictedMean, expectedPredictedMean);
            this.verifyEqualWithAbsTol(actualPredictedCov, expectedPredictedCov);
            this.verifyEqualWithAbsTol(actualPredictedCov, actualPredictedCov'); % should be symmetric
            
            this.verifyEqual(probability, this.inputProbs(1));
            this.verifyEqual(mode, 1);
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.delayFilter.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testPerformPredictionWithModeObservations
        function testPerformPredictionWithModeObservations(this)
            % init the filter and proceed one step
            this.delayFilter.setState(this.stateGaussian);
            this.delayFilter.step(this.delayFilterModel, this.measModel, [], []);
            
            modeObservations = [1 2 2 1];
            modeDelays = [4 2 1 3];
            expectedPreviousTrueMode = 2;
     
            % set the state and perform a prediction
            this.delayFilter.setState(this.stateGaussian);
            this.delayFilter.step(this.delayFilterModel, this.measModel, [], [], modeObservations, modeDelays);
            
            % should reduce to a normal prediction of a Kalman filter
            % (without input uncertainty)
            expectedPredictedMean = this.A  * this.stateGaussianMean ...
                + this.B * this.uncertainInputs(:, expectedPreviousTrueMode);
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.W;                
            
            actualPredictedState = this.delayFilter.getState();
            % the resulting mode distribution should be the a unit vector
            [mode, probability] = this.delayFilter.getPreviousModeEstimate();
            
            this.verifyClass(actualPredictedState, ?Gaussian);
            
            [actualPredictedMean, actualPredictedCov] = actualPredictedState.getMeanAndCov();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualPredictedCov - this.stateGaussianCov), 0);
            
            this.verifyEqualWithAbsTol(actualPredictedMean, expectedPredictedMean);
            this.verifyEqualWithAbsTol(actualPredictedCov, expectedPredictedCov);
            this.verifyEqualWithAbsTol(actualPredictedCov, actualPredictedCov'); % should be symmetric
            
            this.verifyEqual(probability, 1);
            this.verifyEqual(mode, expectedPreviousTrueMode);
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.delayFilter.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
            
            % finally, check if the current mode is correctly predicted
            % the resulting mode distribution should be the stationary one
            [predictedMode, predictedModeProb] = this.delayFilter.getModeEstimate();
            this.verifyEqual(predictedModeProb, this.inputProbs(1));
            this.verifyEqual(predictedMode, 1);
        end
        
        %% testPerformPredictionZeroMaxMeasDelay
        function testPerformPredictionZeroMaxMeasDelay(this)
            % should reduce to a normal prediction of a Kalman filter
            expectedPredictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean;
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.G * this.W_red * this.G' ...
                + this.B * this.inputCov * this.B';
            
            this.zeroDelayFilter.setState(this.stateGaussian);
            this.zeroDelayFilter.step(this.zeroDelayFilterModel, this.measModel, [], []);
            actualPredictedState = this.zeroDelayFilter.getState();
             % the resulting mode distribution should be the stationary one
            [mode, probability] = this.zeroDelayFilter.getPreviousModeEstimate();
            
            this.verifyClass(actualPredictedState, ?Gaussian);
            
            [actualPredictedMean, actualPredictedCov] = actualPredictedState.getMeanAndCov();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualPredictedCov - this.stateGaussianCov), 0);
            
            this.verifyEqualWithAbsTol(actualPredictedMean, expectedPredictedMean);
            this.verifyEqualWithAbsTol(actualPredictedCov, expectedPredictedCov);
            this.verifyEqualWithAbsTol(actualPredictedCov, actualPredictedCov'); % should be symmetric
            
            this.verifyEqual(probability, this.inputProbs(1));
            this.verifyEqual(mode, 1);
            
            % finally check if the update data was stored correctly
            [actualNumUsedMeas, actualNumDiscardedMeas] = this.zeroDelayFilter.getLastUpdateMeasurementData();
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
    end
end

