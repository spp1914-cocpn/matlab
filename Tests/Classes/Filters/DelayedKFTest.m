classdef DelayedKFTest < matlab.unittest.TestCase
    % Test cases for DelayedKF.
    
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
        
        uncertainInputs;
        inputMean;
        inputCov;
        
        sysNoiseMean;
        sysNoise;
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
    end
    
    methods (TestMethodSetup)
         %% initProperties
        function initProperties(this)
            this.maxMeasDelay = 4;
            this.numPossibleInputs = 5;
            % all inputs can be "active" with equal probability
            this.inputProbs = ones(1, this.numPossibleInputs) / this.numPossibleInputs;
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
            this.sysNoiseMean = ones(this.dimX, 1) * 0.75; % sys noise is not zero-mean
            this.sysNoise = Gaussian(this.sysNoiseMean, this.W);
            
            this.V = gallery('minij', this.dimY); 
                            
            this.measModel = LinearMeasurementModel(this.C);
            this.measNoise = Gaussian(zeros(this.dimY, 1), this.V);
            this.measModel.setNoise(this.measNoise);
                   
            this.delayFilterModel = DelayedKFSystemModel(this.A, this.B, this.sysNoise, ...
                this.numPossibleInputs, this.maxMeasDelay, this.inputProbs);
            this.zeroDelayFilterModel = DelayedKFSystemModel(this.A, this.B, this.sysNoise, ...
                this.numPossibleInputs, DelayedKFTest.zeroMeasDelay, this.inputProbs);
            
            this.delayFilterModel.setSystemInput(this.uncertainInputs);
            this.zeroDelayFilterModel.setSystemInput(this.uncertainInputs);
            
            this.delayFilter = DelayedKF(this.maxMeasDelay, this.filtername);
            this.zeroDelayFilter = DelayedKF(DelayedKFTest.zeroMeasDelay, this.zeroDelayFiltername);
        end
    end
    
    methods (Access = private)
        %% initMeasurements
        function initMeasurements(this)
            this.numMeas = 3;
            this.delayedMeasurements = ones(this.dimY, this.numMeas);
            this.delays = [0 1 5];
            % extract the expected measurements, i.e, those with delay <= maxMeasDelay
            this.applicableDelays = this.delays([1 2]);
            this.applicableMeasurements = this.delayedMeasurements(:, [1 2]);
            
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
            augMean = repmat(this.stateGaussianMean, this.maxMeasDelay + 1, 1);
            augCov = repmat(this.stateGaussianCov, this.maxMeasDelay + 1, this.maxMeasDelay + 1);
            % process the delayed measurements one after another
            expectedUpdatedMean = this.stateGaussianMean;
            expectedUpdatedCov = augCov;
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
        %% testDelayedKF
        function testDelayedKFNoName(this)
            filter = DelayedKF(this.maxMeasDelay);
            
            this.verifyEqual(filter.getName(), DelayedKFTest.defaultFiltername);
            this.verifyTrue(filter.getUseAnalyticMeasurementModel());
            
        end
        
        %% testDelayedKFZeroMaxMeasDelay
        function testDelayedKFZeroMaxMeasDelay(this)
            filter = DelayedKF(DelayedKFTest.zeroMeasDelay, this.zeroDelayFiltername);
            
            this.verifyEqual(filter.getMaxMeasurementDelay(), DelayedKFTest.zeroMeasDelay);
            this.verifyEqual(filter.getName(), this.zeroDelayFiltername);
            this.verifyTrue(filter.getUseAnalyticMeasurementModel());
        end
        
        %% testDelayedKF
        function testDelayedKF(this)
            filter = DelayedKF(this.maxMeasDelay, this.filtername);
            
            this.verifyEqual(filter.getMaxMeasurementDelay(), this.maxMeasDelay);
            this.verifyEqual(filter.getName(), this.filtername);
            this.verifyTrue(filter.getUseAnalyticMeasurementModel());
        end
        
        %% testSetStateInvalidState
        function testSetStateInvalidState(this)
            expectedErrId = 'Filter:UnsupportedSystemState'; % issued from super class GaussianFilter
            
            invalidState = this; % not a distribution
            this.verifyError(@() this.delayFilter.setState(invalidState), expectedErrId);
        end
        
        %% testSetStateGaussian
        function testSetStateGaussian(this)
            this.delayFilter.setState(this.stateGaussian);
            actualState = this.delayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
    
            [actualMean, actualCov] = actualState.getMeanAndCovariance();
               
            this.verifyEqualWithAbsTol(actualMean, this.stateGaussianMean);
            this.verifyEqualWithAbsTol(actualCov, this.stateGaussianCov);
            
        end
        
        %% testSetStateGaussianMixture
        function testSetStateGaussianMixture(this)
            % use a Distribution that is not Gaussian
            this.delayFilter.setState(this.stateGaussianMixture);
            actualState = this.delayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
            
            [actualMean, actualCov] = actualState.getMeanAndCovariance();
                        
            this.verifyEqualWithAbsTol(actualMean, this.stateGaussianMixtureMean);
            this.verifyEqualWithAbsTol(actualCov, this.stateGaussianMixtureCov);
            
        end
        
        %% testGetLastUpdateData
        function testGetLastUpdateData(this)
            % this method is not implemented
            expectedErrId = 'Filter:NotImplemented';
            
            this.verifyError(@() this.delayFilter.getLastUpdateData(), expectedErrId);
        end
        
        %% testPerformUpdateInvalidMeasModel
        function testPerformUpdateInvalidMeasModel(this)
            expectedErrId = 'Filter:UnsupportedMeasurementModel';
            
            this.delayFilter.setState(this.stateGaussian);
            % not a LinearMeasurementModel
            invalidMeasModel = this.delayFilterModel;
            this.verifyError(@() this.delayFilter.update(invalidMeasModel, this.delayedMeasurements), expectedErrId);
        end
        
        %% testPerformUpdateNoApplicableMeas
        function testPerformUpdateNoApplicableMeas(this)
                       
            this.delayFilter.setState(this.stateGaussian);
            % try to perform a measurement update with measurements that are too old
            this.delayFilter.update(this.measModel, this.notApplicableMeasurements, this.notApplicableDelays);
            
            % no applicable measurements, thus no update is carried out
            expectedMean = this.stateGaussianMean;
            expectedCov = this.stateGaussianCov;
            
            actualState = this.delayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
            
            [actualMean, actualCov] = actualState.getMeanAndCovariance();
            
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
            this.verifyEqualWithAbsTol(actualCov, actualCov'); % should be symmetric
        end
        
        %% testPerformUpdateZeroMaxMeasDelay
        function testPerformUpdateZeroMaxMeasDelay(this)
            % we expect that only the first measurement is processed
            this.zeroDelayFilter.setState(this.stateGaussian);
            this.zeroDelayFilter.update(this.measModel, this.delayedMeasurements, this.delays);
            
            % update should be a usual Kalman filter update step
            residualCov = this.C * this.stateGaussianCov * this.C' + this.V;
            gain = this.stateGaussianCov * this.C' * inv(residualCov);
            expectedMean = this.stateGaussianMean + gain * (this.zeroDelayMeasurement - this.C * this.stateGaussianMean);
            expectedCov = this.stateGaussianCov - gain * residualCov * gain';
            
            actualState = this.zeroDelayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
            
            [actualMean, actualCov] = actualState.getMeanAndCovariance();

            % first, a sanity check: posterior cov must be smaller or equal than prior
            % cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThanOrEqual(eig(this.stateGaussianMixtureCov - actualCov), 0);
            
            this.verifyEqualWithAbsTol(actualMean, expectedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedCov);
            this.verifyEqualWithAbsTol(actualCov, actualCov'); % should be symmetric
        end
        
        %% testPerformUpdate
        function testPerformUpdate(this)
            % two measurements should be processed
            this.delayFilter.setState(this.stateGaussian);
            this.delayFilter.update(this.measModel, this.delayedMeasurements, this.delays);
            
            [expectedUpdatedMean, expectedUpdatedCov] = this.computeUpdatedMoments();
            
            actualState = this.delayFilter.getState();
            
            this.verifyClass(actualState, ?Gaussian);
            
            [actualMean, actualCov] = actualState.getMeanAndCovariance();
                      
            this.verifyEqualWithAbsTol(actualCov, actualCov'); % should be symmetric
            this.verifyEqualWithAbsTol(actualMean, expectedUpdatedMean);
            this.verifyEqualWithAbsTol(actualCov, expectedUpdatedCov);
        end
  
%%
%%        
        %% testPerformPredictionInvalidSysModel
        function testPerformPredictionInvalidSysModel(this)
            expectedErrId = 'Filter:UnsupportedSystemModel';
            
            this.delayFilter.setState(this.stateGaussian);
            % not a DelayedKFSystemModel
            invalidSysModel = LinearPlant(this.A, this.B, this.W);
            this.verifyError(@() this.delayFilter.predict(invalidSysModel), expectedErrId);
        end
        
        %% testPerformPrediction
        function testPerformPrediction(this)
             % should reduce to a normal prediction of a Kalman filter
            expectedPredictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean + this.sysNoiseMean;
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.W ...
                + this.B * this.inputCov * this.B';
            
            % set the state and perform a prediction
            this.delayFilter.setState(this.stateGaussian);
            this.delayFilter.predict(this.delayFilterModel);
            actualPredictedState = this.delayFilter.getState();
            
            this.verifyClass(actualPredictedState, ?Gaussian);
            
            [actualPredictedMean, actualPredictedCov] = actualPredictedState.getMeanAndCovariance();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualPredictedCov - this.stateGaussianCov), 0);
            
            this.verifyEqualWithAbsTol(actualPredictedMean, expectedPredictedMean);
            this.verifyEqualWithAbsTol(actualPredictedCov, expectedPredictedCov);
            this.verifyEqualWithAbsTol(actualPredictedCov, actualPredictedCov'); % should be symmetric
        end
        
        %% testPerformPredictionZeroMaxMeasDelay
        function testPerformPredictionZeroMaxMeasDelay(this)
            % should reduce to a normal prediction of a Kalman filter
            expectedPredictedMean = this.A  * this.stateGaussianMean + this.B * this.inputMean + this.sysNoiseMean;
            expectedPredictedCov = this.A * this.stateGaussianCov * this.A' + this.W ...
                + this.B * this.inputCov * this.B';
            
            this.zeroDelayFilter.setState(this.stateGaussian);
            this.zeroDelayFilter.predict(this.zeroDelayFilterModel);
            actualPredictedState = this.zeroDelayFilter.getState();
            
            this.verifyClass(actualPredictedState, ?Gaussian);
            
            [actualPredictedMean, actualPredictedCov] = actualPredictedState.getMeanAndCovariance();
            
            % sanity check: resuting cov must be larger
            % than initial cov -> check the eigenvalues of the difference matrix
            this.verifyGreaterThan(eig(actualPredictedCov - this.stateGaussianCov), 0);
            
            this.verifyEqualWithAbsTol(actualPredictedMean, expectedPredictedMean);
            this.verifyEqualWithAbsTol(actualPredictedCov, expectedPredictedCov);
            this.verifyEqualWithAbsTol(actualPredictedCov, actualPredictedCov'); % should be symmetric
        end
    end
end

