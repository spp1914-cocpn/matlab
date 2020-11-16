classdef (Sealed) DelayedKF < DelayedMeasurementsFilter
    % This class represents a Kalman filter that can operate with time delayed control inputs 
    % and/or time delayed measurements. The filter calculates the linear minimum
    % mean square state estimate of a system based on a model of that system,
    % the expected control input applied to the system and possibly
    % time-delayed measurements of the state.
    %
    % This implementation is based on the original one by Jörg Fischer and Maxim Dolgov.
    %
    % Literature: 
    %   Maryam Moayedi, Yung Kuan Foo and Yeng Chai Soha
    %   Filtering for networked control systems with single/multiple measurement packets
    %   subject to multiple-step measurement delays and multiple packet dropouts
    %   International Journal of Systems Science (2011)
    
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
    
    properties (Constant, Access = private)
        modeTransitionMatrixHistoryLength = 5; % educated guess
        % thus, we always store T_{k-1}, T_{k-2}, T_{k-3}, T_{k-4}, T_{k-5}
        % that are used to predict theta_{k}, theta_{k-1}, theta_{k-2},
        % theta_{k-3} and theta_{k-4}, respectively
    end
    
    properties (Access = private)
        augmentedStateMean; 
        augmentedStateCov; % needed as augmented cov not necessarily pd
        modeTransitionMatrix(:,:) double {Validator.validateTransitionMatrix(modeTransitionMatrix)} = eye(1);  % stores T_k
        previousModeEstimate; % estimate of the mode at time k-1, i.e., probability distribution for theta_{k-1} (stored as column vec)
        lastModeObservationDelay; % store the time delay to the last mode observation
        % 3D matrix, where each slice is a mode transition matrix
        % i-th slice indicates the mode transition matrix of the MJLS at
        % time k-i, where k is the current time
        % so the first entry is T_{k-1} and so forth
        modeTransitionMatrixHistory;
    end
    
    properties (GetAccess = private, SetAccess = immutable)
        numModes;
    end
    
    methods (Access = public)
        %% DelayedKF
        function this = DelayedKF(maxMeasDelay, modeTransitionMatrix, name)
            % Class constructor.
            %
            % Parameters:
            %   >> maxMeasDelay (nonnegative integer)
            %      The maximum allowed measurement delay. That is, all
            %      measurements to be processed with a larger delay are discarded by the
            %      filter. If 0 is passed, the filter reduces to a
            %      filter which cannot cope with out-of-sequence and
            %      delayed measurements.
            %
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            %      In particular, the matrix must be square with all
            %      elements in [0,1] and rows summing up to 1.
            %
            %   >> name (Char)
            %      An appropriate filter name / description of the implemented
            %      filter.
            %      Default name: 'Delayed KF'.
            
            % Returns:
            %   << this (DelayedKF)
            %      A new DelayedKF instance.
            if nargin == 2
                name = 'Delayed KF';
            end
            this@DelayedMeasurementsFilter(maxMeasDelay, name);
                       
            this.modeTransitionMatrix = modeTransitionMatrix;
            this.numModes = size(this.modeTransitionMatrix, 1); 
            
            this.previousModeEstimate = [zeros(this.numModes - 1, 1); 1]; % we assume that we know the initial mode
            this.lastModeObservationDelay = 1;
            this.modeTransitionMatrixHistory = repmat(modeTransitionMatrix, 1, 1, DelayedKF.modeTransitionMatrixHistoryLength); % fixed size by default  
        end
       
        %% getState
        function state = getState(this)
            % Get the system state.
            %
            % Returns:
            %   << state (Gaussian)
            %      The system state.
            
            state = Gaussian(this.augmentedStateMean(1:this.dimState), ...
                this.augmentedStateCov(1:this.dimState, 1:this.dimState));
        end
        
        %% getStateMeanAndCov
        function [stateMean, stateCov, stateCovSqrt] = getStateMeanAndCov(this)
            % Get mean and covariance matrix of the system state.
            %
            % Returns:
            %   << stateMean (Column vector)
            %      Mean vector of the system state.
            %
            %   << stateCov (Positive definite matrix)
            %      Covariance matrix of the system state.
            %
            %   << stateCovSqrt (Square matrix, optional)
            %      Lower Cholesky decomposition of the system state covariance matrix.
            
            stateMean = this.augmentedStateMean(1:this.dimState);
            stateCov = this.augmentedStateCov(1:this.dimState, 1:this.dimState);
            if nargout == 3
                stateCovSqrt = chol(stateCov, 'lower');
            end
        end
        
        %% getModeEstimate
        function [mode, probability] = getModeEstimate(this)
            % Get a point estimate of the current system mode, i.e., theta_{k}.
            %
            % Returns:
            %   << mode (Positive Integer)
            %      Estimate of the current system mode which is the maximum value of the current
            %      mode distribution.
            %
            %   << probability (Nonnegative Scalar)
            %      Probability of the maximum mode.
            %
            
            % use T_{k-1} for prediction
            [probability, mode] = max(this.modeTransitionMatrixHistory(:, :, 1)' * this.previousModeEstimate);
        end
        
        %% getPreviousModeEstimate
        function [mode, probability] = getPreviousModeEstimate(this)
            % Get a point estimate of the previous system mode, i.e.,
            % theta_{k-1}.                        
            %
            % Returns:
            %   << mode (Positive Integer)
            %      Estimate of the previous system mode which is the
            %      maximum value of the previous
            %      mode distribution.
            %
            %   << probability (Nonnegative Scalar)
            %      Probability of the maximum mode.
            %            
            
            [probability, mode] = max(this.previousModeEstimate);
        end
        
        %% setModeTransitionMatrix
        function setModeTransitionMatrix(this, modeTransitionMatrix)
            % Set the mode transition matrix used to be used for the prediction from k to k+1,
            % where k is the current time step.
            %
            % Parameters:
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            %      In particular, the matrix must be square with all
            %      elements in [0,1] and rows summing up to 1.
            
            Validator.validateTransitionMatrix(modeTransitionMatrix, this.numModes);
            this.modeTransitionMatrix = modeTransitionMatrix;
        end
        
        %% predict
        function runtime = predict(this, ~)
             this.error('InvalidPredictionStep', ...
                '** Time updates can only be performed in conjunction with a measurement update, so use step() directly **');
        end
        
        %% update
        function runtime = update(this, ~, ~, ~)
             this.error('InvalidUpdateStep', ...
                '** Measurement updates can only be performed in conjunction with a time update, so use step() directly **'); 
        end
        
        %% step
        function runtime = step(this, sysModel, measModel, measurements, delays, modeMeas, modeDelays)
            % Perform a combined time and measurement update, i.e., when
            % called at time k, the estimate is propagated from x^e_{k-1}
            % to x^e_{k}.
            %
            % Parameters:
            %   >> sysModel (DelayedKFSystemModel)
            %      Augmented system model that provides the mapping between the prior system
            %      state and the predicted state (i.e., the system state's temporal evolution).
            %
            %   >> measModel (LinearMeasurementModel)
            %      The linear measurement model.
            %
            %   >> measurements (Matrix, might be empty)
            %      Column-wise arranged measurement vectors, where each column represents an individual
            %      (possibly delayed) measurement. An empty matrix
            %      indicates that no measurements are available.
            %
            %   >> delays (Nonnegative integer or vector of nonnegative integers (one for each measurement))
            %      If a vector is passed here, the i-th element denotes the delay of the i-th measurement.
            %      If a single nonnegative integer is passed instead, its
            %      value is taken as delay for all measurements.
            %
            %   >> modeMeas (Vector, might be empty)
            %      The individual (possibly delayed) mode observations to
            %      be processed. An empty matrix
            %      indicates that no mode observations are available.
            %
            %   >> modeDelays (Vector of nonnegative integers (one for each mode observation))
            %      Nonnegative vector where the i-th element denotes the delay of the i-th mode observation.
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the combined time and measurement update.
            
            if nargin ~= 5 && nargin ~= 7
                this.error('InvalidStep:InvalidNumArgs', ...
                    '** Number of arguments must be 5 or 7, i.e., only <modeMeas> and <modeDelays> are optional **'); 
            end
            measurementDelays = this.checkMeasurementsAndDelays(measurements, delays);
            [applicableMeasurements, applicableDelays] = this.getApplicableMeasurements(measurements, measurementDelays);
            if nargin == 7
                this.checkModeMeasurementsAndDelays(modeMeas, modeDelays);
                [applicableModes, applicableModeDelays] = getApplicableModeMeasurements(this, modeMeas, modeDelays);
                % we only need the most recent mode observation (if any)
                [mostRecentModeDelay, idx] = min(applicableModeDelays);
                mostRecentMode = applicableModes(idx);
            else
                mostRecentModeDelay = [];
                mostRecentMode = [];
            end
                            
            s = tic;
            try
                this.performPrediction(sysModel, mostRecentMode, mostRecentModeDelay);
                
                if ~isempty(applicableMeasurements)
                    this.performUpdate(measModel, applicableMeasurements, applicableDelays);
                end
            catch ex
                % copied from Filter
                if ~strcmp(ex.identifier, 'Filter:IgnoreMeasurement') ...
                    && ~strcmp(ex.identifier, 'Filter:IgnorePrediction')
                    % real error => do not catch it
                    ex.rethrow();
                end
            end
            if nargout == 1
                runtime = toc(s);
            end            
            % save the update data
            numUsedMeas = size(applicableMeasurements, 2);
            this.setLastUpdateMeasurementData(numUsedMeas, size(measurements, 2) - numUsedMeas);
            
            % prepare for the next time step (k+1)
             this.lastModeObservationDelay = this.lastModeObservationDelay + 1;
             this.modeTransitionMatrixHistory = cat(3, this.modeTransitionMatrix, ...
                this.modeTransitionMatrixHistory(:, :, 1:DelayedKF.modeTransitionMatrixHistoryLength-1));
        end
    end

    methods(Access = protected)
        %% performSetState
        function performSetState(this, state)
            [mean, cov] = state.getMeanAndCov();
            
            % augment state
            this.augmentedStateMean = repmat(mean, this.maxMeasurementDelay + 1, 1);
            % this matrix (covariance of augmented state) is not pd
            this.augmentedStateCov = repmat(cov, this.maxMeasurementDelay + 1, this.maxMeasurementDelay + 1);
        end
        
        %% performSetStateMeanAndCov
        function performSetStateMeanAndCov(this, stateMean, stateCov, ~)
            % augment state
            this.augmentedStateMean = repmat(stateMean, this.maxMeasurementDelay + 1, 1);
            % this matrix (covariance of augmented state) is not pd
            this.augmentedStateCov = repmat(stateCov, this.maxMeasurementDelay + 1, this.maxMeasurementDelay + 1);
        end
        
        %% performUpdate
        function performUpdate(this, measModel, applicableMeasurements, applicableDelays)
            if ~Checks.isClass(measModel, 'LinearMeasurementModel')
                this.errorMeasModel('LinearMeasurementModel');
            end
            if ~isempty(applicableDelays)
                % get the measurements per delay as a cell array of matrices
                % grouping by operation (column-wise)
                [groups, ~, timeDelays] = grp2idx(applicableDelays);
                measurementsPerDelay = splitapply(@(group) {group}, applicableMeasurements, groups');
                % now perform the analytical update: incorporate all
                % measurements per delay at once                
                arrayfun(@(i) this.doUpdate(measModel, measurementsPerDelay{i}, timeDelays(i)), 1:numel(timeDelays));
           end
        end
    
        %% performPrediction
        function performPrediction(this, sysModel, mostRecentMode, mostRecentModeDelay)
            if ~Checks.isClass(sysModel, 'DelayedKFSystemModel')
                this.errorSysModel('DelayedKFSystemModel');
            end
            % integrate the mode observation to update the mode distribution theta_{k-1}, if present
            if ~isempty(mostRecentMode) && mostRecentModeDelay < this.lastModeObservationDelay
                % predict from theta_{k-j}, where j = mostRecentModeDelay,
                % to theta_{k-1} using the corresponding transition matrices T_{k-j}, ..., T_{k-2}
                % then, simply extract the corresponding row as mode probabilities are unit vector
                %probMat = this.modeTransitionMatrix ^ (mostRecentModeDelay-1);
                probMat = eye(this.numModes);
                for j=mostRecentModeDelay:-1:2
                    probMat = probMat * this.modeTransitionMatrixHistory(:, :, j);
                end
                this.previousModeEstimate = probMat(mostRecentMode, :)';
                this.lastModeObservationDelay = mostRecentModeDelay;
            else
                % we simply propagate the previous mode estimate theta_{k-2}, using
                % T_{k-2}
                this.previousModeEstimate = this.modeTransitionMatrixHistory(:, :, 2)' * this.previousModeEstimate;
                %this.previousModeEstimate = this.modeTransitionMatrix' * this.previousModeEstimate;
            end
            % these are the new delay weights used for the prediction from x_{k-1} to x_{k}
            sysModel.setDelayWeights(this.previousModeEstimate);
            
            % do not compute square root of state cov as it is not
            % neccessarily positive-definite
            [predictedAugmentedStateMean, ...
                predictedAugmentedStateCov] = sysModel.analyticMoments(this.augmentedStateMean, ...
                    this.augmentedStateCov);
            
            % check if predicted state covariance is valid
            this.checkCovPrediction(predictedAugmentedStateCov(1:this.dimState, 1:this.dimState), ...
                   'Predicted state');

            % all checks successful, store augmented estimates
            this.augmentedStateMean = predictedAugmentedStateMean;
            this.augmentedStateCov = predictedAugmentedStateCov;
        end             
    end
    
    methods (Access = private)
        %% doUpdate
        function doUpdate(this, measModel, measurements, delay)
            [dimMeas, numMeas] = size(measurements);
            baseMeasMatrix = measModel.measMatrix;
           
            [noiseMean, noiseCov] = measModel.noise.getMeanAndCov();
            dimStackedMeas = size(baseMeasMatrix, 1) * numMeas; % what is expected according to the original model
                       
            augmentedMeasMatrix = [zeros(dimMeas, this.dimState * delay), ...
                    baseMeasMatrix, zeros(dimMeas, this.dimState * (this.maxMeasurementDelay - delay))];
            
            % compute required moments
            measMean = repmat(augmentedMeasMatrix * this.augmentedStateMean + noiseMean, numMeas, 1);
            measCov = Utils.blockDiag(noiseCov, numMeas) ...
                + repmat(augmentedMeasMatrix * this.augmentedStateCov * augmentedMeasMatrix', numMeas, numMeas);
            % cross-covariance of original (non-augmented) state and
            % measurements
            stateMeasCrossCov = repmat(this.augmentedStateCov(1:this.dimState, :) * augmentedMeasMatrix', 1, numMeas);
                  
            % Check measurement moments
            if ~Checks.isColVec(measMean, dimStackedMeas) || any(~isfinite(measMean))
                this.error('InvalidMeasurementMean', ...
                          ['** Measurement mean must be a ' ...
                           'real-valued column vector of dimension %d **'], ...
                          dimStackedMeas);
            end
            if ~Checks.isSquareMat(measCov, dimStackedMeas) || any(~isfinite(measCov(:)))
                this.error('InvalidMeasurementCovariance', ...
                          ['** Measurement covariance must be a ' ...
                           'positive definite matrix of dimension %dx%d **'], ...
                          dimStackedMeas, dimStackedMeas);
            end
            if ~Checks.isMat(stateMeasCrossCov, this.dimState, dimStackedMeas) || any(~isfinite(stateMeasCrossCov(:)))
                this.error('InvalidStateMeasurementCrossCovariance', ...
                          ['** State measurement cross-covariance must be a ' ...
                           'matrix of dimension %dx%d **'], ...
                          this.dimState, dimStackedMeas);
            end
            [sqrtMeasCov, isNonPos] = chol(measCov);
            
            if isNonPos
                this.error('InvalidMeasurementCovariance', ...
                      'Measurement covariance matrix is not positive definite.');
            end
                      
            meanGain = (stateMeasCrossCov / sqrtMeasCov) / sqrtMeasCov'; % L
            covGain = [meanGain; ...
                    zeros(this.dimState * this.maxMeasurementDelay, dimStackedMeas)]; % F
            % only update first component of augmented state, i.e. the current state
            % the remaining components remain unaffected by the update
            postMean = this.augmentedStateMean(1:this.dimState) + meanGain * (measurements(:) - measMean);
                        
            % not necessarily positive definite posterior cov of augmented state
            factor = speye(this.dimState * (this.maxMeasurementDelay + 1)) ...
                - covGain * repmat(augmentedMeasMatrix, numMeas, 1);
            augmentedPostCov = factor * this.augmentedStateCov * factor' ...
                + covGain * Utils.blockDiag(noiseCov, numMeas) * covGain';
                        
            this.checkCovUpdate(augmentedPostCov(1:this.dimState, 1:this.dimState), ...
                                                     'Updated state');
            % all checks successful, store augmented estimates
            this.augmentedStateMean(1:this.dimState) = postMean;
            this.augmentedStateCov = augmentedPostCov;
        end        
        
        %% checkMeasurementsAndDelays
        function delays = checkMeasurementsAndDelays(this, measurements, measDelays)
            if ~ismatrix(measurements)
                this.error('InvalidMeasurements', ...
                    '** measurements must be given as a (possibly empty) matrix **');
            end
            numMeas = size(measurements, 2);
            delays = [];            
            if numMeas > 0
                if Checks.isNonNegativeScalar(measDelays) && mod(measDelays, 1) == 0
                    delays = repmat(measDelays, 1, numMeas);
                elseif Checks.isNonNegativeVec(measDelays, numMeas) ...
                    && all(arrayfun(@(measDelay) mod(measDelay, 1) == 0, measDelays))
                    delays = measDelays;
                else
                    this.error('InvalidMeasDelay', ...
                        '** Each measurement delay must be a nonnegative integer **');
                end
            end
        end
        
        %% checkModeMeasurementsAndDelays
        function checkModeMeasurementsAndDelays(this, trueModes, modeDelays)
            if ~isempty(trueModes)
                if ~Checks.isVec(trueModes) || any(arrayfun(@(mode) ~Checks.isScalarIn(mode, 1, this.numModes) ...
                        || mod(mode, 1) ~= 0, trueModes))
                     this.error('InvalidModeObservations', ...
                        '** Mode observations must be a given as a vector of integers from {1, ..., %d} **', ...
                        this.numModes);
                end
                if ~Checks.isVec(modeDelays, numel(trueModes)) ...
                        || any(arrayfun(@(x) ~Checks.isNonNegativeScalar(x) || mod(x, 1) ~= 0, modeDelays))
                     this.error('InvalidModeDelays', ...
                        '** Mode observation delays must be a given as a vector of %d nonnegative integers **', ...
                        numel(trueModes));
                end
            end
        end
        
        %% getApplicableModeMeasurements
        function [applicableModes, applicableModeDelays] = getApplicableModeMeasurements(this, trueModes, modeDelays)
            % find the mode observations with valid delays
            idx = find(modeDelays <= DelayedKF.modeTransitionMatrixHistoryLength);
            applicableModes = trueModes(idx);
            applicableModeDelays = modeDelays(idx);
            if numel(idx) ~= numel(modeDelays)
                this.warning('IgnoringModeObservations', ...
                    '** %s\nIgnoring %d of %d mode observations. **', 'Delay too large', ...
                    numel(modeDelays) - numel(idx), numel(modeDelays));
            end
        end 
    end
end
