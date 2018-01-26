classdef DelayedModeIMMF < DelayedMeasurementsFilter
    % Implementation of an Interacting Multiple Model (IMM) Filter which can handle delayed and out-of-sequence 
    % measurements and mode observations to estimate the plant state of an NCS given
    % in terms of a Markov Jump Linear System (MJLS).
    %
    % Literature: 
    %  	Florian Rosenthal, Benjamin Noack, and Uwe D. Hanebeck,
    %   State Estimation in Networked Control Systems With Delayed And Lossy Acknowledgments,
    %   Proceedings of the 2017 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017),
    %   Daegu, Korea, November 2017.
    
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
    
    properties (SetAccess = immutable, GetAccess = private)
       % the underlying IMM filter instance (IMMF instance)
        immf;
    end
    
    properties (Access = private)
        % cell array (row vector like) to store the last filter states
        % (cell array of Gaussian mixtures)
        % i-th element indicates the posterior at time k-i, where k is
        % the current time
        stateHistory;
        % cell array where each entry is a matrix of measurements (might be empty)
        % i-th element indicates the measurements taken at time k-i for update at k-i, where k is
        % the current time
        measurementHistory;
        % cell array where each entry is either a scalar or an empty matrix
        % i-th element indicates the true mode at time k-i for estimation
        % at k-i, where k is the current time
        trueModeHistory;
        % cell array where each entry is a matrix of inputs applied to the system in each
        % mode (might be empty)
        % inputHistory{i} contains the inputs used at time k-i where k is
        % the current time
        inputHistory;
    end
    
    methods (Access = public)
        %% DelayedModeIMMF
        function this = DelayedModeIMMF(modeFilters, modeTransitionMatrix, maxMeasDelay, name)
            % Class constructor.
            %
            % Parameters:
            %   >> modeFilters (FilterSet containing KF subclasses)
            %      A FilterSet consisting of Kalman filters, one for each
            %      mode of the system.
            %
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            %      In particular, the matrix must be square with all
            %      elements in [0,1] and rows summing up to 1.
            %
            %   >> maxMeasDelay (nonnegative integer)
            %      The maximum allowed measurement delay. That is, all
            %      measurements to be processed with a larger delay are discarded by the
            %      filter. If 0 is passed, the filter reduces to a default
            %      IMM filter which cannot cope with out-of-sequence and
            %      delayed measurements.
            %
            %   >> name (Char)
            %      An appropriate filter name / description of the implemented
            %      filter.
            %      Default name: 'IMM (with Delayed Measurements)'.
            %
            % Returns:
            %   << this (DelayedModeIMMF)
            %      A new DelayedModeIMMF instance.
            
            if nargin == 3
                name = 'IMM (with Delayed Measurements and Modes)';
            end
            this@DelayedMeasurementsFilter(maxMeasDelay, name);
                 
            this.immf = IMMF(modeFilters, modeTransitionMatrix);
            
            this.stateHistory = cell(1, maxMeasDelay + 1);
            this.measurementHistory = cell(1, maxMeasDelay);
            this.inputHistory = cell(1, maxMeasDelay + 1);
            this.trueModeHistory = cell(1, maxMeasDelay + 1);
        end
        
        %% reset
        function reset(this)
            %  Reset the filter by clearing the maintained state,
            %  measurement, input and system mode histories.
            %
            this.stateHistory = cell(1, this.maxMeasurementDelay + 1);
            this.measurementHistory = cell(1, this.maxMeasurementDelay);
            this.inputHistory = cell(1, this.maxMeasurementDelay + 1);
            this.trueModeHistory = cell(1, this.maxMeasurementDelay + 1); 
        end
        
        %% setState
        function setState(this, state)
            % Set the system state.
            %
            % This function is mainly used to set an initial system state, as
            % it is intended that the filter is responsible for modifying the
            % system state by exploiting system and measurement models.
            %
            % Parameters:
            %   >> state (Subclass of Distribution)
            %      The new system state.
            %      If a GaussianMixture is passed, the i-th component is used to as the state
            %      of the filter conditioned on the i-th mode. Likewise, the
            %      weight of the i-th mixture component is taken as probability
            %      of being in the i-th mode.
            %      In case of any other Distribution subclass, mean and
            %      covariance are used for all mode-conditioned filters, and a
            %      uniform distribution for the modes is employed.
            this.immf.setState(state);
            % set history to new state
           [this.stateHistory{:}] = deal(this.getState());
        end
        
        %% getState
        function state = getState(this)
            % Get the current system state.
            %
            % Returns:
            %   << state (Gaussian Mixture)
            %      The current system state.
            state = this.immf.getState();
        end
        
        %% getPointEstimate
        function [pointEstimate, uncertainty] = getPointEstimate(this)
            % Get a point estimate of the current system state.
            %
            % Returns:
            %   << pointEstimate (Column vector)
            %      Point estimate of the current system state which is simply
            %      the mean of the underlying Gaussian mixture.
            %
            %   << uncertainty (Positive definite matrix)
            %      Uncertainty of the current system state point estimate
            %      (covariance matrix of the underlying Gaussian mixture).
            [pointEstimate, uncertainty] = this.immf.getPointEstimate();
        end
        
        %% getModeEstimate
        function [mode, probability] = getModeEstimate(this)
            % Get a point estimate of the current system mode.
            %
            % Returns:
            %   << mode (Positive Integer)
            %      Estimate of the current system mode which is the maximum value of the current
            %      mode distribution.
            %
            %   << probability (Nonnegative Scalar)
            %      Probability of the maximum mode.
            %
            [mode, probability] = this.immf.getModeEstimate();
        end
        
        %% getPreviousModeEstimate
        function [mode, probability] = getPreviousModeEstimate(this)
            % Get a point estimate of the previous system mode.
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
            mixture = this.stateHistory{1};
            [~, ~, modeProbs] = mixture.getComponents();
            [probability, mode] = max(modeProbs);
        end
        
        %% setModeTransitionMatrix
        function setModeTransitionMatrix(this, modeTransitionMatrix)
            % Set the mode transition matrix.
            %
            % Parameters:
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            %      In particular, the matrix must be square with all
            %      elements in [0,1] and rows summing up to 1.
            this.immf.setModeTransitionMatrix(modeTransitionMatrix);
        end
        
        %% getModeTransitionMatrix
        function modeTransitionMatrix = getModeTransitionMatrix(this)
            % Get the mode transition matrix.
            %
            % Parameters:
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            
            modeTransitionMatrix = this.immf.modeTransitionProbs;
        end
        
        %% predict
        function runtime = predict(this, ~)
             this.error('InvalidPredictionStep', ...
                '** Time updates can only be performed in conjunction with a measurement update, so use step() directly');
        end
        
        %% update
        function runtime = update(this, ~, ~, ~)
             this.error('InvalidUpdateStep', ...
                '** Measurement updates can only be performed in conjunction with a time update, so use step() directly **'); 
        end
        
        %% step
        function runtime = step(this, sysModel, measModel, measurements, delays, modeMeas, modeDelays)
            % Perform a combined time and measurement update.
            %
            % Parameters:
            %   >> sysModel (JumpLinearSystemModel)
            %      Jump linear system model that provides the mapping between the prior system
            %      state and the predicted state (i.e., the system state's temporal evolution).
            %
            %   >> measModel (LinearMeasurementModel, or a cell array of LinearMeasurementModel)
            %      The linear measurement model that is used for all modes,
            %      or a cell array consisting of the mode-dependent linear
            %      measurement models.
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
            if nargin == 7
                this.checkModeMeasurementsAndDelays(modeMeas, modeDelays);
                [applicableModes, applicableModeDelays] = this.getApplicableModeMeasurements(modeMeas, modeDelays);
            else
                applicableModes = [];
                applicableModeDelays = [];
            end
            [applicableMeasurements, applicableDelays] = this.getApplicableMeasurements(measurements, measurementDelays);
            
            if nargout == 1
                s = tic;
                this.tryPerformStep(sysModel, measModel, applicableMeasurements, applicableDelays, applicableModes, applicableModeDelays); 
                runtime = toc(s);
            else
                this.tryPerformStep(sysModel, measModel, applicableMeasurements, applicableDelays, applicableModes, applicableModeDelays);
            end
            % save the update data
            numUsedMeas = size(applicableMeasurements, 2);
            this.setLastUpdateMeasurementData(numUsedMeas, size(measurements, 2) - numUsedMeas);
        end
    end
    
    methods (Access = protected)
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
                    && any(arrayfun(@(measDelay) mod(measDelay, 1) == 0, measDelays))
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
                if ~Checks.isVec(trueModes) || any(arrayfun(@(mode) ~Checks.isScalarIn(mode, 1, this.immf.numModes) ...
                        || mod(mode, 1) ~= 0, trueModes))
                     this.error('InvalidModeObservations', ...
                        '** Mode observations must be a given as a vector of integers from {1, ..., %d} **', ...
                        this.immf.numModes);
                end
                if ~Checks.isVec(modeDelays, numel(trueModes)) ...
                        || any(arrayfun(@(x) ~Checks.isNonNegativeScalar(x) || mod(x, 1) ~= 0, modeDelays))
                     this.error('InvalidModeDelays', ...
                        '** Mode observation delays must be a given as a vector of %d nonnegative integers **', ...
                        numel(trueModes));
                end
            end
        end
        
        %% performPrediction
        function performPrediction(this, ~)
            this.error('InvalidPredictionStep', ...
                '** Time updates can only be performed in conjunction with a measurement update, so use step() directly');
        end
        
        %% performUpdate
        function performUpdate(this, ~, ~,~)
             this.error('InvalidUpdateStep', ...
                '** Measurement updates can only be performed in conjunction with a time update, so use step() directly **');
        end
        
        %% performStep
        function performStep(this, sysModel, measModel, applicableMeasurements, applicableDelays, applicableModes, applicableModeDelays)
            if ~Checks.isClass(sysModel, 'JumpLinearSystemModel')
                this.errorSysModel('JumpLinearSystemModel');
            end
            
            maxDelay = max([applicableDelays applicableModeDelays]);
            if isempty(maxDelay)
                % this should only happen in case both
                % <applicableDelays> and <applicableModeDelays> are
                % empty
                maxDelay = 0;
            end
            currentMeasurements = this.incorporateDelayedMeasurements(applicableMeasurements, applicableDelays);
            currentMode = this.incorporateDelayedModeMeasurements(applicableModes, applicableModeDelays);
            % remember the current inputs
            currentInputs = [];%this.getSystemInputs(sysModel);
            this.inputHistory{1} = this.getSystemInputs(sysModel);
            % handle border case that maxDelay is this.maxMeasurementDelay + 1 
            % (maximum allowed delay for a mode observation)
            delays = min(maxDelay, this.maxMeasurementDelay):-1:0;
            for i=delays %i = maxDelay:-1:0 using this expression directly does not seem to work correctly when compiled 
                % check what filter to use, depending on whether true
                % mode i+1 steps before (i.e., the previous mode) is known or not
                trueMode = this.trueModeHistory{i + 1};
                inputs = this.inputHistory{i + 1}; % inputs for prediction

                if ~isempty(trueMode)
                    % we also have to update the posterior (mode probs) at this time
                    [modeStateMeans, modeStateCovs, ~] = this.stateHistory{i + 1}.getComponents();
                    weights = [zeros(1, trueMode - 1), 1 zeros(1, this.immf.numModes - trueMode)];
                    newPosterior = GaussianMixture(modeStateMeans, modeStateCovs, weights);
                    this.stateHistory{i + 1} = newPosterior;
                end
                if i ~= 0
                    measurements = this.measurementHistory{i};
                else
                    measurements = currentMeasurements;
                end
                this.immf.setState(this.stateHistory{i + 1});
                DelayedModeIMMF.applyInputs(sysModel, inputs);
                this.immf.performPrediction(sysModel);
                if ~isempty(measurements)
                    this.immf.performUpdate(measModel, measurements);
                end
                if i ~= 0
                    % save posterior
                    this.stateHistory{i} = this.immf.getState();
                end
            end                   
            % update history for the next step (i.e., proceed to k+1)
            this.updateHistory(currentInputs, currentMeasurements, currentMode);
        end
    end
    
    methods (Access = private)
        %% incorporateDelayedMeasurements
        function nonDelayedMeasurements = incorporateDelayedMeasurements(this, applicableMeasurements, applicableDelays)
            nonDelayedMeasurements = [];
            if ~isempty(applicableDelays)
                % get the available measurement delays
                occuringDelays = unique(applicableDelays); % sorted in ascending order (lowest delay first)
                % get the measurements per delay as a cell array of matrices
                % (lowest delay first); grouping by operation (column-wise)
                measurementsPerDelay = splitapply(@(group) {group}, applicableMeasurements, grp2idx(applicableDelays)');
                % update the history of measurements
                startIndex = 1;
                if occuringDelays(1) == 0
                    % measurements that are not delayed (i.e., from the current time step)
                    nonDelayedMeasurements = measurementsPerDelay{1};
                    startIndex = 2;
                end
                for i = startIndex:numel(occuringDelays)
                    % delay of j time steps-> store at j-th element of the
                    % history
                    this.measurementHistory{occuringDelays(i)} = measurementsPerDelay{i};
                end
            end
        end
        
        %% incorporateDelayedModeMeasurements
        function currentMode = incorporateDelayedModeMeasurements(this, applicableModes, applicableModeDelays)
            currentMode = [];
            if ~isempty(applicableModeDelays)
                [contained, idx] = ismember(0, applicableModeDelays);
                if contained
                    currentMode = applicableModes(idx);
                    indices = [1:idx-1 idx+1:numel(applicableModeDelays)]; 
                else
                    indices = 1:numel(applicableModeDelays);
                end
                for i = indices
                    % delay of j time steps-> store at j-th element of the history
                    if ~isempty(this.trueModeHistory{applicableModeDelays(i)})
                        this.warning('OverwritingModeObservation', ...
                            '** Overwrite existing mode observation: Known mode %d, provided mode %d. **', ...
                            this.trueModeHistory{applicableModeDelays(i)}, applicableModes(i));
                    end
                    this.trueModeHistory{applicableModeDelays(i)} = applicableModes(i);
                   
%                     if isempty(this.trueModeHistory{applicableModeDelays(i)}) ...
%                             || applicableModes(i) < this.trueModeHistory{applicableModeDelays(i)}
%                         this.trueModeHistory{applicableModeDelays(i)} = applicableModes(i);
%                     elseif ~isempty(this.trueModeHistory{applicableModeDelays(i)}) ...
%                             && applicableModes(i) >= this.trueModeHistory{applicableModeDelays(i)}
%                          this.warning('IgnoringModeObservations', ...
%                             '** Ignoring mode observation: Known mode %d, given mode %d. **', ...
%                             this.trueModeHistory{applicableModeDelays(i)}, applicableModes(i));
%                     end
                end
            end
        end
                      
        %% getApplicableModeMeasurements
        function [applicableModes, applicableModeDelays] = getApplicableModeMeasurements(this, trueModes, modeDelays)
            % find the mode observations with valid delays
            idx = find(modeDelays <= this.maxMeasurementDelay + 1);
            applicableModes = trueModes(idx);
            applicableModeDelays = modeDelays(idx);
            if numel(idx) ~= numel(modeDelays)
                this.warning('IgnoringModeObservations', ...
                    '** %s\nIgnoring %d of %d mode observations. **', 'Delay too large', ...
                    numel(modeDelays) - numel(idx), numel(modeDelays));
            end
        end
        
        %% tryPerformStep
        function tryPerformStep(this, sysModel, measModel, applicableMeasurements, applicableDelays, applicableModes, applicableModeDelays)
            try
                this.performStep(sysModel, measModel, applicableMeasurements, applicableDelays, applicableModes, applicableModeDelays);
            catch ex
                 % copied from Filter
                if ~strcmp(ex.identifier, 'Filter:IgnoreMeasurement') ...
                        && ~strcmp(ex.identifier, 'Filter:IgnorePrediction')
                    % real error => do not catch it
                   ex.rethrow();
                end
           end    
        end
        
        %% updateHistory
        function updateHistory(this, currentInputs, currentMeasurements, currentTrueMode)
            if nargin == 2
                currentMeasurements = [];
                currentTrueMode = [];
            elseif nargin == 3
                currentTrueMode = [];
            end
            % update state, input and measurement history for the next step
            this.stateHistory = circshift(this.stateHistory, 1, 2);
            this.stateHistory{1} = this.immf.getState();
            
            this.measurementHistory = circshift(this.measurementHistory, 1, 2);
            this.measurementHistory{1} = currentMeasurements;
            
            this.inputHistory = circshift(this.inputHistory, 1, 2);
            this.inputHistory{1} = currentInputs;
           
            this.trueModeHistory = circshift(this.trueModeHistory, 1, 2);
            this.trueModeHistory{1} = currentTrueMode;
            
            % if called at time k, the history is prepared for k+1
        end
        
        %% getSystemInputs
        function usedInputs = getSystemInputs(this, sysModel)
            % assume that sysModel is a JumpLinearSystemModel
            modeSpecificInputs = sysModel.getSystemInput();
            usedInputs{this.immf.numModes} = [];
            if ~isempty(modeSpecificInputs)
                 usedInputs = modeSpecificInputs;
            end    
        end
              
    end
    
    methods (Access = private, Static)
        %% applyInputs
        function applyInputs(sysModel, inputs, mode)
            % assume that sysModel is a JumpLinearSystemModel and inputs is a cell array
            if nargin == 2
                % i-th input to be applied to i-th system model
                sysModel.setSystemInput(inputs);
            else
                % assume mode is in bounds
                % only apply input for given mode
                if isempty(inputs)
                    input = [];
                else
                    input = inputs(:, mode);
                end
                sysModel.modeSystemModels{mode}.setSystemInput(input);
            end
        end
    end
end

