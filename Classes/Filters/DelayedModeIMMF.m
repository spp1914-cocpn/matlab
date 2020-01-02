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
    %
    %  	Florian Rosenthal, Benjamin Noack, and Uwe D. Hanebeck,
    %   State Estimation in Networked Control Systems with Delayed and Lossy Acknowledgments,
    %   Multisensor Fusion and Integration in the Wake of Big Data, Deep Learning and Cyber Physical System,
    %   Lecture Notes in Electrical Engineering, Volume 501,
    %   Springer, Cham, 2018.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        % inputHistory{i} contains the inputs used at time k-i, where k is
        % the current time
        inputHistory;
        % 3D matrix, where each slice is a mode transition matrix
        % i-th slice indicates the mode transition matrix of the MJLS at
        % time k-i, where k is the current time
        modeTransitionMatrixHistory;
    end
    
    methods (Access = public)
        %% DelayedModeIMMF
        function this = DelayedModeIMMF(modeFilters, modeTransitionMatrix, maxMeasDelay, name)
            % Class constructor.
            %
            % Parameters:
            %   >> modeFilters (Cell array containing LinearGaussianFilter subclasses)
            %      A cell array consisting of Kalman filters, one for each
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
            %   >> name (Char, Optional)
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
            
            this.modeTransitionMatrixHistory = repmat(modeTransitionMatrix, 1, 1, maxMeasDelay + 1);            
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
                
        %% getState
        function state = getState(this)
            % Get the current system state.
            %
            % Returns:
            %   << state (Gaussian Mixture)
            %      The current system state.
            state = this.immf.getState();
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
            
            if nargout == 2
                [stateMean, stateCov] = this.immf.getStateMeanAndCov();
            else
                [stateMean, stateCov, stateCovSqrt] = this.immf.getStateMeanAndCov();
            end
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
        function [mode, probability] = getPreviousModeEstimate(this, afterStep)
            % Get a point estimate of the previous system mode.
            %
            % Parameters:
            %   >> afterStep (Logical Scalar, optional)
            %      A flag to indicate whether step(..) has been called beforehand, that is, 
            %      to indicate whether a time and measurement update has already been carried out.
            %      If left out, the default value true will be utilized.
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
            
            % we assume that function is called after a step (time + measurement update) has been
            % performed, if afterStep is not present
            if nargin == 1
                mixture = this.stateHistory{2};
            elseif Checks.isFlag(afterStep)
                if afterStep
                    mixture = this.stateHistory{2};
                else
                    mixture = this.stateHistory{1};
                end
            else
                this.error('GetPreviousModeEstimate:InvalidFlag', ...
                    '<afterStep> must be a flag (i.e., a logical scalar) or left out ** ');
            end
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
                '** Time updates can only be performed in conjunction with a measurement update, so use step() directly **');
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
                
            s = tic;
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
            if nargout == 1
                runtime = toc(s);
            end            
            % save the update data
            numUsedMeas = size(applicableMeasurements, 2);
            this.setLastUpdateMeasurementData(numUsedMeas, size(measurements, 2) - numUsedMeas);
        end     
    end
    
    methods (Access = protected)
        %% performSetState
        function performSetState(this, state)
            this.immf.setState(state);
            % set history to new state
            for j=1:this.maxMeasurementDelay + 1
                this.stateHistory{j} = this.immf.getState().copy();
            end            
        end
        
        %% performSetStateMeanAndCov
        function performSetStateMeanAndCov(this, stateMean, stateCov, stateCovSqrt)
            this.immf.setStateMeanAndCov(stateMean, stateCov, stateCovSqrt);
            
            % set history to new state
            for j=1:this.maxMeasurementDelay + 1
                this.stateHistory{j} = this.immf.getState().copy();
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
        function performStep(this, sysModel, measModel, applicableMeasurements, applicableDelays, ...
                applicableModes, applicableModeDelays)
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
            currentModeTransitionMatrix = this.immf.modeTransitionProbs; % we have to memorize the current matrix
            currentMeasurements = this.incorporateDelayedMeasurements(applicableMeasurements, applicableDelays);
            currentMode = this.incorporateDelayedModeMeasurements(applicableModes, applicableModeDelays);
            
            this.inputHistory{1} = this.getSystemInputs(sysModel);            
            % handle border case that maxDelay is this.maxMeasurementDelay + 1 
            % (maximum allowed delay for a mode observation)
            delays = min(maxDelay, this.maxMeasurementDelay):-1:0;
            for i=delays %i = maxDelay:-1:0 using this expression directly does not seem to work correctly when compiled 
                % check whether true mode i+1 steps before (i.e., the previous mode) is known or not
                trueMode = this.trueModeHistory{i + 1};
                if ~isempty(trueMode)
                    % we also have to update the posterior (mode probs) at this time
                    [modeStateMeans, modeStateCovs, ~] = this.stateHistory{i + 1}.getComponents();
                    weights = [zeros(1, trueMode - 1), 1 zeros(1, this.immf.numModes - trueMode)];

                    this.stateHistory{i+1}.set(modeStateMeans, modeStateCovs, weights);
                end
                if i ~= 0
                    measurements = this.measurementHistory{i};
                else
                    measurements = currentMeasurements;
                end
                this.immf.setState(this.stateHistory{i + 1});
                this.immf.setModeTransitionMatrix(this.modeTransitionMatrixHistory(:, :, i + 1));
                               
                DelayedModeIMMF.applyInputs(sysModel, this.inputHistory{i + 1}); % inputs for prediction
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
            this.updateHistory(currentMeasurements, currentMode, currentModeTransitionMatrix);            
        end
    end
    
    methods (Access = private)
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
        
        %% updateHistory
        function updateHistory(this, currentMeasurements, currentTrueMode, currentModeTransitionMatrix)
            % update state, input and measurement history for the next step
            this.stateHistory = circshift(this.stateHistory, 1, 2);
            this.stateHistory{1} = this.immf.getState();
            
            this.measurementHistory = circshift(this.measurementHistory, 1, 2);
            this.measurementHistory{1} = currentMeasurements;
            
            this.inputHistory = circshift(this.inputHistory, 1, 2);
            this.inputHistory{1} = []; % not known, as we need estimate for computation of inputs
           
            this.trueModeHistory = circshift(this.trueModeHistory, 1, 2);
            this.trueModeHistory{1} = currentTrueMode;
      
            this.modeTransitionMatrixHistory = cat(3, currentModeTransitionMatrix, ...
                this.modeTransitionMatrixHistory(:, :, 1:this.maxMeasurementDelay));
            
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
        function applyInputs(sysModel, inputs)
            % assume that sysModel is a JumpLinearSystemModel and inputs is
            % an array of inputs
            % so i-th input to be applied to i-th system model
            sysModel.setSystemInput(inputs);
        end
    end
end
