classdef (Sealed) DelayedIMMF < DelayedMeasurementsFilter
    % Implementation of an Interacting Multiple Model (IMM) Filter which can handle delayed and out-of-sequence 
    % measurements to estimate the state of a Markov Jump Linear System (MJLS) in terms of
    % a Gaussian Mixture distribution.
    % The filter was developed by Fischer et al. in the context of sequence-based
    % networked control systems where both control inputs and measurements
    % are transmitted via networks subject to random packet delay and/or loss.
    %
    % Literature: 
    %   Jörg Fischer, Achim Hekler, and Uwe D. Hanebeck,
    %   State Estimation in Networked Control Systems,
    %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
    %   Singapore, July 2012.
    %
    % While this implementation is to some extent more general than the original one of Fischer et al., 
    % which can be found <a href="matlab:
    % web('http://www.cloudrunner.eu/algorithm/60/interacting-multiple-model-imm-filter-for-networked-control-systems-ncs/version/2/')"
    % >here</a>, some parts are directly based thereof.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2016  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (SetAccess = immutable, GetAccess = private)
       % the underlying IMM filter instance (IMMF instance)
        immf;
    end
    
    properties (Access = private)
        % cell array (row vector like) to store the last filter states
        % (cell array of Gaussian mixtures)
        % i-th element indicates the posterior at time k-(i-1), where k is
        % the current time
        % in particular, stateHistory{1} contains the current posterior,
        % i.e., equal to immf.getState()
        stateHistory;
        % cell array where each entry is a matrix of measurements (might be empty)
        % i-th element indicates the measurements taken at time k-(i-1) for update at k-(i-1), where k is
        % the current time
        measurementHistory;
        % cell array where each entry is a cell array of inputs applied to the system in each
        % mode (might be empty)
        % i-th element indicates the inputs used at time k-(i-1) for prediction at k-(i-1), where k is
        % the current time
        inputHistory;
    end
    
    methods (Access = public)
        function this = DelayedIMMF(modeFilters, modeTransitionMatrix, maxMeasDelay, name)
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
            %   >> name (Char)
            %      An appropriate filter name / description of the implemented
            %      filter.
            %      Default name: 'IMM (with Delayed Measurements)'.
            %
            % Returns:
            %   << this (DelayedIMMF)
            %      A new DelayedIMMF instance.
            
            if nargin == 3
                name = 'IMM (with Delayed Measurements)';
            end
            this@DelayedMeasurementsFilter(maxMeasDelay, name);
            
            this.immf = IMMF(modeFilters, modeTransitionMatrix);
            this.stateHistory = cell(1, maxMeasDelay + 1);
            this.measurementHistory = cell(1, maxMeasDelay);
            this.inputHistory = cell(1, maxMeasDelay);
        end
                
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
    end
    
    methods (Access = protected)
        %% performSetState
        function performSetState(this, state)
            this.immf.setState(state);
            % set history to new state
            [this.stateHistory{:}] = deal(this.getState());
        end
        
        %% performSetStateMeanAndCov
        function performSetStateMeanAndCov(this, stateMean, stateCov, stateCovSqrt)
             this.immf.setStateMeanAndCov(stateMean, stateCov, stateCovSqrt);
        end
        
        function performPrediction(this, sysModel)
            this.immf.performPrediction(sysModel);
            % update history for the next step
            this.updateHistory(this.getSystemInputs(sysModel));
        end
        
        function performUpdate(this, ~, ~,~)
            this.error('InvalidUpdateStep', ...
                ['** Measurement updates can only be performed in conjunction with a time update, so use step() directly\n', ...
                'However, performing consecutive time updates is possible by calling predict() successively **']); 
        end
        
        function performStep(this, sysModel, measModels, applicableMeasurements, measDelays)
            % Perform a combined time and measurement update with delayed measurements.
            %
            % Parameters:
            %   >> sysModel (System model instance to be used for all modes 
            %      or JumpLinearSystemModel, which comprises individual models for each mode of the system)
            %      System model that provides the mode-specific mappings between the prior system
            %      state and the predicted state (i.e., the system state's temporal evolution).
            %
            %   >> measModels (Measurement model instance to be used for all modes
            %      or cell array (row-vector like) of measurement models, one for each mode of the system)
            %      Measurement model(s) that provide(s) the mode-specific mappings 
            %      between measurements and the system state.
            %
            %   >> applicableMeasurements (Matrix)
            %      Column-wise arranged measurement vectors, where each column represents an individual
            %      (possibly delayed) measurement.
            %      Empty in case no applicable measurement (i.e.,
            %      measurement with delay <= maxMeasDelay) present.
            %
            %   >> measDelays (Row vector of nonnegative integers (one for each measurement))
            %      The i-th element denotes the delay of the i-th measurement.
            %      Empty in case no applicable measurement (i.e.,
            %      measurement with delay <= maxMeasDelay) present.
            %
            if isempty(measDelays)
                % no measurements are applicable, delays too large, so just predict
                this.performPrediction(sysModel);
                return;
            end
            currentMeasurements = [];
            % remember the current inputs
            currentInputs = this.getSystemInputs(sysModel);
            % get the available delays
            occuringDelays = unique(measDelays); % sorted in ascending order (lowest delay first)
            % get the measurements per delay as a cell array of matrices
            % (lowest delay first); grouping by operation (column-wise)
            measurementsPerDelay = splitapply(@(group) {group}, applicableMeasurements, grp2idx(measDelays)');
            % update the history of measurements
            if occuringDelays(1) == 0
                % measurements that are not delayed (i.e., from the current time step)
                currentMeasurements = measurementsPerDelay{1};
            else
              this.measurementHistory{occuringDelays(1)} = measurementsPerDelay{1};
            end
            for i = 2:numel(occuringDelays)
                this.measurementHistory{occuringDelays(i)} = measurementsPerDelay{i};
            end
            maxDelay = max(occuringDelays);
            if maxDelay ~= 0
                % get the state from the history and reset the underlying IMM filter
                this.immf.setState(this.stateHistory{maxDelay + 1});
                % we have to incorporate delayed measurements
                idx = maxDelay:-1:1;
                for i=idx %i = maxDelay:-1:1 using this expression directly does not seem to work correctly when compiled 
                    this.applyInputs(sysModel, this.inputHistory{i});
                    this.immf.performPrediction(sysModel);
                    if ~isempty(this.measurementHistory{i})
                        this.immf.performUpdate(measModels, this.measurementHistory{i});
                    end
                    % save posterior at that time
                    this.stateHistory{i} = this.immf.getState();
                end
            end
            
            % regular step to estimate the current state
            % don't forget to reset the inputs first
            this.applyInputs(sysModel, currentInputs);
            this.immf.performPrediction(sysModel);
            if ~isempty(currentMeasurements)
                this.immf.performUpdate(measModels, currentMeasurements);
            end
                       
            % update history for the next step
            this.updateHistory(currentInputs, currentMeasurements);
        end
    end
    
    methods (Access = private)
        function applyInputs(this, sysModel, inputs)
            % i-th input to be applied to i-th system model (if jump linear
            % system model)
            if Checks.isClass(sysModel, 'JumpLinearSystemModel')
                sysModel.setSystemInput(cell2mat(inputs));
            elseif Checks.isClass(sysModel, 'LinearSystemModel')
                if numel(inputs) ~= 1
                    this.error('DelayedIMMF:InvalidPastSystemInput', ...
                        '** <inputs> is believed to contain 1 element only, for one linear system model is used for all modes **');
                end
                sysModel.setSystemInput(inputs{1});
            end
        end
        
        function usedInputs = getSystemInputs(this, sysModel)
            if Checks.isClass(sysModel, 'JumpLinearSystemModel')
                modeSpecificInputs = sysModel.getSystemInput();
                usedInputs{this.immf.numModes} = [];
                if ~isempty(modeSpecificInputs)
                    usedInputs = num2cell(modeSpecificInputs, [1 this.immf.numModes]);
                end    
            elseif Checks.isClass(sysModel, 'LinearSystemModel')
                % a single linear model for all modes is used
                usedInputs{1} = sysModels.getSystemInput();
            else
               usedInputs{this.immf.numModes} = [];
            end
        end
        
        function updateHistory(this, usedInputs, usedMeasurements)
            if nargin == 2
                usedMeasurements = [];
            end
            % update state, input and measurement history for the next step
            this.stateHistory = circshift(this.stateHistory, 1, 2);
            % posterior state at current time step, i.e. after
            % step/prediction was called
            this.stateHistory{1} = this.immf.getState();
            
            this.measurementHistory = circshift(this.measurementHistory, 1, 2);
            % measurements used for the update in the current step
            this.measurementHistory{1} = usedMeasurements;
            
            this.inputHistory = circshift(this.inputHistory, 1, 2);
            % inputs used for prediction in the current time step
            this.inputHistory{1} = usedInputs;
        end
    end
end

