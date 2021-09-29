classdef IMMF < Filter
    % Implementation of an Interacting Multiple Model (IMM) Filter to
    % estimate the state of a Markov Jump Linear System (MJLS) in terms of
    % a Gaussian Mixture distribution.
    % The filter is suboptimal since the optimal filter is
    % infinite-dimensional and hence intractable.
    %
    % Literature: 
    %   Henk AP Blom, and Yaakov Bar-Shalom,
    %   The interacting multiple model algorithm for systems with Markovian switching coefficients,
    %   IEEE transactions on Automatic Control 33.8 (1988): 780-783.
    %
    %   X. Rong Li and Vesselin P. Jilkov,
    %   Survey of maneuvering target tracking - Part V: Multiple-Model Methods,
    %   IEEE Transactions on Aerospace and Electronic Systems 41.4 (2005): 1255-1321.
    
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
    
    properties (SetAccess = immutable, GetAccess = private)
        % filter set with mode-conditioned Kalman filters
        modeFilters(1,:) cell;
    end
    
    properties (SetAccess = immutable, GetAccess = public)
         numModes;
    end
    
    properties (Access = private)
        % row vector of dimension <numModes> containing the current mode
        % probabilities
        modeProbabilities;
    end
    
    properties (SetAccess = private, GetAccess = public)
        % <numModes> x <numModes> matrix where the (i,j)-th entry denotes
        % the probabitity of changing from mode i to j (i,j = 1, ..., numModes)
        modeTransitionProbs;
    end
    
    properties (Constant)
        probabilityBound = 1e-50;
    end
    
    methods (Access = public)
        function this = IMMF(modeFilters, modeTransitionMatrix, name)
            % Class constructor.
            %
            % Parameters:
            %   >> modeFilters (Cell array containing LinearGaussianFilter subclasses)
            %      A cell array consisting of Kalman filters, one for each
            %      mode of the system.
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            %      In particular, the matrix must be square with all
            %      elements in [0,1] and rows summing up to 1.
            %   >> name (Char, optional)
            %      An appropriate filter name / description of the implemented
            %      filter.
            %      Default name: 'Interacting Multiple Model Filter'.
            %
            % Returns:
            %   << this (IMMF)
            %      A new IMMF instance.
            
            if nargin == 2
                name = 'Interacting Multiple Model Filter';
            end
            this@Filter(name);

            this.modeFilters = modeFilters;
            this.numModes = numel(modeFilters);
            if any(cellfun(@(filter) ~Checks.isClass(filter, 'LinearGaussianFilter'), this.modeFilters))
                this.error('InvalidModeFilters:InvalidFilterType', ...
                     '** Each mode-conditioned filter must be a LinearGaussianFilter, i.e., a Kalman filter **');
            end
            Validator.validateTransitionMatrix(modeTransitionMatrix, this.numModes);
            this.modeTransitionProbs = this.normalizeModeTransitionMatrix(modeTransitionMatrix);
        end
    end
    
    methods (Access = public)
        
        %% setModeTransitionMatrix
        function setModeTransitionMatrix(this, modeTransitionMatrix)
            % Set the mode transition matrix.
            %
            % Parameters:
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            %      In particular, the matrix must be square with all
            %      elements in [0,1] and rows summing up to 1.
            %
            Validator.validateTransitionMatrix(modeTransitionMatrix, this.numModes);
            this.modeTransitionProbs = this.normalizeModeTransitionMatrix(modeTransitionMatrix);
        end
               
        %% getState        
        function state = getState(this)
            % Get the current system state.
            %
            % Returns:
            %   << state (Gaussian Mixture)
            %      The current system state.
            [modeStateMeans, modeStateCovs] = this.getModeStateMeansAndCovs();
            state = GaussianMixture(modeStateMeans, modeStateCovs, this.modeProbabilities);
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

            [modeStateMeans, modeStateCovs] = this.getModeStateMeansAndCovs();
            [stateMean, stateCov] = Utils.getGMMeanAndCov(modeStateMeans, modeStateCovs, this.modeProbabilities);
            if nargout == 3
                stateCovSqrt = chol(stateCov, 'Lower');
            end
        end
        
        %% getModeEstimate
        function [mode, probability] = getModeEstimate(this)
            % Get a point estimate of the current system mode.
            % This is simply the maximum value of the current mode
            % probability distribution, and its associated probability.
            % In case of multiple maxima, the first one is returned.
            %
            % Returns:
            %   << mode (Positive Integer)
            %      Estimate of the current system mode which is the index of the maximum value of the current
            %      mode distribution.
            %
            %   << probability (Nonnegative Scalar)
            %      Probability of the maximum mode.
            %
            [probability, mode] = max(this.modeProbabilities);
        end
    end
    
    methods (Access = protected)
        %% performSetState
        function performSetState(this, state)
            if Checks.isClass(state, 'GaussianMixture')
                if state.getNumComponents() ~= this.numModes
                    this.error('InvalidStateGaussianMixture' ,...
                        '** Gaussian mixture is expected to have %d components **', this.numModes);
                end
                [modeStateMeans, modeStateCovs, modeProbs] = state.getComponents();
                this.modeProbabilities = this.normalizeModeProbabilities(modeProbs);
                                
                for i=1:this.numModes
                    this.modeFilters{i}.setStateMeanAndCov(modeStateMeans(:, i), modeStateCovs(:, :, i));
                end
            else
                [mean, cov, covSqrt] = state.getMeanAndCov();
                this.performSetStateMeanAndCov(mean, cov, covSqrt);
            end
        end
        
        %% performSetStateMeanAndCov
        function performSetStateMeanAndCov(this, stateMean, stateCov, stateCovSqrt)
            % no prior knowledge on modes available, so assume uniform distribution
            this.modeProbabilities = this.normalizeModeProbabilities(repmat(1 / this.numModes, 1, this.numModes));
            cellfun(@(filter) filter.setStateMeanAndCov(stateMean, stateCov, stateCovSqrt), this.modeFilters);
        end
        
        %% performUpdate
        function performUpdate(this, measModels, measurement)
            % compute the individual model likelihoods
            if (iscell(measModels))
                if numel(measModels) ~= this.numModes
                    this.issueErrorMeasModel();
                end
                arrayfun(@(index) this.modeFilters{index}.update(measModels{index}, measurement), 1:this.numModes);
            else
                cellfun(@(filter) filter.update(measModels, measurement), this.modeFilters);
            end
            this.updateModeProbabilities(measurement);
        end
        
        %% performPrediction
        function performPrediction(this, sysModel)
            if Checks.isClass(sysModel, 'JumpLinearSystemModel')
                modeDependentModels = sysModel.modeSystemModels;
                if numel(modeDependentModels) ~= this.numModes
                    this.issueErrorSysModel();
                end    
                this.performStateEstimateMixing();
                arrayfun(@(index) this.modeFilters{index}.predict(modeDependentModels{index}), 1:this.numModes)
            elseif Checks.isClass(sysModel, 'SystemModel')
                % use the given model for all filters
                this.performStateEstimateMixing();
                cellfun(@(filter) filter.predict(sysModel), this.modeFilters)
            else
                this.issueErrorSysModel();
            end
        end
    end
   
    methods (Access = private)
        %% issueErrorMeasModel
        function issueErrorMeasModel(this)
            this.errorMeasModel('MeasurementModel', ...
                sprintf('Cell array of %d MeasurementModels', this.numModes));
        end
        
        %% issueErrorSysModel
        function issueErrorSysModel(this)
            this.errorSysModel(sprintf('JumpLinearSystemModel (%d modes)', this.numModes), ...
                'SystemModel (for all modes)');
        end
        
        %% performStateEstimateMixing
        function performStateEstimateMixing(this)
            oldModeProbabilities = this.modeProbabilities;
            this.predictModeProbabilities();
            this.mixStateEstimates(oldModeProbabilities);
        end
        
        %% predictModeProbabilities
        function predictModeProbabilities(this)
            %yields a column vector, so transpose the result
            this.modeProbabilities = this.normalizeModeProbabilities(...
                reshape(this.modeTransitionProbs' * this.modeProbabilities', 1, this.numModes));
        end
        
        %% mixStateEstimates
        function mixStateEstimates(this, oldModeProbabilities)
            [modeStateMeans, modeStateCovs] = this.getModeStateMeansAndCovs();
            % multiply (j,i)-th element of transition matrix by old probability of mode j
            mixingWeights = this.modeTransitionProbs .* oldModeProbabilities';
            % divide (j,i)-th element of mixing weights by predicted (new) probability of mode i
            mixingWeights = mixingWeights ./ this.modeProbabilities;
            
            % mix all estimates
            % resulting state (per mode) is a Gaussian approximation of a
            % Gaussian mixture
            for i=1:this.numModes
                [mean, cov] = Utils.getGMMeanAndCov(modeStateMeans, modeStateCovs, mixingWeights(:, i)');
                %initialize filter with mixed estimate
                this.modeFilters{i}.setStateMeanAndCov(mean, cov);
            end
        end
        
        %% getModeStateMeansAndCovs
        function [modeStateMeans, modeStateCovs] = getModeStateMeansAndCovs(this)
            [means, covs] = cellfun(@getStateMeanAndCov, this.modeFilters, 'UniformOutput', false);
            modeStateMeans = cell2mat(means);
            modeStateCovs = reshape(cell2mat(covs), this.dimState, this.dimState, this.numModes);
        end
        
        %% updateModeProbabilities
        function updateModeProbabilities(this, measurement)
            % assume Gaussian likelihoods
            logLikelihoods = zeros(1, this.numModes);
            dimMeas = size(measurement, 1);
            logNormConst = dimMeas * 0.5 * log(2 * pi);
            I = eye(dimMeas);
            for i=1:this.numModes                
                [~, measurementMean, measurementCovariance, ~] = this.modeFilters{i}.getLastUpdateData();     
                covSqrt = chol(measurementCovariance, 'lower');
                
                logLikelihoods(i) = -0.5 * sum((covSqrt \ (measurement - measurementMean)) .^2, 1) ...
                    - logNormConst - sum(log(diag(covSqrt)));  % this is the log det of covSqrt
            end
            
            this.checkLogLikelihoodEvaluations(logLikelihoods, this.numModes);
            % for numerical stability
            maxLogValue = max(logLikelihoods);
            logLikelihoods = logLikelihoods - maxLogValue;
            
            updatedModeProbs = this.modeProbabilities .* exp(logLikelihoods);
            this.modeProbabilities = this.normalizeModeProbabilities(updatedModeProbs / sum(updatedModeProbs));
        end
            
        %% normalizeModeProbabilities
        function normalizedModeProbs = normalizeModeProbabilities(this, modeProbs)
            normalizedModeProbs = Utility.normalizeProbabilities(modeProbs, IMMF.probabilityBound);
        end
        
        %% normalizeModeTransitionMatrix
        function normalizedTransitionMatrix = normalizeModeTransitionMatrix(this, modeTransitionMatrix)
            normalizedTransitionMatrix = Utility.normalizeTransitionMatrix(modeTransitionMatrix, IMMF.probabilityBound);
        end        
    end
end

