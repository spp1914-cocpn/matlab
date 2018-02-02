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
    
    properties (SetAccess = immutable, GetAccess = private)
        % filter set with mode-conditioned Kalman filters
        modeFilters@FilterSet;
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
        probabilityBound = 1e-12;
    end
    
    methods (Access = public)
        function this = IMMF(modeFilters, modeTransitionMatrix, name)
            % Class constructor.
            %
            % Parameters:
            %   >> modeFilters (FilterSet containing KF subclasses)
            %      A FilterSet consisting of Kalman filters, one for each
            %      mode of the system.
            %   >> modeTransitionMatrix (Matrix)
            %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
            %      In particular, the matrix must be square with all
            %      elements in [0,1] and rows summing up to 1.
            %   >> name (Char)
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
            this.numModes = modeFilters.getNumFilters();
            if any(arrayfun(@(i) ~Checks.isClass(modeFilters.get(i), 'KF'), 1:this.numModes))
                this.error('InvalidModeFilters:InvalidFilterType', ...
                     '** Each mode-conditioned filter must be a Kalman filter **');
            end
            Validator.validateTransitionMatrix(modeTransitionMatrix, this.numModes);
            this.modeTransitionProbs = modeTransitionMatrix;
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
            this.modeTransitionProbs = modeTransitionMatrix;
        end
               
        function setState(this, state)
            % Set the system state.
            %
            % This function is mainly used to set an initial system state, as
            % it is intended that the IMM filter is responsible for modifying the
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
            if Checks.isClass(state, 'GaussianMixture')
                if state.getNumComponents() ~= this.numModes
                    this.error('InvalidStateGaussianMixture' ,...
                        '** Gaussian mixture is expected to have %d components **', this.numModes);
                end
                this.dimState = state.getDimension();
                [modeStateMeans, modeStateCovs, modeProbs] = state.getComponents();
                this.modeProbabilities = this.normalizeModeProbabilities(modeProbs);
                for i=1:this.numModes
                    this.modeFilters.get(i).setState(Gaussian(modeStateMeans(:, i), modeStateCovs(:, :, i)));
                end
            elseif Checks.isClass(state, 'Distribution')
                [mean, cov] = state.getMeanAndCovariance();
                this.dimState = numel(mean);
                % no prior knowledge on modes available, so assume uniform distribution
                this.modeProbabilities = this.normalizeModeProbabilities(repmat(1 / this.numModes, 1, this.numModes));
                this.modeFilters.setStates(Gaussian(mean, cov));
            else
                this.error('InvalidState' , ...
                    '** State must either be a Gaussian mixture with %d components or any arbitrary distribution **', ...
                    this.numModes);
            end
        end
                
        function state = getState(this)
            % Get the current system state.
            %
            % Returns:
            %   << state (Gaussian Mixture)
            %      The current system state.
            [modeStateMeans, modeStateCovs] = this.getModeStateMeansAndCovs();
            state = GaussianMixture(modeStateMeans, modeStateCovs, this.modeProbabilities);
        end
                
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
            [modeStateMeans, modeStateCovs] = this.getModeStateMeansAndCovs();
            [pointEstimate, covMeans] = Utils.getMeanAndCov(modeStateMeans, this.modeProbabilities);
            weightedStateCovs = bsxfun(@times, modeStateCovs, reshape(this.modeProbabilities, [1 1 this.numModes]));
            uncertainty = covMeans + sum(weightedStateCovs, 3);
        end
        
        %% getModeEstimate
        function [mode, probability] = getModeEstimate(this)
            % Get a point estimate of the current system mode.
            % This is simply the maximum value of the current mode
            % probability distribution, and its associated probability.
            % in case of multiple maxima, the first one is returned.
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
        function performUpdate(this, measModels, measurements)
            if (iscell(measModels))
                if numel(measModels) ~= this.numModes
                    this.issueErrorMeasModel();
                end
                arrayfun(@(index) this.modeFilters.updateSingle(index, measModels{index}, measurements), 1:this.numModes);
            else
                this.modeFilters.update(measModels, measurements);
            end
            this.updateModeProbabilities();
        end
        
        function performPrediction(this, sysModel)
            if Checks.isClass(sysModel, 'JumpLinearSystemModel')
                modeDependentModels = sysModel.modeSystemModels;
                if numel(modeDependentModels) ~= this.numModes
                    this.issueErrorSysModel();
                end    
                this.performStateEstimateMixing();
                arrayfun(@(index) this.modeFilters.predictSingle(index, modeDependentModels{index}), 1:this.numModes)
            elseif Checks.isClass(sysModel, 'SystemModel')
                % use the given model for all filters
                this.performStateEstimateMixing();
                this.modeFilters.predict(sysModel);
            else
                this.issueErrorSysModel();
            end
        end
    end
   
    methods (Access = private)
        function issueErrorMeasModel(this)
            this.errorMeasModel('MeasurementModel', ...
                        sprintf('Cell array of %d MeasurementModels', this.numModes));
        end
        
        function issueErrorSysModel(this)
            this.errorSysModel(sprintf('JumpLinearSystemModel (%d modes)', this.numModes), ...
                    'SystemModel (for all modes)');
        end
        
        function performStateEstimateMixing(this)
            oldModeProbabilities = this.modeProbabilities;
            this.predictModeProbabilities();
            this.mixStateEstimates(oldModeProbabilities);
        end
        
        function predictModeProbabilities(this)
            %yields a column vector, so transpose the result
            this.modeProbabilities = reshape(this.modeTransitionProbs' * this.modeProbabilities', 1, this.numModes);
            this.modeProbabilities = this.normalizeModeProbabilities(this.modeProbabilities);
        end
        
        function mixStateEstimates(this, oldModeProbabilities)
            [modeStateMeans, modeStateCovs] = this.getModeStateMeansAndCovs();
            % multiply (j,i)-th element of transition matrix by old probability of mode j
            mixingWeights = this.modeTransitionProbs .* oldModeProbabilities';
            % divide (j,i)-th element of mixing weights by predicted (new) probability of mode i
            mixingWeights = mixingWeights ./ this.modeProbabilities;
            % mix all estimates
            % resulting state (per mode) is a Gaussian approximation of a
            % Gaussian mixture
            [mixedMeans, mixedMeansCovs] = arrayfun(@(col) Utils.getMeanAndCov(modeStateMeans, mixingWeights(:, col)'), ...
                1:this.numModes, 'UniformOutput', false);
                
            weightedCovs = arrayfun(@(col) bsxfun(@times, modeStateCovs, reshape(mixingWeights(:, col)', [1 1 this.numModes])), ...
                1:this.numModes, 'UniformOutput', false);
            mixedCovs = cell2mat(mixedMeansCovs) + sum(cell2mat(weightedCovs), 3);

            for i=1:this.numModes
                covIdx = this.dimState * (i - 1) + 1;
                %initialize filter with mixed estimate
                this.modeFilters.get(i).setState(Gaussian(mixedMeans{i}, ...
                    mixedCovs(:, covIdx:(covIdx + this.dimState - 1))));
            end
        end
        
        function [modeStateMeans, modeStateCovs] = getModeStateMeansAndCovs(this)
            [means, covs] = cellfun(@(dist) dist.getMeanAndCovariance(), this.modeFilters.getStates(), ...
                'UniformOutput', false);
            modeStateMeans = cell2mat(means);
            modeStateCovs = reshape(cell2mat(covs), [this.dimState, this.dimState, this.numModes]);
        end
        
        function updateModeProbabilities(this)
            logLikelihoods = zeros(1, this.numModes);
            for i=1:this.numModes
                [measurement, measurementMean, measurementCovariance, ~, ~] = this.modeFilters.get(i).getLastUpdateData();
                % assume Gaussian likelihood
                dist = Gaussian(measurementMean, measurementCovariance);
                logLikelihoods(i) = dist.logPdf(measurement);
            end
            this.checkLogLikelihoodEvaluations(logLikelihoods, this.numModes);
            % for numerical stability
            maxLogValue = max(logLikelihoods);
            logLikelihoods = logLikelihoods - maxLogValue;
            
            likelihoods = exp(logLikelihoods);
            updatedModeProbs = this.modeProbabilities .* likelihoods;
            normalizationConstant = sum(updatedModeProbs);
            this.modeProbabilities = this.normalizeModeProbabilities(updatedModeProbs / normalizationConstant);
        end
        
        function normalizedModeProbs = normalizeModeProbabilities(this, modeProbs)
            idx = find(modeProbs <= IMMF.probabilityBound);
            normalizedModeProbs = modeProbs;
            if ~isempty(idx)
                %this.warning('NormalizingModeProbabilities', ...
                %    '** %d of %d mode probablities too small. Perform normalization **', numel(idx), this.numModes); 
                normalizedModeProbs(idx) = IMMF.probabilityBound;
                normalizedModeProbs = normalizedModeProbs / sum(normalizedModeProbs);
            end
        end
    end
end

