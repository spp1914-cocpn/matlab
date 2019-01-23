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
    %    Copyright (C) 2017-2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Access = private)
        augmentedStateMean; 
        augmentedStateCov; % needed as augmented cov not necessarily pd
    end
   
    methods (Access = public)
        %% DelayedKF
        function this = DelayedKF(maxMeasDelay, name)
            % Class constructor.
            %
            % Parameters:
            %   >> maxMeasDelay (nonnegative integer)
            %      The maximum allowed measurement delay. That is, all
            %      measurements to be processed with a larger delay are discarded by the
            %      filter. If 0 is passed, the filter reduces to a
            %      filter which cannot cope with out-of-sequence and
            %      delayed measurements.
            %   >> name (Char)
            %      An appropriate filter name / description of the implemented
            %      filter.
            %      Default name: 'Delayed KF'.
            
            % Returns:
            %   << this (DelayedKF)
            %      A new DelayedKF instance.
            if nargin == 1
                name = 'Delayed KF';
            end
            this@DelayedMeasurementsFilter(maxMeasDelay, name);
        end % function DelayedKF
       
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
        function performPrediction(this, sysModel)
             if ~Checks.isClass(sysModel, 'DelayedKFSystemModel')
                this.errorSysModel('DelayedKFSystemModel');
             end
            this.doPrediction(sysModel);
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
        
        %% doPrediction
        function doPrediction(this, augmentedSystemModel)
            % do not compute square root of state cov as it is not
            % neccessarily positive-definite
           [predictedAugmentedStateMean, ...
                predictedAugmentedStateCov] = augmentedSystemModel.analyticMoments(this.augmentedStateMean, ...
                                                                    this.augmentedStateCov);
            
            % check if predicted state covariance is valid
            this.checkCovPrediction(predictedAugmentedStateCov(1:this.dimState, 1:this.dimState), ...
                                                           'Predicted state');
                     

            % all checks successful, store augmented estimates
            this.augmentedStateMean = predictedAugmentedStateMean;
            this.augmentedStateCov = predictedAugmentedStateCov;

        end       
    end
end % classdef
