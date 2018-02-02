classdef (Sealed) DelayedKF < AnalyticKF & DelayedMeasurementsFilter
    % This class represents a Kalman filter that can operate with time delayed control inputs 
    % and/or time delayed measurements. The filter calculates the linear minimum
    % mean square state estimate of a system based on a model of that system,
    % the expected control input applied to the system and possibly
    % time-delayed measurements of the state.
    %
    % Literature: 
    %   Maryam Moayedi, Yung Kuan Foo and Yeng Chai Soha
    %   Filtering for networked control systems with single/multiple measurement packets
    %   subject to multiple-step measurement delays and multiple packet dropouts
    %   International Journal of Systems Science (2011)
    %
    %
    % AUTHOR:       Jörg Fischer
    % LAST UPDATE:  Maxim Dolgov, 26.06.2013
    %               Florian Rosenthal, 12.01.2017
    
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
    
    properties (Access = private)
        augmentedStateMean; 
        augmentedStateCov; % needed as augmented cov not necessarily pd
    end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %% DelayedKF
        function s = DelayedKF(maxMeasDelay, name)
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
            s@DelayedMeasurementsFilter(maxMeasDelay, name);
            s@AnalyticKF(name);
            s.setUseAnalyticMeasurementModel(true);    
        end % function DelayedKF

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
            
            setState@GaussianFilter(this, state);
            [mean, cov] = state.getMeanAndCovariance();
            %augment state
            this.augmentedStateMean = repmat(mean, this.maxMeasurementDelay + 1, 1);
            % this matrix (covariance of augmented state) is not necessarily pd
            this.augmentedStateCov = repmat(cov, this.maxMeasurementDelay + 1, this.maxMeasurementDelay + 1);
        end
        
        %% getLastUpdateData
        function [measurement, ...
                  measMean, ...
                  measCov, ...
                  stateMeasCrossCov, ...
                  numIterations] = getLastUpdateData(this)
              this.error('NotImplemented', 'NotImplemented');
        end
        
    end % methods public

    methods(Access = protected)
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
                % now perform the analytical update: incoprorate all
                % measurements per delay at once
                arrayfun(@(i) this.updateAnalytic(measModel, measurementsPerDelay{i}, timeDelays(i)), 1:numel(timeDelays));
           end
        end
    
        %% performPrediction
        function performPrediction(this, sysModel)
             if ~Checks.isClass(sysModel, 'DelayedKFSystemModel')
                this.errorSysModel('DelayedKFSystemModel');
             end
            this.predictAnalytic(sysModel);
        end
        
        %% updateAnalytic
        function updateAnalytic(this, measModel, measurements, delay)
            [dimMeas, numMeas] = size(measurements);
            baseMeasMatrix = measModel.measMatrix;
           
            [noiseMean, noiseCov] = measModel.noise.getMeanAndCovariance();
            dimStackedMeas = size(baseMeasMatrix, 1) * numMeas; % what is expected according to the original model
                       
            augmentedMeasMatrix = [zeros(dimMeas, this.dimState * delay), ...
                    baseMeasMatrix, zeros(dimMeas, this.dimState * (this.maxMeasurementDelay - delay))];
            
            % compute required moments
            measMean = repmat(augmentedMeasMatrix * this.augmentedStateMean + noiseMean, numMeas, 1);
            measCov = Utils.baseBlockDiag(augmentedMeasMatrix * this.augmentedStateCov * augmentedMeasMatrix', ...
                noiseCov, numMeas);
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
                        
            this.checkAndSaveUpdate(postMean, augmentedPostCov(1:this.dimState, 1:this.dimState));
            % all checks successful, store augmented estimates
            this.augmentedStateMean(1:this.dimState) = postMean;
            this.augmentedStateCov = augmentedPostCov;
        end
        
        %% predictAnalytic
        function predictAnalytic(this, augmentedSystemModel)
           [predictedAugmentedStateMean, ...
                predictedAugmentedStateCov] = augmentedSystemModel.analyticPredictedMoments(this.augmentedStateMean, ...
                                                                    this.augmentedStateCov);
            
            this.checkAndSavePrediction(predictedAugmentedStateMean(1:this.dimState), ...
                predictedAugmentedStateCov(1:this.dimState, 1:this.dimState));
            % all checks successful, store augmented estimates
            this.augmentedStateMean = predictedAugmentedStateMean;
            this.augmentedStateCov = predictedAugmentedStateCov;
        end
    end
end % classdef
