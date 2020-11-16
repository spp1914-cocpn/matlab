 classdef (Abstract) DelayedMeasurementsFilter < Filter
    % Abstract base class for filters that can perform measurement updates
    % with delayed and out-of-sequence measurements.
    
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
    
    properties (SetAccess = immutable, GetAccess = protected)
        % the maximum allowed measurement delay (i.e., measurements with longer delay are ignored)
        maxMeasurementDelay(1,1) double {Validator.validateMaxPacketDelay(maxMeasurementDelay, 0)} = 1;
    end
    
    properties (SetAccess = public, GetAccess = private)
        issueWarningIgnoreDelayedMeas(1,1) logical = false; % default false
    end
    
    properties (Access = private)
        lastNumUsedMeas;
        lastNumDiscardedMeas;
    end
            
    methods(Access = protected)
        %% DelayedMeasurementsFilter
        function this = DelayedMeasurementsFilter(maxMeasDelay, name)
            % Class constructor.
            %
            % Parameters:
            %   >> maxMeasDelay (Nonnegative integer)
            %      The maximum allowed measurement delay. That is, all
            %      measurements to be processed with a larger delay are discarded by the
            %      filter. If 0 is passed, the filter reduces to a
            %      filter which cannot cope with out-of-sequence and
            %      delayed measurements.
            %   >> name (Char)
            %      An appropriate filter name / description of the implemented
            %      filter.
            %
            % Returns:
            %   << this (DelayedMeasurementsFilter)
            %      A new DelayedMeasurementsFilter instance.
            this@Filter(name);
            
            this.maxMeasurementDelay = maxMeasDelay;
            this.lastNumUsedMeas = 0;
            this.lastNumDiscardedMeas = 0;
        end
        
        %% getApplicableMeasurements
        function [applicableMeas, applicableDelays] = getApplicableMeasurements(this, measurements, measDelays)
            % find the measurements with valid delays
            idx = find(measDelays <= this.maxMeasurementDelay);
            applicableMeas = measurements(:, idx);
            applicableDelays = measDelays(idx);
            if numel(idx) ~= numel(measDelays)
                this.warnIgnoreDelayedMeas(numel(idx), numel(measDelays));
            end
        end
    
        %% warnIgnoreDelayedMeas
        function warnIgnoreDelayedMeas(this, numApplicableMeas, numMeas)
            if this.issueWarningIgnoreDelayedMeas
                this.warning('IgnoringMeasurements', ...
                    '** %s\nIgnoring %d of %d measurements. **', 'Delay too large', ...
                    numMeas - numApplicableMeas, numMeas);
            end
        end
               
        %% performStep
        function performStep(this, sysModel, measModel, applicableMeasurements, measDelays)
            this.performPrediction(sysModel);
            this.performUpdate(measModel, applicableMeasurements, measDelays);
        end
        
        %% setLastUpdateMeasurementData
        function setLastUpdateMeasurementData(this, numUsedMeas, numDiscardedMeas)
            % Store information on the last performed measurement update (due to calls of step or update).
            %
            % Parameters:
            %   >> numUsedMeas (Nonnegative integer)
            %      The number of measurements used by the last measurement update.
            %
            %   >> numDiscardedMeas (Nonnegative integer)
            %      The number of measurements discarded during the last
            %      measurement update due to their delays.
            
            this.lastNumUsedMeas = numUsedMeas;
            this.lastNumDiscardedMeas = numDiscardedMeas;
        end
    end
    
    methods (Abstract, Access = protected)
        performUpdate(this, measModel, applicableMeasurements, measDelays);
    end
    
    methods (Access = public)
        %% update
        function runtime = update(this, measModel, measurements, delays)
            % Perform a measurement update (filter step) using the given measurement(s).
            %
            % Parameters:
            %   >> measModel (Arbitrary class (filter dependent))
            %      Measurement model that provides the mapping between measurements and the system state.
            %
            %   >> measurements (Matrix)
            %      Column-wise arranged measurement vectors, where each column represents an individual
            %      (possibly delayed) measurement. 
            %
            %   >> delays (Vector of nonnegative integers (one for each measurement))
            %      Vector of nonnegative integers, the i-th element denotes the delay of the i-th measurement.
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the measurement update.
            
            if isempty(measurements) || ~ismatrix(measurements) || any(~isfinite(measurements(:)))
                this.error('InvalidMeasurements', ...
                    '** Measurements must be column-wise arranged and real-valued **');
            end
            numMeas = size(measurements, 2);
            
            if ~Checks.isNonNegativeVec(delays, numMeas) ...
                    || any(arrayfun(@(measDelay) mod(measDelay, 1) ~= 0, delays))
                this.error('InvalidMeasDelay', ...
                    '** Each measurement delay must be a nonnegative integer **');
            end

            [applicableMeasurements, applicableDelays] = this.getApplicableMeasurements(measurements, delays);
            if nargout == 1
                s = tic;
                this.tryPerformUpdate(measModel, applicableMeasurements, applicableDelays);
                runtime = toc(s);
            else
                this.tryPerformUpdate(measModel, applicableMeasurements, applicableDelays);
            end
            numUsedMeas = size(applicableMeasurements, 2);
            this.setLastUpdateMeasurementData(numUsedMeas, size(measurements, 2) - numUsedMeas);
        end
        
        %% step
        function runtime = step(this, sysModel, measModel, measurements, delays)
            % Perform a combined time and measurement update.
            %
            % Parameters:
            %   >> sysModel (Arbitrary class (filter dependent))
            %      System model that provides the mapping between the prior system
            %      state and the predicted state (i.e., the system state's temporal evolution).
            %
            %   >> measModel (Arbitrary class (filter dependent))
            %      Measurement model that provides the mapping between measurements and the system state.
            %
            %   >> measurements (Matrix, can be empty)
            %      Column-wise arranged measurement vectors, where each column represents an individual
            %      (possibly delayed) measurement. 
            %      If no measurements to be processed, pass the empty
            %      matrix here.
            %
            %   >> delays (Nonnegative integer or vector of nonnegative integers (one for each measurement))
            %      If a vector is passed here, the i-th element denotes the delay of the i-th measurement.
            %      If no measurements to be processed, pass the empty
            %      matrix here.
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the combined time and measurement update.

            if ~ismatrix(measurements) || any(~isfinite(measurements(:)))
                this.error('InvalidMeasurements', ...
                    '** Measurements must be column-wise arranged and real-valued **');
            end
            numMeas = size(measurements, 2);
            
            if ~Checks.isNonNegativeVec(delays, numMeas) ...
                    || any(arrayfun(@(measDelay) mod(measDelay, 1) ~= 0, delays))
                this.error('InvalidMeasDelay', ...
                    '** Each measurement delay must be a nonnegative integer **');
            end

            [applicableMeasurements, applicableDelays] = this.getApplicableMeasurements(measurements, delays);
            if nargout == 1
                s = tic;
                this.tryPerformStep(sysModel, measModel, applicableMeasurements, applicableDelays); 
                runtime = toc(s);
            else
                this.tryPerformStep(sysModel, measModel, applicableMeasurements, applicableDelays);
            end
            numUsedMeas = size(applicableMeasurements, 2);
            this.setLastUpdateMeasurementData(numUsedMeas, size(measurements, 2) - numUsedMeas);
        end
        
        %% getLastUpdateMeasurementData
        function [numUsedMeas, numDiscardedMeas] = getLastUpdateMeasurementData(this)
            % Get information on the last performed measurement update (due to calls of step or update).
            %
            % Returns:
            %   << numUsedMeas (Nonnegative integer)
            %      The number of measurements used by the last measurement update.
            %
            %   << numDiscardedMeas (Nonnegative integer)
            %      The number of measurements discarded during the last
            %      measurement update due to their delays.
            
            numUsedMeas = this.lastNumUsedMeas;
            numDiscardedMeas = this.lastNumDiscardedMeas;
        end
        
        %% getMaxMeasurementDelay
        function maxMeasDelay = getMaxMeasurementDelay(this)
            % Get the maximum allowed measurement delay.
            %
            % Returns:
            %   << maxMeasDelay (Nonnegative integer)
            %      The maximum allowed measurement delay (in time steps).
            %   
            maxMeasDelay = this.maxMeasurementDelay;
            % explicit getter mainly for compability with base class Filter
            % from NonlinearEstimationToolbox
        end
    end
    
    methods (Access = private)
        %% tryPerformUpdate
        function tryPerformUpdate(this, measModel, applicableMeasurements, applicableDelays)
            try
                this.performUpdate(measModel, applicableMeasurements, applicableDelays);            
            catch ex
                % copied from Filter
                if ~strcmp(ex.identifier, 'Filter:IgnoreMeasurement')
                    % real error => do not catch it
                    ex.rethrow();
                end
            end
        end
        
        %% tryPerformStep
        function tryPerformStep(this, sysModel, measModel, applicableMeasurements, applicableDelays)
            try
                this.performStep(sysModel, measModel, applicableMeasurements, applicableDelays);
            catch ex
                 % copied from Filter
                if ~strcmp(ex.identifier, 'Filter:IgnoreMeasurement') ...
                        && ~strcmp(ex.identifier, 'Filter:IgnorePrediction')
                    % real error => do not catch it
                   ex.rethrow();
                end
           end    
        end
    end
end

