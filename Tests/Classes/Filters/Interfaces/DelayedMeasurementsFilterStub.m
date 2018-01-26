classdef (Sealed) DelayedMeasurementsFilterStub < DelayedMeasurementsFilter
    % This class is a dummy implementation of a filter which is used for testing
    % purposes.
       
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
   
    methods (Access = public)
        %% DelayedMeasurementsFilterStub
        function this = DelayedMeasurementsFilterStub(maxMeasDelay, name)
            this@DelayedMeasurementsFilter(maxMeasDelay, name);
        end
        
        %% getState
        function state = getState(~)
            state = Gaussian(0, 1);
        end
        
        %% setState
        function setState(~, ~)
            % this dummy implementation does nothing
        end
        
        %% getPointEstimate
        function [pointEstimate, uncertainty] = getPointEstimate(~)
            pointEstimate = 0;
            uncertainty = 1;
        end
        
        %% obtainApplicableMeasurements
        function [applicableMeas, applicableDelays] = obtainApplicableMeasurements(this, measurements, measDelays)
            % just call the protected base class method
            [applicableMeas, applicableDelays] = this.getApplicableMeasurements(measurements, measDelays);
        end
        
        %% setMeasurementUpdateData
        function setMeasurementUpdateData(this, numUsedMeas, numDiscardedMeas)
            this.setLastUpdateMeasurementData(numUsedMeas, numDiscardedMeas);
        end
    end
    
    methods (Access = protected)
        %% performPrediction
        function performPrediction(~,~)
             % this dummy implementation does nothing
        end
        
        %% performUpdate
        function performUpdate(~, ~, ~, ~)
            % this dummy implementation does nothing
        end
    end
end

