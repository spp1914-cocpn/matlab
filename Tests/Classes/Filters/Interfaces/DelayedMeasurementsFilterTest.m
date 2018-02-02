classdef DelayedMeasurementsFilterTest < matlab.unittest.TestCase
    % Test cases for DelayedMeasurementsFilter.
    
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
    
    properties (Constant, Access = private)
        zeroMeasDelay = 0;
        emptyChar = '';
    end
    
    properties (Access = private)
        dimX;
        dimY;
       
        filterName;
        
        maxMeasDelay;
        numMeas;
        delayedMeasurements;
        applicableMeasurements;
        delays;
        applicableDelays;
         
        filterUnderTest;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.filterName = 'TestFilter';
            this.maxMeasDelay = 4;
            this.numMeas = 6;
            
            this.dimX = 3;
            this.dimY = 2;
            
            this.delayedMeasurements = ones(this.dimY, this.numMeas);
            this.delays = [1 5 2 4 3 6];
            % extract the expected measurements, i.e, those with delay <= maxMeasDelay
            applicableIdx = [1 3 4 5];
            this.applicableDelays = this.delays(applicableIdx);
            this.applicableMeasurements = this.delayedMeasurements(:, applicableIdx);
             
            this.filterUnderTest = DelayedMeasurementsFilterStub(this.maxMeasDelay, this.filterName);
        end
    end
    
    methods (Test)
        %% testDelayedMeasurementsFilterInvalidMaxMeasDelay
        function testDelayedMeasurementsFilterInvalidMaxMeasDelay(this)
            expectedErrId = 'Filter:InvalidMaxMeasDelay';
            
            invalidMaxMeasDelay = [1 2]; % not a scalar
            this.verifyError(@() DelayedMeasurementsFilterStub(invalidMaxMeasDelay, this.filterName), expectedErrId);
            
            invalidMaxMeasDelay = -1; % negative scalar
            this.verifyError(@() DelayedMeasurementsFilterStub(invalidMaxMeasDelay, this.filterName), expectedErrId);
            
            invalidMaxMeasDelay = 1.5; % not an integer
            this.verifyError(@() DelayedMeasurementsFilterStub(invalidMaxMeasDelay, this.filterName), expectedErrId);
            
            invalidMaxMeasDelay = inf; % not finite
            this.verifyError(@() DelayedMeasurementsFilterStub(invalidMaxMeasDelay, this.filterName), expectedErrId);
        end
        
        %% testDelayedMeasurementsFilterInvalidName
        function testDelayedMeasurementsFilterInvalidName(this)
            expectedErrId = 'Filter:InvalidFilterName';
            
            invalidName = {'name'}; % not a char array
            this.verifyError(@() DelayedMeasurementsFilterStub(this.maxMeasDelay, invalidName), expectedErrId);
            
            invalidName = eye(3); % not a char array
            this.verifyError(@() DelayedMeasurementsFilterStub(this.maxMeasDelay, invalidName), expectedErrId);
        end
        
        %% testDelayedMeasurementsFilter
        function testDelayedMeasurementsFilter(this)
            filter = DelayedMeasurementsFilterStub(this.maxMeasDelay, this.filterName);
            
            this.verifyEqual(filter.getName(), this.filterName);
            this.verifyEqual(filter.getMaxMeasurementDelay(), this.maxMeasDelay);
           
            % zero meas delay and empty char are allowed
            filter = DelayedMeasurementsFilterStub(DelayedMeasurementsFilterTest.zeroMeasDelay, ...
                DelayedMeasurementsFilterTest.emptyChar);
            
            this.verifyEmpty(filter.getName());
            this.verifyEqual(filter.getMaxMeasurementDelay(), DelayedMeasurementsFilterTest.zeroMeasDelay);
        end
        
        %% testGetApplicableMeasurementsNoWarning
        function testGetApplicableMeasurementsNoWarning(this)
            % discard two measurements, due to a large delay
            % no warning should be issued, however
            
            [actualApplicableMeas, actualApplicableDelays] = ...
                this.verifyWarningFree(...
                    @() this.filterUnderTest.obtainApplicableMeasurements(this.delayedMeasurements, this.delays));
             
            this.verifySize(actualApplicableMeas, [this.dimY, this.maxMeasDelay]);
            this.verifySize(actualApplicableDelays, [1 this.maxMeasDelay]);
            
            this.verifyEqual(actualApplicableMeas, this.applicableMeasurements); 
            this.verifyEqual(actualApplicableDelays, this.applicableDelays);
            
        end
         
        %% testGetApplicableMeasurementsWarning
        function testGetApplicableMeasurementsWarning(this)
            % discard two measurements, due to a large delay
            % warning should be issued
            this.filterUnderTest.issueWarningIgnoreDelayedMeas = true;
            expectedWarningId = 'Filter:IgnoringMeasurements';
            
            [actualApplicableMeas, actualApplicableDelays] = ...
                this.verifyWarning(...
                    @() this.filterUnderTest.obtainApplicableMeasurements(this.delayedMeasurements, this.delays), ...
                    expectedWarningId);
            
            this.verifySize(actualApplicableMeas, [this.dimY, this.maxMeasDelay]);
            this.verifySize(actualApplicableDelays, [1 this.maxMeasDelay]);
            
            this.verifyEqual(actualApplicableMeas, this.applicableMeasurements); 
            this.verifyEqual(actualApplicableDelays, this.applicableDelays);
            
        end
         
        %% testGetApplicableMeasurementsNothingToDiscard
        function testGetApplicableMeasurementsNothingToDiscard(this)
            % no measurements to discard now
            
            [actualApplicableMeas, actualApplicableDelays] = ...
                this.verifyWarningFree(...
                    @() this.filterUnderTest.obtainApplicableMeasurements(this.applicableMeasurements, this.applicableDelays));
            
            this.verifySize(actualApplicableMeas, [this.dimY, this.maxMeasDelay]);
            this.verifySize(actualApplicableDelays, [1 this.maxMeasDelay]);
            
            this.verifyEqual(actualApplicableMeas, this.applicableMeasurements); 
            this.verifyEqual(actualApplicableDelays, this.applicableDelays);
            
        end
        
        %% testUpdateInvalidMeasurementsNoDelays
        function testUpdateInvalidMeasurementsNoDelays(this)
            expectedErrId = 'Filter:InvalidMeasurements';
            
            % measurements must not be empty
            invalidMeasurements = [];
            % perform an "update" (model not required), do not pass
            % delays
            this.verifyError(@() this.filterUnderTest.update([], invalidMeasurements), expectedErrId);
            
            % measurements must be a matrix or vector
            invalidMeasurements = ones(3, 3, 3);
            % perform an "update" (model not required), do not pass
            % delays
            this.verifyError(@() this.filterUnderTest.update([], invalidMeasurements), expectedErrId);
        end
     
         %% testUpdateInvalidMeasurements
        function testUpdateInvalidMeasurements(this)
            expectedErrId = 'Filter:InvalidMeasurements';
            
            % measurements must not be empty
            invalidMeasurements = [];
             % perform an "update" (model not required)
            this.verifyError(@() this.filterUnderTest.update([], invalidMeasurements, this.delays), expectedErrId);
            
            % measurements must be a matrix or vector
            invalidMeasurements = ones(3, 3, 3);
             % perform an "update" (model not required)
            this.verifyError(@() this.filterUnderTest.update([], invalidMeasurements, this.delays), expectedErrId);
        end
        
         %% testUpdateInvalidDelays
        function testUpdateInvalidDelays(this)
            expectedErrId = 'Filter:InvalidMeasDelay';
            
            % delays must not be a negative scalar
            invalidDelay = -1;
             % perform an "update" (model not required)
            this.verifyError(@() this.filterUnderTest.update([], this.delayedMeasurements, invalidDelay), ...
                expectedErrId);
            
            % if delays are given as nonnegative vector, must consist of
            % <numMeas> elements
            measurements = this.delayedMeasurements;
            invalidDelays = this.delays(1:end-1); % vector is too small
            % perform an "update" (model not required)
            this.verifyError(@() this.filterUnderTest.update([], measurements, invalidDelays), ...
                expectedErrId);
            
            % if delays are given as nonnegative vector, must consist of
            % <numMeas> elements
            measurements = this.delayedMeasurements(:, 1);
            invalidDelays = this.delays; % vector is too large
             % perform an "update" (model not required)
            this.verifyError(@() this.filterUnderTest.update([], measurements, invalidDelays), ...
                expectedErrId);
            
            % if delays are given as vector of appropriate dimension, all
            % entries must be nonnegative
            measurements = this.delayedMeasurements;
            invalidDelays = [this.delays(1:end-1), -1];
             % perform an "update" (model not required)
            this.verifyError(@() this.filterUnderTest.update([], measurements, invalidDelays), ...
                expectedErrId);
            
            % if delays are given as vector of appropriate dimension, all
            % entries must be integers
            measurements = this.delayedMeasurements;
            invalidDelays = [this.delays(1:end-1), 1.5];
             % perform an "update" (model not required)
            this.verifyError(@() this.filterUnderTest.update([], measurements, invalidDelays), ...
                expectedErrId);
        end
        
        %% testUpdateNoDelays
        function testUpdateNoDelays(this)
             % perform an "update" where no measurements are to be discarded
             % since delays are assumed to be zero
            measurements = this.delayedMeasurements;
            
            this.filterUnderTest.update([], measurements);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            expectedNumDiscardedMeas = 0;
            expectedNumUsedMeas = size(measurements, 2);
            
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
        end
        
        %% testUpdateDiscardMeasurements
        function testUpdateDiscardMeasurements(this)
            % perform an "update" where measurements are to be discarded
                      
            this.filterUnderTest.update([], this.delayedMeasurements, this.delays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            expectedNumDiscardedMeas = numel(this.delays) - numel(this.applicableDelays);
            expectedNumUsedMeas = numel(this.delays) - expectedNumDiscardedMeas;
            
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
        end
        
        %% testUpdateNothingToDiscard
        function testUpdateNothingToDiscard(this)
            % perform an "update" where where no measurements are to be discarded
                      
            this.filterUnderTest.update([], this.applicableMeasurements, this.applicableDelays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            expectedNumDiscardedMeas = 0;
            expectedNumUsedMeas = numel(this.applicableDelays);
            
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
        end
%%
%%
        %% testStepInvalidMeasurementsNoDelays
        function testStepInvalidMeasurementsNoDelays(this)
            expectedErrId = 'Filter:InvalidMeasurements';
            
            % measurements must not be empty
            invalidMeasurements = [];
            % perform a "step" (models are not required), do not pass
            % delays
            this.verifyError(@() this.filterUnderTest.step([], [], invalidMeasurements), expectedErrId);
            
            % measurements must be a matrix or vector
            invalidMeasurements = ones(3, 3, 3);
            % perform a "step" (models are not required), do not pass
            % delays
            this.verifyError(@() this.filterUnderTest.step([], [], invalidMeasurements), expectedErrId);
        end
     
         %% testStepInvalidMeasurements
        function testStepInvalidMeasurements(this)
            expectedErrId = 'Filter:InvalidMeasurements';
            
            % measurements must not be empty
            invalidMeasurements = [];
            % perform a "step" (models are not required)
            this.verifyError(@() this.filterUnderTest.step([], [], invalidMeasurements, this.delays), expectedErrId);
            
            % measurements must be a matrix or vector
            invalidMeasurements = ones(3, 3, 3);
            % perform a "step" (models are not required)
            this.verifyError(@() this.filterUnderTest.step([], [], invalidMeasurements, this.delays), expectedErrId);
        end
         
        %% testStepInvalidDelays
        function testStepInvalidDelays(this)
            expectedErrId = 'Filter:InvalidMeasDelay';
            
            % delays must not be a negative scalar
            invalidDelay = -1;
            % perform a "step" (models are not required)
            this.verifyError(@() this.filterUnderTest.step([], [], this.delayedMeasurements, invalidDelay), ...
                expectedErrId);
            
            % if delays are given as nonnegative vector, must consist of
            % <numMeas> elements
            measurements = this.delayedMeasurements;
            invalidDelays = this.delays(1:end-1); % vector is too small
            % perform a "step" (models are not required)
            this.verifyError(@() this.filterUnderTest.step([], [], measurements, invalidDelays), ...
                expectedErrId);
            
            % if delays are given as nonnegative vector, must consist of
            % <numMeas> elements
            measurements = this.delayedMeasurements(:, 1);
            invalidDelays = this.delays; % vector is too large
            % perform a "step" (models are not required)
            this.verifyError(@() this.filterUnderTest.step([], [], measurements, invalidDelays), ...
                expectedErrId);
            
            % if delays are given as vector of appropriate dimension, all
            % entries must be nonnegative
            measurements = this.delayedMeasurements;
            invalidDelays = [this.delays(1:end-1), -1];
            % perform a "step" (models are not required)
            this.verifyError(@() this.filterUnderTest.step([], [], measurements, invalidDelays), ...
                expectedErrId);
            
            % if delays are given as vector of appropriate dimension, all
            % entries must be integers
            measurements = this.delayedMeasurements;
            invalidDelays = [this.delays(1:end-1), 1.5];
            % perform a "step" (models are not required)
            this.verifyError(@() this.filterUnderTest.step([], [], measurements, invalidDelays), ...
                expectedErrId);
        end
        
        %% testStepNoDelays
        function testStepNoDelays(this)
             % perform a "step" where no measurements are to be discarded
             % since delays are assumed to be zero
            measurements = this.delayedMeasurements;
            
            this.filterUnderTest.step([], [], measurements);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            expectedNumDiscardedMeas = 0;
            expectedNumUsedMeas = size(measurements, 2);
            
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
        end
        
        %% testStepDiscardMeasurements
        function testStepDiscardMeasurements(this)
            % perform a "step where measurements are to be discarded
                      
            this.filterUnderTest.step([], [], this.delayedMeasurements, this.delays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            expectedNumDiscardedMeas = numel(this.delays) - numel(this.applicableDelays);
            expectedNumUsedMeas = numel(this.delays) - expectedNumDiscardedMeas;
            
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
        end
        
        %% testStepNothingToDiscard
        function testStepNothingToDiscard(this)
            % perform a "step where where no measurements are to be discarded
                      
            this.filterUnderTest.step([], [], this.applicableMeasurements, this.applicableDelays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            expectedNumDiscardedMeas = 0;
            expectedNumUsedMeas = numel(this.applicableDelays);
            
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
        end
        
        %% testGetLastUpdateMeasurement
        function testGetLastUpdateMeasurementData(this)
             [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            % should be both zero initially
            expectedNumDiscardedMeas = 0;
            expectedNumUsedMeas = 0;
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
            
            % now call set function
            expectedNumUsedMeas = 2;
            expectedNumDiscardedMeas = 3;
            
            this.filterUnderTest.setMeasurementUpdateData(expectedNumUsedMeas, expectedNumDiscardedMeas);
            [actualNumUsedMeas, actualNumDiscardedMeas] = ...
                this.filterUnderTest.getLastUpdateMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, expectedNumUsedMeas);
            this.verifyEqual(actualNumDiscardedMeas, expectedNumDiscardedMeas);
        end
        
        %% testGetMaxMeasurementDelay
        function testGetMaxMeasurementDelay(this)
            actualMaxMeasDelay = this.filterUnderTest.getMaxMeasurementDelay();
            
            this.verifyEqual(actualMaxMeasDelay, this.maxMeasDelay);
        end
    end
end

