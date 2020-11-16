classdef SequenceBasedTrackingControllerTest < matlab.unittest.TestCase
    % Test cases for SequenceBasedTrackingController.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Access = protected)
        dimX;
        dimU;
        sequenceLength;
        needsFilter;
        stateTrajectory;
        inputTrajectory;
        plantState;
        Z;
        refTrajectory;
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.dimX = 3;
            this.dimU = 2;
            
            this.sequenceLength = 6;
            this.needsFilter = false;
            this.stateTrajectory = ones(this.dimX, 10);
            this.inputTrajectory = ones(this.dimU, 10);
        
            this.plantState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            
            this.Z = [1 1 1];
            this.refTrajectory = ones(1, this.sequenceLength + 1);
            this.controllerUnderTest = SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, this.Z, this.refTrajectory, this.sequenceLength + 1);
        end
    end
    
    methods (Test)
        %% testSequenceBasedTrackingControllerInvalidZMatrix
        function testSequenceBasedTrackingControllerInvalidZMatrix(this)
            expectedErrId = 'SequenceBasedTrackingController:InvalidZMatrix';
            
            invalidZ = ones(1, this.dimU); % invalid number of cols
            this.verifyError(@() SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, invalidZ, this.refTrajectory, this.sequenceLength + 1), expectedErrId);
               
            invalidZ = this.Z; 
            invalidZ(2) = nan; % must be finite
            this.verifyError(@() SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, invalidZ, this.refTrajectory, this.sequenceLength + 1), expectedErrId);
            
        end
        
        %% testSequenceBasedTrackingControllerInvalidRefTrajectory
        function testSequenceBasedTrackingControllerInvalidRefTrajectory(this)
            expectedErrId = 'SequenceBasedTrackingController:InvalidReferenceTrajectory';
               
            invalidRefTrajectory = this.refTrajectory(:, 1:end-1); % too short
            this.verifyError(@() SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, this.Z, invalidRefTrajectory, this.sequenceLength + 1), expectedErrId);
            
            invalidRefTrajectory = [this.refTrajectory; this.refTrajectory]; % wrong dimension
            this.verifyError(@() SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, this.Z, invalidRefTrajectory, this.sequenceLength + 1), expectedErrId);
            
            % now we do not check the length
            invalidRefTrajectory = [this.refTrajectory; this.refTrajectory]; % wrong dimension
            this.verifyError(@() SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, this.Z, invalidRefTrajectory, []), expectedErrId);
        end
        
        %% testSequenceBasedTrackingController
        function testSequenceBasedTrackingController(this)
            % successful call of ctor
            controller = SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, this.Z, this.refTrajectory, this.sequenceLength + 1);
            
            [actualZ, actualDimRef] = controller.getProperties();
            
            this.verifyEqual(controller.refTrajectory, this.refTrajectory);
            this.verifyEqual(actualZ, this.Z);
            this.verifyEqual(actualDimRef, 1);
            this.verifyEqual(controller.requiresExternalStateEstimate, this.needsFilter);
                                   
            newRefTrajectory = [this.refTrajectory this.refTrajectory];
            % now a successfull call whith longer trajectory but no check
            controller = SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, this.Z, newRefTrajectory, []);
            
            [actualZ, actualDimRef] = controller.getProperties();
            
            this.verifyEqual(controller.refTrajectory, newRefTrajectory);
            this.verifyEqual(actualZ, this.Z);
            this.verifyEqual(actualDimRef, 1);
            this.verifyEqual(controller.requiresExternalStateEstimate, this.needsFilter);
            
            % now again a successful call with a longer trajectory and check
            newRefTrajectory = [this.refTrajectory this.refTrajectory];
            controller = SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, this.Z, newRefTrajectory, this.sequenceLength);
            
            [actualZ, actualDimRef] = controller.getProperties();
            
            this.verifyEqual(controller.refTrajectory, newRefTrajectory);
            this.verifyEqual(actualZ, this.Z);
            this.verifyEqual(actualDimRef, 1);
            this.verifyEqual(controller.requiresExternalStateEstimate, this.needsFilter);
            
            % finally, a successful call without Z
            controller = SequenceBasedTrackingControllerStub(this.dimX, this.dimU, this.sequenceLength, ...
                this.needsFilter, [], this.refTrajectory, this.sequenceLength + 1);
            
            [actualZ, actualDimRef] = controller.getProperties();
            this.verifyEmpty(controller.refTrajectory);
            this.verifyEmpty(actualZ);
            this.verifyEqual(actualDimRef, this.dimX);
            this.verifyEqual(controller.requiresExternalStateEstimate, this.needsFilter);
        end
        
        %% testGetDeviationFromRefForStateInvalidState
        function testGetDeviationFromRefForStateInvalidState(this)
            expectedErrId = 'SequenceBasedTrackingController:GetDeviationFromRefForState:InvalidTrueState';
                        
            invalidState = this; % not a vector
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(invalidState, 1), ...
                expectedErrId);
            
            invalidState = ones(this.dimU, 1); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(invalidState, 1), ...
                expectedErrId);
        end
        
        %% testGetDeviationFromRefForState
        function testGetDeviationFromRefForState(this)
            state = zeros(this.dimX, 1);
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(state, 1);
            
            % dummy always returns 0
            this.verifyEqual(actualDeviation, 0);
        end
    end
end

