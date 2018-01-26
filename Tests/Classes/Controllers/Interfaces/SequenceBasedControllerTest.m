%classdef SequenceBasedControllerTest < matlab.mock.TestCase % requires R2017a or later
classdef SequenceBasedControllerTest < matlab.unittest.TestCase     
    % Test cases for SequenceBasedController.
    
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
        dimX;
        dimU;
        sequenceLength;
        stateTrajectory;
        inputTrajectory;
        plantState;
        
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.dimX = 3;
            this.dimU = 2;
            
            this.sequenceLength = 6;
            this.stateTrajectory = ones(this.dimX, 10);
            this.inputTrajectory = ones(this.dimU, 10);
        
            this.plantState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            
            this.controllerUnderTest = SequenceBasedControllerStub(this.dimX, this.dimU, this.sequenceLength);
        end
    end

    methods (Test)
        
        %% testSequenceBasedControllerInvalidSequenceLength
        function testSequenceBasedControllerInvalidSequenceLength(this)
            expectedErrId = 'Validator:ValidateSequenceLength:InvalidSequenceLength';
            
            % must be positive
            invalidSequenceLength = 0;
            this.verifyError(@() SequenceBasedControllerStub(this.dimX, this.dimU, invalidSequenceLength), ...
                expectedErrId);
            
            % must be integer
            invalidSequenceLength = 4.5;
            this.verifyError(@() SequenceBasedControllerStub(this.dimX, this.dimU, invalidSequenceLength), ...
                expectedErrId);
            
            % must be finite
            invalidSequenceLength = inf;
            this.verifyError(@() SequenceBasedControllerStub(this.dimX, this.dimU, invalidSequenceLength), ...
                expectedErrId);
            
            % must be finite
            invalidSequenceLength = nan;
            this.verifyError(@() SequenceBasedControllerStub(this.dimX, this.dimU, invalidSequenceLength), ...
                expectedErrId);
        end
        
        %% testSequenceBasedController
        function testSequenceBasedController(this)
            controller = SequenceBasedControllerStub(this.dimX, this.dimU, this.sequenceLength);
            [actualDimX, actualDimU, actualSeqLength] = controller.getProperties();
            
            this.verifyEqual(actualDimX, this.dimX);
            this.verifyEqual(actualDimU, this.dimU);
            this.verifyEqual(actualSeqLength, this.sequenceLength);
        end
        
        %% testComputeCostsInvalidParams
        function testComputeCostsInvalidParams(this)
            expectedErrId = 'SequenceBasedController:ComputeCosts';
            
            % first, invalid state trajectory: matrix of wrong dimension
            invalidStateTrajectory = this.inputTrajectory;
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, this.inputTrajectory), ...
                expectedErrId);
            
            % now, invalid input trajectory; wrong dimensions
            invalidInputTrajectory = this.stateTrajectory;
            this.verifyError(@() this.controllerUnderTest.computeCosts(this.stateTrajectory, invalidInputTrajectory), ...
                expectedErrId);
        end
        
        %% testComputeCosts
        function testComputeCosts(this)
            % the controller stub always returns costs 0
            this.verifyEqual(this.controllerUnderTest.computeCosts(this.stateTrajectory, this.inputTrajectory), 0);
        end
        
        %% testComputeControlSequenceInvalidState
        function testComputeControlSequenceInvalidState(this)
            expectedErrId = 'SequenceBasedController:ComputeControlSequence';
            
            invalidPlantState = this.sequenceLength; % invalid type
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(invalidPlantState), ...
                expectedErrId);
            
            invalidPlantState = Gaussian(zeros(this.dimX + 1, 1), eye(this.dimX + 1)); % Distribution of invalid dim
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(invalidPlantState), ...
                expectedErrId);
        end
        
        %% testComputeControlSequence
        function testComputeControlSequence(this)
            expectedSequence = ones(this.dimU * this.sequenceLength, 1);
            actualSequence = this.controllerUnderTest.computeControlSequence(this.plantState);
            
            this.verifyEqual(actualSequence, expectedSequence);
        end
    end
end




