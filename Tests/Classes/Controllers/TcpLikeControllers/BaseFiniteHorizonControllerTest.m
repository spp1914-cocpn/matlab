classdef (Abstract) BaseFiniteHorizonControllerTest < BaseTcpLikeControllerTest
    % This class contains basic test cases and additional functions to facilitate
    % testing of finite horizon, sequence-based controllers.
    
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
    
    properties (SetAccess = private, GetAccess = protected)
        horizonLength;
    end
         
    methods (Access=protected)
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            this.horizonLength = 1;
        end
    end
    
    methods (Test)
        %% testDoControlSequenceComputationInvalidTimestep
        function testDoControlSequenceComputationInvalidTimestep(this)
            % call function without a timestep given
            expectedErrId = 'MATLAB:minrhs';
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, this.numModes), ...
                expectedErrId);
            
            expectedErrId = [class(this.controllerUnderTest) ':DoControlSequenceComputation:InvalidTimestep'];
                       
            invalidTimestep = this; % not a scalar
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, this.numModes, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = -1; % scalar, but negative
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, this.numModes, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = this.horizonLength + 1; % scalar, but out of bounds
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, this.numModes, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = 0.5; % scalar, but not an integer
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, this.numModes, invalidTimestep), ...
                expectedErrId);
        end
        
        %% testDoControlSequenceComputationInvalidMode
        function testDoControlSequenceComputationInvalidMode(this)
            % call function without a mode given
            expectedErrId = 'MATLAB:minrhs';
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution), ...
                expectedErrId);
            
            expectedErrId = [class(this.controllerUnderTest) ':DoControlSequenceComputation:InvalidMode'];
                       
            invalidMode = this; % not a scalar
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode, this.horizonLength), ...
                expectedErrId);
            
            invalidMode = 0; % scalar, but non-positive
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode, this.horizonLength), ...
                expectedErrId);
            
            invalidMode = this.numModes + 1; % scalar, but out of bounds
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode, this.horizonLength), ...
                expectedErrId);
            
            invalidMode = 0.5; % scalar, but not an integer
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(this.stateDistribution, invalidMode, this.horizonLength), ...
                expectedErrId);
        end
        
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            appliedInputs = ones(this.dimU, this.horizonLength);
            expectedErrId = [class(this.controllerUnderTest) ':DoCostsComputation:InvalidStateTrajectory'];
            
            invalidStateTrajectory = ones(this.dimX, this.horizonLength + 10); % state trajectory too long
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
            
            invalidStateTrajectory = ones(this.dimX, this.horizonLength); % state trajectory too short
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
        end
        
         %% testDoCostsComputationInvalidInputTrajectory
        function testDoCostsComputationInvalidInputTrajectory(this)
            states = ones(this.dimX, this.horizonLength + 1);
            expectedErrId = [class(this.controllerUnderTest) ':DoCostsComputation:InvalidInputTrajectory'];
            
            invalidInputTrajectory = ones(this.dimU, this.horizonLength + 1); %  % trajectory too long
            this.verifyError(@() this.controllerUnderTest.computeCosts(states, invalidInputTrajectory), ...
                expectedErrId);
        end
  
    end
end

