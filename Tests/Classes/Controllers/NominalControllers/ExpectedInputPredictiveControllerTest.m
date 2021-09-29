classdef ExpectedInputPredictiveControllerTest < NominalPredictiveControllerTest
    % Test cases for ExpectedInputPredictiveController.
    
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

    
    properties (Access = private)        
        caDelayProbs;  
    end
    
    methods (Access = protected)
        %% initClosedLoopSystem
        function closedLoopSystem = initClosedLoopSystem(this)
            plant = ss(this.A, this.B, eye(this.dimX), [], []);
            plant.InputDelay = 2;
            % we have a fixed input delay of 2 time steps
            closedLoopSystem = feedback(plant, ss(this.expectedL), +1); % discrete-time model
        end 
        
        %% initControllerUnderTest
        function controller = initControllerUnderTest(this)
            controller = ExpectedInputPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength, this.caDelayProbs);
        end
        
        %% initAdditionalProperties
        function initAdditionalProperties(this)            
            this.caDelayProbs = [0 0 1 0 0]; % we have a fixed delay of 2 time steps           

            this.sequenceLength = 3;
            
            this.horizonLength = this.sequenceLength;            
        end
        
        %% callBaseCtor
        function controller = callBaseCtor(this, varargin)
            if numel(varargin) == 6
                % a setpoint is passed
                controller = ExpectedInputPredictiveController(varargin{1:5}, this.caDelayProbs, varargin{end});
            else
                controller = ExpectedInputPredictiveController(varargin{:}, this.caDelayProbs);
            end
        end
    end
    
    methods (Test)     
        %% testExpectedInputPredictiveControllerInvalidCaDelayProbs
        function testExpectedInputPredictiveControllerInvalidCaDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() ExpectedInputPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, invalidDelayProbs), expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() ExpectedInputPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, invalidDelayProbs, this.setpoint), expectedErrId);
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() ExpectedInputPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, invalidDelayProbs, this.setpoint), expectedErrId);
        end
%%        
%%
        %% testChangeCaDelayProbs
        function testChangeCaDelayProbs(this)
            this.assertTrue(isa(this.controllerUnderTest, 'CaDelayProbsChangeable'));
            
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            
            newCaDelayProbs = [0 1 0 0 0]; % we have a fixed delay of 1 time step now
            this.controllerUnderTest.changeCaDelayProbs(newCaDelayProbs);
            
            % now check if computation of input sequence is affected            
            plant = ss(this.A, this.B, eye(this.dimX), [], []);
            plant.InputDelay = 1;
            % we have a fixed input delay of 1 time steps
            newClosedLoopSystem = feedback(plant, ss(this.expectedL), +1); % discrete-time model          
            
            [~, newInputTrajectory] = this.computeTrajectories(newClosedLoopSystem);
                        
            actualSequence = this.controllerUnderTest.computeControlSequence(state);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            expectedSequence = newInputTrajectory(:);
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence, expectedSequence, 'AbsTol', 9*1e-4);            
        end
        
        %% testChangeSequenceLength
        function testChangeSequenceLength(this)
            expectedErrId = 'ExpectedInputPredictiveController:ChangeSequenceLength:NotSupported';
            
            this.assertEqual(this.controllerUnderTest.sequenceLength, this.sequenceLength);
            
            newSeqLength = this.sequenceLength + 2;
            % the operation is not supported
            this.verifyError(@() this.controllerUnderTest.changeSequenceLength(newSeqLength), ...
                expectedErrId);
        end
        
        %% testDoControlSequenceComputationZeroState
        function testDoControlSequenceComputationZeroState(this)
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));

            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            this.verifyEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));            
        end
        
         %% testDoControlSequenceComputationZeroStateWithSetpoint
        function testDoControlSequenceComputationZeroStateWithSetpoint(this)
            % assert that the set point is changed
            this.controllerUnderTest.changeSetPoint(this.setpoint);
            this.assertEqual(this.controllerUnderTest.setpoint, this.setpoint);
            
            % perform a sanity check: given state is origin, so first
            % element of the sequence should be the feedforward
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));

            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];
            
            this.verifySize(actualSequence, expectedSize);
            this.verifyEqual(actualSequence(1:this.dimU), this.feedforward, ...
                'AbsTol', 9*1e-4);           
        end
        
    end
end

