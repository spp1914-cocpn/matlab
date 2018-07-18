classdef (Abstract) BaseTcpLikeControllerTest < matlab.unittest.TestCase
    % This class contains basic test cases and additional functions to facilitate
    % testing of sequence-based TCP-like controllers.
    %
    
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
    
    properties (SetAccess = private, GetAccess = protected)
        dimX;
        dimU;
        A;
        B;
        Q;
        R;
        
        delayProbs;
        sequenceLength;
        numModes;
        
        transitionMatrix;
        
        state;
        stateDistribution;
        zeroState;
        zeroStateDistribution;
        
        controllerUnderTest;
     end
    
     methods (Access = protected, Abstract)
        controller = initControllerUnderTest(this);
        initAdditionalProperties(this);
        expectedCosts = computeExpectedCosts(this, states, inputs);
        expectedStageCosts = computeExpectedStageCosts(this, state, input, timestep);        
    end
     
    methods (Sealed, TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.dimX = 3;
            this.dimU = 2;
            
            %this.horizonLength = 1;
            this.sequenceLength = 1; % so we have two modes
            this.numModes = this.sequenceLength + 1;
            this.delayProbs = ones(1, 5) / 5;
            
            % from the delay probs and the 2 modes we get the following
            % transition matrix
            this.transitionMatrix = [1/5 4/5; 1/5 4/5];
                
            this.A = eye(this.dimX);
            this.B = ones(this.dimX, this.dimU);
            this.Q = eye(this.dimX);
            this.R = eye(this.dimU);
              
            this.state = 2 * ones(this.dimX, 1);
            this.stateDistribution = Gaussian(this.state, eye(this.dimX));
            this.zeroState = zeros(this.dimX, 1);
            this.zeroStateDistribution = Gaussian(this.zeroState, eye(this.dimX));
            
            this.initAdditionalProperties();
            
            this.controllerUnderTest = this.initControllerUnderTest();
        end
    end
    
     methods (Test)
         
         %% testDoStageCostsComputation
         function testDoStageCostsComputation(this)
             input = ones(this.dimU, 1);
             timestep = 1;
             
             expectedStageCosts = this.computeExpectedStageCosts(this.state, input, timestep);
             actualStageCosts = this.controllerUnderTest.computeStageCosts(this.state, input, timestep);
             this.verifyEqual(actualStageCosts, expectedStageCosts);    
         end
         
        %% testDoCostsComputation
        function testDoCostsComputation(this)
            states = ones(this.dimX, 2);
            inputs = ones(this.dimU, 1);
           
            expectedCosts = this.computeExpectedCosts(states, inputs);
            actualCosts = this.controllerUnderTest.computeCosts(states, inputs);
            this.verifyEqual(actualCosts, expectedCosts);    
        end
    end
end

