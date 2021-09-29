classdef NonlinearPlantTest < matlab.unittest.TestCase
    % Test cases for NonlinearPlant.
    
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
    
    properties (Access = private)
        dimX;
        dimU;
        sysEquation;
        W; % sys noise cov
        plantNoise;
        sysInput;
        plantState;
        
        plantUnderTest;
    end
    
    methods (Access = private, Static)
        function newState = sysFunction(stateSamples, input, ~)
            % deterministic system equation
            if isempty(input)
                newState = stateSamples .^ 2;
            else
                newState = stateSamples .^ 2 + [input; input(1, :)];
            end
        end
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.dimX = TestPlantModel.dimX;
            this.dimU = TestPlantModel.dimU;
             
            this.W = gallery('moler', this.dimX);
            this.plantNoise = Gaussian(zeros(this.dimX, 1), this.W);
       
            this.sysInput = [4 2]';
            this.plantState = [2 2 2]';
            
           this.plantUnderTest = TestPlantModel();
        end
    end
    
    methods (Test)        
        %% testNonlinearPlant
        function testNonlinearPlant(this)
            plant = TestPlantModel(this.dimX, this.dimU);
            
            this.verifyEmpty(plant.getSystemInput());
            this.verifyEqual(plant.dimState, this.dimX);
            this.verifyEqual(plant.dimInput, this.dimU);
        end
           
        
        %% testSetSystemInputInvalidInput
        function testSetSystemInputInvalidInput(this)
            expectedErrId = 'NonlinearPlant:SetSystemInput:InvalidSystemInput';
            
            invalidInput = this.sysInput'; % not a column vector
            this.verifyError(@() this.plantUnderTest.setSystemInput(invalidInput), expectedErrId);
            
            invalidInput = [1 nan]'; % not finite
            this.verifyError(@() this.plantUnderTest.setSystemInput(invalidInput), expectedErrId);
        end
        
        %% testSetSystemInput
        function testSetSystemInput(this)
           % should be error free
           this.plantUnderTest.setSystemInput(this.sysInput);
           this.plantUnderTest.setSystemInput([]);
        end
        
         %% testGetSystemInput
        function testGetSystemInput(this)
           % set and immediately retrieve input
           this.plantUnderTest.setSystemInput(this.sysInput);
           this.verifyEqual(this.plantUnderTest.getSystemInput(), this.sysInput);
           
           this.plantUnderTest.setSystemInput([]);
           this.verifyEmpty(this.plantUnderTest.getSystemInput());
        end
        
        %% testDerivative
        function testDerivative(this)
            this.plantUnderTest.setSystemInput(this.sysInput);
            % compute the derivatives and compare with the analytic ones
            nominalState = this.plantState;
            nominalNoise = ones(this.dimX, 1);
            [A, B, G, AA, BB, GG] = this.plantUnderTest.derivative(nominalState, nominalNoise);
            
            expectedA = diag(2 * nominalState);
            expectedB = [1 0 ; 0 1 ; 1 0 ];
            expectedG = zeros(this.dimX);
                        
            this.verifyEqual(A, expectedA);
            this.verifyEqual(B, expectedB);
            this.verifyEqual(G, expectedG);
            
            expectedAA = zeros(this.dimX, this.dimX, this.dimX);
            expectedAA(:, :, 1) = diag([2 0 0]);
            expectedAA(:, :, 2) = diag([0 2 0]);
            expectedAA(:, :, 3) = diag([0 0 2]);
            
            this.verifyEqual(AA, expectedAA);
            this.verifyEqual(BB, zeros(this.dimU, this.dimU, this.dimX));
            this.verifyEqual(GG, zeros(this.dimX, this.dimX, this.dimX));
        end
        
        %% testDerivativeNoInput
        function testDerivativeNoInput(this)
            % compute the derivatives and compare with the analytic ones
            nominalState = this.plantState;
            nominalNoise = ones(this.dimX, 1);
            [A, B, G, AA, BB, GG] = this.plantUnderTest.derivative(nominalState, nominalNoise);
            
            expectedA = diag(2 * nominalState);
            expectedB = [1 0 ; 0 1 ; 1 0 ];
            expectedG = zeros(this.dimX);            
            
            this.verifyEqual(A, expectedA);
            this.verifyEqual(B, expectedB);
            this.verifyEqual(G, expectedG);
            
            expectedAA = zeros(this.dimX, this.dimX, this.dimX);
            expectedAA(:, :, 1) = diag([2 0 0]);
            expectedAA(:, :, 2) = diag([0 2 0]);
            expectedAA(:, :, 3) = diag([0 0 2]);
            
            this.verifyEqual(AA, expectedAA);
            this.verifyEqual(BB, zeros(this.dimU, this.dimU, this.dimX));
            this.verifyEqual(GG, zeros(this.dimX, this.dimX, this.dimX));
        end
        
         %% testSimulateInvalidState
        function testSimulateInvalidState(this)
            expectedErrId = 'NonlinearPlant:Simulate:InvalidSystemState';
            
            invalidState = this.plantState'; % not a column vector
            this.verifyError(@() this.plantUnderTest.simulate(invalidState), expectedErrId);
            
            invalidState = [1 inf 3]'; % not finite
            this.verifyError(@() this.plantUnderTest.simulate(invalidState), expectedErrId);
        end
        
        %% testSimulateNoNoise
        function testSimulateNoNoise(this)
            % first, noise free
            this.assertEmpty(this.plantUnderTest.noise);
            
            % no input set
            expectedState = this.plantUnderTest.simulateForInput(this.plantState, []);
            
            actualState = this.plantUnderTest.simulate(this.plantState);
            this.verifyEqual(actualState, expectedState);
            
            % input set
            this.plantUnderTest.setSystemInput(this.sysInput);
            expectedState = this.plantUnderTest.simulateForInput(this.plantState, this.sysInput);
            
            actualState = this.plantUnderTest.simulate(this.plantState);
            this.verifyEqual(actualState, expectedState);
        end
        
        %% testSimulateWithNoise
        function testSimulateWithNoise(this)
            % first, noise free
            this.assertEmpty(this.plantUnderTest.noise);
            % now set noise, should not affect the results, as system is
            % deterministic
            this.plantUnderTest.setNoise(this.plantNoise);
            this.assertNotEmpty(this.plantUnderTest.noise);
            
            % no input set
            expectedState = this.plantUnderTest.simulateForInput(this.plantState, []);
            
            actualState = this.plantUnderTest.simulate(this.plantState);
            this.verifyEqual(actualState, expectedState);
            
            % input set
            this.plantUnderTest.setSystemInput(this.sysInput);
            expectedState = this.plantUnderTest.simulateForInput(this.plantState, this.sysInput);
            
            actualState = this.plantUnderTest.simulate(this.plantState);
            this.verifyEqual(actualState, expectedState);
        end
        
        %% testSystemEquationNoNoise
        function testSystemEquationNoNoise(this)
                        
            stateSamples = repmat(this.plantState, 1, 4);
            noiseSamples = this.plantNoise.drawRndSamples(4);
            inputSamples = repmat(this.sysInput, 1, 4);
            zeroInputSamples = zeros(this.dimU, 4);
            
            % no input set
            % noise should not affect the result
            expectedStates = this.plantUnderTest.simulateForInput(stateSamples, []);
            
            actualStates = this.plantUnderTest.systemEquation(stateSamples, noiseSamples);
            this.verifyEqual(actualStates, expectedStates);
            
            % input set
            % noise should not affect the result
            this.plantUnderTest.setSystemInput(this.sysInput);
            expectedStates = this.plantUnderTest.simulateForInput(stateSamples, this.sysInput);
            
            actualStates = this.plantUnderTest.systemEquation(stateSamples, noiseSamples);
            this.verifyEqual(actualStates, expectedStates);
        end
    end
end

