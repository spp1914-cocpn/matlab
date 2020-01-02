classdef InvertedPendulumTest < matlab.unittest.TestCase
    % Test cases for InvertedPendulum.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        massCart; % mass of cart
        massPendulum; % mass of pendulum
        friction; % friction of the cart
        inertia; % moment of inertia of the pendulum
    
        length; % length of pendulum
        samplingInterval;
        
        pendulumUnderTest;
        
        upwardEquilibrium;
        downwardEquilibrium;
        
        processNoise;
    end
    
     methods (TestMethodSetup)
        function initProperties(this)
            this.massCart = 0.5;
            this.massPendulum = 0.5;
            this.friction = 0.1;
            this.inertia = 0.015;
            this.length = 0.3;
            
            this.samplingInterval = 0.01;
            
            this.pendulumUnderTest = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval);
            
            this.upwardEquilibrium = [0 0 pi 0]';
            this.downwardEquilibrium = [2 0 0 0]'; % different position though
                        
            plantNoiseCov = eye(4) * 0.001;
            this.processNoise = Gaussian(zeros(4, 1), plantNoiseCov);            
        end
     end
    
    methods (Test)
        %% testInvertedPendulum
        function testInvertedPendulum(this)
            pendulum = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval);
            % check if moment of inertia is computed as expected (J = ml²/3)
            this.verifyEqual(pendulum.inertia, this.inertia);
        end       
        
        %% testLinearizeAroundUpwardEquilibriumCont
        function testLinearizeAroundUpwardEquilibriumCont(this)
            q = this.inertia * (this.massPendulum + this.massCart) ...
                + this.massCart * this.massPendulum * this.length^2; % denominator for the A and B matrices
            A = [0 1 0 0;
                0 -(this.inertia + this.massPendulum * this.length ^ 2) * this.friction / q (this.massPendulum ^ 2 * 9.81 * this.length ^ 2)/ q 0;
                0 0 0 1;
                0 -(this.massPendulum * this.length * this.friction) / q this.massPendulum * 9.81 * this.length * (this.massPendulum + this.massCart) / q 0];
            B = [0;
                (this.inertia + this.massPendulum * this.length ^ 2) / q;
                0;
                this.massPendulum * this.length / q];
            C = [1 0 0 0;
                 0 0 1 0];
            
            expectedSys = ss(A, B, C, []);
            [A_cont, B_cont, C_cont] = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont();
            actualSys = ss(A_cont, B_cont, C_cont, []);
            this.verifyEqual(actualSys, expectedSys);
        end
        
         %% testLinearizeAroundUpwardEquilibrium
        function testLinearizeAroundUpwardEquilibrium(this)
            % use c2d to obtain a discretized dynamics
            q = this.inertia * (this.massPendulum + this.massCart) ...
                + this.massCart * this.massPendulum * this.length^2; % denominator for the A and B matrices
            A_cont = [0 1 0 0;
                0 -(this.inertia + this.massPendulum * this.length ^ 2) * this.friction / q (this.massPendulum ^ 2 * 9.81 * this.length ^ 2)/ q 0;
                0 0 0 1;
                0 -(this.massPendulum * this.length * this.friction) / q this.massPendulum * 9.81 * this.length * (this.massPendulum + this.massCart) / q 0];
            B_cont = [0;
                (this.inertia + this.massPendulum * this.length ^ 2) / q;
                0;
                this.massPendulum * this.length / q];
            C_cont = [1 0 0 0;
                 0 0 1 0];
            
            W_cont = 2 * eye(4);
            expectedSys = c2d(ss(A_cont, B_cont, C_cont, []), this.samplingInterval, 'zoh');
            expectedW = integral(@(x) expm(A_cont*x) * W_cont * expm(A_cont'*x), ...
                0, this.samplingInterval, 'ArrayValued', true);
            
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(W_cont);
            actualSys = ss(A, B, C, [], this.samplingInterval);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W, expectedW, 'AbsTol', 1e-12);
        end
        
         
        %% testIsValidStateError
        function testIsValidStateError(this)
            expectedErrId = 'InvertedPendulum:IsValidState:InvalidSystemState';
            
            errState = this; % not a vector
            this.verifyError(@() this.pendulumUnderTest.isValidState(errState), expectedErrId);
            
            errState = [1 1 1]'; % wrong dimension
            this.verifyError(@() this.pendulumUnderTest.isValidState(errState), expectedErrId);
            
            errState = [1 1 1 1]; % not a col vector
            this.verifyError(@() this.pendulumUnderTest.isValidState(errState), expectedErrId);
            
            errState = [1 1 inf 1]'; % not finite
            this.verifyError(@() this.pendulumUnderTest.isValidState(errState), expectedErrId);
        end
        
        %% testIsValidState
        function testIsValidState(this)
            validState = [0 0 pi 0]';
            validState2 = [-5 0 deg2rad(180-45) 0]';
            validState3 = [0 0 deg2rad(180+45) 0]';
            invalidState = [0 0 0 0]';
            invalidState2 = [-5.1 0 deg2rad(180-30) 0]';
            invalidState3 = [5.5 0 deg2rad(180+2) 0]';
            this.verifyTrue(this.pendulumUnderTest.isValidState(validState));
            this.verifyTrue(this.pendulumUnderTest.isValidState(validState2));
            this.verifyTrue(this.pendulumUnderTest.isValidState(validState3));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState2));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState3));
        end
        
        %% testNonlinearDynamicsNoiseFree
        function testNonlinearDynamicsNoiseFree(this)
            % noise-free
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyEqual(actualState, this.upwardEquilibrium, 'AbsTol', 1e-8);
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyEqual(actualState, this.downwardEquilibrium, 'AbsTol', 1e-8);
            
        end
        
         %% testNonlinearDynamicsWithNoise
        function testNonlinearDynamicsWithNoise(this)
            % add process noise
            this.pendulumUnderTest.setNoise(this.processNoise);
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyNotEqual(actualState, this.upwardEquilibrium);
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyNotEqual(actualState, this.downwardEquilibrium);            
        end
    end
end

