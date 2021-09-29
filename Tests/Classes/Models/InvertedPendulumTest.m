classdef InvertedPendulumTest < matlab.unittest.TestCase
    % Test cases for InvertedPendulum.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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

    end
    
    methods (Access = private)
        %% integrateDynamics
        function finalState = integrateDynamics(this, initialState, input, noisePend, noiseAct, tmax)
            if nargin < 6 || isempty(tmax)
                tspan = [0 this.samplingInterval];
            else
                tspan = [0 tmax];
            end
            % integrate dynamics using Runge-Kutta method
            [t,x]=ode45(@fun, tspan, initialState);
            finalState = x(end, :)';            
            
            function dxdt = fun(~, x)
                % the nonlinear pendulum dynamics is time-invariant                
                dxdt = zeros(4,1);
                % first state: cart position
                dxdt(1) = x(2);
                % third state: angle of pendulum rod
                dxdt(3) = x(4);
                
                dxdt(2) = (this.inertia + this.massPendulum * this.length ^ 2) * (input + noiseAct + this.massPendulum * this.length * sin(x(3)) * x(4) ^ 2 - this.friction * x(1));
                dxdt(2) = dxdt(2) + this.massPendulum * this.length ^ 2 * cos(x(3)) * (this.massPendulum * 9.81 * sin(x(3)) - noisePend);
                dxdt(2) = dxdt(2) / ((this.inertia + this.massPendulum * this.length ^ 2) * (this.massCart + this.massPendulum)-this.massPendulum ^2 * this.length ^2 * cos(x(3))^2);                
                              
                dxdt(4) = (this.massCart + this.massPendulum) * this.length * (noisePend - this.massPendulum * 9.81 * sin(x(3)));
                dxdt(4) = dxdt(4) - this.massPendulum * this.length * cos(x(3)) * (input + noiseAct + this.massPendulum * this.length * sin(x(3)) * x(4) ^ 2 - this.friction * x(1));
                dxdt(4) = dxdt(4) / ((this.inertia + this.massPendulum * this.length ^ 2) * (this.massCart + this.massPendulum)-this.massPendulum ^2 * this.length ^2 * cos(x(3))^2);                
            end            
        end
    end
    
    methods (TestMethodSetup)
         %% initProperties
        function initProperties(this)
            this.massCart = 0.5;
            this.massPendulum = 0.5;
            this.friction = 0.1;
            this.inertia = 0.015;
            this.length = 0.3;
                       
            this.samplingInterval = 0.0001; % 10 kHz
            
            this.pendulumUnderTest = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval);
            
            this.upwardEquilibrium = [0 0 pi 0]';
            this.downwardEquilibrium = [2 0 0 0]'; % different position though         
        end
     end
    
    methods (Test)
        %% testInvertedPendulum
        function testInvertedPendulum(this)
            pendulum = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval);
            % check if moment of inertia is computed as expected (J = ml²/3)
            this.verifyEqual(pendulum.inertia, this.inertia);
            this.verifyEqual(pendulum.samplingInterval, this.samplingInterval);
            
            pendulum = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction);
            % check if moment of inertia is computed as expected (J = ml²/3)
            this.verifyEqual(pendulum.inertia, this.inertia);
            this.verifyEqual(pendulum.samplingInterval, 0.001);
            
            % by default, no noise is used
            this.verifyEqual(pendulum.varDisturbanceForcePendulum, 0);
            this.verifyEqual(pendulum.varDisturbanceForceActuator, 0);
            this.verifyEqual(pendulum.varDisturbanceForcePendulumContLin, 0);
            this.verifyEqual(pendulum.varDisturbanceForceActuatorContLin, 0);
            
            % by default, mex is used
            this.verifyTrue(pendulum.useMexImplementation);
        end       
%%
%%
        %% testSetNoiseInvalidNoise
        function testSetNoiseInvalidNoise(this)
            expectedErrId = 'InvertedPendulum:SetNoise:InvalidNoise';
            
            invalidNoise = this; % not a distribution
            this.verifyError(@() this.pendulumUnderTest.setNoise(invalidNoise), expectedErrId);
            
            invalidNoise = Gaussian(0, 1); % invalid dimension
            this.verifyError(@() this.pendulumUnderTest.setNoise(invalidNoise), expectedErrId);
        end
        
        %% testSetNoise
        function testSetNoise(this)
            this.assertEqual(this.pendulumUnderTest.varDisturbanceForcePendulum, 0);
            this.assertEqual(this.pendulumUnderTest.varDisturbanceForceActuator, 0);
            
            varNoise = 0.01;
            noise = Gaussian([0 0]', blkdiag(varNoise, varNoise));
            this.pendulumUnderTest.setNoise(noise);
            
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForcePendulum, varNoise);
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForceActuator, varNoise);
            
            % the noise for the linearization remains unaffected
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForcePendulumContLin, 0);
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForceActuatorContLin, 0);
        end
%%
%%
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
            G1 = [0;
                 (this.massPendulum * this.length ^ 2) / q; 
                 0; 
                 this.length * (this.massPendulum + this.massCart) / q];
            G = [G1, B];
            expectedSys = ss(A, [B G], C, []); % with noise as additional input
            
            % no noise given
            [A_cont, B_cont, C_cont, G_cont, W_cont] = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont();
            actualSys = ss(A_cont, [B_cont G_cont], C_cont, []);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W_cont, zeros(2));
            
            % also check with only disturbance acting on pendulum rod
            varDisturbancePend = 0.01;
            [A_cont, B_cont, C_cont, G_cont, W_cont] = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont(varDisturbancePend);
            actualSys = ss(A_cont, [B_cont G_cont], C_cont, []);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W_cont, blkdiag(varDisturbancePend, 0));
            
            % finally, check with both disturbance forces present
            varDisturbancePend = 0.012;
            varDisturbanceAct = 0.1;
                        
            [A_cont, B_cont, C_cont, G_cont, W_cont] ...
                = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont(varDisturbancePend, varDisturbanceAct);
            actualSys = ss(A_cont, [B_cont G_cont], C_cont, []);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W_cont, blkdiag(varDisturbancePend, varDisturbanceAct));
        end
%%
%%        
        
        %% testLinearizeAroundUpwardEquilibrium
        function testLinearizeAroundUpwardEquilibrium(this)            
            newSamplingInterval = this.samplingInterval / 4;
                        
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
            G1 = [0;
                 (this.massPendulum * this.length ^ 2) / q; 
                 0; 
                 this.length * (this.massPendulum + this.massCart) / q];
            G_cont = [G1, B_cont];
            
            expectedSys = c2d(ss(A_cont, B_cont, C_cont, []), newSamplingInterval, 'zoh');
            
            % without any noise to be used for the linearization
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(newSamplingInterval);
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W, zeros(4));
            
            % also check with only disturbance acting on pendulum rod
            varDisturbancePend = 0.01;            
            expectedW = integral(@(x) expm(A_cont*x) * G_cont * blkdiag(varDisturbancePend, 0) * G_cont' * expm(A_cont'*x), ...
                0, newSamplingInterval, 'ArrayValued', true);
            
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(newSamplingInterval, varDisturbancePend);
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W, expectedW, 'AbsTol', 1e-12);
            
            % finally, check with both disturbance forces present
            varDisturbancePend = 0.012;
            varDisturbanceAct = 0.1;                        
            expectedW = integral(@(x) expm(A_cont*x) * G_cont * blkdiag(varDisturbancePend, varDisturbanceAct) * G_cont' * expm(A_cont'*x), ...
                0, newSamplingInterval, 'ArrayValued', true);
            
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(newSamplingInterval, varDisturbancePend, varDisturbanceAct);
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W, expectedW, 'AbsTol', 1e-12);            
        end
        
        %% testIsValidState
        function testIsValidState(this)
            validState = [0 0 pi 0]';
            validState2 = [-5 0 deg2rad(180-45) 0];
            validState3 = [0 0 deg2rad(180+45) 0]';
            invalidState = [0 0 0 0];
            invalidState2 = [-5.1 0 deg2rad(180-30) 0]';
            invalidState3 = [5.5 0 deg2rad(180+2) 0]';
            invalidState4 = [0 0 pi -inf];
            invalidState5 = [5 nan pi 1e10];
            this.verifyTrue(this.pendulumUnderTest.isValidState(validState));
            this.verifyTrue(this.pendulumUnderTest.isValidState(validState2));
            this.verifyTrue(this.pendulumUnderTest.isValidState(validState3));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState2));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState3));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState4));
            this.verifyFalse(this.pendulumUnderTest.isValidState(invalidState5));
        end
%%
%%
        %% testNonlinearDynamicsNoiseFreeMex
        function testNonlinearDynamicsNoiseFreeMex(this)
            this.assertTrue(this.pendulumUnderTest.useMexImplementation);

            % noise-free
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyEqual(actualState, this.upwardEquilibrium, 'AbsTol', 1e-8);
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyEqual(actualState, this.downwardEquilibrium, 'AbsTol', 1e-8);
            
            % now test with an initial disturbance of 2 degrees, but without noise
            % and an input of 0.1 N
            input = 0.1;
            initialState = this.upwardEquilibrium + [0 0 deg2rad(2) 0]';
            this.pendulumUnderTest.setSystemInput(input);
            actualState = this.pendulumUnderTest.simulate(initialState);
          
            expectedState = this.integrateDynamics(initialState, input, 0, 0);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
        
        %% testNonlinearDynamicsNoiseFreeNoMex
        function testNonlinearDynamicsNoiseFreeNoMex(this)            
            pendulum = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval, false);
            
            this.assertEqual(pendulum.samplingInterval, this.samplingInterval);
            this.assertFalse(pendulum.useMexImplementation);
            
            % noise-free
            actualState = pendulum.simulate(this.upwardEquilibrium);
            this.verifyEqual(actualState, this.upwardEquilibrium, 'AbsTol', 1e-8);
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyEqual(actualState, this.downwardEquilibrium, 'AbsTol', 1e-8);
            
            % now test with an initial disturbance of 2 degrees, but without noise
            % and an input of 0.1 N
            input = 0.1;
            initialState = this.upwardEquilibrium + [0 0 deg2rad(2) 0]';
            pendulum.setSystemInput(input);
            actualState = pendulum.simulate(initialState);
          
            expectedState = this.integrateDynamics(initialState, input, 0, 0);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
%%
%%
        %% testNonlinearDynamicsWithNoiseMex
        function testNonlinearDynamicsWithNoiseMex(this)
            this.assertTrue(this.pendulumUnderTest.useMexImplementation);
            
            % sanity check first: add process noise (discrete-time noise)
            this.pendulumUnderTest.setNoise(Gaussian([0 0]', eye(2)));
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyNotEqual(actualState, this.upwardEquilibrium); % due to noise
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyNotEqual(actualState, this.downwardEquilibrium); 
            
            % now use a "deterministic noise": Dirac mixture with only one component
            mixture = DiracMixture([0.01, -0.12]', 1);
            this.pendulumUnderTest.setNoise(mixture);
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            % compute the expected state, driven by noise
            expectedState = this.integrateDynamics(this.upwardEquilibrium, 0, 0.01, -0.12);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
        
        %% testNonlinearDynamicsWithNoiseNoMex
        function testNonlinearDynamicsWithNoiseNoMex(this)
             pendulum = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval, false);
            
            this.assertEqual(pendulum.samplingInterval, this.samplingInterval);
            this.assertFalse(pendulum.useMexImplementation);
            
            % sanity check first: add process noise (discrete-time noise)
            pendulum.setNoise(Gaussian([0 0]', eye(2)));
            actualState = pendulum.simulate(this.upwardEquilibrium);
            this.verifyNotEqual(actualState, this.upwardEquilibrium); % due to noise
            
            actualState = pendulum.simulate(this.downwardEquilibrium);
            this.verifyNotEqual(actualState, this.downwardEquilibrium); 
            
            % now use a "deterministic noise": Dirac mixture with only one component
            mixture = DiracMixture([0.01, -0.12]', 1);
            pendulum.setNoise(mixture);
            actualState = pendulum.simulate(this.upwardEquilibrium);
            % compute the expected state, driven by noise
            expectedState = this.integrateDynamics(this.upwardEquilibrium, 0, 0.01, -0.12);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);           
        end
        
        %% testNonlinearDynamicsMultipleStepsDisturbance
        function testNonlinearDynamicsMultipleStepsDisturbance(this)
            this.assertTrue(this.pendulumUnderTest.useMexImplementation);
            
            % disturbance modeled as "deterministic noise": Dirac mixture with only one component
            mixture = DiracMixture([0.00001, -0.0002]', 1);
            this.pendulumUnderTest.setNoise(mixture);
            
            numTimeSteps = 0.1 / this.samplingInterval; % 0.1 seconds
            actualState = this.upwardEquilibrium;
            for k=1:numTimeSteps
               actualState = this.pendulumUnderTest.simulate(actualState); 
            end
            expectedState = this.integrateDynamics(this.upwardEquilibrium, 0, 0.00001, -0.0002, 0.1);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-6);
        end
        
        %% testNonlinearDynamicsMultipleStepsDisturbanceNoMex
        function testNonlinearDynamicsMultipleStepsDisturbanceNoMex(this)
            pendulum = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval, false);
            
            this.assertEqual(pendulum.samplingInterval, this.samplingInterval);
            this.assertFalse(pendulum.useMexImplementation);
            
            % disturbance modeled as "deterministic noise": Dirac mixture with only one component
            mixture = DiracMixture([0.00001, -0.0002]', 1);
            pendulum.setNoise(mixture);
            
            numTimeSteps = 0.1 / this.samplingInterval; % 0.1 seconds
            actualState = this.upwardEquilibrium;
            for k=1:numTimeSteps
               actualState = pendulum.simulate(actualState); 
            end
            expectedState = this.integrateDynamics(this.upwardEquilibrium, 0, 0.00001, -0.0002, 0.1);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-6);
        end
    end
end

