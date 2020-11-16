classdef InvertedPendulumTest < matlab.unittest.TestCase
    % Test cases for InvertedPendulum.
    
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
        plantNoiseCov;
        noisePsd;
    end
    
    methods (Access = private)
        %% integrateDynamics
        function finalState = integrateDynamics(this, initialState, input)            
            tspan = [0 this.samplingInterval];
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
                
                dxdt(2) = (this.inertia + this.massPendulum * this.length^2)*(input-this.friction*x(2)) + this.massPendulum * this.length * sin(x(3));
                dxdt(2) = dxdt(2) * (this.inertia * x(4)^2+ this.massPendulum *this.length *(this.length * x(4)^2 + 9.81*cos(x(3))));
                dxdt(2) = dxdt(2) / (this.inertia * (this.massCart + this.massPendulum) + this.massPendulum * this.length^2 *(this.massCart + this.massPendulum * (1-cos(x(3))^2)));
                
                dxdt(4) = (this.massCart + this.massPendulum) * (-this.massPendulum * 9.81 * this.length * sin(x(3)) - this.massPendulum * this.length * cos(x(3)) * (input-this.friction*x(2)+this.massPendulum * this.length * x(4) ^ 2*sin(x(3))));
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
            
            q = this.inertia * (this.massPendulum + this.massCart) ...
                + this.massCart * this.massPendulum * this.length^2;
            A_cont = [0 1 0 0;
                0 -(this.inertia + this.massPendulum * this.length ^ 2) * this.friction / q (this.massPendulum ^ 2 * 9.81 * this.length ^ 2)/ q 0;
                0 0 0 1;
                0 -(this.massPendulum * this.length * this.friction) / q this.massPendulum * 9.81 * this.length * (this.massPendulum + this.massCart) / q 0];
            
            this.samplingInterval = 0.0001; % 10 kHz
            
            this.pendulumUnderTest = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval);
            
            this.upwardEquilibrium = [0 0 pi 0]';
            this.downwardEquilibrium = [2 0 0 0]'; % different position though
            
            this.noisePsd = eye(4) * 0.001;
            % corresponding discrete-time noise
            this.plantNoiseCov = integral(@(x) expm(A_cont*x) * this.noisePsd * expm(A_cont'*x), ...
                    0, this.samplingInterval, 'ArrayValued', true);
            this.processNoise = Gaussian(zeros(4, 1), this.plantNoiseCov);            
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
        end       
        
        %% testLinearizeAroundUpwardEquilibriumCont
        function testLinearizeAroundUpwardEquilibriumCont(this)
            this.assertEqual(this.pendulumUnderTest.W_cont, zeros(4));
            
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
            
            % also check with 4 return values
            [A_cont, B_cont, C_cont, psd] = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont();
            actualSys = ss(A_cont, B_cont, C_cont, []);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(psd, zeros(4));
            
            % this time with noise present
            psd = 2 * gallery('moler', 4, 4);            
            this.pendulumUnderTest.W_cont = psd;
            this.assertEqual(this.pendulumUnderTest.W_cont, psd);
            
            [A_cont, B_cont, C_cont, W_cont] = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont();
            actualSys = ss(A_cont, B_cont, C_cont, []);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W_cont, psd);
        end
        
        %% testLinearizeAroundUpwardEquilibriumOneArg
        function testLinearizeAroundUpwardEquilibriumOneArg(this)
            newSamplingInterval = this.samplingInterval / 2;
            % change the sampling interval first
            this.pendulumUnderTest.samplingInterval = newSamplingInterval;
            
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
            
            % first a test without noise acting on the continuous-time
            % linearization
            this.assertEqual(this.pendulumUnderTest.W_cont, zeros(4));
            expectedSys = c2d(ss(A_cont, B_cont, C_cont, []), newSamplingInterval, 'zoh');
                        
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium();
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEmpty(W);
            
            % now add some noise acting on the continuous-time
            % linearization
            W_cont = 2 * gallery('moler', 4, 4);            
            this.pendulumUnderTest.W_cont = W_cont;
            this.assertEqual(this.pendulumUnderTest.W_cont, W_cont);
            expectedSys = c2d(ss(A_cont, B_cont, C_cont, []), newSamplingInterval, 'zoh');
            expectedW = integral(@(x) expm(A_cont*x) * W_cont * expm(A_cont'*x), ...
                0, newSamplingInterval, 'ArrayValued', true);
            
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium();
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEqual(W, expectedW, 'AbsTol', 1e-12);
        end
        
        %% testLinearizeAroundUpwardEquilibriumTwoArgs
        function testLinearizeAroundUpwardEquilibriumTwoArgs(this)            
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
            
            % first a test without noise acting on the continuous-time
            % linearization
            this.assertEqual(this.pendulumUnderTest.W_cont, zeros(4));
            expectedSys = c2d(ss(A_cont, B_cont, C_cont, []), this.samplingInterval, 'zoh');
                        
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(this.samplingInterval);
            actualSys = ss(A, B, C, [], this.samplingInterval);
            this.verifyEqual(actualSys, expectedSys);
            this.verifyEmpty(W);
            
            % now add some noise acting on the continuous-time
            % linearization
            W_cont = 2 * gallery('moler', 4, 4);            
            this.pendulumUnderTest.W_cont = W_cont;
            this.assertEqual(this.pendulumUnderTest.W_cont, W_cont);
            expectedSys = c2d(ss(A_cont, B_cont, C_cont, []), this.samplingInterval, 'zoh');
            expectedW = integral(@(x) expm(A_cont*x) * W_cont * expm(A_cont'*x), ...
                0, this.samplingInterval, 'ArrayValued', true);
            
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(this.samplingInterval);
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
            
            % now test with an initial disturbance of 2 degrees, but without noise
            initialState = this.upwardEquilibrium + [0 0 deg2rad(2) 0]';
            actualState = this.pendulumUnderTest.simulate(initialState);
          
            expectedState = this.integrateDynamics(initialState, 0);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
        
         %% testNonlinearDynamicsWithNoise
        function testNonlinearDynamicsWithNoise(this)
            % add process noise (discrete-time noise)
            this.pendulumUnderTest.setNoise(this.processNoise);
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyNotEqual(actualState, this.upwardEquilibrium); % due to noise
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyNotEqual(actualState, this.downwardEquilibrium);            
        end
    end
end

