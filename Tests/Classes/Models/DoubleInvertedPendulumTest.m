classdef DoubleInvertedPendulumTest < matlab.unittest.TestCase
    % Test cases for DoubleInvertedPendulum.
    
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
        massPendulum1; % mass of first second pendulum
        massPendulum2; % mass of second pendulum
        friction; % friction of the cart
        friction1; % friction of the first joint
        friction2; % friction of the second joint
        
        length1; % length of first pendulum
        length2; % length of second pendulum
        samplingInterval;
        
        pendulumUnderTest;
        
        upwardEquilibrium;
        downwardEquilibrium;
    end
    
    methods (Access = private)
        
        %% integrateDynamics
        function finalState = integrateDynamics(this, initialState, input, noisePend1, noisePend2, noiseAct, tmax)
            if nargin < 7 || isempty(tmax)
                tspan = [0 this.samplingInterval];
            else
                tspan = [0 tmax];
            end
            % integrate dynamics using Runge-Kutta method
            [~,x]=ode45(@fun, tspan, initialState);
            finalState = x(end, :)';            
            
            function dxdt = fun(~, x)
                invM = this.computeInverseMassMatrixForState(x);
                % states in x are ordered as follows: xc, theta1, theta2, then
%               % the corresponding velocities
                sinDiff = sin(x(2) -x(3));
                dxdt = zeros(6, 1);
                dxdt(1:3) = x(4:6); 
                % Coriolis terms
                dxdt(4) = -this.friction * x(4) + this.length1 * (this.massPendulum1 + this.massPendulum2) * sin(x(2)) * (x(5)^2) + this.massPendulum2 * this.length2 * sin(x(3)) * (x(6)^2);
                dxdt(5) = -this.friction1 * x(5) -this.length1 * this.length2 * this.massPendulum2 * sinDiff * (x(6)^2);
                dxdt(6) = -this.friction2 * x(6) + this.length1 * this.length2 * this.massPendulum2 * sinDiff * (x(5) ^2);
                
                % gravity terms
                dxdt(5) = dxdt(5) + 9.81 * this.length1 * (this.massPendulum1 + this.massPendulum2) .* sin(x(2));
                dxdt(6) = dxdt(6) + 9.81 * this.length2 * this.massPendulum2 .* sin(x(3));
                
                % input and noise                
                dxdt(4) = dxdt(4) + input + noiseAct;                
                dxdt(5:6) = dxdt(5:6) + [noisePend1; noisePend2];
                
                dxdt(4:6) = invM * dxdt(4:6);
            end            
        end
        
        %% computeInverseMassMatrixForState
        function invM = computeInverseMassMatrixForState(this, state)
            massPend = this.massPendulum1 + this.massPendulum2;
            invM = zeros(3,3);
            % mass matrix M is symmetric and 3-by-3, so we can compute the inverse in
            % closed form using the determinant, which has also a closed
            % form expression (rule of Sarrus)
            if nargin == 1
                % inverse of mass matrix M for q=0
                % sin(0) = 0
                detM = (this.length1 * this.length2) ^ 2 * this.massPendulum2 * ...
                    (this.massCart * this.massPendulum1);

                % entries of the symmetric mass matrix M (state-dependent)
                M11 = this.massCart + this.massPendulum1 + this.massPendulum2;
                M12 = this.length1 * massPend; %cos(0) = 1
                M13 = this.massPendulum2 * this.length2; %cos(0) = 1
                M22 = this.length1^2 * massPend;
                M23 = this.length1 * this.length2 * this.massPendulum2; %cos(0-0)=1
                M33 = this.length2^2 * this.massPendulum2;
                
                invM(1,1) = (M22 * M33-M23 ^ 2) / detM;
                invM(1,2) = (M13 * M23-M12 * M33) / detM;
                invM(1,3) = (M12 * M23 - M13 * M22) / detM;
                invM(2,1) = invM(1,2);
                invM(2,2) = (M11 * M33 - M13 ^ 2) / detM;
                invM(2,3) = (M13 * M12 - M11 * M23) / detM;
                invM(3,1) = invM(1,3);
                invM(3,2) = invM(2,3);
                invM(3,3) = (M11 * M22 - M12 ^ 2) / detM;
                return
            end            
            
            detM = (this.length1 * this.length2) ^ 2 * this.massPendulum2 * ...
                (this.massCart * this.massPendulum1 + this.massPendulum1^2 * (sin(state(2, :)) ^ 2) ...
                + this.massCart * this.massPendulum2 * (sin(state(2, :)-state(3, :)) ^ 2) + this.massPendulum1 * this.massPendulum2 * (sin(state(2, :)) ^ 2));
            
            % entries of the symmetric mass matrix M (state-dependent)
            M11 = this.massCart + this.massPendulum1 + this.massPendulum2;
            M12 = this.length1 * massPend * cos(state(2, :));
            M13 = this.massPendulum2 * this.length2 * cos(state(3, :));
            M22 = this.length1^2 * massPend;
            M23 = this.length1 * this.length2 * this.massPendulum2 * cos(state(2, :)-state(3, :));
            M33 = this.length2^2 * this.massPendulum2;
         
            invM(1,1) = (M22 * M33-M23 ^ 2) / detM;
            invM(1,2) = (M13 * M23-M12 * M33) / detM;
            invM(1,3) = (M12 * M23 - M13 * M22) / detM;
            invM(2,1) = invM(1,2, :);
            invM(2,2) = (M11 * M33 - M13 ^ 2) / detM;
            invM(2,3) = (M13 * M12 - M11 * M23) / detM;
            invM(3,1) = invM(1,3, :);
            invM(3,2) = invM(2,3, :);
            invM(3,3) = (M11 * M22 - M12 ^ 2) / detM;
        end
        
        %% createExpectedContLinearization
        function contSys = createExpectedContLinearization(this, useNoise)            
            % system dynamics is in standard form M(q)q_dotdot = C(q,q_dot)q_dot + G(q) + H*u
            % H = [1 0 0]'
            % C contains Coriolis terms and damping factors
            % G contains gravity terms
            % linearization in state space form is given by
            % x_dot = [0, I; M(0)^-1*dG(0)/dx, M(0)^-1*dC(0)/dx]*x + [0; M(0)^-1*H]*u
            
            H = [1 0 0]';
            % inverse of mass matrix M for q=0
            invM = this.computeInverseMassMatrixForState();
            
            % derivative of G w.r.t. q, evaluated at 0
            % sin'=cos, cos(0)=1
            dG = zeros(3);
            dG(2,2) = 9.81 * (this.massPendulum1 + this.massPendulum2) * this.length1; 
            dG(3,3) = 9.81 * this.length2 * this.massPendulum2;
                       
            % derivative of C w.r.t. q_dot, evaluated at 0
            % linear friction/damping
            % derivatives w.r.t. q do not matter since evaluated at q_dot=0, q=0
            dC = zeros(3);
            dC(1,1) = -this.friction;
            dC(2,2) = -this.friction1;
            dC(3,3) = -this.friction2;
            
            A = zeros(6);
            A(1:3, 4:6) = eye(3);
            A(4:6, 1:3) = invM * dG;
            A(4:6, 4:6) = invM * dC;
            
            B = [zeros(3,1); invM*H];
            
            C = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 1 -1 0 0 0];
            
            if useNoise
                % noise matrix
                % we must consider that first noise variable w1 acts and lower pendulum
                % w2 acts on upper pendulum and w3 acts on cart
                % so last column of G equals B
                G = [zeros(3,3); invM * [0 0 1; 0 1 0; 1 0 0]];
                contSys = ss(A, [B G], C, []); % with noise as additional input
            else
                contSys = ss(A, B, C, []); % with noise as additional input  
            end
        end
    end
    
    methods (TestMethodSetup)
         %% initProperties
        function initProperties(this)
            this.massCart = 1.7;
            this.massPendulum1 = 1.5;
            this.massPendulum2 = 1.3;
            this.friction = 1.2;
            this.friction1 = 0.7;
            this.friction2 = 0.3;
            
            this.length1 = 0.75;
            this.length2 = 0.35;
                       
            this.samplingInterval = 0.0001; % 10 kHz
            
            this.pendulumUnderTest = DoubleInvertedPendulum(this.massCart, this.massPendulum1, ...
                this.massPendulum2, this.length1, this.length2, this.friction, ...
                this.friction1, this.friction2, this.samplingInterval);
            
            this.upwardEquilibrium = [0 0 0 0 0 0]'; % unstable
            this.downwardEquilibrium = [2 pi pi 0 0 0]'; % stable, different position though         
        end
    end
    
    methods (Test)
        %% testDoubleInvertedPendulum
        function testDoubleInvertedPendulum(this)
            pendulum = DoubleInvertedPendulum(this.massCart, this.massPendulum1, ...
                this.massPendulum2, this.length1, this.length2, this.friction, ...
                this.friction1, this.friction2, this.samplingInterval);
            
            % check if ode options are set correctly
            this.verifyNotEmpty(pendulum.odeOpts);
            this.verifyFalse(pendulum.odeOpts.MassSingular);
            this.verifyEqual(pendulum.odeOpts.MStateDependence, 'strong');
            
            % compute the mass matrix for the stable equilibrium using the function stored in option            
            massMatrix = pendulum.odeOpts.Mass(0, this.downwardEquilibrium);
            % and validate
            expectedMassMatrix = zeros(6);
            expectedMassMatrix(1,1) = 1;
            expectedMassMatrix(2,2) = 1;
            expectedMassMatrix(3,3) = 1;
            expectedMassMatrix(4,4) = this.massCart + this.massPendulum1 + this.massPendulum2;
            expectedMassMatrix(4,5) = -this.length1 * (this.massPendulum1 + this.massPendulum2); % cos(pi) = -1
            expectedMassMatrix(4,6) = -this.massPendulum2 * this.length2; % cos(pi) = 1
            expectedMassMatrix(5,4) = expectedMassMatrix(4,5);
            expectedMassMatrix(5,5) = this.length1^2 * (this.massPendulum1 + this.massPendulum2);
            expectedMassMatrix(5,6) = this.length1 * this.length2 * this.massPendulum2; % cos(pi-pi)=cos(0)=1
            expectedMassMatrix(6,4) = expectedMassMatrix(4,6);
            expectedMassMatrix(6,5) = expectedMassMatrix(5,6);
            expectedMassMatrix(6,6) = this.length2^2*this.massPendulum2;            
            
            this.verifyTrue(issymmetric(massMatrix));
            this.verifyEqual(massMatrix, expectedMassMatrix);
            
            % by default, no noise is used
            this.verifyEqual(pendulum.varDisturbanceForcePendulum1, 0);
            this.verifyEqual(pendulum.varDisturbanceForcePendulum2, 0);
            this.verifyEqual(pendulum.varDisturbanceForceActuator, 0);
            this.verifyEqual(pendulum.varDisturbanceForcePendulum1ContLin, 0);
            this.verifyEqual(pendulum.varDisturbanceForcePendulum2ContLin, 0);
            this.verifyEqual(pendulum.varDisturbanceForceActuatorContLin, 0);
            
            pendulum = DoubleInvertedPendulum(this.massCart, this.massPendulum1, ...
                this.massPendulum2, this.length1, this.length2, this.friction, ...
                this.friction1, this.friction2);
            
            this.verifyEqual(pendulum.samplingInterval, 1/1000);
            
            % by default, mex is used
            this.verifyTrue(pendulum.useMexImplementation);
        end       
%%
%%
        %% testSetNoiseInvalidNoise
        function testSetNoiseInvalidNoise(this)
            expectedErrId = 'DoubleInvertedPendulum:SetNoise:InvalidNoise';
            
            invalidNoise = this; % not a distribution
            this.verifyError(@() this.pendulumUnderTest.setNoise(invalidNoise), expectedErrId);
            
            invalidNoise = Gaussian(0, 1); % invalid dimension
            this.verifyError(@() this.pendulumUnderTest.setNoise(invalidNoise), expectedErrId);
        end
        
        %% testSetNoise
        function testSetNoise(this)
            this.assertEqual(this.pendulumUnderTest.varDisturbanceForcePendulum1, 0);
            this.assertEqual(this.pendulumUnderTest.varDisturbanceForcePendulum2, 0);
            this.assertEqual(this.pendulumUnderTest.varDisturbanceForceActuator, 0);
            
            varNoise = 0.01;
            noise = Gaussian([0 0 0]', blkdiag(varNoise, varNoise, varNoise));
            this.pendulumUnderTest.setNoise(noise);
            
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForcePendulum1, varNoise);
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForcePendulum2, varNoise);
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForceActuator, varNoise);
            
            % the noise for the linearization remains unaffected
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForcePendulum1ContLin, 0);
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForcePendulum2ContLin, 0);
            this.verifyEqual(this.pendulumUnderTest.varDisturbanceForceActuatorContLin, 0);
        end
%%
%%
        %% testLinearizeAroundUpwardEquilibriumCont
        function testLinearizeAroundUpwardEquilibriumCont(this)
            
            expectedSys = this.createExpectedContLinearization(true);
            
            % no noise given
            [A_cont, B_cont, C_cont, G_cont, W_cont] = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont();      
            actualSys = ss(A_cont, [B_cont G_cont], C_cont, []);
            this.verifyEqual(actualSys.A, expectedSys.A, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.B, expectedSys.B, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.C, expectedSys.C, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.D, expectedSys.D, 'AbsTol', 1e-10);
            this.verifyEqual(W_cont, zeros(3));
            
            % also check with only disturbance acting on first pendulum rod
            varDisturbancePend1 = 0.01;
            [A_cont, B_cont, C_cont, G_cont, W_cont] = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont(varDisturbancePend1);
            actualSys = ss(A_cont, [B_cont G_cont], C_cont, []);
            this.verifyEqual(actualSys.A, expectedSys.A, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.B, expectedSys.B, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.C, expectedSys.C, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.D, expectedSys.D, 'AbsTol', 1e-10);
            this.verifyEqual(W_cont, blkdiag(varDisturbancePend1, 0, 0));
                        
            % finally, check with all disturbance forces present            
            varDisturbancePend2 = 0.014;
            varDisturbanceAct = 0.1;
                        
            [A_cont, B_cont, C_cont, G_cont, W_cont] ...
                = this.pendulumUnderTest.linearizeAroundUpwardEquilibriumCont(varDisturbancePend1, varDisturbancePend2, varDisturbanceAct);
            actualSys = ss(A_cont, [B_cont G_cont], C_cont, []);
            this.verifyEqual(actualSys.A, expectedSys.A, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.B, expectedSys.B, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.C, expectedSys.C, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.D, expectedSys.D, 'AbsTol', 1e-10);
            this.verifyEqual(W_cont, blkdiag(varDisturbancePend1, varDisturbancePend2, varDisturbanceAct));
        end
%%
%%
        %% testLinearizeAroundUpwardEquilibrium
        function testLinearizeAroundUpwardEquilibrium(this)            
            newSamplingInterval = this.samplingInterval / 4;
            
            contSysWithNoise = this.createExpectedContLinearization(true);
            contSys = this.createExpectedContLinearization(false);                     
            expectedSys = c2d(contSys, newSamplingInterval, 'zoh');
            
            % without any noise to be used for the linearization
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(newSamplingInterval);
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys.A, expectedSys.A, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.B, expectedSys.B, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.C, expectedSys.C, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.D, expectedSys.D, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.Ts, expectedSys.Ts);
            this.verifyEqual(W, zeros(6));
        
            A_cont = contSys.A;
            G_cont = contSysWithNoise.B(:, 2:end);
            
            % also check with only disturbance acting on lower pendulum rod
            varDisturbancePend1 = 0.01;            
            expectedW = integral(@(x) expm(A_cont*x) * G_cont * blkdiag(varDisturbancePend1, 0, 0) * G_cont' * expm(A_cont'*x), ...
                0, newSamplingInterval, 'ArrayValued', true);
            
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(newSamplingInterval, varDisturbancePend1);
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys.A, expectedSys.A, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.B, expectedSys.B, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.C, expectedSys.C, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.D, expectedSys.D, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.Ts, expectedSys.Ts);
            this.verifyTrue(issymmetric(W));
            this.verifyEqual(W, expectedW, 'AbsTol', 1e-12);
            
            % finally, check with all disturbance forces present
            varDisturbancePend2 = 0.012;
            varDisturbanceAct = 0.1;                        
            expectedW = integral(@(x) expm(A_cont*x) * G_cont * blkdiag(varDisturbancePend1, varDisturbancePend2, varDisturbanceAct) * G_cont' * expm(A_cont'*x), ...
                0, newSamplingInterval, 'ArrayValued', true);
            
            [A, B, C, W] = this.pendulumUnderTest.linearizeAroundUpwardEquilibrium(newSamplingInterval, varDisturbancePend1, varDisturbancePend2, varDisturbanceAct);
            actualSys = ss(A, B, C, [], newSamplingInterval);
            this.verifyEqual(actualSys.A, expectedSys.A, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.B, expectedSys.B, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.C, expectedSys.C, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.D, expectedSys.D, 'AbsTol', 1e-10);
            this.verifyEqual(actualSys.Ts, expectedSys.Ts);
            this.verifyTrue(issymmetric(W));
            this.verifyEqual(W, expectedW, 'AbsTol', 1e-12);            
        end
%%
%%
        %% testIsValidState
        function testIsValidState(this)
            validState = [0 0 0 0 0 0]';
            validState2 = [-5 deg2rad(45) 0 0 0 0];
            validState3 = [0 deg2rad(45) 0 0 0 0]';
            invalidState = [0 0 deg2rad(100) 0 0 0];
            invalidState2 = [-5.1 0 deg2rad(45) 0 0 0]';
            invalidState3 = [5.5 deg2rad(180+2) 0 0 0 0]';
            invalidState4 = [0 0 0 0 0 -inf];
            invalidState5 = [5 nan deg2rad(45) 1e10 1e12 1e13];
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
            this.assertEqual(this.pendulumUnderTest.samplingInterval, this.samplingInterval);            
            this.assertTrue(this.pendulumUnderTest.useMexImplementation);
            
            % noise-free
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyEqual(actualState, this.upwardEquilibrium, 'AbsTol', 1e-8);
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyEqual(actualState, this.downwardEquilibrium, 'AbsTol', 1e-8);
            
            % now test with an initial disturbance of 2 for the upper pendulum degrees, but without noise
            % and an input of 0.1 N
            input = 0.1;
            initialState = this.upwardEquilibrium + [0 0 deg2rad(2) 0 0 0]';
            this.pendulumUnderTest.setSystemInput(input);
            actualState = this.pendulumUnderTest.simulate(initialState);
          
            expectedState = this.integrateDynamics(initialState, input, 0, 0, 0);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
            
            % finally, for a given initial state != 0
            input = 0.1;
            initialState = [0.1 deg2rad(-20) deg2rad(12) -0.25 deg2rad(2) deg2rad(-15)]';
            this.pendulumUnderTest.setSystemInput(input);
            actualState = this.pendulumUnderTest.simulate(initialState);
            
            expectedState = this.integrateDynamics(initialState, input, 0, 0, 0);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
        
        %% testNonlinearDynamicsNoiseFreeNoMex
        function testNonlinearDynamicsNoiseFreeNoMex(this)
            pendulum = DoubleInvertedPendulum(this.massCart, this.massPendulum1, ...
                this.massPendulum2, this.length1, this.length2, this.friction, ...
                this.friction1, this.friction2, this.samplingInterval, false);
            
            this.assertEqual(pendulum.samplingInterval, this.samplingInterval);            
            this.assertFalse(pendulum.useMexImplementation);
            
            % noise-free
            actualState = pendulum.simulate(this.upwardEquilibrium);
            this.verifyEqual(actualState, this.upwardEquilibrium, 'AbsTol', 1e-8);
            
            actualState = pendulum.simulate(this.downwardEquilibrium);
            this.verifyEqual(actualState, this.downwardEquilibrium, 'AbsTol', 1e-8);
            
            % now test with an initial disturbance of 2 for the upper pendulum degrees, but without noise
            % and an input of 0.1 N
            input = 0.1;
            initialState = this.upwardEquilibrium + [0 0 deg2rad(2) 0 0 0]';
            pendulum.setSystemInput(input);
            actualState = pendulum.simulate(initialState);
          
            expectedState = this.integrateDynamics(initialState, input, 0, 0, 0);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
            
            % finally, for a given initial state != 0
            input = 0.1;
            initialState = [0.1 deg2rad(-20) deg2rad(12) -0.25 deg2rad(2) deg2rad(-15)]';
            pendulum.setSystemInput(input);
            actualState = pendulum.simulate(initialState);
            
            expectedState = this.integrateDynamics(initialState, input, 0, 0, 0);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
%%
%%
        %% testNonlinearDynamicsWithNoiseMex
        function testNonlinearDynamicsWithNoiseMex(this)
            this.assertEqual(this.pendulumUnderTest.samplingInterval, this.samplingInterval);  
            this.assertTrue(this.pendulumUnderTest.useMexImplementation);
            
            % sanity check first: add process noise
            this.pendulumUnderTest.setNoise(Gaussian([0 0 0]', eye(3)));
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyNotEqual(actualState, this.upwardEquilibrium); % due to noise
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyNotEqual(actualState, this.downwardEquilibrium); 
            
            % now use a "deterministic noise": Dirac mixture with only one component
            mixture = DiracMixture([0.01, -0.12, 0.02]', 1);
            this.pendulumUnderTest.setNoise(mixture);
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            % compute the expected state, driven by noise
            expectedState = this.integrateDynamics(this.upwardEquilibrium, 0, 0.01, -0.12, 0.02);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
        
        %% testNonlinearDynamicsWithNoiseNoMex
        function testNonlinearDynamicsWithNoiseNoMex(this)
             pendulum = DoubleInvertedPendulum(this.massCart, this.massPendulum1, ...
                this.massPendulum2, this.length1, this.length2, this.friction, ...
                this.friction1, this.friction2, this.samplingInterval, false);
            
            this.assertEqual(pendulum.samplingInterval, this.samplingInterval);            
            this.assertFalse(pendulum.useMexImplementation);
            
            % sanity check first: add process noise
            pendulum.setNoise(Gaussian([0 0 0]', eye(3)));
            actualState = pendulum.simulate(this.upwardEquilibrium);
            this.verifyNotEqual(actualState, this.upwardEquilibrium); % due to noise
            
            actualState = pendulum.simulate(this.downwardEquilibrium);
            this.verifyNotEqual(actualState, this.downwardEquilibrium); 
            
            % now use a "deterministic noise": Dirac mixture with only one component
            mixture = DiracMixture([0.01, -0.12, 0.02]', 1);
            pendulum.setNoise(mixture);
            actualState = pendulum.simulate(this.upwardEquilibrium);
            % compute the expected state, driven by noise
            expectedState = this.integrateDynamics(this.upwardEquilibrium, 0, 0.01, -0.12, 0.02);
            this.verifyEqual(actualState, expectedState, 'AbsTol', 1e-8);
        end
    end
end

