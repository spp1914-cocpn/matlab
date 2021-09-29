classdef DoubleInvertedPendulum < NonlinearPlant
    % Implementation of a double inverted pendulum on a cart with  
    % nonlinear dynamics according to M(x)x_dot = f(x,u,w), 
    % where M(x) is a nonsingular, state-dependent mass matrix,
    % x=[q;q_dot] is the 6-dimensional state vector consisting of
    % the generalized coordinates q and their time derivatives.
    % q consists of cart position (in meters) and the two pendulum angles 
    % (upwards, clockwise direction corresponds to positive values, in radians).
    % u is the control input (force, in newtons) applied to move the cart,
    % and w is an optional process noise, modeling generalized disturbance forces (in
    % newton-meters and newtons, respectively) acting in directions of each of the generalized coordinates.
    % Note that the state was chosen such that the unstable upward
    % equilibrium corresponds to pendulum angles of 0 rad.
    % Moreover, the pendulum rods are assumed to be massless, with masses located on either top.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Constant, Access = private)
        graviationalAcceleration = 9.81; % in N/kg
        defaultInternalSamplingRate = 0.001; % 1 kHz
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        massCart;
        massPendulum1;
        massPendulum2;
        lengthPendulum1;
        lengthPendulum2;
        frictionCart;
        frictionPendulum1;
        frictionPendulum2;
    end
    
    properties (SetAccess = immutable, GetAccess = ?DoubleInvertedPendulumTest)
        odeOpts;
    end
    
    properties (SetAccess = private, GetAccess = public)
        %
        % variance of disturbance force driving the rod of pendulum 1 (entering at the pendulum tip)
        % default zero
        varDisturbanceForcePendulum1 (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        % variance of disturbance force driving the rod of pendulum 2 (entering at the pendulum tip)
        % default zero
        varDisturbanceForcePendulum2 (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        % variance of disturbance force driving the cart (entering at the actuator, actuation noise)
        % default zero
        varDisturbanceForceActuator (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        %        
    end
    
    properties (Access = public)
        samplingInterval(1,1) double {mustBePositive, mustBeFinite} = DoubleInvertedPendulum.defaultInternalSamplingRate;
        
        % we need noise terms for the continuous-time linearization of the model
        % can be different from the ones used to simulate the
        % discrete-time nonlinear dynamics
        % the controllers use discrete-time linearizations, the noise in
        % these models is dependent on the sampling rate
        % the terms are assumed to be independent of each other
        %
        varDisturbanceForcePendulum1ContLin (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        varDisturbanceForcePendulum2ContLin (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        varDisturbanceForceActuatorContLin (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
    end
    
    methods (Access = public)
        %% DoubleInvertedPendulum
        function this = DoubleInvertedPendulum(massCart, massPendulum1, massPendulum2, lengthPendulum1, lengthPendulum2, ...
                frictionCart, frictionPendulum1, frictionPendulum2, samplingInterval)
            % Class constructor.
            %
            % Parameters:
            %   >> massCart (Positive Scalar)
            %      A positive scalar denoting the mass (in kg) of the cart.
            %
            %   >> massPendulum1 (Positive scalar)
            %      A positive scalar denoting the mass of the bob (in kg) located at the top of the lower pendulum.
            %
            %   >> massPendulum2 (Positive scalar)
            %      A positive scalar denoting the mass of the bob (in kg) located at the top of the upper pendulum.
            %
            %   >> lengthPendulum1 (Positive scalar)
            %      A positive scalar denoting the length the lower pendulum rod (in m).           
            %
            %   >> lengthPendulum2 (Positive scalar)
            %      A positive scalar denoting the length the upper pendulum rod (in m).
            %
            %   >> frictionCart (Positive scalar)
            %      A positive scalar denoting the coefficient of viscous friction for the cart (in Ns/m).
            %
            %   >> frictionPendulum1(Positive scalar)
            %      A nonnegative scalar denoting the coefficient of viscous friction for the joint at the lower pendulum (in N*s*m).
            %
            %   >> frictionPendulum2(Positive scalar)
            %      A nonnegative scalar denoting the coefficient of viscous friction for the joint at the upper pendulum (in N*s*m).
            %
            %   >> samplingInterval (Positive scalar, optional)
            %      The sampling interval (in seconds) that shall be used for the numerical simulation of the
            %      underlying ODE (the continuous-time plant model).
            %      If left out, the default vale 0.001 is used, which
            %      corresponds to a sampling rate of 1 kHz.         
            %
            % Returns:
            %   << this (DoubleInvertedPendulum)
            %      A new DoubleInvertedPendulum instance.
            
            this@NonlinearPlant(6, 1);
            
            this.massCart = massCart;
            this.massPendulum1 = massPendulum1;
            this.massPendulum2 = massPendulum2;
            this.lengthPendulum1 = lengthPendulum1;
            this.lengthPendulum2 = lengthPendulum2;
            this.frictionCart = frictionCart;
            this.frictionPendulum1 = frictionPendulum1;
            this.frictionPendulum2 = frictionPendulum2;
            
            if nargin > 8
                this.samplingInterval = samplingInterval;
            end
            
            this.odeOpts = odeset('Mass', @(t, y) blkdiag(eye(3), this.createMassMatrix(y)));
            this.odeOpts.MStateDependence = 'strong';
            this.odeOpts.MassSingular = false;
        end
        
        %% setNoise
        function setNoise(this, noise)
            % Set the additional noise driving the pendulum.
            % The noise driving the pendulum consists of two, typically
            % independent, components: A disturbance force driving the
            % pendulum rod, which enters at both the pendulum tips, and an
            % actuation noise acting on the force applied to the cart.
            % Thus, the 3-dimensional joint distribution is passed here
            % with first two components the marginal distribution of the
            % disturbance force (more precisely moments) acting on the pendulum tips. 
            %
            % Parameters:
            %   >> noise (Subclass of Distribution)
            %      The 3-dimensional system noise.
            
            assert(Checks.isClass(noise, 'Distribution') && noise.getDim() == 3, ...
                'DoubleInvertedPendulum:SetNoise:InvalidNoise', ...
                '** <noise> must be a 3-dimensional distribution **');                
            
            setNoise@SystemModel(this, noise);
                
            [~, noiseCov] = noise.getMeanAndCov();
            % store for convenience
            this.varDisturbanceForceActuator = noiseCov(3, 3);
            this.varDisturbanceForcePendulum1 = noiseCov(1, 1);
            this.varDisturbanceForcePendulum2 = noiseCov(2, 2);
        end
        
        %% linearizeAroundUpwardEquilibriumCont
        function [A_cont, B_cont, C_cont, G_cont, W_cont] = linearizeAroundUpwardEquilibriumCont(this, linContDisturbanceForcePend1, linContDisturbanceForcePend2, linContDisturbanceForceAct)
            % Obtain a continuous-time linearization of the pendulum dynamics around the unstable upward equilibrium 
            % which corresponds to the state [0 0 0 0 0 0]'.
            % The linearization is a differential equation of the form x_dot = A_cont * x + B_cont * u + G_cont * w 
            % and y = C_cont * x is the output.
            % The state x is given by [q, q_dot]' 
            % with q = [position of the cart; deviation of lower pendulum from upward equilibrium; deviation of upper pendulum from upward equilibrium].            
            % The input u is the force (in newtons) applied to move the cart.
            % The noise w=[w1 w2 w3] is composed of disturbances (torques) acting at the
            % tips of both pendulum rods (w1 and w2) and a disturbance force (w3) acting on the
            % cart (additive on u, actuation noise). All three are assumed white
            % and mutually independent, the returned noise covariance W_cont is thus diagonal.
            %
            % Parameters:
            %   >> linContDisturbanceForcePend1 (Nonnegative scalar, optional)
            %      The variance of the disturbance (torque), in newton-meters
            %      squared, entering at the tip of the lower pendulum rod.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForcePendulum1ContLin> property is used.
            %      If no disturbance torque is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            %   >> linContDisturbanceForcePend2 (Nonnegative scalar, optional)            
            %      The variance of the disturbance (torque), in newton-meters
            %      squared, entering at the tip of the upper pendulum rod.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForcePendulum2ContLin> property is used.
            %      If no disturbance torque is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            %   >> linContDisturbanceForceAct (Nonnegative scalar, optional)
            %      The variance of the disturbance (force), in newtons squared,
            %      addtionally entering at the actuator to drive the cart.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForceActuatorContLin> property is used.
            %      If no disturbance force is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            % Returns:
            %   << A_cont (6-by-6 matrix)
            %      The system matrix of the continuous-time linearization.
            %
            %   << B_cont (6-dimensional column vector)
            %      The input matrix of the continuous-time linearization.
            %
            %   << C_cont (3-by-6 matrix)
            %      The measurement matrix, assuming that x_1 (position of
            %      the cart) and x_3 (deviation of pendulum angle from
            %      upward equilibrium) are measured.
            %
            %   << G_cont (6-by-3 matrix)
            %      The system noise matrix, as described above, of the continuous-time linearization.
            %
            %   << W_cont (3-by-3 matrix)
            %      The covariance of the white noise acting on the
            %      continuous-time linearized dynamics as described above.
            %      Note: If either of the noise components is left out (variance = 0), this matrix is singular.
            
            arguments
                this
                linContDisturbanceForcePend1 (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForcePendulum1ContLin;
                linContDisturbanceForcePend2 (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForcePendulum2ContLin;
                linContDisturbanceForceAct (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceActuatorContLin;
            end
            
            % we have M*y''+C*y'+g*y=u+w
            % hence, as M is invertible, y''+M_inv*Cy'+M_inv*g*y=M_inv*u+M_inv*w
            massPend = this.massPendulum1 + this.massPendulum2;            
       
            % compute based on linearized continuous-time dynamics
            A_cont = diag([1 1 1], 3);
            A_cont(4, 2) = -DoubleInvertedPendulum.graviationalAcceleration * massPend / this.massCart;
            A_cont(4, 4) = -this.frictionCart / this.massCart;
            A_cont(4, 5) =  this.frictionPendulum1/ (this.massCart * this.lengthPendulum1);
            A_cont(5, 2) = DoubleInvertedPendulum.graviationalAcceleration * massPend * (this.massCart + this.massPendulum1) / (this.massCart * this.massPendulum1 * this.lengthPendulum1);
            A_cont(5, 3) = -DoubleInvertedPendulum.graviationalAcceleration * this.massPendulum2 / (this.massPendulum1 * this.lengthPendulum1);
            A_cont(5, 4) = this.frictionCart / (this.massCart * this.lengthPendulum1);
            A_cont(5, 5) = -this.frictionPendulum1 * (this.massCart + this.massPendulum1) / (this.massCart * this.massPendulum1 * this.lengthPendulum1 ^ 2);
            A_cont(5, 6) = this.frictionPendulum2 / (this.massPendulum1 * this.lengthPendulum1 * this.lengthPendulum2);
            A_cont(6, 2) = -DoubleInvertedPendulum.graviationalAcceleration * massPend / (this.massPendulum1 * this.lengthPendulum2);
            A_cont(6, 3) = DoubleInvertedPendulum.graviationalAcceleration * massPend / (this.massPendulum1 * this.lengthPendulum2);
            A_cont(6, 5) = this.frictionPendulum1 / (this.massPendulum1 * this.lengthPendulum1 * this.lengthPendulum2);
            A_cont(6, 6) = -this.frictionPendulum2 * massPend / (this.massPendulum1 * this.massPendulum2 * this.lengthPendulum2 ^ 2);
            
            B_cont = [0; 0; 0; 1/this.massCart;  -1/(this.massCart * this.lengthPendulum1); 0];            
      
            C_cont = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 1 -1 0 0 0];
            
            % we must consider that first noise variable w1 acts on lower pendulum
            % w2 acts on upper pendulum and w3 acts on cart
            % so last column of G must equal B
            G_cont = zeros(6, 3);
            G_cont(4, 3) = 1/this.massCart;
            G_cont(4, 2) = -1/(this.massCart * this.lengthPendulum1);
            G_cont(5, 3) = G_cont(4,2);
            G_cont(5, 2) = (this.massCart + this.massPendulum1)  / (this.massCart * this.massPendulum1 * this.lengthPendulum1  ^ 2);
            G_cont(5, 1) = - 1/(this.massPendulum1 * this.lengthPendulum1 * this.lengthPendulum2);
            G_cont(6, 2) = G_cont(5, 1);
            G_cont(6, 1) = massPend / (this.massPendulum1 * this.massPendulum2 * this.lengthPendulum2  ^ 2);

            W_cont = diag([linContDisturbanceForcePend1; linContDisturbanceForcePend2; linContDisturbanceForceAct]); % semidefinite if one variance is 0
        end
        
        %% linearizeAroundUpwardEquilibrium
        function [A, B, C, W] = linearizeAroundUpwardEquilibrium(this, samplingInterval, contDisturbanceForcePend1, contDisturbanceForcePend2, contDisturbanceForceAct)
            % Obtain a discrete-time linearization of the dynamics around the unstable upward equilibrium 
            % which corresponds to the state [0 0 0 0 0 0]'.
            % The linearization is a difference equation of the form x_k = A * x_k-1 + B * u_k-1 + w_k-1 
            % and y_k = C * x_k is the output.
            % and y = C_cont * x is the output.
            % The state x is given by [q, q_dot]' 
            % with q = [position of the cart; deviation of lower pendulum from upward equilibrium; deviation of upper pendulum from upward equilibrium].
            % The input u_k is the force (in newtons) applied to move the cart.
            % The noise w_k is the discrete-time equivalent of the continuous-time noise G_cont * w and has covariance
            % W = G_cont * W_cont * G_cont' with w, G_cont and W_cont as described and computed by linearizeAroundUpwardEquilibriumCont().
            % Note: If either of the noise components is left out (variance = 0), this matrix is singular.
            %
            % Parameters:
            %   >> samplingInterval (Positive scalar, optional)
            %      The sampling interval (in seconds) to be used for the
            %      discretization of the linearized continuous-time
            %      dynamics. If none is provided, the value given by this
            %      instance's <samplingInterval> property is used.
            %
            % Parameters:
            %   >> linContDisturbanceForcePend1 (Nonnegative scalar, optional)
            %      The variance of the disturbance (torque), in newton-meters
            %      squared, entering at the tip of the lower pendulum rod.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForcePendulum1ContLin> property is used.
            %      If no disturbance torque is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            %   >> linContDisturbanceForcePend2 (Nonnegative scalar, optional)
            %      The variance of the disturbance (torque), in newton-meters
            %      squared, entering at the tip of the upper pendulum rod.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForcePendulum2ContLin> property is used.
            %      If no disturbance torque is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            %   >> linContDisturbanceForceAct (Nonnegative scalar, optional)
            %      The variance of the disturbance force (in newtons squared)
            %      additionally entering at the actuator to drive the cart in the continuous-time linearization.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForceActuatorContLin> property is used.
            %      If no disturbance force is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            % Returns:
            %   << A (6-by-6 matrix)
            %      The system matrix of the discrete-time linearization.
            %
            %   << B (6-dimensional column vector)
            %      The input matrix of the discrete-time linearization.
            %
            %   << C (3-by-6 matrix)
            %      The measurement matrix, assuming that x_1 (position of
            %      the cart) and x_3 (deviation of pendulum angle from
            %      upward equilibrium) are measured.
            %
            %   << W (6-by-6 matrix)
            %      The covariance matrix of the noise, as described above.
            
            arguments
                this                
                samplingInterval(1,1) double {mustBePositive} = this.samplingInterval;
                contDisturbanceForcePend1 (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForcePendulum1ContLin;
                contDisturbanceForcePend2 (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForcePendulum2ContLin;
                contDisturbanceForceAct (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceActuatorContLin;                
            end
            
            % compute based on linearized continuous-time dynamics
            [A_cont, B_cont, C_cont, G_cont, W_cont] = this.linearizeAroundUpwardEquilibriumCont(contDisturbanceForcePend1, ...
                contDisturbanceForcePend2, contDisturbanceForceAct);
            
            % discretize A and B in one shot, assuming zero order hold for the
            % input (cf. Raymond DeCarlo, Linear Systems: A State Variable
            % Approach with Numerical Implementation, page 215)
            tmp = expm([A_cont B_cont; zeros(1,7)] * samplingInterval);
            A = tmp(1:6, 1:6);
            B = tmp(1:6, 7:end);
            C = C_cont;        
            % assume that noise is zero mean, white      
            % use the following trick due to van Loan to discretize the noise            
            F = [-A_cont G_cont * W_cont * G_cont'; zeros(6) A_cont'] * samplingInterval;
            G = expm(F);
            W = A * G(1:6, 7:12); 
            % symmetrize W
            W = (W + W') / 2;
        end
        
        %% isValidState
        function isValid = isValidState(~, state)
            % Function to check whether a given plant state is valid in the
            % sense that the pendulum rod cannot be considered fallen over and the displacement of the cart is small.
            %
            % Parameters:
            %   >> state (6-dimensional vector)
            %      The system state to check.
            %
            % Returns:
            %   << isValid (Flag, i.e., boolean)
            %      Flag to indicate whether the given state state is valid.
            %      True is returned in case the absolute deviation of the rod from the upward
            %      equilibrium is within [0°, 45°] and the position of the cart is within [-5m, 5m], false otherwise.
            %      False is also returned if any of the state variables is +/-inf or +/-nan.
            
            arguments
                ~
                state(6,1) double {mustBeReal} % if row vector is passed, it is converted automatically
            end            
            isValid = all(isfinite(state)) && abs(rad2deg(state(2))) <= 55 ...
                && abs(rad2deg(state(3))) <= 90 ...
                && abs(state(1)) <= 5;
        end
    end
    
    methods (Access = protected)
        %% nonlinearDynamics
        function predictedStates = nonlinearDynamics(this, stateSamples, inputSamples, noiseSamples)
            % dynamics is of the form M(q)q_dotdot = f(q,q_dot, u, w) (1)
            % with q the generalized coordinates cart position (in m),
            % angular deviation from upright for both pendulums (in rad)
            % q_dot and q_dotdot are the corresponding time derivatives
            % (i.e., generalized velocities and accelerations)
            % u is the input force (in N) applied to move the cart
            % w are the disturbance forces acting on the cart and both pendulum rods
            % M is the nonsingular mass matrix for all possible configurations q
            % f contains coriolis and gravity terms
            %
            % (1) can be rewritten in ODE form:
            % [I 0; 0 M]x_dot = [q_dot; f(x, u, w)] (2)
            % by introducing the state variable x = [q, q_dot]', hence x_dot = [q_dot, q_dotdot]
            % assuming that u and w are constant over the sampling interval, we can solve this ODE numerically
            massPend = this.massPendulum1 + this.massPendulum2;            
            
            tspan = [0 this.samplingInterval];            
            
            if isempty(inputSamples)
                inputSamples = zeros(1, size(stateSamples, 2));
            end
            
            if isempty(noiseSamples)
                noiseSamples = zeros(3, size(stateSamples, 2));
            end
            
            predictedStates = zeros(6, size(stateSamples, 2));
            for j=1:size(stateSamples, 2)                
                [~,x]=ode113(@(t, x) fun(t, x, inputSamples(:, j), noiseSamples(:, j)), tspan, stateSamples(:, j), this.odeOpts);
                predictedStates(:, j) = x(end, :)';
            end
        
            function dxdt = fun(~, x, u, w)
                % states in x are ordered as follows: xc, theta1, theta2, then
%               % the corresponding velocities
                sinDiff = sin(x(2) - x(3));
                dxdt = zeros(6, 1);
                dxdt(1:3, :) = x(4:6); 
                % Coriolis terms
                dxdt(4) = -this.frictionCart * x(4) + this.lengthPendulum1 * massPend * sin(x(2)) * x(5)^2 + this.massPendulum2 * this.lengthPendulum2 * sin(x(3)) * x(6)^2;
                dxdt(5) = -this.frictionPendulum1 * x(5) - this.lengthPendulum1 * this.lengthPendulum2 * this.massPendulum2 * sinDiff * x(6)^2;
                dxdt(6) = -this.frictionPendulum2 * x(6) + this.lengthPendulum1 * this.lengthPendulum2 * this.massPendulum2 * sinDiff * x(5)^2;
                
                % gravity terms
                dxdt(5) = dxdt(5) + DoubleInvertedPendulum.graviationalAcceleration * this.lengthPendulum1 * massPend * sin(x(2));
                dxdt(6) = dxdt(6) + DoubleInvertedPendulum.graviationalAcceleration * this.lengthPendulum2 * this.massPendulum2 * sin(x(3));
                
                % input and noise                
                dxdt(4) = dxdt(4) + u + w(3);
                dxdt(5:6) = dxdt(5:6) + w(1:2);
            end
        end
    end
    
    methods (Access = private)
        
        %% getInverseMassMatrix
        function invM = getInverseMassMatrix(this, y)
            massPend = this.massPendulum1 + this.massPendulum2;
            numStateSamples = size(y, 2);
            
            detM = (this.lengthPendulum1 * this.lengthPendulum2) ^ 2 * this.massPendulum2 * ...
                (this.massCart * this.massPendulum1 + this.massPendulum1^2 .* (sin(y(2, :)) .^ 2) ...
                + this.massCart * this.massPendulum2 .* (sin(y(2, :)-y(3, :)) .^ 2) + this.massPendulum1 * this.massPendulum2 .* (sin(y(2, :)) .^ 2));
            
            % entries of the symmetric mass matrix M (state-dependent)
            M11 = this.massCart + this.massPendulum1 + this.massPendulum2;
            M12 = this.lengthPendulum1 * massPend .* cos(y(2, :));
            M13 = this.massPendulum2 * this.lengthPendulum2 .* cos(y(3, :));
            M22 = this.lengthPendulum1^2 * massPend;
            M23 = this.lengthPendulum1 * this.lengthPendulum2 * this.massPendulum2 .* cos(y(2, :)-y(3, :));
            M33 = this.lengthPendulum2^2 * this.massPendulum2;
            
            invM = zeros(3,3, numStateSamples);
            invM(1,1, :) = (M22 * M33-M23 .^ 2) ./ detM;
            invM(1,2, :) = (M13 .* M23-M12 .* M33) ./ detM;
            invM(1,3, :) = (M12 .* M23 - M13 .* M22) / detM;
            invM(2,1, :) = invM(1,2, :);
            invM(2,2, :) = (M11 * M33 - M13 .^ 2) ./ detM;
            invM(2,3, :) = (M13 .* M12 - M11 .* M23) ./ detM;
            invM(3,1, :) = invM(1,3, :);
            invM(3,2, :) = invM(2,3, :);
            invM(3,3, :) = (M11 * M22 - M12 .^ 2) ./ detM;
        end
        
        %% createMassMatrix
        function M = createMassMatrix(this, y)
            massPend = this.massPendulum1 + this.massPendulum2;
            % states in y are ordered as follows: q=[xc, theta1, theta2], then
            % the corresponding velocities q_dot
            M = zeros(3,3); % M affects only the verlocities q_dot
            M(1,1) = this.massCart + this.massPendulum1 + this.massPendulum2;
            M(1,2) = this.lengthPendulum1 * massPend * cos(y(2));
            M(1,3) = this.massPendulum2 * this.lengthPendulum2 * cos(y(3));
            M(2,1) = M(1,2);
            M(2,2) = this.lengthPendulum1^2 * massPend;
            M(2,3) = this.lengthPendulum1 * this.lengthPendulum2 * this.massPendulum2 * cos(y(2)-y(3));
            M(3,1) = M(1,3);
            M(3,2) = M(2,3);
            M(3,3) = this.lengthPendulum2^2 * this.massPendulum2;
        end
    
    end
end

