classdef InvertedPendulum < NonlinearPlant
    % Implementation of an inverted pendulum on a cart with discrete-time 
    % nonlinear dynamics according to x_{k+1} = f(x_{k}, u_{k}) + w_{k}, 
    % where x_{k} is the 4-dimensional pendulum state consisting of
    % cart position, cart velocity, angle of the pendulum (downwards, in radians), and
    % angular velocity of the pendulum, and where w_{k} is an optional
    % process noise.
    % Note that the state was chosen such that the unstable upward
    % equilibrium corresponds to an angle of pi.
    % Moreover, the pendulum is assumed to be a uniform rod with center of
    % mass located in the middle.
    
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
        massPendulum;
        lengthPendulum; % effective length of pendulum
        inertia;
        friction;
        
        useMexImplementation(1,1) logical = true; 
        % by default, we use the C++ (mex) implementation for computation
        % of the new pendulum state
        % this is usually faster, but can produce slightly different results
    end
    
    properties (SetAccess = private, GetAccess = public)
        %
        % variance of disturbance force driving the pendulum rod (entering at the pendulum tip)
        % default zero
        varDisturbanceForcePendulum (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        % variance of disturbance force driving the cart (entering at the actuator, actuation noise)
        % default zero
        varDisturbanceForceActuator (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        %        
    end
    
    properties (Access = public)
        samplingInterval(1,1) double {mustBePositive, mustBeFinite} = InvertedPendulum.defaultInternalSamplingRate;
        
        % we need noise terms for the continuous-time linearization of the model
        % can be different from the ones used to simulate the
        % discrete-time nonlinear dynamics
        % the controllers use discrete-time linearizations, the noise in
        % these models is dependent on the sampling rate
        % the terms are assumed to be independent of each other
        %
        varDisturbanceForcePendulumContLin (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        varDisturbanceForceActuatorContLin (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
    end
    
    methods (Access = public)
        %% InvertedPendulum
        function this = InvertedPendulum(massCart, massPendulum, lengthPendulum, friction, samplingInterval, useMexImplementation)
            % Class constructor.
            %
            % Parameters:
            %   >> massCart (Positive Scalar)
            %      A positive scalar denoting the mass (in kg) of the cart.
            %
            %   >> massPendulum (Positive scalar)
            %      A positive scalar denoting the mass (in kg) of the pendulum.
            %
            %   >> lengthPendulum (Positive scalar)
            %      A positive scalar denoting the length to the center of mass of the pendulum rod (in m).
            %      The center of mass of the rod is assumed to be in the
            %      middle, so this value is the effective pendulum length
            %      and half the overall length.
            %
            %   >> friction (Positive scalar)
            %      A positive scalar denoting the coefficient of friction for the cart (in Ns/m).
            %
            %   >> samplingInterval (Positive scalar, optional)
            %      The sampling interval (in seconds) that was used for the discretization of the
            %      underlying continuous-time plant model.
            %      If left out, the default vale 0.001 is used, which
            %      corresponds to a sampling rate of 1 kHz.
            %
            %   >> useMexImplementation (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for computation of the pendulum dynamics 
            %      which is typically faster than the Matlab implementation. 
            %      If left out, the default value true is used.
            %
            % Returns:
            %   << this (InvertedPendulum)
            %      A new InvertedPendulum instance.
            
            this@NonlinearPlant(4, 1);
            
            this.massCart = massCart;
            this.massPendulum = massPendulum;
            this.lengthPendulum = lengthPendulum;
            this.inertia = massPendulum * lengthPendulum^2 / 3; % moment of inertia of the pendulum, assuming a uniform rod (center of mass in the middle)
            this.friction = friction;
            if nargin > 4
                this.samplingInterval = samplingInterval;
            end
            if nargin > 5
                this.useMexImplementation = useMexImplementation;
            end
        end 
        
        %% setNoise
        function setNoise(this, noise)
            % Set the additional noise driving the pendulum.
            % The noise driving the pendulum consists of two, typically
            % independent, components: A disturbance force driving the
            % pendulum rod, which enters at the pendulum tip, and an
            % actuation noise acting on the force applied to the cart.
            % Thus, the 2-dimensional joint distribution is passed here
            % with first component the marginal distribution of the force 
            % acting on the pendulum tip. 
            %
            % Parameters:
            %   >> noise (Subclass of Distribution)
            %      The 2-dimensional system noise.
            
            assert(Checks.isClass(noise, 'Distribution') && noise.getDim() == 2, ...
                'InvertedPendulum:SetNoise:InvalidNoise', ...
                '** <noise> bust be a 2-dimensional distribution **');                
            
            setNoise@SystemModel(this, noise);
                
            [~, noiseCov] = noise.getMeanAndCov();
            % store for convenience
            this.varDisturbanceForceActuator = noiseCov(2, 2);
            this.varDisturbanceForcePendulum = noiseCov(1, 1);
        end
        
        %% linearizeAroundUpwardEquilibriumCont
        function [A_cont, B_cont, C_cont, G_cont, W_cont] = linearizeAroundUpwardEquilibriumCont(this, linContDisturbanceForcePend, linContDisturbanceForceAct)
            % Obtain a continuous-time linearization of the pendulum dynamics around the unstable upward equilibrium 
            % which corresponds to the state [0 0 pi 0]'.
            % The linearization is a differential equation of the form x_dot = A_cont * x + B_cont * u + G_cont * w 
            % and y = C_cont * x is the output.
            % The state x is given by: 
            % [position of the cart; velocity of the cart; deviation of pendulum angle from upward equilibrium; angular velocity]
            % The input u is the force (in newtons) applied to move the cart.
            % The noise w is composed of a disturbance force acting at the
            % tip of the pendulum rod and a disturbance force acting on the
            % cart (additive on u, actuation noise). Both are assumed white
            % and mutually independent, the returned noise covariance W_cont is thus diagonal.
            %
            % Parameters:
            %   >> linContDisturbanceForcePend (Nonnegative scalar, optional)
            %      The variance of the disturbance force (in newtons
            %      squared) entering at the tip of the pendulum rod.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForcePendulumContLin> property is used.
            %      If no disturbance force is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            %   >> linContDisturbanceForceAct (Nonnegative scalar, optional)
            %      The variance of the disturbance force (in newtons
            %      squared) addtionally entering at the actuator to drive the cart.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForceActuatorContLin> property is used.
            %      If no disturbance force is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            % Returns:
            %   << A_cont (4-by-4 matrix)
            %      The system matrix of the continuous-time linearization.
            %
            %   << B_cont (4-dimensional column vector)
            %      The input matrix of the continuous-time linearization.
            %
            %   << C_cont (2-by-4 matrix)
            %      The measurement matrix, assuming that x_1 (position of
            %      the cart) and x_3 (deviation of pendulum angle from
            %      upward equilibrium) are measured.
            %
            %   << G_cont (4-by-2 matrix)
            %      The system noise matrix, as described above, of the continuous-time linearization.
            %
            %   << W_cont (2-by-2 matrix)
            %      The covariance of the white noise acting on the
            %      continuous-time linearized dynamics as described above.
            %      Note: If either of the noise components is left out (variance = 0), this matrix is singular.
            
            arguments
                this
                linContDisturbanceForcePend (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForcePendulumContLin;
                linContDisturbanceForceAct (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceActuatorContLin;
            end
            
            % compute based on linearized continuous-time dynamics
            denom = this.inertia * (this.massCart + this.massPendulum) + this.massPendulum * this.lengthPendulum ^ 2 * this.massCart;
            A_cont = [0 1 0 0;
                 0  -(this.inertia+this.massPendulum*this.lengthPendulum^2)*this.friction/denom ((this.massPendulum * this.lengthPendulum)^2*InvertedPendulum.graviationalAcceleration)/denom 0;
                 0 0 0 1;
                 0 -(this.massPendulum*this.lengthPendulum* this.friction)/denom (this.massPendulum*this.lengthPendulum*(this.massPendulum + this.massCart)*InvertedPendulum.graviationalAcceleration)/denom 0];
     
            B_cont = [0;  (this.inertia+this.massPendulum*this.lengthPendulum^2)/denom; 0; (this.massPendulum*this.lengthPendulum)/denom];
            
            C_cont = [1 0 0 0; 0 0 1 0];
                        
            G_cont = [[0; this.massPendulum*this.lengthPendulum^2/denom; 0; this.lengthPendulum*(this.massPendulum + this.massCart)/denom], B_cont];
            W_cont = diag([linContDisturbanceForcePend; linContDisturbanceForceAct]); % semidefinite if one variance is 0
        end
        
        %% linearizeAroundUpwardEquilibrium
        function [A, B, C, W] = linearizeAroundUpwardEquilibrium(this, samplingInterval, contDisturbanceForcePend, contDisturbanceForceAct)
            % Obtain a discrete-time linearization of the pendulum dynamics around the unstable upward equilibrium 
            % which corresponds to the state [0 0 pi 0]'.
            % The linearization is a difference equation of the form x_k = A * x_k-1 + B * u_k-1 + w_k-1 
            % and y_k = C * x_k is the output.
            % The state x_k is given by: 
            % [position of the cart; velocity of the cart; deviation of pendulum angle from upward equilibrium; angular velocity]
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
            %   >> linContDisturbanceForcePend (Nonnegative scalar, optional)
            %      The variance of the disturbance force (in newtons squared)
            %      entering at the tip of the pendulum rod in the continuous-time linearization.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForcePendulumContLin> property is used.
            %      If no disturbance force is to be used, either set
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
            %   << A (4-by-4 matrix)
            %      The system matrix of the discrete-time linearization.
            %
            %   << B (4-dimensional column vector)
            %      The input matrix of the discrete-time linearization.
            %
            %   << C (2-by-4 matrix)
            %      The measurement matrix, assuming that x_1 (position of
            %      the cart) and x_3 (deviation of pendulum angle from
            %      upward equilibrium) are measured.
            %
            %   << W (4-by-4 matrix)
            %      The covariance matrix of the noise, as described above.
            
            arguments
                this                
                samplingInterval(1,1) double {mustBePositive} = this.samplingInterval;
                contDisturbanceForcePend (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForcePendulumContLin;
                contDisturbanceForceAct (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceActuatorContLin;                
            end
            
            % compute based on linearized continuous-time dynamics
            [A_cont, B_cont, C_cont, G_cont, W_cont] = this.linearizeAroundUpwardEquilibriumCont(contDisturbanceForcePend, contDisturbanceForceAct);
            
            % discretize A and B in one shot, assuming zero order hold for the
            % input (cf. Raymond DeCarlo, Linear Systems: A State Variable
            % Approach with Numerical Implementation, page 215)
            tmp = expm([A_cont B_cont; zeros(1,5)] * samplingInterval);
            A = tmp(1:4, 1:4);
            B = tmp(1:4, 5:end);
            C = C_cont;           
            % assume that noise is zero mean, white
            % discretize the noise by integration
%             W = integral(@(x) expm(A_cont*x) * G_cont * W_cont * G_cont' * expm(A_cont'*x), ...
%                 0, samplingInterval, 'ArrayValued', true);
            % use the following trick due to van Loan instead to discretize the noise
            F = [-A_cont G_cont * W_cont * G_cont'; zeros(4) A_cont'] * samplingInterval;
            G = expm(F);
            W = A * G(1:4, 5:8); 
        end
        
        %% isValidState
        function isValid = isValidState(~, state)
            % Function to check whether a given plant state is valid in the
            % sense that the pendulum rod cannot be considered fallen over and the displacement of the cart is small.
            %
            % Parameters:
            %   >> state (4-dimensional vector)
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
                state(4,1) double {mustBeReal} % if row vector is passed, it is converted automatically
            end
            
            isValid = all(isfinite(state)) && ...
                135 <= rad2deg(state(3)) && rad2deg(state(3)) <= 225 && -5 <= state(1) && state(1) <= 5;
        end
    end
    
    methods (Access = protected)
        %% nonlinearDynamics
        function predictedStates = nonlinearDynamics(this, stateSamples, inputSamples, noiseSamples)
            if this.useMexImplementation            
                predictedStates = mex_PendulumNonlinearDynamics(this.samplingInterval, this.inertia, this.massPendulum, this.massCart,...
                    this.lengthPendulum, this.friction, stateSamples, inputSamples, noiseSamples); 
            else
                % new cart position (x_1) and pendulum angle (x_3)
                predictedStates(1, :) = this.samplingInterval * stateSamples(2, :) + stateSamples(1, :);
                predictedStates(3, :) = this.samplingInterval * stateSamples(4, :) + stateSamples(3, :);

                if isempty(inputSamples)
                    inputSamples = zeros(1, size(stateSamples, 2));
                end

                % compute new cart velocity (x_2) and pendulum angular velocity (x_4)
                sins = sin(stateSamples(3, :));
                coss = cos(stateSamples(3, :));

                massLength = this.massPendulum * this.lengthPendulum; % m*l
                massLengthCos = massLength .* coss; % m* l * cos(x_3)
                inertiaSum = this.inertia + this.lengthPendulum * massLength; % I+m*l^2
                sumMasses = this.massPendulum + this.massCart; % M + m
                denominator = sumMasses * inertiaSum - massLengthCos .^2; % (M+m)*(I+m*l^2)-m^2*l^2*cos(x_3)^2
                % u + m*l*sin(x_3)*(x_4)^2-b*x_2 + w_u
                % noisy input w_u here (actuation forces enters the same way as control input)
                forces = inputSamples + massLength * sins .* (stateSamples(4, :) .^ 2) - this.friction * stateSamples(2, :);
                % m * g * sin(x_3) - w_p
                % disturbance force w_p here (disturbance force acting on the pendulum tip)
                forces2 = this.massPendulum * InvertedPendulum.graviationalAcceleration * sins; 
                if ~isempty(noiseSamples)
                    % add the noise
                    forces = forces + noiseSamples(2, :);
                    forces2 = forces2 - noiseSamples(1, :); 
                end           

                predictedStates(2, :) = this.samplingInterval * (inertiaSum * forces + this.lengthPendulum * massLengthCos .* forces2) ...
                    ./ denominator + stateSamples(2, :);

                predictedStates(4, :) = this.samplingInterval * (sumMasses * this.lengthPendulum * (-forces2) - massLengthCos .* forces) ...
                    ./ denominator + stateSamples(4, :);
            end          
        end
    end
end

