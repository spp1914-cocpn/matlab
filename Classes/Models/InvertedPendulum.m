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
    %    Copyright (C) 2018-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        massCart;
        massPendulum;
        lengthPendulum;
        inertia;
        friction;  
    end
    
    properties (Access = public)
        samplingInterval;
    end
    
    methods
        function set.samplingInterval(this, newSamplingInterval)
            assert(Checks.isPosScalar(newSamplingInterval), ...
                'InvertedPendulum:InvalidSamplingInterval', ...
                '** <samplingInterval> must be a positive scalar. **');
            this.samplingInterval = newSamplingInterval;
        end
    end
    
    methods (Access = public)
        %% InvertedPendulum
        function this = InvertedPendulum(massCart, massPendulum, lengthPendulum, friction, samplingInterval)
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
            %      A positive scalar denoting the length of the pendulum rod (in m).
            %
            %   >> friction (Positive scalar)
            %      A positive scalar denoting the coefficient of friction for the cart (in Ns/m).
            %
            %   >> samplingInterval (Positive scalar)
            %      The sampling interval (in seconds) that was used for the discretization of the
            %      underlying continuous-time plant model.
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
            this.samplingInterval = samplingInterval;
        end 
        
        %% linearizeAroundUpwardEquilibriumCont
        function [A_cont, B_cont, C_cont] = linearizeAroundUpwardEquilibriumCont(this)
            % Obtain a continuous-time linearization of the pendulum dynamics around the unstable upward equilibrium 
            % which corresponds to the state [0 0 pi 0]'.
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
            
            % compute based on linearized continuous-time dynamics
            denom = this.inertia * (this.massCart + this.massPendulum) + this.massPendulum * this.lengthPendulum ^ 2 * this.massCart;
            A_cont = [0 1 0 0;
                 0  -(this.inertia+this.massPendulum*this.lengthPendulum^2)*this.friction/denom ((this.massPendulum * this.lengthPendulum)^2*InvertedPendulum.graviationalAcceleration)/denom 0;
                 0 0 0 1;
                 0 -(this.massPendulum*this.lengthPendulum* this.friction)/denom (this.massPendulum*this.lengthPendulum*(this.massPendulum + this.massCart)*InvertedPendulum.graviationalAcceleration)/denom 0];
     
            B_cont = [0;  (this.inertia+this.massPendulum*this.lengthPendulum^2)/denom; 0; (this.massPendulum*this.lengthPendulum)/denom];
            
            C_cont = [1 0 0 0; 0 0 1 0];
        end
        
        %% linearizeAroundUpwardEquilibrium
        function [A, B, C, W] = linearizeAroundUpwardEquilibrium(this, W_cont)
            % Obtain a discrete-time linearization of the pendulum dynamics around the unstable upward equilibrium 
            % which corresponds to the state [0 0 pi 0]'. 
            % To that end, the continuous-time dynamics is first linearized
            % and then discretized.
            %
            % Parameters:
            %   >> W_cont (Positive semidefinite 4-by-4 matrix, optional)
            %      A 4-by-4 positive semidefinite matrix denoting the intensity (PSD)/covariance of the noise acting on continuous-time linearized dynamics.
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
            %      If W_cont was passed, the covariance of the
            %      corresponding (equivalent) discrete-time noise process is returned,
            %      otherwise the empty matrix.
            
            % compute based on linearized continuous-time dynamics
            [A_cont, B_cont, C_cont] = this.linearizeAroundUpwardEquilibriumCont();
            
            % discretize A and B in one shot, assuming zero order hold for the
            % input (cf. Raymond DeCarlo, Linear Systems: A State Variable
            % Approach with Numerical Implementation, page 215)
            tmp = expm([A_cont B_cont; zeros(1,5)] * this.samplingInterval);
            A = tmp(1:4, 1:4);
            B = tmp(1:4, 5:end);
            C = C_cont; 
            
            W = [];
            if nargin == 2 && ~isempty(W_cont)
                % assume that noise is zero mean, white
                % discretize the noise by integration
                % W = integral(@(x) expm(A_cont*x) * plantNoiseCov * expm(A_cont'*x), ...
                    %0, samplingInterval, 'ArrayValued', true);
                % use the following trick due to van Loan instead to discretize the noise
                F = [-A_cont W_cont; zeros(4) A_cont'] * this.samplingInterval;
                G = expm(F);
                W = A * G(1:4, 5:8);
            end
        end
        
        %% isValidState
        function isValid = isValidState(~, state)
            % Function to check whether a given plant state is valid in the
            % sense that the pendulum rod cannot be considered fallen over and the displacement of the cart is small.
            %
            % Parameters:
            %   >> state (4-dimensional column vector)
            %      The system state to check.
            %
            % Returns:
            %   << isValid (Flag, i.e., boolean)
            %      Flag to indicate whether the given state state is valid.
            %      True is returned in case the absolute deviation of the rod from the upward
            %      equilibrium is within [0°, 45°] and the position of the cart is within [-5m, 5m], false otherwise.
            
             assert(Checks.isColVec(state, 4) && all(isfinite(state)), ...
                'InvertedPendulum:IsValidState:InvalidSystemState', ...
                '** <state> must be a real-valued 4-dimensional column vector  **');
            
            isValid = Checks.isScalarIn(rad2deg(state(3)), 135, 225) ...
                && Checks.isScalarIn(state(1), -5, 5);
        end
    end
    
    methods (Access = protected)
        %% nonlinearDynamics
        function predictedStates = nonlinearDynamics(this, stateSamples, inputSamples, noiseSamples)
            % new cart position (x_1) and pendulum angle (x_3)
            predictedStates(1, :) = this.samplingInterval * stateSamples(2, :) + stateSamples(1, :);
            predictedStates(3, :) = this.samplingInterval * stateSamples(4, :) + stateSamples(3, :);

            if isempty(inputSamples)
                inputSamples = zeros(1, size(stateSamples, 2));
            end

            % compute new cart velocity (x_2) and pendulum angular velocity (x_4)
            sins = sin(stateSamples(3, :));
            coss = cos(stateSamples(3, :));
            part1 = this.inertia + this.massPendulum * this.lengthPendulum ^ 2;
            part2 = (this.lengthPendulum * this.massPendulum) ^ 2;
            part3 = this.massPendulum * this.lengthPendulum * sins .* (stateSamples(4, :) .^2) ...
                - this.friction * stateSamples(2, :) + inputSamples;

            enumerator = part1 * part3 + part2 * InvertedPendulum.graviationalAcceleration * sins .* coss;
            denominator = part1 * (this.massCart + this.massPendulum) - part2 * coss .^ 2;

            predictedStates(2, :) = this.samplingInterval * enumerator ./ denominator + stateSamples(2, :);

            enumerator = (this.massCart + this.massPendulum) * (-this.massPendulum * InvertedPendulum.graviationalAcceleration * this.lengthPendulum * sins) ...
                - this.massPendulum * this.lengthPendulum * coss .* part3;

            predictedStates(4, :) = this.samplingInterval * enumerator ./ denominator + stateSamples(4, :);
            
            if ~isempty(noiseSamples)
                predictedStates = predictedStates + noiseSamples;
            end
        end
    end
end

