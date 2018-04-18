classdef InvertedPendulum < NonlinearPlant
    % Implementation of an inverted pendulum on a cart with discrete-time 
    % nonlinear noise-free dynamics according to x_{k+1} = f(x_{k}, u_{k}), 
    % where x_{k} is the 4-dimensional pendulum state consisting of
    % cart position, cart velocity, angle of the pendulum (downwards, in radians) and
    % angular velocity of the pendulum.
    % Note that the state was chosen such that the unstable upward
    % equilibrium corresponds to an angle of pi.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Constant, Access = private)
        graviationalAcceleration = 9.81; % in N/kg
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        massCart;
        massPendulum;
        lengthPendulum;
        inertia;
        friction;
        
        samplingInterval;
    end
    
    methods (Access = public)
        %% InvertedPendulum
        function this = InvertedPendulum(massCart, massPendulum, lengthPendulum, inertia, friction, samplingInterval)
            this@NonlinearPlant(4, 1);
            
            this.massCart = massCart;
            this.massPendulum = massPendulum;
            this.lengthPendulum = lengthPendulum;
            this.inertia = inertia;
            this.friction = friction;
            this.samplingInterval = samplingInterval;
        end
    end
    
    methods (Access = protected)
        %% nonlinearDynamics
        function predictedStates = nonlinearDynamics(this, stateSamples, inputSamples, ~)
            predictedStates(1, :) = this.samplingInterval * stateSamples(2, :) + stateSamples(1, :);
            predictedStates(3, :) = this.samplingInterval * stateSamples(4, :) + stateSamples(3, :);

            if isempty(inputSamples)
                inputSamples = zeros(1, size(stateSamples, 2));
            end

            % compute new velocity, x_2
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
        end
    end
  
end

