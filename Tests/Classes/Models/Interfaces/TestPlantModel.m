classdef TestPlantModel < NonlinearPlant
    % This class represents a simple nonlinear plant which is used for testing
    % purposes.
    
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
    
    properties (Access = public, Constant)
        dimX = 3;
        dimU = 2;
    end
    
    methods (Access = public)
        function this = TestPlantModel(dimX, dimU)
            if nargin == 0
                dimX = TestPlantModel.dimX;
                dimU = TestPlantModel.dimU;
            end
            this@NonlinearPlant(dimX, dimU);
        end
        
        function states = simulateForInput(this, stateSamples, input)
            states = this.nonlinearDynamics(stateSamples, repmat(input, 1, size(stateSamples, 2)), []);
        end
        
        function isValid = isValidState(this, state)
            isValid = true;
        end
    end
    
    methods (Access = protected)
        function predictedStates = nonlinearDynamics(this, stateSamples, inputSamples, noiseSamples)
            % deterministic system equation
            if isempty(inputSamples)
                predictedStates = stateSamples .^ 2;
            else
                predictedStates = stateSamples .^ 2 + [inputSamples; inputSamples(1, :)];
            end
        end
    end
    
end

