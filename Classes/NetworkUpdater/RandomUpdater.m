classdef RandomUpdater < Updater
    % Implements an Updater which just returns a new probabilitz matrix
    % with 'random' values
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2019  Joanna Mueller <joanna.mueller@student.kit.edu>
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
     methods (Access = public)                
        function newTransitionMatrix = updateTransitionProbabilityMatrix(...
                    this, measurement, previousTransitionMatrix)
            % Updates transition probability matrix by a new randomly
            % distributed matrix. The random variables are drawn from
            % an uniform distribution.
            %
            % Parameters:
            %   >> measurements: two-dimensional vector of measurements
            %   >> previousTransitionMatrix: Previous transition matrix 
            % Returns:
            %   << newTransitionMatrix: two-dimentional transition
            %   probability matrix for the next time step    
            dim = size(previousTransitionMatrix, 1);
            
            caDelayProb = randn([1, dim]);
            caDelayProb = abs(caDelayProb);
            caDelayProb = caDelayProb / sum(caDelayProb);
            
            newTransitionMatrix = Utility.calculateDelayTransitionMatrix(caDelayProb);
        end
    end
end
