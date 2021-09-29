classdef (Abstract) ModeTransitionMatrixChangeable < handle
    % Superclass enabling filers to change the transition matrix
    % of the Markov chain theta_k to that governs the switching between the modes
    % of the augmented dynamical system.
        
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
    
    methods (Abstract, Access = public)
        % Set the mode transition matrix used to be used for the prediction from k to k+1,
        % where k is the current time step, and for all future predictions until
        % this method is called again.
        %
        % Parameters:
        %   >> modeTransitionMatrix (Matrix)
        %      A stochastic matrix whose (i,j)-th entry defines the probability of switching from mode i to j.
        %      In particular, the matrix must be square with all
        %      elements in [0,1] and rows summing up to 1.
        setModeTransitionMatrix(this, newModeTransitionMatrix)
    end
end

