classdef Updater < handle & matlab.mixin.Copyable
    % Structure class for every network updater module to change the 
    % probability transition matrix (TPM) in the DelayedModeIMMF.
    %
    % Updater provides functions to update delay probabilities and the TPM
    % as well as compatibily functions to ensure that a user or developer
    % can freely choose an updater without changing the filter code.
    
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
        function [newTransitionProbabilities, newTransitionMatrix] ...
                = updateTransitionProbabilities(...
                this, previousTransitionProbabilities)
            % Updates new probabilities by old probabilities.
            %
            % Will be moved to private access possibly.
            %
            % Parameters:
            %   >> previousTransitionProbabilities: Distribution of the 
            % transition probabilities from previous time step
            
            newTransitionProbabilities = previousTransitionProbabilities;            
            newTransitionMatrix = Utility.calculateDelayTransitionMatrix(newTransitionProbabilities);
        end
        
        function newTransitionMatrix = updateTransitionProbabilityMatrix(...
                this, measurement, previousTransitionMatrix)
            % Updates transition probability by the previous transition
            % matrix. In result, the transition matrix stays the same.
            %
            % Parameters:
            %   >> previousTransitionMatrix: Previous TPM
            %   >> measurements: Only for compatibility reasons
            newTransitionMatrix = previousTransitionMatrix;
        end
        
        function incrementk(this)
            % Exists only for compatibility and does nothing in Updater.
            % Function to decrease distance between a threshould and the 
            % current time step. This ensures that the convex optimization
            % update gathers enough data and calculates the first update
            % reasonably.
        end
        
        function initializeUpdateStep(this, stateMean, stateCovariance, ...
                modeProbabilities, inputs)
            % Exists only for compatibility and does nothing in Updater.
            % Updates attributes which are mandatory for calculating the
            % update of the TPM.
        end
    end
    
     %% Static functions
    methods (Static)
        %% Bound probabilities
        function boundedProbabilities = boundProbabilities(probabilities)
            % Bounds values of the probabilities to achieve numerical
            % stability.
            %
            % Parameters:
            %   >> probabilities: one-dimensional probability vector  
            idx = find(probabilities <= 1e-12);
            normalizedModeProbs = probabilities;
            if ~isempty(idx)
                normalizedModeProbs(idx) = 1e-12;
                normalizedModeProbs = normalizedModeProbs / sum(normalizedModeProbs);
            end            
            boundedProbabilities = normalizedModeProbs;
        end
    end
end