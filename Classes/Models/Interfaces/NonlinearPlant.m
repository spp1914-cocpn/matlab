classdef (Abstract) NonlinearPlant < SystemModel
    % Implementation of a general, nonlinear
    % plant evolving according to x_{k+1} = f(x_{k}, u_{k}, w_{k}), 
    % where f is an arbitrary system function, 
    % and w_{k} (optional) is arbitrarily distributed process
    % noise.
        
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
    
    properties (GetAccess = public, SetAccess = immutable)
        dimInput;
        dimState;
    end
    
    properties (Access = private)
        sysInput;
    end
    
    methods (Access = public)
        function this = NonlinearPlant(dimState, dimInput)
            % Class constructor.
            %
            % Parameters:
            %
            %   >> dimState (Positive integer)
            %      The dimension of the plant's state.
            %
            %   >> dimInput (Positive integer)
            %      The dimension of the inputs applied to the plant.
            %
            % Returns:
            %   << this (NonlinearPlant)
            %      A new NonlinearPlant instance.
            
            assert(Checks.isNonNegativeScalar(dimState) && mod(dimState, 1) == 0, ...
                'NonlinearPlant:InvalidStateDim', ...
                '** <dimState> must be a nonnegative integer **');
            assert(Checks.isNonNegativeScalar(dimInput) && mod(dimInput, 1) == 0, ...
                'NonlinearPlant:InvalidInputDim', ...
                '** <dimInput> must be a nonnegative integer **');
            
            this.dimState = dimState;
            this.dimInput = dimInput;
            this.sysInput = [];
            % by default, no noise is set
        end
       
        function setSystemInput(this, sysInput)
            % Set the system input vector.
            %
            % By default, the system input is an empty matrix.
            %
            % Parameters:
            %   >> sysInput (Column vector or empty matrix)
            %      The new system input vector.
            %      An empty matrix means no input vector.
            
            assert(isempty(sysInput) || (Checks.isColVec(sysInput, this.dimInput) && all(isfinite(sysInput))), ...
                'NonlinearPlant:SetSystemInput:InvalidSystemInput', ...
                '** <sysInput> must be a real-valued %d-dimensional column vector  **', this.dimInput);
            
            this.sysInput = sysInput;
        end
        
        function input = getSystemInput(this)
            % Get the system input vector.
            %
            % Returns:
            %   << input (Column vector or empty matrix)
            %      The currently applied input vector.
            %      An empty matrix means no input vector.
            
            input = this.sysInput;
        end
        
        function [stateJacobian, inputJacobian, noiseJacobian, ...
                  stateHessians, inputHessians, noiseHessians] = derivative(this, nominalState, nominalNoise)
            % Compute the first-order and second-order derivatives of the
            % implemented system equation by using first- and second-order
            % Taylor expansion around the provided nominal state and noise
            % values, and the current system input. If none is set, a
            % nominal input of zero is assumed.
            %
            % By default, the derivatives are computed using difference quotients.
            %
            % Parameters:
            %   >> nominalState (Column vector)
            %      The nominal system state vector.
            %
            %   >> nominalNoise (Column vector)
            %      The nominal system noise vector.
            %
            % Returns:
            %   << stateJacobian (Square matrix)
            %      The Jacobian of the state variables.
            %
            %   << inputJacobian (Square matrix)
            %      The Jacobian of the input variables.
            %
            %   << noiseJacobian (Matrix)
            %      The Jacobian of the noise variables.
            %
            %   << stateHessians (3D matrix)
            %      The Hessians of the state variables.
            %
            %   << inputHessians (3D matrix)
            %      The Hessians of the input variables.
            %
            %   << noiseHessians (3D matrix)
            %      The Hessians of the noise variables.
            
            % linearize input around the one currently set input, or zero
            if isempty(this.sysInput)
                nominalInput = zeros(this.dimInput, 1);
            else
                nominalInput = this.sysInput;
            end
            if nargout == 3
                [stateJacobian, inputJacobian, noiseJacobian] = Utils.diffQuotientStateInputAndNoise(@this.nonlinearDynamics, ...
                                                                                 nominalState, nominalInput, nominalNoise);
            else
                [stateJacobian, inputJacobian, noiseJacobian, ...
                 stateHessians, inputHessians, noiseHessians] = Utils.diffQuotientStateInputAndNoise(@this.nonlinearDynamics, ...
                                                                                 nominalState, nominalInput, nominalNoise);
            end
        end
        
        function predictedState = simulate(this, state)
            % Simulate the temporal evolution for a given system state.
            %
            % Parameters:
            %   >> state (Column vector)
            %      The system state to predict.
            %
            % Returns:
            %   << predictedState (Column vector)
            %      The simulated temporal system state evolution.
            
            assert(Checks.isColVec(state, this.dimState) && all(isfinite(state)), ...
                'NonlinearPlant:Simulate:InvalidSystemState', ...
                '** <state> must be a real-valued %d-dimensional column vector  **', this.dimState);
                        
            if isempty(this.noise)
                noiseSample = [];
            else
                noiseSample = this.noise.drawRndSamples(1);
            end
            
            predictedState = this.nonlinearDynamics(state, this.sysInput, noiseSample);
        end
        
        function predictedStates = systemEquation(this, stateSamples, inputSamples, noiseSamples)
            % The system equation.
            %
            % Parameters:
            %   >> stateSamples (Matrix)
            %      L column-wise arranged state samples.
            %
            %   >> inputSamples (Matrix)
            %      L column-wise arranged input samples.
            %
            %   >> noiseSamples (Matrix)
            %      L column-wise arranged system noise samples.
            %
            % Returns:
            %   << predictedStates (Matrix)
            %      L column-wise arranged predicted state samples.
            
            % simple use function to call true system function
            predictedStates = this.nonlinearDynamics(stateSamples, inputSamples, noiseSamples);
        end
    end
    
    methods (Access = public, Abstract)
        % Function to check whether a given plant state is valid (e.g.,
        % does not violate constraints).
        %
        % Parameters:
        %   >> state (Column vector)
        %      The system state to check.
        %
        % Returns:
        %   << isValid (Flag, i.e., boolean)
        %      Flag to indicate whether the given state is valid (e.g.,
        %      admissible with regards to contraints).
        
        isValid = isValidState(this, state);
    end
    
    methods (Access = protected, Abstract)
        % The nonlinear plant dynamics, i.e., the system equation.
        %
        % Parameters:
        %   >> stateSamples (Matrix)
        %      L column-wise arranged state samples.
        %
        %   >> inputSamples (Matrix)
        %      L column-wise arranged input samples.
        %
        %   >> noiseSamples (Matrix)
        %      L column-wise arranged system noise samples.
        %
        % Returns:
        %   << predictedStates (Matrix)
        %      L column-wise arranged predicted state samples.
        predictedStates = nonlinearDynamics(this, stateSamples, inputSamples, noiseSamples);
    end
end

