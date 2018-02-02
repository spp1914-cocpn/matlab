classdef (Abstract) SequenceBasedController < handle
    % Abstract base class for sequence-based controllers for (linear)
    % networked control systems.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (SetAccess = immutable, GetAccess = protected)
        dimPlantState;
        dimPlantInput;
    end
    
    properties (Access = protected)
        debugging = [];
    end
    
    properties (SetAccess = protected, GetAccess = public)
        % setter is protected, so that subclass can change this property
        % after creation, if supported
        sequenceLength;
    end
    
    methods 
        function set.sequenceLength(this, seqLength)
            Validator.validateSequenceLength(seqLength);
            this.sequenceLength = seqLength;
        end
    end
    
    methods (Access = protected)
        function this = SequenceBasedController(dimPlantState, dimPlantInput, sequenceLength)
            % Class constructor.
            %
            % Parameters:
            %   >> dimPlantState (Positive integer)
            %      The dimension of the plant's state.
            %
            %   >> dimPlantInput (Positive integer)
            %      The dimension of the inputs applied to the plant.
            %
            %  >> sequenceLength (Positive integer)
            %     The length of the input sequence (i.e., the number of
            %     control inputs) to be computed by the controller.
            %
            % Returns:
            %   << obj (SequenceBasedController)
            %      A new SequenceBasedController instance.
            
            this.dimPlantState = dimPlantState;
            this.dimPlantInput = dimPlantInput;
            this.sequenceLength = sequenceLength;
        end
    end
    
    methods (Access = protected, Abstract)
        inputSequence = doControlSequenceComputation(this, plantState, varargin);
        costs = doCostsComputation(this, stateTrajectory, appliedInputs);
    end
    
    methods (Access = public)
        %% computeCosts
        function costs = computeCosts(this, stateTrajectory, appliedInputs)
            % Compute accrued costs for the given state and input
            % trajectory according to this controller's underlying cost functional.
            %
            % Parameters:
            %   >> stateTrajectory (Matrix of dimension dimPlantState-by-n)
            %      A matrix representing a state trajectory, i.e., adjacent
            %      columns contain succesive plant states.
            %
            %   >> appliedInputs (Matrix of dimension dimPlantInput-by-m)
            %      A matrix representing an input trajectory, i.e., adjacent
            %      columns contain succesive control inputs.
            %
            % Returns:
            %   << costs (Nonnegative scalar)
            %      The accrued costs according to this controller's underlying cost functional.
            %
            if ~Checks.isFixedRowMat(stateTrajectory, this.dimPlantState) ...
                    || ~Checks.isFixedRowMat(appliedInputs, this.dimPlantInput)
                error('SequenceBasedController:ComputeCosts', ...
                    ['** <stateTrajectory> must be a real matrix with %d rows ' ...
                    'and <appliedInputs> must be a real matrix with %d rows'], this.dimPlantState, this.dimPlantInput);
            end
            costs = this.doCostsComputation(stateTrajectory, appliedInputs);
        end
        
        %% computeControlSequence
        function inputSequence = computeControlSequence(this, plantState, varargin)
            % Computes control input sequence based on the given plant state
            % (true state or estimate).
            %
            % Parameters:
            %   >> plantState (Subclass of Distribution)
            %      The current state of the plant, expressed as an
            %      estimate in terms of a distribution.
            %      In case of state-feedback with exactly known state, a
            %      DiracMixture with just on component should be employed.
            %
            %   >> varargin (Optional arguments)
            %      Any optional arguments, such as current time step (e.g. if the controller's horizon is not infinite) or 
            %      the current (estimated) mode of the controllor-actuator-plant subsystem.
            %
            % Returns:
            %   << inputSequence (Column Vector)
            %      A vector containing the stacked control inputs.
            %
            if ~Checks.isClass(plantState, 'Distribution') || plantState.getDimension() ~= this.dimPlantState
               error('SequenceBasedController:ComputeControlSequence', ...
                   ['** Input parameter <state> (state estimate) must be ' ...
                   'a %d-dimensional Distribution **'], this.dimPlantState);
            end
            inputSequence = this.doControlSequenceComputation(plantState, varargin{:});
        end
    end
    
    methods (Access = public, Abstract)
        % Resets the controller to its initial state, in case the
        % controller is a dynamical system itself.
        reset(this);
    end
    
end

