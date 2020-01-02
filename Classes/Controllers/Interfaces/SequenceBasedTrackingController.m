classdef (Abstract) SequenceBasedTrackingController < SequenceBasedController
    % Abstract base class for sequence-based tracking controllers for (linear)
    % networked control systems.
    
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
    
    properties (SetAccess = immutable, GetAccess = protected)
        Z;
        dimRef;
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        refTrajectory;
    end
    
    methods (Access = protected)
        %% SequenceBasedTrackingController
        function this = SequenceBasedTrackingController(dimPlantState, dimPlantInput, sequenceLength, ...
                needsStateEstimates, Z, refTrajectory, expectedTrajectoryLength)
            % Class constructor.
            %
            % Parameters:
            %   >> dimPlantState (Positive integer)
            %      The dimension of the plant's state.
            %
            %   >> dimPlantInput (Positive integer)
            %      The dimension of the inputs applied to the plant.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %  >> needsStateEstimates (Logical scalar, i.e, a flag)
            %     Flag to indicate whether the tracking controller needs an external
            %     filter to supply it with state estimates.
            %
            %   >> Z (Matrix, n-by-dimPlantState)
            %      The time-invariant plant output (performance) matrix, 
            %      i.e., z_k = Z*x_k
            %
            %   >> refTrajectory (Matrix, n-by-horizonLength+1)
            %      The reference trajectory to track, given as a matrix
            %      with the reference plant outputs column-wise arranged.
            %
            %   >> expectedTrajectoryLength (Positive integer)
            %      The expected length reference trajectory.
            %      Pass the empty matrix here, if no checks shall be done.
            %
            % Returns:
            %   << this (SequenceBasedTrackingController)
            %      A new SequenceBasedTrackingController instance.
            
            this = this@SequenceBasedController(dimPlantState, dimPlantInput, sequenceLength, needsStateEstimates);
            if ~isempty(Z)
                assert(Checks.isFixedColMat(Z, dimPlantState) && all(isfinite(Z(:))), ...
                    'SequenceBasedTrackingController:InvalidZMatrix', ...
                    '** Input parameter <Z> (Plant output/performance maxtrix) must be a real-valued matrix with %d cols **', ...
                    dimPlantState); 
                                
                this.Z = Z;
                this.dimRef = size(Z, 1);
                if isempty(expectedTrajectoryLength)
                    this.validateReferenceTrajectory(refTrajectory);
                else
                    this.validateReferenceTrajectory(refTrajectory, expectedTrajectoryLength);
                end
                this.refTrajectory = refTrajectory;
            else
                this.Z = [];
                this.dimRef = dimPlantState;
                this.refTrajectory = [];
            end
        end
    end
    
    methods (Access = public)
        %% getDeviationFromRefForState
        function deviation = getDeviationFromRefForState(this, trueState, timestep)
            % Compute deviation of the performance output for a given true
            % state from the reference at a given time step.
            %
            % Parameters:
            %   >> trueState (Vector of dimension dimPlantState)
            %      The plant true state at the given time step
            %
            %   >> timestep (Positive integer)
            %      The time step for which to retrieve the deviation from the
            %      reference trajectory.
            %
            % Returns:
            %   << deviation (Vector)
            %      The deviation of the performance output of the given
            %      true state and the reference, i.e., Z * x_true -z_ref.
            %
            assert(Checks.isVec(trueState, this.dimPlantState), ...
                'SequenceBasedTrackingController:GetDeviationFromRefForState:InvalidTrueState', ...
                '** Input parameter <trueState>  must be a %d-dimensional vector **', this.dimPlantState);
            % pass state as a column vector
            deviation = this.doGetDeviationFromRefForState(trueState(:), timestep);
        end
    end
    
    methods (Abstract, Access = protected)
        deviation = doGetDeviationFromRefForState(this, state, timestep);
    end
    
    methods (Access = private)
        %% validateReferenceTrajectory
        function validateReferenceTrajectory(this, refTrajectory, expectedTrajectoryLength)
            if nargin == 3
                assert(Checks.isMat(refTrajectory, this.dimRef, expectedTrajectoryLength) && all(isfinite(refTrajectory(:))), ...
                    'SequenceBasedTrackingController:InvalidReferenceTrajectory', ...
                    '** Reference trajectory <referenceTrajectory> must be a real-valued %d-by-%d matrix **',...
                    this.dimRef, expectedTrajectoryLength);                    
            else
                assert(Checks.isFixedRowMat(refTrajectory, this.dimRef) && all(isfinite(refTrajectory(:))), ...
                    'SequenceBasedTrackingController:InvalidReferenceTrajectory', ...
                    '** Reference trajectory <referenceTrajectory> must be a real-valued matrix with %d rows**',...
                    this.dimRef);
            end           
        end
    end
end

