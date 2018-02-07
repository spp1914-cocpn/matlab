classdef (Abstract) SequenceBasedTrackingController < SequenceBasedController
    % Abstract base class for sequence-based tracking controllers for (linear)
    % networked control systems.
    
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
    
    properties (SetAccess = immutable, GetAccess = protected)
        Z;
        dimRef;
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        refTrajectory;
    end
    
    methods (Access = protected)
        function this = SequenceBasedTrackingController(dimPlantState, dimPlantInput, sequenceLength, ...
                Z, refTrajectory, expectedTrajectoryLength)
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
            %   >> Z (Matrix, n-by-dimPlantState)
            %      The time-invariant plant output (performance) matrix, 
            %      i.e., z_k = Z*x_k
            %
            %   >> refTrajectory (Matrix, n-by-horizonLength+1)
            %      The reference trajectory to track, given as a matrix
            %      with the reference plant outputs column-wise arranged.
            %
            %  >> expectedTrajectoryLength (Positive integer)
            %     The expected length reference trajectory.
            %     Pass the empty matrix here, if no checks shall be done.
            %
            % Returns:
            %   << this (SequenceBasedTrackingController)
            %      A new SequenceBasedTrackingController instance.
            
            this = this@SequenceBasedController(dimPlantState, dimPlantInput, sequenceLength);
            if ~isempty(Z) 
                if ~Checks.isFixedColMat(Z, dimPlantState) || any(~isfinite(Z(:)))
                    error('SequenceBasedTrackingController:InvalidZMatrix', ...
                        '** Input parameter <Z> (Plant output/performance maxtrix) must be a real-valued matrix with %d cols **', ...
                        dimPlantState); 
                end
                                
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
            %   >> trueState (Vector dimension dimPlantState)
            %      The plant true state at the given time step
            %
            %   >> timestep (Positive integer)
            %      The time step for which to retrieve the deviation from the
            %      reference trajectory.
            %
            % Returns:
            %   << deviation (Nonnegative scalar)
            %      The deviation of the performance output of the given
            %      true state and the reference, i.e., Z * x_true -z_ref.
            %
            if ~Checks.isVec(trueState, this.dimPlantState)
                 error('SequenceBasedTrackingController:GetDeviationFromRefForState:InvalidTrueState', ...
                   ['** Input parameter <trueState>  must be ' ...
                   'a %d-dimensional vector **'], this.dimPlantState);
            end
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
                if ~Checks.isMat(refTrajectory, this.dimRef, expectedTrajectoryLength) ...
                        || any(~isfinite(refTrajectory(:)))
                    error('SequenceBasedTrackingController:InvalidReferenceTrajectory', ...
                        '** Reference trajectory <referenceTrajectory> must be a real-valued %d-by-%d matrix **',...
                        this.dimRef, expectedTrajectoryLength);
                end
            elseif ~Checks.isFixedRowMat(refTrajectory, this.dimRef) || any(~isfinite(refTrajectory(:)))
                error('SequenceBasedTrackingController:InvalidReferenceTrajectory', ...
                    '** Reference trajectory <referenceTrajectory> must be a real-valued matrix with %d rows**',...
                    this.dimRef);
            end
        end
    end
end

