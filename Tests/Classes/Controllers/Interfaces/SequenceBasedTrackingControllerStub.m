classdef SequenceBasedTrackingControllerStub < SequenceBasedTrackingController
    % This class represents a dummy sequence-based tracking controller which is used for testing
    % purposes.
    
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
    
    properties
    end
    
    methods (Access = public)
         %% SequenceBasedTrackingControllerStub
         function this = SequenceBasedTrackingControllerStub(dimX, dimU, sequenceLength, ...
                 needsFilter, Z, refTrajectory, expectedTrajectoryLength)
             this@SequenceBasedTrackingController(dimX, dimU, sequenceLength, needsFilter, ...
                 Z, refTrajectory, expectedTrajectoryLength);
         end
         
         function reset(~)
         end
         
         %% getProperties
         function [Z, dimRef] = getProperties(this)
             Z = this.Z;
             dimRef = this.dimRef;
         end
    end
    
    methods (Access = protected)
          %% doControlSequenceComputation
         function inputSequence = doControlSequenceComputation(this, ~, ~)
             inputSequence = ones(this.dimPlantInput * this.sequenceLength, 1);
         end
         
         %% doStageCostsComputation
         function stageCosts = doStageCostsComputation(~, ~, ~, ~)
             stageCosts = 0;
         end
         
         %% doCostsComputation
         function costs = doCostsComputation(~, ~, ~)
             costs = 0;
         end
         
         %% doGetDeviationFromRefForState
         function deviation = doGetDeviationFromRefForState(this, state, timestep)
             deviation = 0;
         end
    end
end

