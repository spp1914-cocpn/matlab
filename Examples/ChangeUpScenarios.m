classdef ChangeUpScenarios < SetUpScenarios
    % Class for updating controllers, estimators, networks, and plant in 
    % case of changing properties of the network.
    % In this class, there are default setups for each component which can
    % be used in example code.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
    %                        Joanna Mueller <joanna.mueller@student.kit.edu>
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
    methods(Access = public)
        function this = ChangeUpScenarios(simTime, delayNumber, estimatorType)
            % Constructor of class which is used to change some attributes
            % of the components of a networked control system. This is used
            % to simulate changes of the network.
            %
            % parameters:
            %   >> simTime: duration of the simulation in number of time 
            %   steps 
            %   >> delayNumber: determines the new delay  distribution
            %   >> estimatorType: type determining if and how the
            %   transition probability matrix is updated
            % Returns:
            %   << this: object of type SetUpScenarios
            this = this@SetUpScenarios(simTime, estimatorType);
            
            % In the SimulatedNetworks.mat is a set of different example
            % distributions which are used in this example setup.
            NetworkDelayDistribution = load('SimulatedNetworks.mat');
            NetworkDelayDistribution = NetworkDelayDistribution.NetworkDelayDistribution;
            NetworkDelayDistribution = NetworkDelayDistribution(delayNumber,:);

            this.caDelayProb = NetworkDelayDistribution;
            this.scDelayProb = NetworkDelayDistribution;
            this.acDelayProb = NetworkDelayDistribution;
            
            this.caPacketDelays = DiscreteSample(this.caDelayProb,this.simTime)-1;
            this.scPacketDelays = DiscreteSample(this.scDelayProb,this.simTime)-1;
            this.acPacketDelays = DiscreteSample(this.acDelayProb,this.simTime)-1;
            
            this.caDelayProb = Utility.truncateDiscreteProbabilityDistribution(this.caDelayProb, ...
                    ( this.numModes ));

            this.transitionMatrix = Utility.calculateDelayTransitionMatrix(this.caDelayProb);
        end
    end
end