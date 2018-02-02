classdef CommunicationNetwork < handle
    % Defines a simple communication network which allos to send and
    % receive data packets.
    %
    %
    % AUTHOR:       Jörg Fischer
    % LAST UPDATE:  Maxim Dolgov, 30.04.2013
    %               Florian Rosenthal, 20.02.2017   
    
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
      
    properties (SetAccess = private)
        % cell array with data packets to store;
        % (array of DataPacket)
        nwData = {};
        % number of packets that leave the network to current time step;
        % Index is the time step; (array of positive integer)
        numPacketsOut = [];
       
    end % properties
    
    properties (GetAccess = private, SetAccess = immutable)
        % used protocol (has no relevance to behaviour of the network. Only
        % needed for other components as actuator and controller to know if age
        % of buffer (mode) is available
        % (scalar of 1 (TCP) or 0 (UDP)
        usesTcp = 1;
        % length of simulation time in time steps; (positive integer)
        simTime = -1;
        % maximal time delay;
        % The network does not output a data packet that has a longer time 
        % delay; (positive integer)
        maxDelay = -1;
    end
  
    methods (Access = public)
        %% CommunicationNetwork
        function this = CommunicationNetwork(simTime, maxDelay, protocol)
            % Class constructor.
            %
            % Parameters:
            %   >> simTime (Positive integer)
            %      The total number of simulation time steps.
            %
            %   >> maxDelay (Nonnegative integer)
            %      The maximum delay (in time step) a single DataPacket may experience in the network
            %      before it gets lost.
            %
            %   >> protocol (Char array, either TCP or UDP, optional)
            %      Char array to indicate whther this network is TCP-like or UDP-like.
            %      If value is not provided here, TCP-like is used by
            %      default.
            %
            % Returns:
            %   << this (CommunicationNetwork)
            %      A new CommunicationNetwork instance.
           
            if nargin == 2
              this.usesTcp = 1;
            elseif ~ischar(protocol)
              error('CommunicationNetwork:InvalidProtocol', ...
                  '** Input parameter <protocol> must be string with value "TCP" or "UDP" **');
            elseif strcmpi(protocol, 'UDP')
              this.usesTcp = 0;
            elseif strcmpi(protocol, 'TCP')
              this.usesTcp = 1;
            else
              error('CommunicationNetwork:InvalidProtocolChar', ...
                  '** Input parameter <protocol> must be string with value "TCP" or "UDP" **');
            end
            
            if ~Checks.isPosScalar(simTime) || mod(simTime, 1) ~= 0
                error('CommunicationNetwork:InvalidSimTime', ...
                    ['** Input parameter <simTime> (simulation time) must be',...
                    ' a positive integer **']);
            end
            this.simTime = simTime;

            % maxDelay
            if ~Checks.isNonNegativeScalar(maxDelay)  || mod(maxDelay, 1) ~= 0
              error('CommunicationNetwork:InvalidMaxDelay', ...
                  ['** Input parameter <maxDelay> (maximum possible packet',...
                     ' delay) must be a nonnegative integer **'])
            end
            this.maxDelay = maxDelay;

            % init network data
            this.numPacketsOut = zeros(simTime, 1);
            this.nwData = {};
        end 
  
        %% sendPacket
        function sendPacket(this, packet, packetDelay, timeStep)
            if ~(Checks.isScalarIn(timeStep, 0, this.simTime) && mod(timeStep, 1) == 0)
                error ('CommunicationNetwork:SendPacket:InvalidTimeStep', ...
                    '** Time step must be an integer within [0, %d] **', this.simTime);    
            end
            
            if ~Checks.isClass(packet, 'DataPacket')
                 error('CommunicationNetwork:SendPacket:InvalidPacket', ...
                     '** Input parameter <packet> (sent data packet) is not of class DataPacket **');
            end
            % check packet and delay
            if ~(Checks.isNonNegativeScalar(packetDelay) && mod(packetDelay, 1) == 0)
                 error('CommunicationNetwork:SendPacket:InvalidPacketDelay', ...
                     '** Input parameter <packetDelay> (time delay to experience by the packet) must be a nonnegative integer **');
            end
            if packet.timeStamp ~= timeStep
                % error if packet time stamp does not correspond to transmission time    
                error('CommunicationNetwork:SendPacket:InvalidPacketTimeStamp', ...
                    '** Packet time stamp does not correspond to transmission time (time step) **');
            elseif packetDelay <= this.maxDelay && packetDelay <= this.simTime - timeStep
                % handle packet as its delay is not too large and will leav
                % the network before simulation ends
                packet.packetDelay = packetDelay;
                i = this.numPacketsOut(timeStep + packetDelay) + 1;
                this.numPacketsOut(timeStep + packetDelay) = i;
                this.nwData{timeStep + packetDelay, i} = packet;
            end
        end
        
        %% receivePackets
        function [numPackets, packets] = receivePackets(this, timeStep)
            if ~(Checks.isScalarIn(timeStep, 0, this.simTime) && mod(timeStep, 1) == 0)
                error ('CommunicationNetwork:ReceivePackets:InvalidTimeStep', ...
                    '** Time step must be an integer within [0, %d] **', this.simTime);    
            end
            % first output parameter is number of packets that will leave the network
            numPackets = this.numPacketsOut(timeStep);
            % the second output parameter is(are) the packet(s) that will leave
            % the network at current time step
            if numPackets == 0
                % no packets
                packets = [];
            else
                packets(1) = this.nwData{timeStep, 1};
                for i=2:numPackets
                    packets(i) = this.nwData{timeStep, i}; %#ok
                end
            end
        end
        
        %% reset
        function reset(this)
            % init network data
            this.numPacketsOut = zeros(this.simTime, 1);
            % To keep track of packets that will leave the network to the same
            % time step the network data array is two-dimensional
            this.nwData = {};
        end
    end % methods public
end % classdef

