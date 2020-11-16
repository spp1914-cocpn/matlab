classdef CommunicationNetwork < handle
    % Defines a simple communication network which allows to send and
    % receive data packets, based on the original implementation by Jörg
    % Fischer and Maxim Dolgov.
    %
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
      
    properties (SetAccess = private)
        % map with data packets to store;
        % (timestep) -> (array of DataPacket)
        nwData;
        % number of packets that leave the network to current time step;
        % Index is the time step; (array of positive integer)
        numPacketsOut = [];       
    end
    
    properties (GetAccess = private, SetAccess = immutable)
        % length of simulation time in time steps; (positive integer)
        simTime(1,1) double {mustBePositive, mustBeInteger} = 1;
        % maximal time delay;
        % The network does not output a data packet that has a longer time 
        % delay; (nonnegative integer)
        maxDelay(1,1) double {Validator.validateMaxPacketDelay(maxDelay, 0)} = 1;
    end
  
    methods (Access = public)
        %% CommunicationNetwork
        function this = CommunicationNetwork(simTime, maxDelay)
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
            % Returns:
            %   << this (CommunicationNetwork)
            %      A new CommunicationNetwork instance.

            this.simTime = simTime;
            this.maxDelay = maxDelay;

            % init network data
            this.numPacketsOut = zeros(simTime, 1);
            this.nwData = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
        end 
  
        %% sendPacket
        function sendPacket(this, packet, packetDelay, timeStep)
            % Send a data packet that shall be received with a given delay.
            % If the specified delay is larger than the maximum delay of the network, it is dropped.
            % The packet is also dropped if the time of reception is out of
            % bounds, i.e., after the specified simulation time.
            %
            % Parameters:
            %   >> packet (DataPacket)
            %      The DataPacket to be transmitted.
            %
            %   >> packetDelay (Nonnegative integer)
            %      A nonnegative integer denoting the transmission delay (in time steps)
            %      this packet is going to have.
            %
            %   >> timeStep (Nonnegative integer)
            %      A positive integer denoting the time step the packet is sent.
              
            assert(Checks.isScalarIn(timeStep, 0, this.simTime) && mod(timeStep, 1) == 0, ...
                'CommunicationNetwork:SendPacket:InvalidTimeStep', ...
                '** Time step must be an integer within [0, %d] **', this.simTime);    
            assert(Checks.isClass(packet, 'DataPacket'), ...            
                'CommunicationNetwork:SendPacket:InvalidPacket', ...
                '** Input parameter <packet> (sent data packet) is not of class DataPacket **');
            assert(Checks.isNonNegativeScalar(packetDelay) && mod(packetDelay, 1) == 0, ...
                'CommunicationNetwork:SendPacket:InvalidPacketDelay', ...
                '** Input parameter <packetDelay> (time delay to experience by the packet) must be a nonnegative integer **');
            % error if packet time stamp does not correspond to transmission time    
            assert(packet.timeStamp == timeStep, ...    
                'CommunicationNetwork:SendPacket:InvalidPacketTimeStamp', ...
                '** Packet time stamp does not correspond to transmission time (time step) **');
            
            if packetDelay <= this.maxDelay && packetDelay <= this.simTime - timeStep
                % handle packet as its delay is not too large and will leave
                % the network before simulation ends
                packet.packetDelay = packetDelay;
                if ~this.nwData.isKey(timeStep + packetDelay)
                    this.nwData(timeStep + packetDelay) = packet;
                else
                    this.nwData(timeStep + packetDelay) = [this.nwData(timeStep + packetDelay), packet];
                end
                this.numPacketsOut(timeStep + packetDelay) = this.numPacketsOut(timeStep + packetDelay) + 1;
            end
        end
        
        %% receivePackets
        function [numPackets, packets] = receivePackets(this, timeStep)
            % Get all packets that are received at a given time step.
            %
            % Parameters:
            %   >> timeStep (Positive integer)
            %      A positive integer denoting the time step.
            %
            % Returns:
            %   << numPackets (Nonnegative integer)
            %      The number of packets that are received at the given time step.
            %
            %   << packets (Row vector of DataPackets)
            %      The received packets, column-wise arranged (no particular order).
            
            assert(Checks.isScalarIn(timeStep, 0, this.simTime) && mod(timeStep, 1) == 0, ...
                'CommunicationNetwork:ReceivePackets:InvalidTimeStep', ...
                '** Time step must be an integer within [0, %d] **', this.simTime);    
            
            % first output parameter is number of packets that will leave the network
            numPackets = this.numPacketsOut(timeStep);
            % the second output parameter is(are) the packet(s) that will leave
            % the network at current time step
            if numPackets == 0
                % no packets
                packets = [];
            else
                packets = this.nwData(timeStep);
                this.nwData.remove(timeStep); % clear the map, remove the key       
            end
        end
        
        %% reset
        function reset(this)
            % Reset the communication network by dropping all packets that
            % are still to be transmitted.
            
            this.numPacketsOut = zeros(this.simTime, 1);
            this.nwData.remove(this.nwData.keys());
        end
        
        %% getNumBufferedPackets
        function bufferSize = getNumBufferedPackets(this, timeStep)
            % Get the number of packets that are still buffered, i.e., still to be transmitted at the given time step or later.
            %
            % Parameters:
            %   >> timeStep (Nonnegative integer)
            %      A positive integer denoting a time step.
            %
            % Returns:
            %   << bufferSize (Nonnegative integer)
            %      The number of packets still to be transmitted at the given time step or later.
            
            assert(Checks.isScalarIn(timeStep, 0, this.simTime) && mod(timeStep, 1) == 0, ...
                'CommunicationNetwork:GetNumBufferedPackets:InvalidTimeStep', ...
                '** Time step must be an integer within [0, %d] **', this.simTime); 
            
            bufferSize = sum(this.numPacketsOut(timeStep:this.simTime));
        end
        
        %% getPacketRate
        function packetRate = getPacketRate(this, timeStep, windowSize)
            % Get the packet rate, computed as a trailing sum of number of packets transmitted in the last
            % windowSize time steps (including the given time step).
            %
            % Parameters:
            %   >> timeStep (Nonnegative integer)
            %      A positive integer denoting a time step.
            %
            %   >> windowSize (Nonnegative integer)
            %      A nonnegative integer specifying the time window to be used for the computation of the packet rate.
            %      If windowSize specifies the number of time steps per
            %      second, the returned value will thus denote the packet
            %      rate given in packets/s.
            %
            % Returns:
            %   << packetRate (Nonnegative integer)
            %      The number of packets transmitted within the specified time window.
            
            assert(Checks.isScalarIn(timeStep, 0, this.simTime) && mod(timeStep, 1) == 0, ...
                'CommunicationNetwork:GetPacketRate:InvalidTimeStep', ...
                '** Time step must be an integer within [0, %d] **', this.simTime); 
            
            assert(Checks.isNonNegativeScalar(windowSize) && mod(windowSize, 1) == 0, ...
                'CommunicationNetwork:GetPacketRate:InvalidWindowSize', ...
                '** window size must be a nonnegative integer **'); 
            
            beginIdx = max(1, timeStep - windowSize + 1); % trailing sum
            packetRate = sum(this.numPacketsOut(beginIdx:timeStep));
        end
    end
end

