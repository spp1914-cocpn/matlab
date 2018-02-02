classdef DataPacket < handle
    % Defines a data packet used to transmit data over an abstract/virtual
    % communication network.
    %
    %
    % AUTHOR:       Jörg Fischer
    % LAST UPDATE:  Maxim Dolgov, 12.11.2013
    %               Florian Rosenthal, 27.10.2017    
    
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
    
    properties (Access = public)
        % time delay of the data packet, measured in time steps
        % (positive integer or NaN)
        packetDelay = -1;
        sourceAddress;
        destinationAddress;
        isAck@logical = false;
    end % properties
    
    properties (GetAccess = public, SetAccess = immutable)
        %% input parameters
        % time stamp that indicates time (step) when data packet was generated
        % (positive integer)
        timeStamp = 1;
        % value of data packet. Can be any data (application dependend data 
        % type) 
        payload = [];
        id; % packet id, positive integer
    end
    
    methods
        function set.packetDelay(this, delay)
            if ~Checks.isNonNegativeScalar(delay) || mod(delay, 1) ~= 0
                 error('DataPacket:InvalidPacketDelay', ...
                    '** Packet delay must be a nonnegative integer **');
            end
            this.packetDelay = delay;
        end
        
        function set.sourceAddress(this, address)
            if ~Checks.isPosScalar(address) || mod(address, 1) ~= 0
                error('DataPacket:InvalidSourceAddress', ...
                    '** Source address must be a positive integer **');
            end
            this.sourceAddress = address;
        end
        
        function set.destinationAddress(this, address)
            if ~Checks.isPosScalar(address) || mod(address, 1) ~= 0
                error('DataPacket:InvalidDestinationAddress', ...
                    '** Destination address must be a positive integer **');
            end
            this.destinationAddress = address;
        end
    end
   
    methods (Access = public)
        %% DataPacket
        function this = DataPacket(payload, timeStamp, id)
            % Class constructor.
            %
            % Parameters:
            %   >> payload
            %      The payload of the new DataPacket.
            %
            %   >> timestamp (Positive Integer)
            %      A positive integer indicating the time this packet is
            %      created and sent.
            %   
            %   >> id (Positive integer, optional)
            %      A positive integer ro identify the packet, which,
            %      however, may not be unique.
            %      If none is provided, a valid id is obtained by DataPacket.getNextId();
            %
            % Returns:
            %   << this (DataPacket)
            %      A new DataPacket instance.
            %
            if nargin == 2
                % generate id by calling getNextId()
                id = DataPacket.getNextId();
            elseif nargin == 3 && ~(Checks.isPosScalar(id) && mod(id, 1) == 0)
            	error('DataPacket:InvalidId', ...
                    '** Id of packet must be a positive integer **');   
            end
            if ~Checks.isPosScalar(timeStamp) || mod(timeStamp, 1) ~= 0
                error('DataPacket:InvalidTimeStamp', ...
                    '** Time stamp of packet must be a positive integer **');
            end
            this.payload = payload;
            this.timeStamp = timeStamp;
            this.id = id;
        end 
        
        %% isNewerThan
        function result = isNewerThan(this, otherPacket)
            if ~Checks.isClass(otherPacket, 'DataPacket')
                error('DataPacket:InvalidType', ...
                    '** <otherPacket> must be of type DataPacket **');
            end
            result = this.timeStamp > otherPacket.timeStamp;
        end
        
        %% createAckForDataPacket
        function ackPacket = createAckForDataPacket(this, timeStamp, ackPayload)
            % Create an acknowledgment packet (ACK) for this packet which
            % is intended to be returned to the sender.
            %
            % Parameters:
            %   >> timestamp (Positive Integer)
            %      A positive integer indicating the time the ACK is created and sent.
            %
            %   >> ackPayload (Default: empty matrix)
            %      Additional payload of the ACK packet to be piggy-backed
            %      to the sender.
            %
            % Returns:
            %   << this (DataPacket)
            %      A DataPacket with destinationAddress set to the sourceAddress of this packet,
            %      sourceAddress set to the destinationAddress of this
            %      packet, and time stamp of this packet as payload.
            
            %ack payload: piggy-back information if desired
            if nargin == 3
                pload = {this.timeStamp, ackPayload};
            else
                pload = this.timeStamp; % time stamp of the ACK'ed packet
            end
            
            ackPacket = DataPacket(pload, timeStamp);
            ackPacket.sourceAddress = this.destinationAddress;
            ackPacket.destinationAddress = this.sourceAddress;
            ackPacket.isAck = true;
        end
        
    end
    
    methods(Static, Access = public)
        function id = getNextId()
            persistent nextId;
            if isempty(nextId)
                nextId = 1;
            end
            id = nextId;
            nextId = nextId + 1;
        end
    end
end % classdef