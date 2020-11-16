classdef DataPacketTest < matlab.unittest.TestCase
    % Test cases for DataPacket.
    
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
    
    properties
        srcAddress;
        dstAddress;
        payload;
        timestamp;
        id;
        delay;
        
        packetUnderTest;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.srcAddress = 1;
            this.dstAddress = 5;
            this.timestamp = 10;
            this.id = 42;
            this.delay = 3;

            this.payload = gallery('moler', 12);
            
            this.packetUnderTest = DataPacket(this.payload, this.timestamp, this.id);
        end
    end
    
    methods (Test)
        %% testDataPacketInvalidId
        function testDataPacketInvalidId(this)
            invalidId = -1; % ids must be positive
            this.verifyError(@() DataPacket(this.payload, this.timestamp, invalidId), 'MATLAB:validators:mustBePositive');
            
            invalidId = 2.5; % fractional ids are not permitted
            this.verifyError(@() DataPacket(this.payload, this.timestamp, invalidId), 'MATLAB:validators:mustBeInteger');
        end
        
        %% testDataPacketInvalidTimestamp
        function testDataPacketInvalidTimestamp(this)
            invalidTimestamp = 0; % time stamp must be positive
            this.verifyError(@() DataPacket(this.payload, invalidTimestamp, this.id), 'MATLAB:validators:mustBePositive');
            
            invalidTimestamp = 2.5; % fractional time stamps are not permitted
            this.verifyError(@() DataPacket(this.payload, invalidTimestamp, this.id), 'MATLAB:validators:mustBeInteger');
        end
        
        %% testDataPacket
        function testDataPacket(this)
            packet = DataPacket(this.payload, this.timestamp, this.id);
            
            this.verifyFalse(packet.isAck);
            this.verifyEqual(packet.timeStamp, this.timestamp);
            this.verifyEqual(packet.id, this.id);
            this.verifyEqual(packet.payload, this.payload);
        end
        
        %% testSetPacketDelayInvalidDelay
        function testSetPacketDelayInvalidDelay(this)
            invalidDelay = -1; % delay must not be negative
            this.verifyError(@() this.setPacketDelay(invalidDelay), 'MATLAB:validators:mustBeNonnegative');
            
            invalidDelay = 2.5;  % delay in time stamps, thus not fractional
            this.verifyError(@() this.setPacketDelay(invalidDelay), 'MATLAB:validators:mustBeInteger');
        end
        
        %% testSetPacketDelay
        function testSetPacketDelay(this)
            this.setPacketDelay(this.delay);
            
            this.verifyEqual(this.packetUnderTest.packetDelay, this.delay);
        end
        
        %% testSetSourceAddressInvalidAddress
        function testSetSourceAddressInvalidAddress(this)
            invalidAddress = 0; % addresses must be > 0
            this.verifyError(@() this.setSourceAddress(invalidAddress), 'MATLAB:validators:mustBePositive');
            
            invalidAddress = 2.5;  % not fractional
            this.verifyError(@() this.setSourceAddress(invalidAddress), 'MATLAB:validators:mustBeInteger');
        end
        
        %% testSetSourceAddress
        function testSetSourceAddress(this)
            this.setSourceAddress(this.srcAddress);
            
            this.verifyEqual(this.packetUnderTest.sourceAddress, this.srcAddress);
        end
        
        %% testSetDestinationAddressInvalidAddress
        function testSetDestinationAddressInvalidAddress(this)
            invalidAddress = 0; % addresses must be > 0
            this.verifyError(@() this.setDestinationAddress(invalidAddress), 'MATLAB:validators:mustBePositive');
            
            invalidAddress = 2.5;  % not fractional
            this.verifyError(@() this.setDestinationAddress(invalidAddress), 'MATLAB:validators:mustBeInteger');
        end
        
        %% testSetDestinationAddress
        function testSetDestinationAddress(this)
            this.setDestinationAddress(this.dstAddress);
            
            this.verifyEqual(this.packetUnderTest.destinationAddress, this.dstAddress);
        end       
        
        %% testIsNewerThan
        function testIsNewerThan(this)
            otherTimestamp = this.timestamp + 1;
            otherId = this.id + 1;
            otherPacket = DataPacket(this.payload, otherTimestamp, otherId);
            
            this.verifyFalse(this.packetUnderTest.isNewerThan(otherPacket));
            this.verifyTrue(otherPacket.isNewerThan(this.packetUnderTest));
        end
        
        %% testCreateAckForDataPacket
        function testCreateAckForDataPacket(this)
            ackTimestamp = this.timestamp + 2;
            % prepare the data packet
            this.setSourceAddress(this.srcAddress);
            this.setDestinationAddress(this.dstAddress);
            
            ackPacket = this.packetUnderTest.createAckForDataPacket(ackTimestamp);
            actualAckPayload = ackPacket.payload;
            % check if the addresses are set correctly
            this.verifyEqual(ackPacket.sourceAddress, this.dstAddress);
            this.verifyEqual(ackPacket.destinationAddress, this.srcAddress);
            % is ack?
            this.verifyTrue(ackPacket.isAck);
            % check the payload of the ack: should be only the time stamp
            % of the ack'ed packet
            this.verifyClass(actualAckPayload, ?double);
            this.verifyEqual(actualAckPayload, this.timestamp);
            
            % now piggy-back payload
            ackPayload = this.payload;
            
            ackPacket = this.packetUnderTest.createAckForDataPacket(ackTimestamp, ackPayload);
            actualAckPayload = ackPacket.payload;
            % check if the addresses are set correctly
            this.verifyEqual(ackPacket.sourceAddress, this.dstAddress);
            this.verifyEqual(ackPacket.destinationAddress, this.srcAddress);
            % is ack?
            this.verifyTrue(ackPacket.isAck);
            % now the ack payload should be a 1x2 cell array
            this.verifyClass(actualAckPayload, ?cell);
            this.verifySize(actualAckPayload, [1 2]);
            this.verifyEqual(actualAckPayload{1}, this.timestamp);
            this.verifyEqual(actualAckPayload{2}, ackPayload);
        end
        
        %% testGetNextId
        function testGetNextId(this)
            clear DataPacket;
            
            this.verifyEqual(DataPacket.getNextId(), 1);
            this.verifyEqual(DataPacket.getNextId(), 2);
        end
    end
    
    methods (Access = private)
        function setPacketDelay(this, delay)
            this.packetUnderTest.packetDelay = delay;
        end
        
        function setSourceAddress(this, srcAddr)
            this.packetUnderTest.sourceAddress = srcAddr;
        end
        
        function setDestinationAddress(this, dstAddr)
            this.packetUnderTest.destinationAddress = dstAddr;
        end
    end
end

