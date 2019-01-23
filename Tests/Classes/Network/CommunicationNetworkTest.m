classdef CommunicationNetworkTest < matlab.unittest.TestCase
   % Test cases for CommunicationNetwork.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        maxDelay;
        simTime;
        
        networkUnderTest;
        
        srcAddress;
        dstAddress;
        payload;
        timestamp
        id;
        delay;
        dataPacket;
        dataPacket2;
        dataPacket3;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.maxDelay = 10;
            this.simTime = 100;
            
            this.networkUnderTest = CommunicationNetwork(this.simTime, this.maxDelay);
            
            this.srcAddress = 1;
            this.dstAddress = 5;
            this.timestamp = 10;
            this.id = 42;
            this.delay = 3;
            this.payload = gallery('moler', 12);
            this.dataPacket = DataPacket(this.payload, this.timestamp, this.id);
            this.dataPacket.sourceAddress = this.srcAddress;
            this.dataPacket.destinationAddress = this.dstAddress;
            
            % prepare packet to be sent one time step later
            this.dataPacket2 = DataPacket(this.payload, this.timestamp + 1, this.id + 1);
            this.dataPacket2.sourceAddress = this.srcAddress;
            this.dataPacket2.destinationAddress = this.dstAddress;
            
            % prepare packet to be sent two time steps later
            this.dataPacket3 = DataPacket(this.payload, this.timestamp + 2, this.id + 2);
            this.dataPacket3.sourceAddress = this.srcAddress;
            this.dataPacket3.destinationAddress = this.dstAddress;

        end
    end
    
    methods (Test)
        %% testCommunicationNetworkInvalidSimTime
        function testCommunicationNetworkInvalidSimTime(this)
            expectedErrId = 'CommunicationNetwork:InvalidSimTime';
            
            invalidSimTime = 0; % must be > 0
            this.verifyError(@() CommunicationNetwork(invalidSimTime, this.maxDelay), expectedErrId);
            
            invalidSimTime = 1000.5; % must be integer
            this.verifyError(@() CommunicationNetwork(invalidSimTime, this.maxDelay), expectedErrId);
            
            invalidSimTime = inf; % must be finite
            this.verifyError(@() CommunicationNetwork(invalidSimTime, this.maxDelay), expectedErrId);
        end
        
        %% testCommunicationNetworkInvalidMaxDelay
        function testCommunicationNetworkInvalidMaxDelay(this)
            expectedErrId = 'CommunicationNetwork:InvalidMaxDelay';
            
            invalidMaxDelay = -10; % must be >= 0
            this.verifyError(@() CommunicationNetwork(this.simTime, invalidMaxDelay), expectedErrId);
            
            invalidMaxDelay = 1000.5; % must be integer
            this.verifyError(@() CommunicationNetwork(this.simTime, invalidMaxDelay), expectedErrId);
            
            invalidMaxDelay = inf; % must be finite
            this.verifyError(@() CommunicationNetwork(this.simTime, invalidMaxDelay), expectedErrId);
        end
        
        %% testCommunicationNetworkInvalidProtocol
        function testCommunicationNetworkInvalidProtocol(this)
            expectedErrId = 'CommunicationNetwork:InvalidProtocol';
            
            invalidProtocol = -10; % must be a char array
            this.verifyError(@() CommunicationNetwork(this.simTime, this.maxDelay, invalidProtocol), expectedErrId);
            
            expectedErrId = 'CommunicationNetwork:InvalidProtocolChar';
            invalidProtocol = 'TPC'; % neither 'TCP', nor 'UDP'
            this.verifyError(@() CommunicationNetwork(this.simTime, this.maxDelay, invalidProtocol), expectedErrId);
        end
        
        %% testCommunicationNetwork
        function testCommunicationNetwork(this)
            network = CommunicationNetwork(this.simTime, this.maxDelay);
            expectedNumPacketsOut = zeros(this.simTime, 1);
            
            this.verifyClass(network.nwData, ?containers.Map);
            this.verifyEmpty(network.nwData);
            
            this.verifyEqual(network.numPacketsOut, expectedNumPacketsOut);
        end
        
        %% testSendPacketInvalidTimeStep
        function testSendPacketInvalidTimeStep(this)
            expectedErrId = 'CommunicationNetwork:SendPacket:InvalidTimeStep';
            
            invalidTimestep = -1; % must not be negative
            this.verifyError(@() this.networkUnderTest.sendPacket(this.dataPacket, this.delay, invalidTimestep),...
                expectedErrId);
            
            invalidTimestep = this.simTime + 1; % must be within [0, simTime]
            this.verifyError(@() this.networkUnderTest.sendPacket(this.dataPacket, this.delay, invalidTimestep),...
                expectedErrId);
            
            invalidTimestep = 2.5; % must be integer
            this.verifyError(@() this.networkUnderTest.sendPacket(this.dataPacket, this.delay, invalidTimestep),...
                expectedErrId);
        end
        
        %% testSendPacketInvalidPacket
        function testSendPacketInvalidPacket(this)
            expectedErrId = 'CommunicationNetwork:SendPacket:InvalidPacket';
            timestep = this.timestamp;
            
            invalidPacket = zeros(3,3,3); % must be DataPacket instance
            this.verifyError(@() this.networkUnderTest.sendPacket(invalidPacket, this.delay, timestep),...
                expectedErrId);
            
            % use valid packet instance, but time mismatch
            expectedErrId = 'CommunicationNetwork:SendPacket:InvalidPacketTimeStamp';
            timestep = this.timestamp + 1;
            this.verifyError(@() this.networkUnderTest.sendPacket(this.dataPacket, this.delay, timestep),...
                expectedErrId);
        end
        
        %% testSendPacketInvalidInvalidPacketDelay
        function testSendPacketInvalidInvalidPacketDelay(this)
            expectedErrId = 'CommunicationNetwork:SendPacket:InvalidPacketDelay';
            
            invalidPacketDelay = -1; % must not be negative
            this.verifyError(@() this.networkUnderTest.sendPacket(this.dataPacket, invalidPacketDelay, this.timestamp),...
                expectedErrId);
            
            invalidPacketDelay = 3.5; % must be an integer
            this.verifyError(@() this.networkUnderTest.sendPacket(this.dataPacket, invalidPacketDelay, this.timestamp),...
                expectedErrId);
        end
        
        %% testSendPacket
        function testSendPacket(this)
            timestep = this.timestamp;
            
            this.networkUnderTest.sendPacket(this.dataPacket, this.delay, timestep);
 
            actualBufferSize = size(this.networkUnderTest.nwData);
            expectedBufferSize = [1 1];
           
            % now check the internal buffers
            this.verifyEqual(this.dataPacket.packetDelay, this.delay);
            this.verifyEqual(this.networkUnderTest.numPacketsOut(timestep + this.delay), 1);
            this.verifyEqual(sum(this.networkUnderTest.numPacketsOut), 1);
            this.verifyEqual(actualBufferSize, expectedBufferSize);
            this.verifyClass(this.networkUnderTest.nwData(timestep + this.delay), ?DataPacket);
            this.verifyEqual(this.networkUnderTest.nwData(timestep + this.delay), this.dataPacket);
            
             % now send another packet, one time step later, with the same delay
            this.networkUnderTest.sendPacket(this.dataPacket2, this.delay, timestep + 1);
             
            this.verifyEqual(this.dataPacket2.packetDelay, this.delay);
            this.verifyEqual(this.networkUnderTest.numPacketsOut(timestep + this.delay), 1);
            this.verifyEqual(this.networkUnderTest.numPacketsOut(timestep + 1 + this.delay), 1);
             % two packets in total now
            this.verifyEqual(sum(this.networkUnderTest.numPacketsOut), 2);
            
            actualBufferSize = size(this.networkUnderTest.nwData);
            expectedBufferSize = [2 1];
            
            this.verifyEqual(actualBufferSize, expectedBufferSize);
            this.verifyClass(this.networkUnderTest.nwData(timestep + this.delay), ?DataPacket);
            this.verifyEqual(this.networkUnderTest.nwData(timestep + this.delay), this.dataPacket);
            this.verifyClass(this.networkUnderTest.nwData(timestep + this.delay + 1), ?DataPacket);
            this.verifyEqual(this.networkUnderTest.nwData(timestep + this.delay + 1), this.dataPacket2);
%             
            % finally, send a third packet, but with a lower delay
            % should arrive at the same time as the second packet
            this.networkUnderTest.sendPacket(this.dataPacket3, this.delay - 1, timestep + 2);
             
            this.verifyEqual(this.dataPacket3.packetDelay, this.delay - 1);
            % no changes expected here
            this.verifyEqual(this.networkUnderTest.numPacketsOut(timestep + this.delay), 1);
            this.verifyClass(this.networkUnderTest.nwData(timestep + this.delay), ?DataPacket);
            this.verifyEqual(this.networkUnderTest.nwData(timestep + this.delay), this.dataPacket);
            % three packets in total now
            this.verifyEqual(this.networkUnderTest.numPacketsOut(timestep + 1 + this.delay), 2);
            this.verifyEqual(sum(this.networkUnderTest.numPacketsOut), 3);
            
            actualBufferSize = size(this.networkUnderTest.nwData);
            expectedBufferSize = [2 1]; % packet arrives at the same time

            this.verifyEqual(actualBufferSize, expectedBufferSize);
            packets = this.networkUnderTest.nwData(timestep + this.delay + 1);
            this.verifySize(packets, [1 2]);
            this.verifyClass(packets(1), ?DataPacket);
            this.verifyEqual(packets(1), this.dataPacket2);
            this.verifyClass(packets(2), ?DataPacket);
            this.verifyEqual(packets(2), this.dataPacket3);
        end
        
        %% testReceivePacketsInvalidTimeStep
        function testReceivePacketsInvalidTimeStep(this)
            expectedErrId = 'CommunicationNetwork:ReceivePackets:InvalidTimeStep';
            
            invalidTimestep = -1; % must not be negative
            this.verifyError(@() this.networkUnderTest.receivePackets(invalidTimestep), expectedErrId);
            
            invalidTimestep = this.simTime + 1; % must be within [0, simTime]
            this.verifyError(@() this.networkUnderTest.receivePackets(invalidTimestep), expectedErrId);
            
            invalidTimestep = 2.5; % must be integer
            this.verifyError(@() this.networkUnderTest.receivePackets(invalidTimestep), expectedErrId);
        end
        
        %% testReceivePackets
        function testReceivePackets(this)
            timestep = this.timestamp;
            
            % send three packets: second and third should be received at the
            % same time step
            this.networkUnderTest.sendPacket(this.dataPacket, this.delay, timestep);
            this.networkUnderTest.sendPacket(this.dataPacket2, this.delay, timestep + 1);
            this.networkUnderTest.sendPacket(this.dataPacket3, this.delay - 1, timestep + 2);
            
            this.assertEqual(size(this.networkUnderTest.nwData), [2 1]);
            
            % now try to receive the first packet
            [actualNumPackets, actualPackets] = this.networkUnderTest.receivePackets(timestep + this.delay);
            
            this.verifyEqual(actualNumPackets, 1);
            this.verifySize(actualPackets, [1 1]);
            this.verifyEqual(actualPackets(1), this.dataPacket);
            
            % check the buffer size
            this.verifyEqual(size(this.networkUnderTest.nwData), [1 1]);
            
            % now try to receive the remaining two packets
            [actualNumPackets, actualPackets] = this.networkUnderTest.receivePackets(timestep + this.delay + 1);
            
            this.verifyEqual(actualNumPackets, 2);
            this.verifySize(actualPackets, [1 2]);
            % packets should be received in the same order as they were sent
            this.verifyEqual(actualPackets(1), this.dataPacket2);
            this.verifyEqual(actualPackets(2), this.dataPacket3);
            
            % check the buffer size
            this.verifyEmpty(this.networkUnderTest.nwData);
            
            % there are now more packets to receive
            [actualNumPackets, actualPackets] = this.networkUnderTest.receivePackets(timestep + this.delay + 2);
            this.verifyEqual(actualNumPackets, 0);
            this.verifyEmpty(actualPackets);
        end
        
        %% testReset
        function testReset(this)
            timestep = this.timestamp;
            % first, send some packets
            this.networkUnderTest.sendPacket(this.dataPacket, this.delay, timestep);
            this.networkUnderTest.sendPacket(this.dataPacket2, this.delay, timestep + 1);
            this.networkUnderTest.sendPacket(this.dataPacket3, this.delay - 1, timestep + 2);
            
            % now reset the network
            this.networkUnderTest.reset();
            % the buffers should be empty again now
            this.verifyClass(this.networkUnderTest.nwData, ?containers.Map);
            this.verifyEmpty(this.networkUnderTest.nwData);
            
            this.verifyEqual(this.networkUnderTest.numPacketsOut, zeros(this.simTime, 1));
        end
        
        %% testGetNumBufferedPacketsInvalidTimestep
        function testGetNumBufferedPacketsInvalidTimestep(this)
            expectedErrId = 'CommunicationNetwork:GetNumBufferedPackets:InvalidTimeStep';
            
            invalidTimestep = -1; % must not be negative
            this.verifyError(@() this.networkUnderTest.getNumBufferedPackets(invalidTimestep), expectedErrId);
            
            invalidTimestep = this.simTime + 1; % must be within [0, simTime]
            this.verifyError(@() this.networkUnderTest.getNumBufferedPackets(invalidTimestep), expectedErrId);
            
            invalidTimestep = 2.5; % must be integer
            this.verifyError(@() this.networkUnderTest.getNumBufferedPackets(invalidTimestep), expectedErrId);
        end
        
        %% testGetNumBufferedPacketsEmptyBuffer
        function testGetNumBufferedPacketsEmptyBuffer(this)
            for k= 1:this.simTime
                this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(k), 0);
            end            
        end
        
        %% testGetNumBufferedPackets
        function testGetNumBufferedPackets(this)
            timestep = this.timestamp;
            
            this.networkUnderTest.sendPacket(this.dataPacket, this.delay, timestep);
            
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep), 1); 
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + this.delay), 1);
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + this.delay + 1), 0);
            
            % now send another packet, one time step later, with the same delay
            this.networkUnderTest.sendPacket(this.dataPacket2, this.delay, timestep + 1);
            
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep), 2);
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + 1), 2);
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + 1 + this.delay), 1);
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + 2 + this.delay), 0);
            %             
            % finally, send a third packet, but with a lower delay
            % should arrive at the same time as the second packet
            this.networkUnderTest.sendPacket(this.dataPacket3, this.delay - 1, timestep + 2);             
            
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep), 3);
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + 1), 3);
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + 1 + this.delay), 2);
            this.verifyEqual(this.networkUnderTest.getNumBufferedPackets(timestep + 2 + this.delay), 0);
        end
        
        %% testGetPacketRateInvalidTimestep
        function testGetPacketRateInvalidTimestep(this)
            expectedErrId = 'CommunicationNetwork:GetPacketRate:InvalidTimeStep';
            
            invalidTimestep = -1; % must not be negative
            this.verifyError(@() this.networkUnderTest.getPacketRate(invalidTimestep, 1), expectedErrId);
            
            invalidTimestep = this.simTime + 1; % must be within [0, simTime]
            this.verifyError(@() this.networkUnderTest.getPacketRate(invalidTimestep, 1), expectedErrId);
            
            invalidTimestep = 2.5; % must be integer
            this.verifyError(@() this.networkUnderTest.getPacketRate(invalidTimestep, 1), expectedErrId);
        end
         
        %% testGetPacketRateInvalidWindowSize
        function testGetPacketRateInvalidWindowSize(this)
            expectedErrId = 'CommunicationNetwork:GetPacketRate:InvalidWindowSize';
            
            invalidWindowSize = -1; % must not be negative
            this.verifyError(@() this.networkUnderTest.getPacketRate(1, invalidWindowSize), expectedErrId);
                        
            invalidWindowSize = 2.5; % must be integer
            this.verifyError(@() this.networkUnderTest.getPacketRate(1, invalidWindowSize), expectedErrId);
        end
        
        %% testGetPacketRateEmpty
        function testGetPacketRateEmpty(this)
            windowSize = 1;
            for k= 1:this.simTime
                this.verifyEqual(this.networkUnderTest.getPacketRate(k, windowSize), 0);
            end 
        end
        
         %% testGetPacketRate
        function testGetPacketRate(this)
            timestep = this.timestamp;
            windowSize = 1;
            
            this.networkUnderTest.sendPacket(this.dataPacket, this.delay, timestep);
            
            for k=1:timestep+this.delay -1
                this.verifyEqual(this.networkUnderTest.getPacketRate(k, windowSize), 0); 
            end
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay, windowSize), 1); 
            
            % now send another packet, one time step later, with the same delay
            this.networkUnderTest.sendPacket(this.dataPacket2, this.delay, timestep + 1);
            
            for k=1:timestep+this.delay -1
                this.verifyEqual(this.networkUnderTest.getPacketRate(k, windowSize), 0); 
            end
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay, windowSize), 1); 
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay + 1, windowSize), 1);
            % 2 packets in the two time steps expected
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay + 1, 2 * windowSize), 2);            
           
            % finally, send a third packet, but with a lower delay
            % should arrive at the same time as the second packet
            this.networkUnderTest.sendPacket(this.dataPacket3, this.delay - 1, timestep + 2);             
            
            for k=1:timestep+this.delay-1
                this.verifyEqual(this.networkUnderTest.getPacketRate(k, windowSize), 0); 
            end
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay, windowSize), 1); 
            % 2 packets in this time step expected
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay + 1, windowSize), 2);
            % 3 packets in the two time steps expected
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay + 1, 2 * windowSize), 3); 
            % sanity check: three packrets in total
            this.verifyEqual(this.networkUnderTest.getPacketRate(timestep + this.delay + 1, timestep + this.delay + 1), 3);             
        end
    end
end

