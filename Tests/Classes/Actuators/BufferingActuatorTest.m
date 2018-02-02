classdef BufferingActuatorTest < matlab.unittest.TestCase
    % Test cases for BufferingActuator.
    
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
    
    properties
        dimU;
        sequenceLength;
        maxSequenceDelay;
        defaultU;
        dimSequence;
        
        srcAddress;
        dstAddress;
        payload;
        timestamp
        id;
        dataPacket;
        dataPacket2;
        dataPacket3;
        
        controlSequence;
        controlSequence2;
        controlSequence3;
        
        actuatorUnderTest;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.dimU = 4;
            this.sequenceLength = 5;
            this.maxSequenceDelay = 7;
            this.timestamp = 10;
            
            this.dimSequence = this.dimU * this.sequenceLength;
            this.defaultU = ones(this.dimU, 1);
            
            maxValue = max(this.dimU, this.sequenceLength);
            tmp = gallery('moler', maxValue);
            this.controlSequence = tmp(1:this.dimU, 1:this.sequenceLength);
            tmp = gallery('minij', maxValue);
            this.controlSequence2 = tmp(1:this.dimU, 1:this.sequenceLength);
            tmp = gallery('binomial', maxValue);
            this.controlSequence3 = tmp(1:this.dimU, 1:this.sequenceLength);
                        
            this.actuatorUnderTest = BufferingActuator(this.sequenceLength, this.maxSequenceDelay, this.defaultU);
            
            this.initControllerPackets();
        end
        
    end
    
    methods (Access = private)
        function initControllerPackets(this)
            this.srcAddress = 1;
            this.dstAddress = 5;
            
            this.id = 1;
              
            this.dataPacket = DataPacket(this.controlSequence, this.timestamp -1, this.id);
            this.dataPacket.sourceAddress = this.srcAddress;
            this.dataPacket.destinationAddress = this.dstAddress;
            
            % prepare packet with smallest delay
            this.dataPacket2 = DataPacket(this.controlSequence2, this.timestamp, this.id + 1);
            this.dataPacket2.sourceAddress = this.srcAddress;
            this.dataPacket2.destinationAddress = this.dstAddress;
            
            % prepare packet 
            this.dataPacket3 = DataPacket(this.controlSequence3, this.timestamp -2, this.id + 2);
            this.dataPacket3.sourceAddress = this.srcAddress;
            this.dataPacket3.destinationAddress = this.dstAddress;
            
            this.dataPacket.packetDelay = this.timestamp - this.dataPacket.timeStamp; % 1 
            this.dataPacket2.packetDelay = this.timestamp - this.dataPacket2.timeStamp; % 0
            this.dataPacket3.packetDelay = this.timestamp - this.dataPacket3.timeStamp; % 2
        end
    end
    
    methods (Test)
        %% testBufferingActuatorInvalidSequenceLength
        function testBufferingActuatorInvalidSequenceLength(this)
            expectedErrId = 'Validator:ValidateSequenceLength:InvalidSequenceLength';
            
            invalidSequenceLength = 0; % must be positive
            this.verifyError(@() BufferingActuator(invalidSequenceLength, this.maxSequenceDelay, this.defaultU), ...
                expectedErrId);
            
            invalidSequenceLength = 4.5; % must not be fractional
            this.verifyError(@() BufferingActuator(invalidSequenceLength, this.maxSequenceDelay, this.defaultU), ...
                expectedErrId);
            
            invalidSequenceLength = inf; % must not be inf
            this.verifyError(@() BufferingActuator(invalidSequenceLength, this.maxSequenceDelay, this.defaultU), ...
                expectedErrId);
        end
        
        %% testBufferingActuatorInvalidMaxPacketDelay
        function testBufferingActuatorInvalidMaxPacketDelay(this)
            expectedErrId = 'BufferingActuator:InvalidMaxPacketDelay';
            
            invalidMaxPacketDelay = -1; % must not be negative
            this.verifyError(@() BufferingActuator(this.sequenceLength, invalidMaxPacketDelay, this.defaultU), ...
                expectedErrId);
            
            invalidMaxPacketDelay = 4.5; % must not be fractional
            this.verifyError(@() BufferingActuator(this.sequenceLength, invalidMaxPacketDelay, this.defaultU), ...
                expectedErrId);
        end
        
        %% testBufferingActuatorInvalidDefaultU
        function testBufferingActuatorInvalidDefaultU(this)
            expectedErrId = 'BufferingActuator:InvalidDefaultU';
            
            invalidDefaultU = ones(this.dimU, this.dimU); % must be a vector
            this.verifyError(@() BufferingActuator(this.sequenceLength, this.maxSequenceDelay, invalidDefaultU), ...
                expectedErrId);
            
            invalidDefaultU = [4 nan]; % must not contain nan
            this.verifyError(@() BufferingActuator(this.sequenceLength, this.maxSequenceDelay, invalidDefaultU), ...
                expectedErrId);
            
            invalidDefaultU = [inf 4]; % must not contain inf
            this.verifyError(@() BufferingActuator(this.sequenceLength, this.maxSequenceDelay, invalidDefaultU), ...
                expectedErrId);
        end
        
        %% testBufferingActuator
        function testBufferingActuator(this)
            % inf is allowed for maxPacketDelay
            actuator = BufferingActuator(this.sequenceLength, inf, this.defaultU);
            this.verifyEqual(actuator.dimU, this.dimU);
            this.verifyEqual(actuator.defaultInput, this.defaultU);
            this.verifyEmpty(actuator.bufferedPacket);
            
            % 0 is allowed for maxPacketDelay
            actuator = BufferingActuator(this.sequenceLength, 0, this.defaultU);
            this.verifyEqual(actuator.dimU, this.dimU);
            this.verifyEqual(actuator.defaultInput, this.defaultU);
            this.verifyEmpty(actuator.bufferedPacket);
            
            % row vector is allowed for defaultU
            actuator = BufferingActuator(this.sequenceLength, this.maxSequenceDelay, this.defaultU');
            this.verifyEqual(actuator.dimU, this.dimU);
            % check if internally stored as column vector
            this.verifyEqual(actuator.defaultInput, this.defaultU);
            this.verifyEmpty(actuator.bufferedPacket);
        end
        
        %% testGetCurrentInputInvalidTimestep
        function testGetCurrentInputInvalidTimestep(this)
            expectedErrId = 'BufferingActuator:GetCurrentInput:InvalidTimeStep';
            
            invalidTimestep = -1; % must not be negative
            this.verifyError(@() this.actuatorUnderTest.getCurrentInput(invalidTimestep), expectedErrId);
            
            invalidTimestep = 2.5; % must be an integer
            this.verifyError(@() this.actuatorUnderTest.getCurrentInput(invalidTimestep), expectedErrId);
            
            invalidTimestep = inf; % must be finite
            this.verifyError(@() this.actuatorUnderTest.getCurrentInput(invalidTimestep), expectedErrId);
        end
        
        %% testGetCurrentInput
        function testGetCurrentInput(this)
            timestep = this.timestamp;
            % first, try to get an input if no sequence is buffered
            [actualInput, actualMode] = this.actuatorUnderTest.getCurrentInput(timestep);
            
            expectedMode = this.sequenceLength + 1; % default input is applied
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
            
            this.actuatorUnderTest.processControllerPackets([this.dataPacket; this.dataPacket2; this.dataPacket3]);
            % dataPacket2 is buffered as it has smallest delay
             
            % now get the input for the the current timestep
            for j = 0:this.sequenceLength - 1
                [actualInput, actualMode] = this.actuatorUnderTest.getCurrentInput(timestep + j);
                expectedMode = j + 1;
                expectedInput = this.controlSequence2(:, expectedMode);
                this.verifyEqual(actualMode, expectedMode);
                this.verifyEqual(actualInput, expectedInput);
            end
            % now try to get the input for a future time steps, which is not
            % part of the sequence any more
            [actualInput, actualMode] = this.actuatorUnderTest.getCurrentInput(timestep + this.sequenceLength);
            expectedMode = this.sequenceLength + 1; % default input is applied
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
            
            [actualInput, actualMode] = this.actuatorUnderTest.getCurrentInput(timestep + this.sequenceLength + 1);
            expectedMode = this.sequenceLength + 1; % default input is applied
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
        end
        
        %% testProcessControllerPacketsInvalidControllerPackets
        function testProcessControllerPacketsInvalidControllerPackets(this)
            expectedErrId = 'BufferingActuator:ProcessControllerPackets:InvalidControllerPackets';
            
            % no cell arrays of DataPackets allowed
            invalidPackets = {this.dataPacket; this.dataPacket2; this.dataPacket3};
            this.verifyError(@() this.actuatorUnderTest.processControllerPackets(invalidPackets), expectedErrId);
            
            % no arrays/matrices of other type allowed
            invalidPackets = eye(4);
            this.verifyError(@() this.actuatorUnderTest.processControllerPackets(invalidPackets), expectedErrId);
        end
        
        %% testProcessControllerPacketsInvalidControlSequences
        function testProcessControllerPacketsInvalidControlSequences(this)
            expectedErrId = 'BufferingActuator:ProcessControllerPackets:InvalidControlSequences';
            
            % send a packet with an invalid payload: not the correct input
            % dimension
            invalidSequence = ones(this.dimU + 1, this.sequenceLength);
            invalidPacket = DataPacket(invalidSequence, this.timestamp, this.id + 3);
            packets = [this.dataPacket; this.dataPacket2; invalidPacket];
            this.verifyError(@() this.actuatorUnderTest.processControllerPackets(packets), expectedErrId);
            
            % send a packet with an invalid payload: not the correct sequence length
            invalidSequence = ones(this.dimU, this.sequenceLength + 1);
            invalidPacket = DataPacket(invalidSequence, this.timestamp, this.id + 4);
            packets = [this.dataPacket; invalidPacket; this.dataPacket3];
            this.verifyError(@() this.actuatorUnderTest.processControllerPackets(packets), expectedErrId);
        end
        
        %% testProcessControllerPackets
        function testProcessControllerPackets(this)
            % first, process two packets; those with delay 0 and 2
            packets = [this.dataPacket2; this.dataPacket3];
            ackPacket = this.actuatorUnderTest.processControllerPackets(packets);
            actualBufferedPacket = this.actuatorUnderTest.bufferedPacket;
            actualAckTimestamp = ackPacket.timeStamp;
            
            this.verifyNotEmpty(actualBufferedPacket);
            this.verifyEqual(actualBufferedPacket, this.dataPacket2);
            this.verifySameHandle(actualBufferedPacket, this.dataPacket2);
            this.verifyEqual(actualAckTimestamp, this.timestamp);
                       
            % now process another packet: that with delay 1
            ackPacket = this.actuatorUnderTest.processControllerPackets(this.dataPacket);
            actualBufferedPacket = this.actuatorUnderTest.bufferedPacket;
            
            % still the same packet and no ack expected
            this.verifyEmpty(ackPacket);
            this.verifyEqual(actualBufferedPacket, this.dataPacket2);
            this.verifySameHandle(actualBufferedPacket, this.dataPacket2);
             
        end
        
        %% testProcessControllerPacketsDelayTooLarge
        function testProcessControllerPacketsDelayTooLarge(this)
            % receive a packet whose delay is way too large to be accepted,
            % although the buffer is empty
            delayedPacket = DataPacket(this.controlSequence, this.timestamp - this.maxSequenceDelay -1, this.id);
            delayedPacket.packetDelay = this.timestamp - delayedPacket.timeStamp;
            ackPacket = this.actuatorUnderTest.processControllerPackets(delayedPacket);
            
            % still the same packet and no ack expected
            this.verifyEmpty(ackPacket);
            this.verifyEmpty(this.actuatorUnderTest.bufferedPacket);
        end
        
        %% testReset
        function testReset(this)
            % first, process two packets; those with delay 0 and 2
            % thus, a packet is buffered
            packets = [this.dataPacket2; this.dataPacket3];
            this.actuatorUnderTest.processControllerPackets(packets);
            
            this.assertNotEmpty(this.actuatorUnderTest.bufferedPacket);
            
            % now call reset
            this.actuatorUnderTest.reset();
            % buffer should be empty now
            this.verifyEmpty(this.actuatorUnderTest.bufferedPacket);
        end
    end
end

