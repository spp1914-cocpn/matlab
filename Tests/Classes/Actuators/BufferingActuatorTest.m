classdef BufferingActuatorTest < matlab.unittest.TestCase
    % Test cases for BufferingActuator.
    
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
        dimU;
        sequenceLength;
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
        %% initProperties
        function initProperties(this)
            this.dimU = 4;
            this.sequenceLength = 5;
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
                        
            this.actuatorUnderTest = BufferingActuator(this.sequenceLength, this.defaultU);
            
            this.initControllerPackets();
        end        
    end
    
    methods (Access = private)
        %% initControllerPackets
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
            this.verifyError(@() BufferingActuator(invalidSequenceLength, this.defaultU), expectedErrId);
            
            invalidSequenceLength = 4.5; % must not be fractional
            this.verifyError(@() BufferingActuator(invalidSequenceLength, this.defaultU), expectedErrId);
            
            invalidSequenceLength = inf; % must not be inf
            this.verifyError(@() BufferingActuator(invalidSequenceLength, this.defaultU), expectedErrId);
        end

        %% testBufferingActuatorInvalidDefaultU
        function testBufferingActuatorInvalidDefaultU(this)
            if verLessThan('matlab', '9.8')
                % Matlab R2018 or R2019
                expectedErrId = 'MATLAB:type:InvalidInputSize';
            else
                expectedErrId = 'MATLAB:validation:IncompatibleSize';
            end
            
            invalidDefaultU = ones(this.dimU, this.dimU); % must be a vector
            this.verifyError(@() BufferingActuator(this.sequenceLength, invalidDefaultU), expectedErrId);
            
            expectedErrId = 'MATLAB:validators:mustBeFinite';
            
            invalidDefaultU = [4 nan]; % must not contain nan
            this.verifyError(@() BufferingActuator(this.sequenceLength, invalidDefaultU), expectedErrId);
            
            invalidDefaultU = [inf 4]; % must not contain inf
            this.verifyError(@() BufferingActuator(this.sequenceLength, invalidDefaultU), expectedErrId);
        end
        
        %% testBufferingActuator
        function testBufferingActuator(this)
            % inf is allowed for maxPacketDelay
            actuator = BufferingActuator(this.sequenceLength, this.defaultU);
            this.verifyEqual(actuator.controlSequenceLength, this.sequenceLength);
            this.verifyEqual(actuator.maxPacketDelay, this.sequenceLength - 1);
            this.verifyEqual(actuator.dimU, this.dimU);
            this.verifyEqual(actuator.defaultInput, this.defaultU);
            this.verifyEmpty(actuator.bufferedPacket);
            
            % 0 is allowed for maxPacketDelay
            actuator = BufferingActuator(this.sequenceLength, this.defaultU);
            this.verifyEqual(actuator.controlSequenceLength, this.sequenceLength);
            this.verifyEqual(actuator.maxPacketDelay, this.sequenceLength - 1);
            this.verifyEqual(actuator.dimU, this.dimU);
            this.verifyEqual(actuator.defaultInput, this.defaultU);
            this.verifyEmpty(actuator.bufferedPacket);
            
            % row vector is allowed for defaultU
            actuator = BufferingActuator(this.sequenceLength, this.defaultU');
            this.verifyEqual(actuator.controlSequenceLength, this.sequenceLength);
            this.verifyEqual(actuator.maxPacketDelay, this.sequenceLength - 1);
            this.verifyEqual(actuator.dimU, this.dimU);
            % check if internally stored as column vector
            this.verifyEqual(actuator.defaultInput, this.defaultU);
            this.verifyEmpty(actuator.bufferedPacket);
        end
        
        %% testStepInvalidTimestep
        function testStepInvalidTimestep(this)
            expectedErrId = 'BufferingActuator:Step:InvalidTimeStep';
            
            invalidTimestep = -1; % must not be negative
            this.verifyError(@() this.actuatorUnderTest.step(invalidTimestep, []), expectedErrId);
            
            invalidTimestep = 2.5; % must be an integer
            this.verifyError(@() this.actuatorUnderTest.step(invalidTimestep, []), expectedErrId);
            
            invalidTimestep = inf; % must be finite
            this.verifyError(@() this.actuatorUnderTest.step(invalidTimestep, []), expectedErrId);
        end
        
        %% testStepInvalidControllerPackets
        function testStepInvalidControllerPackets(this)
            expectedErrId = 'BufferingActuator:Step:InvalidControllerPackets';
            timestep = this.timestamp;  
            
            % no cell arrays of DataPackets allowed
            invalidPackets = {this.dataPacket; this.dataPacket2; this.dataPacket3};
            this.verifyError(@() this.actuatorUnderTest.step(timestep, invalidPackets), expectedErrId);
            
            % no arrays/matrices of other type allowed
            invalidPackets = eye(4);
            this.verifyError(@() this.actuatorUnderTest.step(timestep, invalidPackets), expectedErrId);
        end
        
        %% testStepInvalidControlSequences
        function testStepInvalidControlSequences(this)
            expectedErrId = 'BufferingActuator:Step:InvalidControlSequences';
            timestep = this.timestamp;  
            
            % send a packet with an invalid payload: not the correct input
            % dimension
            invalidSequence = ones(this.dimU + 1, this.sequenceLength);
            invalidPacket = DataPacket(invalidSequence, this.timestamp, this.id + 3);
            packets = [this.dataPacket; this.dataPacket2; invalidPacket];
            this.verifyError(@() this.actuatorUnderTest.step(timestep, packets), expectedErrId);
            
            % send a packet with an invalid payload: not the correct sequence length
            invalidSequence = ones(this.dimU, this.sequenceLength + 1);
            invalidPacket = DataPacket(invalidSequence, this.timestamp, this.id + 4);
            packets = [this.dataPacket; invalidPacket; this.dataPacket3];
            this.verifyError(@() this.actuatorUnderTest.step(timestep, packets), expectedErrId);
        end
        
        %% testStepNoPackets
        function testStepNoPackets(this)
            expectedMode = this.sequenceLength + 1; % default input is applied
            
            timestep = this.timestamp;            
            [actualInput, actualMode, acks] = this.actuatorUnderTest.step(timestep, []);
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
            this.verifyEmpty(acks);
        end
        
        %% testStepMultipleTimesteps
        function testStepMultipleTimesteps(this)
            timestep = this.timestamp;
            % dataPacket2 has the smallest delay
            packets = [this.dataPacket; this.dataPacket2; this.dataPacket3];            
            [actualInput, actualMode, acks] = this.actuatorUnderTest.step(timestep, packets);
            
            this.verifyEqual(actualMode, 1);
            this.verifyEqual(actualInput, this.controlSequence2(:, 1));
            this.verifySize(acks, [3 1]);
            
            % check the payload of the 3 acks
            ackPayload = acks(1).payload;
            this.verifyNotEmpty(ackPayload);
            this.verifyTrue(iscell(ackPayload));
            this.verifySize(ackPayload, [1 2]);
            this.verifyEqual(ackPayload{1}, this.dataPacket.timeStamp); % the timestamp of the ACK'ed packet
            this.verifyTrue(isstruct(ackPayload{2})); % the information about actuator mode and packet delay (tau and theta) 
            this.verifyTrue(isfield(ackPayload{2}, 'timeStep'));
            this.verifyTrue(isfield(ackPayload{2}, 'theta'));
            this.verifyTrue(isfield(ackPayload{2}, 'tau')); % the packet delay
            this.verifyEqual(ackPayload{2}.timeStep, this.dataPacket.timeStamp);
            this.verifyEqual(ackPayload{2}.theta, this.sequenceLength + 1); % max mode
            this.verifyEqual(ackPayload{2}.tau, this.dataPacket.packetDelay); 
                        
            ackPayload = acks(2).payload;
            this.verifyNotEmpty(ackPayload);
            this.verifyTrue(iscell(ackPayload));
            this.verifySize(ackPayload, [1 2]);
            this.verifyEqual(ackPayload{1}, this.dataPacket2.timeStamp); % the timestamp of the ACK'ed packet
            this.verifyTrue(isstruct(ackPayload{2})); % the information about actuator mode and packet delay (tau and theta) 
            this.verifyTrue(isfield(ackPayload{2}, 'timeStep'));
            this.verifyTrue(isfield(ackPayload{2}, 'theta'));
            this.verifyTrue(isfield(ackPayload{2}, 'tau')); % the packet delay
            this.verifyEqual(ackPayload{2}.timeStep, this.dataPacket2.timeStamp);
            this.verifyEqual(ackPayload{2}.theta, 1); % first mode 
            this.verifyEqual(ackPayload{2}.tau, 0); % delay is zero
            
            ackPayload = acks(3).payload;
            this.verifyNotEmpty(ackPayload);
            this.verifyTrue(iscell(ackPayload));
            this.verifySize(ackPayload, [1 2]);
            this.verifyEqual(ackPayload{1}, this.dataPacket3.timeStamp); % the timestamp of the ACK'ed packet
            this.verifyTrue(isstruct(ackPayload{2})); % the information about actuator mode and packet delay (tau and theta) 
            this.verifyTrue(isfield(ackPayload{2}, 'timeStep'));
            this.verifyTrue(isfield(ackPayload{2}, 'theta'));
            this.verifyTrue(isfield(ackPayload{2}, 'tau')); % the packet delay
            this.verifyEqual(ackPayload{2}.timeStep, this.dataPacket3.timeStamp);
            this.verifyEqual(ackPayload{2}.theta, this.sequenceLength + 1); % max mode
            this.verifyEqual(ackPayload{2}.tau, this.dataPacket3.packetDelay); 
            
            % now get the input for the next timesteps
            for j = 1:this.sequenceLength - 1
                [actualInput, actualMode, acks] = this.actuatorUnderTest.step(timestep + j, []);
                expectedMode = j + 1;
                expectedInput = this.controlSequence2(:, expectedMode);
                this.verifyEqual(actualMode, expectedMode);
                this.verifyEqual(actualInput, expectedInput);
                this.verifyEmpty(acks);
            end
            
            % now try to get the input for a future time steps, which is not
            % part of the sequence any more
            [actualInput, actualMode, acks] = this.actuatorUnderTest.step(timestep + this.sequenceLength, []);
            expectedMode = this.sequenceLength + 1; % default input is applied
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
            this.verifyEmpty(acks);
        end       
        
        %% testStepDelayTooLarge
        function testStepDelayTooLarge(this)
            timestep = this.timestamp;
            % receive a packet whose delay is way too large to be accepted,
            % although the buffer is empty
            delayedPacket = DataPacket(this.controlSequence, this.timestamp - this.sequenceLength, this.id);
            delayedPacket.packetDelay = this.timestamp - delayedPacket.timeStamp;
            [actualInput, actualMode, ackPacket] = this.actuatorUnderTest.step(timestep, delayedPacket);
            
            % still the same packet           
            this.verifyEmpty(this.actuatorUnderTest.bufferedPacket);
            this.verifyEqual(actualInput, this.defaultU);
            this.verifyEqual(actualMode, this.sequenceLength + 1); % max mode
            
            this.verifyNotEmpty(ackPacket);
            ackPayload = ackPacket.payload;
            % check the payload of the ack
            this.verifyNotEmpty(ackPayload);
            this.verifyTrue(iscell(ackPayload));
            this.verifySize(ackPayload, [1 2]);
            this.verifyEqual(ackPayload{1}, delayedPacket.timeStamp); % the timestamp of the ACK'ed packet
            this.verifyTrue(isstruct(ackPayload{2})); % the information about actuator mode and packet delay (tau and theta) 
            this.verifyTrue(isfield(ackPayload{2}, 'timeStep'));
            this.verifyTrue(isfield(ackPayload{2}, 'theta'));
            this.verifyTrue(isfield(ackPayload{2}, 'tau')); % the packet delay
            this.verifyEqual(ackPayload{2}.timeStep, delayedPacket.timeStamp);
            this.verifyEqual(ackPayload{2}.theta, this.sequenceLength + 1); % max mode
            this.verifyEqual(ackPayload{2}.tau, delayedPacket.packetDelay); 
        end
        
        %% testReset
        function testReset(this)
            % first, process two packets; those with delay 0 and 2
            % thus, a packet is buffered
            timestep = this.timestamp;
            packets = [this.dataPacket2; this.dataPacket3];
            this.actuatorUnderTest.step(timestep, packets);
            
            this.assertNotEmpty(this.actuatorUnderTest.bufferedPacket);
            
            % now call reset
            this.actuatorUnderTest.reset();
            % buffer should be empty now
            this.verifyEmpty(this.actuatorUnderTest.bufferedPacket);
        end
        
        %% testChangeControlSequenceLengthInvalidSequenceLength
        function testChangeControlSequenceLengthInvalidSequenceLength(this)
            expectedErrId = 'Validator:ValidateSequenceLength:InvalidSequenceLength';
            
            invalidSequenceLength = 0; % must be positive
            this.verifyError(@() this.actuatorUnderTest.changeControlSequenceLength(invalidSequenceLength), ...
                expectedErrId);
            
            invalidSequenceLength = 4.5; % must not be fractional
            this.verifyError(@() this.actuatorUnderTest.changeControlSequenceLength(invalidSequenceLength), ...
                expectedErrId);
            
            invalidSequenceLength = inf; % must not be inf
            this.verifyError(@() this.actuatorUnderTest.changeControlSequenceLength(invalidSequenceLength), ...
                expectedErrId);
        end
        
        %% testChangeControlSequenceLengthInvalidControlSequence
        function testChangeControlSequenceLengthInvalidControlSequence(this)
            expectedErrId = 'BufferingActuator:Step:InvalidControlSequences';
            
            % change the sequence length
            this.actuatorUnderTest.changeControlSequenceLength(this.sequenceLength + 1);
            
            % send a packet with an invalid payload: not the correct sequence length
            timestep = this.timestamp;
            invalidSequence = ones(this.dimU, this.sequenceLength);
            invalidPacket = DataPacket(invalidSequence, this.timestamp, this.id + 4);            
            packets = [this.dataPacket; invalidPacket; this.dataPacket3];
            this.verifyError(@() this.actuatorUnderTest.step(timestep, packets), expectedErrId);
        end
        
        %% testChangeControlSequenceShorterSequence
        function testChangeControlSequenceShorterSequence(this)
            newSeqLength = 1;
            timestep = this.timestamp;
            
            % first, process two packets; those with delay 0 and 2
            % thus, a packet is buffered
            packets = [this.dataPacket2; this.dataPacket3];
            [actualInput, actualMode] = this.actuatorUnderTest.step(timestep, packets);            
        
            this.actuatorUnderTest.changeControlSequenceLength(newSeqLength);
            this.assertNotEmpty(this.actuatorUnderTest.bufferedPacket); % a packet should still be buffered
            
            % the first element of the buffered sequence should still be
            % available            
            expectedMode = newSeqLength;
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.controlSequence2(:, newSeqLength));
            
            % now try to get the input for a future time step, which is not
            % part of the buffered, truncated sequence any more
            [actualInput, actualMode] = this.actuatorUnderTest.step(timestep + 1, []);             
            expectedMode = newSeqLength + 1; % default input is applied
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
        end
        
        %% testChangeControlSequenceLongerSequence
        function testChangeControlSequenceLongerSequence(this)
            newSeqLength = this.sequenceLength + 1;
            timestep = this.timestamp;
            
            [actualInput, actualMode] = this.actuatorUnderTest.step(timestep, [this.dataPacket; this.dataPacket2; this.dataPacket3]);
            % dataPacket2 is buffered as it has smallest delay
            
            this.actuatorUnderTest.changeControlSequenceLength(newSeqLength);
            this.assertNotEmpty(this.actuatorUnderTest.bufferedPacket); % a packet should still be buffered
            this.verifyEqual(actualInput, this.controlSequence2(:, 1));
            this.verifyEqual(actualMode, 1);
            
            % now get the input for the the current timestep
            for j = 1:this.sequenceLength-1
                [actualInput, actualMode] = this.actuatorUnderTest.step(timestep + j, []);
                expectedMode = j + 1;
                expectedInput = this.controlSequence2(:, expectedMode);
                this.verifyEqual(actualMode, expectedMode);
                this.verifyEqual(actualInput, expectedInput);
            end
            
            % try to get the input for a future time step, which is now part of the buffered, extended sequence
            [actualInput, actualMode] = this.actuatorUnderTest.step(timestep + newSeqLength - 1, []);
            expectedMode = newSeqLength; 
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU); % default input is applied as this is appended to the sequence
        end
        
        %% testChangeControlSequenceSameLength
        function testChangeControlSequenceSameLength(this)
            newSeqLength = this.sequenceLength;
            timestep = this.timestamp;
            
            % first, process two packets; those with delay 0 and 2
            % thus, a packet is buffered
            packets = [this.dataPacket2; this.dataPacket3];
            [actualInput, actualMode] = this.actuatorUnderTest.step(timestep, packets);            
        
            this.actuatorUnderTest.changeControlSequenceLength(newSeqLength);
            this.assertNotEmpty(this.actuatorUnderTest.bufferedPacket); % a packet should still be buffered
            this.verifyEqual(actualInput, this.controlSequence2(:, 1));
            this.verifyEqual(actualMode, 1);
            
            % now get the inputs
            for j = 1:newSeqLength - 1
                [actualInput, actualMode] = this.actuatorUnderTest.step(timestep + j, []);
                expectedMode = j + 1;
                expectedInput = this.controlSequence2(:, expectedMode);
                this.verifyEqual(actualMode, expectedMode);
                this.verifyEqual(actualInput, expectedInput);
            end
            % now try to get the input for a future time steps, which is not
            % part of the sequence any more
            [actualInput, actualMode] = this.actuatorUnderTest.step(timestep + newSeqLength, []);
            expectedMode = newSeqLength + 1; % default input is applied
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
            
            [actualInput, actualMode] = this.actuatorUnderTest.step(timestep + newSeqLength + 1, []);
            expectedMode = newSeqLength + 1; % default input is applied
            this.verifyEqual(actualMode, expectedMode);
            this.verifyEqual(actualInput, this.defaultU);
        end
    end
end

