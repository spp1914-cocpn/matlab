classdef BufferingActuator < handle
    % This class simulates an actuator in a packet-based 
    % NCS configuration. The input of the actuator are data packets that
    % consist of time-stamped sequences of control inputs. Output is a 
    % single control input that is applied to the plant. The actuator
    % is equipped with a buffer where the control input sequence with the
    % most recent time stamp is stored. The actuator applies the control 
    % input of that buffered sequence that corresponds to the current time 
    % step. If there is no control input that corresponds to the current 
    % time step, a configurable default value is applied.
    %
    % This implementation is based on the original one by Jörg Fischer.
    %
    % Literature: 
    %   Jörg Fischer, Achim Hekler and Uwe D. Hanebeck,
    %   State Estimation in Networked Control Systems,
    %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
    %   Singapore, July 2012.
    %
    %   Jörg Fischer, Achim Hekler, Maxim Dolgov and Uwe D. Hanebeck,
    %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject
    %   to Random Transmission Delays and Packet Losses,
    %   Proceedings of the 2013 American Control Conference (ACC 2013), 
    %   Washington D. C., USA, June 2013.
    %
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2016-2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (SetAccess = immutable, GetAccess = public)
        % maximum allowed packet delay (packets with control input sequences)
        maxPacketDelay = -1;
        % dimension of control input; (positive integer)
        dimU = -1;
        % output value in case the buffer runs empty; (column vector)
        defaultInput= 0;
    end
    
    properties (SetAccess = private, GetAccess = public)
        % length of applicable control input sequences
        % one sequence per received packet; (positive integer)
        controlSequenceLength = -1;
        % buffer to store packet with active (i.e., most recent) control input sequence; (DataPacket)
        bufferedPacket = [];
    end

    methods (Access = public)
        %% BufferingActuator
        function this = BufferingActuator(controlSequenceLength, maxPacketDelay, defaultU)
            % Class constructor.
            %
            % Parameters:
            %   >> controlSequenceLength (Positive integer)
            %      The length of a valid control input sequence.
            %
            %   >> maxPacketDelay (Nonnegative integer of Inf)
            %      The maximum allowed delay experienced by received control sequences, older ones will be discarded.
            %      May be Inf (nothing will be discarded) or 0 (any delayed sequence will be discarded).
            %
            %   >> defaultU (Vector or Scalar)
            %      The default input that is to be applied in case the
            %      internal buffer runs empty.
            %          
            % Returns:
            %   << this (BufferingActuator)
            %      A new BufferingActuator instance.
            %
            Validator.validateSequenceLength(controlSequenceLength);
            this.controlSequenceLength = controlSequenceLength;
           
            assert(Checks.isNonNegativeScalar(maxPacketDelay) ...
                    && mod(maxPacketDelay, 1) == 0 || isinf(maxPacketDelay), ...
                'BufferingActuator:InvalidMaxPacketDelay', ...
                '** Input parameter <maxPacketDelay> (maximum allowed packet delay) must be a nonnegative integer or inf **');
           
            this.maxPacketDelay = maxPacketDelay;
         
            assert(Checks.isVec(defaultU) && all(isfinite(defaultU)), ...
                'BufferingActuator:InvalidDefaultU', ...
                '** Input parameters <defaultU> (default input) must be a real-valued vector or a real number  **');

            this.defaultInput = defaultU(:);
            this.dimU = size(this.defaultInput, 1);
        end

        function [controlInput, mode] = getCurrentInput(this, timeStep)
            % Get the input to be applied at the given time step.
            %
            % Parameters:
            %   >> timeStep (Nonnegative integer)
            %      The time step for which the control input is sought.
            %
            % Returns:
            %   << controlInput (Column vector or Scalar)
            %      The buffered control input for the given time step, if
            %      present. If not, then the default input is returned.
            %   << mode (Positive integer)
            %      The mode of the underlying MJLS that corresponds the the
            %      returned input (i.e., the age of the buffered sequence).
            %
            assert(Checks.isNonNegativeScalar(timeStep) && mod(timeStep, 1) == 0, ...
                'BufferingActuator:GetCurrentInput:InvalidTimeStep', ...
                '** Time step must be a nonnegative integer **');    
            
            if ~isempty(this.bufferedPacket)...
                    && this.bufferedPacket.timeStamp + this.controlSequenceLength > timeStep
                mode = timeStep - this.bufferedPacket.timeStamp + 1;
                controlInput = this.bufferedPacket.payload(:, mode);
            else
                % buffer is empty or content is obsolete
                controlInput = this.defaultInput;
                mode = this.controlSequenceLength + 1;
            end
        end
  
        function ackPacket = processControllerPackets(this, controllerPackets)
            % Process a given set of packets received from the controller,
            % each of which containing a control input sequence.
            % From this set, that one with the smallest delay (i.e., the
            % one which was sent most recently) is buffered if it is more
            % recent than the one currently stored, all others are
            % discarded (past packets rejection logic).
            %
            % Parameters:
            %   >> controllerPackets (Array of DataPackets)
            %      An array containing DataPackets with a control input sequence.
            %
            % Returns:
            %   << ackPacket (Empty matrix or DataPacket)
            %      The ACK for the DataPacket within the given set that has
            %      become active. An empty matrix is returned in case non
            %      became active.
            %
            if nargout == 1
                ackPacket = [];
            end
            numSeq = numel(controllerPackets);
            if numSeq ~= 0
                assert(Checks.isClass(controllerPackets, 'DataPacket', numSeq), ...
                    'BufferingActuator:ProcessControllerPackets:InvalidControllerPackets', ...
                    '** Input parameter (received control input packets) must consist of %d DataPacket(s) **', numSeq);
                
                [inputSequences{1:numSeq}] = controllerPackets(:).payload;
                assert(sum(cell2mat(cellfun(@(x) Checks.isMat(x, this.dimU, this.controlSequenceLength), ...
                        inputSequences, 'UniformOutput', false))) == numSeq, ...
                    'BufferingActuator:ProcessControllerPackets:InvalidControlSequences', ...
                     '** Each control input sequence must be given as real-valued a %d-by-%d matrix **', ...
                     this.dimU, this.controlSequenceLength);
                
                [delays{1:numSeq}] = controllerPackets(:).packetDelay;
                % get the newest sequence, i.e., that one with smallest delay
                % (should be unique)
                [minDelay, idx] = min(cell2mat(delays));
                if minDelay <= this.maxPacketDelay && this.controlSequenceLength - minDelay > 0 ...
                    && (isempty(this.bufferedPacket) || controllerPackets(idx).isNewerThan(this.bufferedPacket))
                    this.bufferedPacket = controllerPackets(idx);
                     % packet rejection logic: only ACK packet if it becomes
                     % active, all others are discarded
                     if nargout == 1
                         % compute the time stamp
                         k = this.bufferedPacket.timeStamp + minDelay;
                         ackPacket = this.bufferedPacket.createAckForDataPacket(k);
                     end
                end
            end
        end
        
        %% reset
        function reset(this)
            % Reset the actuator by removing the currently buffered
            % control sequence.
            %
            this.bufferedPacket = [];
        end 
        
        %% changeControlSequenceLength
        function changeControlSequenceLength(this, newSequenceLength)
            % Change the length of the control sequence expected to be
            % received from the controller.
            % If an old sequence is currently buffered, its size is adapted
            % as follows. In case the new sequence length is greather than
            % the old one, the specified default input is appended. In case
            % the new sequence length is less than the old one, the
            % last, now superfluous entries are discarded.
            %
            % Parameters:
            %   >> newSequenceLength (Positive integer)
            %      The new sequence length to used.
            
            Validator.validateSequenceLength(newSequenceLength);
            if ~isempty(this.bufferedPacket) && newSequenceLength ~= this.controlSequenceLength
                % fit the size of the buffered sequence
                % that is, trim buffered sequence or append default input
                % copy of DataPacket is required
                        
                srcAddr = this.bufferedPacket.sourceAddress;
                dstAddr = this.bufferedPacket.destinationAddress;
                delay = this.bufferedPacket.packetDelay;
                
                bufferedSequence = [this.bufferedPacket.payload(:, 1:min(this.controlSequenceLength, newSequenceLength)) ...
                    repmat(this.defaultInput, 1, newSequenceLength - this.controlSequenceLength)];
                
                % don't forget to set properties of packet
                this.bufferedPacket = DataPacket(bufferedSequence, this.bufferedPacket.timeStamp, this.bufferedPacket.id);
                this.bufferedPacket.packetDelay = delay;
                this.bufferedPacket.sourceAddress = srcAddr;
                this.bufferedPacket.destinationAddress = dstAddr;
            end
            this.controlSequenceLength = newSequenceLength;
        end
    end
end
