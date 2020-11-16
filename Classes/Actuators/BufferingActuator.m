classdef BufferingActuator < Actuator
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
    %  	Florian Rosenthal, Benjamin Noack, and Uwe D. Hanebeck,
    %   State Estimation in Networked Control Systems with Delayed and Lossy Acknowledgments,
    %   Multisensor Fusion and Integration in the Wake of Big Data, Deep Learning and Cyber Physical System,
    %   Lecture Notes in Electrical Engineering, Volume 501,
    %   Springer, Cham, 2018.
    %
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2016-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    
    properties (SetAccess = private, GetAccess = public)
        % length of applicable control input sequences
        % one sequence per received packet; (positive integer)
        controlSequenceLength(1,1) double = 1;
    end
    
    properties (Access = private)
        % the history of previous true modes, used by ACKs
        modeHistory = [];
    end
    
    methods (Access = public)
        %% BufferingActuator
        function this = BufferingActuator(controlSequenceLength,  defaultU)
            % Class constructor.
            %
            % Parameters:
            %   >> controlSequenceLength (Positive integer)
            %      The length of a valid control input sequence.
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
            this = this@Actuator(defaultU, controlSequenceLength - 1);
            
            this.controlSequenceLength = controlSequenceLength;  
            
            % initially, the buffer is empty and thus: theta_0 = N+1
            % where N+1 is the control sequence length
            % initialize the mode history, but remember that we enumerate
            % the modes from [1, N+2] instead of [0, N+1] as in the
            % publications
            this.modeHistory = zeros(controlSequenceLength, 1) + controlSequenceLength + 1;
        end

        %% step
        function [controlInput, mode, ackPackets] = step(this, timeStep, controllerPackets)
            % Process a given set of packets received from the controller,
            % each of which containing a control input sequence and get
            % the input to be applied at the given time step.
            %
            % Parameters:
            %   >> timeStep (Nonnegative integer)
            %      The time step for which the control input is sought.
            %
            %   >> controllerPackets (Array of DataPackets, might be empty)
            %      An array containing DataPackets with a control input sequence.
            %      From the of received control sequences, the one with the smallest delay
            %      (i.e., the one which was sent most recently) is buffered if it is more
            %      recent than the one currently stored,
            %      and all others are discarded (past packets rejection logic).
            %
            % Returns:
            %   << controlInput (Column vector or Scalar)
            %      The control input for the given time step, if
            %      present. If not, then the default input is returned.
            %
            %   << mode (Positive integer)
            %      The mode (theta_k) of the underlying MJLS that corresponds to the
            %      returned input (i.e., the age of the buffered sequence).
            %      In constrast to the publications, here we always have
            %      theta_k in [1,2, ...,N+2] instead of [0, 1,...,N+1]
            %      with N+1 the control sequence length.
            %
            %   << ackPackets (Array of DataPackets, might be empty (optional))
            %      An array of acknowledgment packets, one for each
            %      received packet from the controller. An empty matrix is
            %      returned in case n packets were received.

            assert(Checks.isNonNegativeScalar(timeStep) && mod(timeStep, 1) == 0, ...
                'BufferingActuator:Step:InvalidTimeStep', ...
                '** Time step must be a nonnegative integer **');    
            
            numSeq = numel(controllerPackets);
            if numSeq ~= 0
                assert(Checks.isClass(controllerPackets, 'DataPacket', numSeq), ...
                    'BufferingActuator:Step:InvalidControllerPackets', ...
                    '** Input parameter (received control input packets) must consist of %d DataPacket(s) **', numSeq);
                
                [inputSequences{1:numSeq}] = controllerPackets(:).payload;
                assert(sum(cell2mat(cellfun(@(x) Checks.isMat(x, this.dimU, this.controlSequenceLength), ...
                        inputSequences, 'UniformOutput', false))) == numSeq, ...
                    'BufferingActuator:Step:InvalidControlSequences', ...
                     '** Each control input sequence must be given as real-valued a %d-by-%d matrix **', ...
                     this.dimU, this.controlSequenceLength);
                
                [delays{1:numSeq}] = controllerPackets(:).packetDelay;
                                
                % get the newest sequence, i.e., that one with smallest delay (should be unique)
                % packet rejection logic: most recent packet will be stored, all others are discarded 
                [minDelay, idx] = min(cell2mat(delays));                
                if minDelay <= this.maxPacketDelay ...
                    && (isempty(this.bufferedPacket) || controllerPackets(idx).isNewerThan(this.bufferedPacket))
                    this.bufferedPacket = controllerPackets(idx);
                end
            end
            
            if ~isempty(this.bufferedPacket)...
                    && this.bufferedPacket.timeStamp + this.controlSequenceLength > timeStep
                mode = timeStep - this.bufferedPacket.timeStamp + 1; % first mode has index 1!
                controlInput = this.bufferedPacket.payload(:, mode);
            else
                % buffer is empty or content is obsolete
                controlInput = this.defaultInput;
                mode = this.controlSequenceLength + 1; % N+2 (last mode)
            end
            
            %ackPackets
            % prepare ACKs for all the received packets
             if nargout == 3
                ackPackets = [];
                for i=1:numSeq
                    info.timeStep = controllerPackets(i).timeStamp; % the time step the sequence was sent
                    info.tau = min(controllerPackets(i).packetDelay, this.maxPacketDelay + 1); % delay of the packet, or considered a loss
                    if info.timeStep == timeStep % packet delay was zero
                        info.theta = mode; % the current mode (starting from 1)
                    else
                        info.theta = this.modeHistory(info.tau); % extract a mode from the past
                    end
                    ackPackets = [ackPackets; controllerPackets(i).createAckForDataPacket(timeStep, info)]; %#ok
                end
            end
            % update the history
            this.modeHistory = [mode; this.modeHistory(1:end-1)];
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
