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
    % AUTHOR:       Jörg Fischer
    % LAST UPDATE:  Jörg Fischer, 06.11.2013
    %               Florian Rosenthal, 06.12.2016
    
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
        %% input properties
        % length of applicable control input sequences
        % one sequence per received packet; (positive integer)
        controlSequenceLength = -1;
        % maximum allowed packet delay (packets with control input sequences)
        maxPacketDelay = -1;
        % dimension of control input; (positive integer)
        dimU = -1;
        % output value in case the buffer runs empty; (column vector)
        defaultInput= 0;
    end % properties

    properties (SetAccess = private, GetAccess = public)
        %% derived properties
        % buffer to store packet with active (i.e., most recent) control input sequence; (DataPacket)
        bufferedPacket = [];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %% BufferingActuator
        function s = BufferingActuator(controlSequenceLength, maxPacketDelay, defaultU)
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
            s.controlSequenceLength = controlSequenceLength;
           
            if Checks.isNonNegativeScalar(maxPacketDelay) && ...
                    (mod(maxPacketDelay, 1) == 0 || isinf(maxPacketDelay))
                s.maxPacketDelay = maxPacketDelay;
            else
                error('BufferingActuator:InvalidMaxPacketDelay', ...
                    ['** Input parameter <maxPacketDelay> (maximum allowed ',...
                       'packet delay) must be a nonnegative integer or inf **']);
            end

            if ~Checks.isVec(defaultU) || any(~isfinite(defaultU))
                error('BufferingActuator:InvalidDefaultU', ...
                    ['** Input parameters <defaultU> (default input) must ',...
                     'be a real-valued vector or a real number  **'])
            end
            s.defaultInput = defaultU(:);
            s.dimU = size(s.defaultInput, 1);
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
            if ~Checks.isNonNegativeScalar(timeStep) || mod(timeStep, 1) ~= 0
                error ('BufferingActuator:GetCurrentInput:InvalidTimeStep', ...
                    '** Time step must be a nonnegative integer **');    
            end
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
                if ~Checks.isClass(controllerPackets, 'DataPacket', numSeq)
                        error('BufferingActuator:ProcessControllerPackets:InvalidControllerPackets', ...
                            ['** Input parameter (received control input packets) must', ...
                             ' consist of %d DataPacket(s) **'], numSeq);
                end
                [inputSequences{1:numSeq}] = controllerPackets(:).payload;
                if sum(cell2mat(cellfun(@(x) Checks.isMat(x, this.dimU, this.controlSequenceLength), ...
                        inputSequences, 'UniformOutput', false))) ~= numSeq
                     error('BufferingActuator:ProcessControllerPackets:InvalidControlSequences', ...
                         '** Each control input sequence must be given as real-valued a %d-by-%d matrix **', ...
                         this.dimU, this.controlSequenceLength);
                end

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
            % Resets the actuator by removing the currently buffered
            % control sequence.
            %
            this.bufferedPacket = [];
        end 
    end
end