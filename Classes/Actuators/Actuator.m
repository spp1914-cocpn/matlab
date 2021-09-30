classdef Actuator < handle
    % This class simulates an actuator in a packet-based 
    % NCS configuration. The input of the actuator are time-stamped data packets that
    % carry single control inputs. Output is a control input that is applied to the plant. 
    % The actuator is equipped with a buffer where most recent (indicated
    % by the time stamp) control input is stored. The actuator applies the control 
    % input until it gets outdated, as configured by the user, or a newer packet arrives. 
    % If no control input has been received yet, or the currently kept one has become outdated, 
    % a configurable default value is applied.
    %
    % Literature: 
    %  	Florian Rosenthal, Markus Jung, Martina Zitterbart, and Uwe D. Hanebeck,
    %   CoCPN - Towards Flexible and Adaptive Cyber-Physical Systems Through Cooperation,
    %   Proceedings of the 2019 16th IEEE Annual Consumer Communications & Networking Conference,
    %   Las Vegas, Nevada, USA, January 2019.
    %      
    %   Markus Jung, Florian Rosenthal, and Martina Zitterbart,
    %   CoCPN-Sim: An Integrated Simulation Environment for Cyber-Physical Systems,
    %   Proceedings of the 2018 IEEE/ACM Third International Conference on Internet-of-Things Design and Implementation (IoTDI), 
    %   Orlando, FL, USA, April 2018.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2016-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        % dimension of control input; (positive integer)
        dimU(1,1) double {mustBePositive, mustBeInteger, mustBeFinite} = 1;
        % output value in case the buffer runs empty; (column vector)
        % mainly used at the beginning of the simulation, until first
        % packet from controller arrives
        defaultInput(:,1) double {mustBeFinite} = 0;
        % ca be inf
        maxPacketDelay(1,1)  {Validator.validateMaxPacketDelay(maxPacketDelay, 1)} = 1;
    end
    
    properties (SetAccess = protected, GetAccess = public)       
        % buffer to store packet with active (i.e., most recent) control input; (DataPacket)
        bufferedPacket = [];
    end
    
    methods (Access = public)
        %% Actuator
        function this = Actuator(defaultInput, maxPacketDelay)  
            % Class constructor.
            %
            % Parameters:
            %   >> defaultInput (Vector or Scalar)
            %      The default input that is to be applied in case the
            %      buffred control input has become outdated, or none is
            %      present yet.
            %
            %   >> maxPacketDelay (Positive integer or inf)
            %      The maximum delay (in time steps) a received control
            %      input can be have to be applicable before becoming
            %      outdated. That is, if <maxPacketDelay> = n, an input
            %      computed at time k is only applicable at times k, k+1,
            %      ..., k+n. It becomes outdated at time k+n+1, and is no
            %      longer applied even if no fresher one has yet been
            %      received. If <maxPacketDelay> = inf, a received control
            %      input will be held until a fresher one arrives.
            %
            % Returns:
            %   << this (Actuator)
            %      A new Actuator instance.
            %
            this.defaultInput = defaultInput;
            this.dimU = size(this.defaultInput, 1);
            
            this.maxPacketDelay = maxPacketDelay;      
        end
        
        %% step
        function [controlInput, mode, ackPackets] = step(this, timeStep, controllerPackets)
            % Process a given set of packets received from the controller,
            % each of which containing a control input and get
            % the input to be applied at the given time step.
            %
            % Parameters:
            %   >> timeStep (Nonnegative integer)
            %      The time step for which the control input is sought.
            %
            %   >> controllerPackets (Array of DataPackets, might be empty)
            %      An array containing DataPackets with a control inputs.
            %      From the of received control inputs, the one with the smallest delay
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
            %      returned input (i.e., the age of the buffered packet).            
            %
            %   << ackPackets (Array of DataPackets, might be empty (optional))
            %      An array of acknowledgment packets, one for each
            %      received packet from the controller. An empty matrix is
            %      returned in case n packets were received.

            assert(Checks.isNonNegativeScalar(timeStep) && mod(timeStep, 1) == 0, ...
                'Actuator:Step:InvalidTimeStep', ...
                '** Time step must be a nonnegative integer **');    
            
            numPackets = numel(controllerPackets);
            if numPackets ~= 0
                assert(Checks.isClass(controllerPackets, 'DataPacket', numPackets), ...
                    'Actuator:Step:InvalidControllerPackets', ...
                    '** Input parameter (received control input packets) must consist of %d DataPacket(s) **', numPackets);
                
                [controlInputs{1:numPackets}] = controllerPackets(:).payload;
                assert(sum(cell2mat(cellfun(@(x) Checks.isVec(x, this.dimU), ...
                        controlInputs, 'UniformOutput', false))) == numPackets, ...
                    'Actuator:Step:InvalidControlInput', ...
                     '** Each control input must be given as real-valued %d-dimensional vector **', this.dimU);
                
                [delays{1:numPackets}] = controllerPackets(:).packetDelay;                                
                % get the newest input, i.e., that one with smallest delay (should be unique)
                % packet rejection logic: all others are discarded 
                [~, idx] = min(cell2mat(delays));
                this.bufferedPacket = controllerPackets(idx);
            end
            
            if ~isempty(this.bufferedPacket)
                mode = min(timeStep - this.bufferedPacket.timeStamp + 1,this.maxPacketDelay + 1); % first mode has index 1!
                if mode == this.maxPacketDelay + 1
                    % packet is too old, so use default input instead
                    controlInput = this.defaultInput;
                else                    
                    controlInput = this.bufferedPacket.payload;
                end
            else
                % buffer is empty, no packet has yet arrived
                mode = this.maxPacketDelay + 1;
                controlInput = this.defaultInput;
            end        
            
            %ackPackets
            % prepare ACKs for all the received packets
             if nargout == 3
                ackPackets = [];                
             end 
        end
        
        %% reset
        function reset(this)
            % Reset the actuator by removing the currently buffered
            % control sequence.
            %
            this.bufferedPacket = [];
        end 
    end
end

