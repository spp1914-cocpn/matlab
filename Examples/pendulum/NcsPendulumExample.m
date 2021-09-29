classdef NcsPendulumExample < handle
    % This class is an example how to simulate a complete networked control
    % system (NCS) (in this case an inverted pendulum shall be stabilized) in Matlab
    % by using the API-functions provided by libncs_matlab.
    properties
        ncs;
        caNetwork;
        scNetwork;
        acNetwork;
        numTimesteps;
        samplingInterval;
        simTimeSeconds;
    end
    
    methods (Access = public)
        %% NcsPendulumExample
        function this = NcsPendulumExample(config, simTime)
            this.ncs = ncs_initialize(ConvertToPicoseconds(simTime), 'Pendulum NCS', config, 'libncs_matlab/matlab/config/inverted_pendulum_short.mat');
            
            % define the network for measurements, control inputs and ACKs
            % from the actuator, respectively
            this.caNetwork = CommunicationNetwork(config.numTimeSteps, config.maxControlSequenceDelay);
            this.scNetwork = CommunicationNetwork(config.numTimeSteps,config.maxMeasDelay);
            this.acNetwork = CommunicationNetwork(config.numTimeSteps,config.maxMeasDelay + 1);
            
            this.numTimesteps = config.numTimeSteps; 
            this.samplingInterval = config.samplingInterval;
            this.simTimeSeconds = simTime;
        end
                      
        %% simulate
        function [statsPerStep, stats, costs] = simulate(this, caPacketDelays, scPacketDelays, acPacketDelays)
            statsPerStep = [];
            for k = 1:this.numTimesteps
                timestamp = uint64(ConvertToPicoseconds(this.samplingInterval * k));
                ncs_doPlantStep(this.ncs, timestamp);
                %receive all packets first
                % receive measurements and mode observations
                [~, arrivedAcPackets] = this.acNetwork.receivePackets(k);
                [~, arrivedScPackets] = this.scNetwork.receivePackets(k);
                [~, arrivedCaPackets] = this.caNetwork.receivePackets(k);  

                for packet = arrivedAcPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                for packet = arrivedScPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                for packet = arrivedCaPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                
                [pktsOut, qocOut, res] = ncs_doLoopStep(this.ncs, timestamp, []);
                statsPerStep = [statsPerStep res];
                for i = 1:numel(pktsOut)
                    packet = pktsOut{i};
                    switch packet.destinationAddress
                        case 1                           
                            this.caNetwork.sendPacket(packet, caPacketDelays(k), k);                       
                        case 2
                            switch packet.sourceAddress
                                case 1
                                    % expected to be an ACK packet, so check
                                    assert(packet.isAck, ...
                                        'ncs_doLoopStep:InvalidACKPacket', ...
                                        '** Packet from 1 (actuator) to 2 (controller) should be an ACK **');

                                    this.acNetwork.sendPacket(packet, acPacketDelays(k), k);
                                case 3
                                    this.scNetwork.sendPacket(packet, scPacketDelays(k), k);
                                otherwise                                    
                            end
                        case 3
                             % so far, do nothing
                             %csPackets = [csPackets, packet];
                        otherwise
                    end
                end
            end
            [costs, stats] = ncs_finalize(this.ncs);
        end
        
        %% simulateWithLosses
        function [statsPerStep, stats, costs, numDroppedPackets] = simulateWithLosses(this, caPacketDelays, scPacketDelays, acPacketDelays, maxNumBufferedPackets)
            statsPerStep = [];
            numDroppedPackets = 0;
            for k = 1:this.numTimesteps
                timestamp = uint64(ConvertToPicoseconds(this.samplingInterval * k));
                ncs_doPlantStep(this.ncs, timestamp);
                %receive all packets first
                % receive measurements and mode observations
                [~, arrivedAcPackets] = this.acNetwork.receivePackets(k);
                [~, arrivedScPackets] = this.scNetwork.receivePackets(k);
                [~, arrivedCaPackets] = this.caNetwork.receivePackets(k);  
                
                numBufferedPackets = this.caNetwork.getNumBufferedPackets(k);
                
                for packet = arrivedAcPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                for packet = arrivedScPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                for packet = arrivedCaPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                
                [pktsOut, qocOut, res] = ncs_doLoopStep(this.ncs, timestamp, []);
                statsPerStep = [statsPerStep res];
                for i = 1:numel(pktsOut)
                    packet = pktsOut{i};
                    switch packet.destinationAddress
                        case 1
                            if numBufferedPackets < maxNumBufferedPackets%maxPacketRate
                                this.caNetwork.sendPacket(packet, caPacketDelays(k), k);                                
                            else
                                numDroppedPackets = numDroppedPackets + 1;
                            end                            
                        case 2
                            switch packet.sourceAddress
                                case 1
                                    % expected to be an ACK packet, so check
                                    assert(packet.isAck, ...
                                        'ncs_doLoopStep:InvalidACKPacket', ...
                                        '** Packet from 1 (actuator) to 2 (controller) should be an ACK **');

                                    this.acNetwork.sendPacket(packet, acPacketDelays(k), k);
                                case 3
                                    this.scNetwork.sendPacket(packet, scPacketDelays(k), k);
                                otherwise                                    
                            end
                        case 3
                             % so far, do nothing
                             %csPackets = [csPackets, packet];
                        otherwise
                    end
                end
            end
            [costs, stats] = ncs_finalize(this.ncs);
        end

        
        
        %% simulateWithNetworkFeedback
        function [statsPerStep, stats, costs] = simulateWithNetworkFeedback(this, caPacketDelays, scPacketDelays, acPacketDelays, networkFeedback)
            statsPerStep = [];

            for k = 1:this.numTimesteps
                timestamp = uint64(ConvertToPicoseconds(this.samplingInterval * k));
                ncs_doPlantStep(this.ncs, timestamp);
                %receive all packets first
                % receive measurements and mode observations
                [~, arrivedAcPackets] = this.acNetwork.receivePackets(k);
                [~, arrivedScPackets] = this.scNetwork.receivePackets(k);
                [~, arrivedCaPackets] = this.caNetwork.receivePackets(k);                
               
                for packet = arrivedAcPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                for packet = arrivedScPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                for packet = arrivedCaPackets
                    ncs_doHandlePacket(this.ncs, timestamp, packet);
                end
                
                if ~isempty(networkFeedback{k})
                    config.caDelayProbs = networkFeedback{k};
                else
                    config = [];
                end
                
                [pktsOut, qocOut, res] = ncs_doLoopStep(this.ncs, timestamp, config);
                statsPerStep = [statsPerStep res];
                for i = 1:numel(pktsOut)
                    packet = pktsOut{i};
                    switch packet.destinationAddress
                        case 1
                            
                                this.caNetwork.sendPacket(packet, caPacketDelays(k), k);                                
       
                        case 2
                            switch packet.sourceAddress
                                case 1
                                    % expected to be an ACK packet, so check
                                    assert(packet.isAck, ...
                                        'ncs_doLoopStep:InvalidACKPacket', ...
                                        '** Packet from 1 (actuator) to 2 (controller) should be an ACK **');

                                    this.acNetwork.sendPacket(packet, acPacketDelays(k), k);
                                case 3
                                    this.scNetwork.sendPacket(packet, scPacketDelays(k), k);
                                otherwise                                    
                            end
                        case 3
                             % so far, do nothing
                             %csPackets = [csPackets, packet];
                        otherwise
                    end
                end
            end
            [costs, stats] = ncs_finalize(this.ncs);
        end
    end
    
end

