classdef RecedingHorizonUdpLikeController < SequenceBasedController
    % Implementation of a receding horizon linear sequence-based controller for
    % NCS with UDP-like networks connecting the controller and actuator, 
    % and sensor and controller, respectively, where application layer ACKs
    % are sent out by the actuator upon reception of applicable control
    % inputs.
    %
    % Literature: 
    %   Florian Rosenthal, Maxim Dolgov, and Uwe D. Hanebeck,
    %   Sequence-Based Receding Horizon Control over Networks With Delays and Data Losses,
    %   Proceedings of the 2019 American Control Conference (ACC 2019),
    %   Philadelphia, Pennsylvania, USA, July 2019.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
     
    properties (GetAccess = public, Constant)
        defaultMaxNumIterations = 1000;%1000;%5000;
    end
         
    properties (Access = private, Constant)
        convergenceDiff = 1e-20;
    end
    
    properties (SetAccess = immutable, GetAccess = private)
        dimMeas;
        dimAugmentedMeas;
        measBufferLength;
        Q;
        R;
        
        dimEta;
        %etaState;
        
        dimAugState; % dimension of the controller's augmented state       
        
        augA;
        augB;
          
        JRJ; % only nonzero for the first mode
                
        augmentedNoiseCov;
        augW;
        augV;
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        horizonLength double {Validator.validateHorizonLength(horizonLength)} = 1;
                        
        useMexImplementation(1,1) logical = true; 
        % by default, we use the C++ (mex) implementation for computation of controller gains
        % this is usually faster, but can produce slightly different results
    end   
        
    properties (SetAccess = private, GetAccess = public)
        lastNumIterations = 0;        
    end
    
    properties (Access = public)
        % max num number of iterations carried out during computation of
        % gains
        % can be adapted (decreased) to reduce computation time
        maxNumIterations(1,1) double = RecedingHorizonUdpLikeController.defaultMaxNumIterations;
    end
    
    properties (Access = private)
        S;
        augmentedMeasMatrix; % sparse!
        
        augQ;
        terminalAugQ; % the same for all modes
        
        transitionMatrixCa;
        transitionMatrixSc_part; % only a third of the matrix
        
        lastNumUsedMeas = 0;
        lastNumDiscardedMeas = 0;
        
        modeSc; % this is observed (stored as integer)
        modeCa; % this one is only partially observed (stored as probability (column) vec)
        lastModeObservationDelay = inf; % store the time delay to the last mode observation
        
        L;
        M;
        K;
        augState; % the controller's augmented plant state
        augStateSecondMoment; % internal estimate of second moment per ca mode
        augStateCov; % internal estimate of covariance per ca mode
    end   
      
    properties (GetAccess = public, Dependent)
        maxMeasurementDelay;
    end
     
    methods
        function maxMeasDelay = get.maxMeasurementDelay(this)
            maxMeasDelay = this.measBufferLength - 1;
        end
    end
    
    methods (Access = public)
        %% RecedingHorizonUdpLikeController
        function this = RecedingHorizonUdpLikeController(A, B, C, Q, R, caModeTransitionMatrix, scDelayProb, ...
                sequenceLength, maxMeasDelay, W, V, horizonLength, x0, x0Cov, useMexImplementation)
            % Class constructor.
            %
            % Parameters:
            %   >> A (Square Matrix)
            %      The system matrix of the plant.
            %
            %   >> B (Matrix)
            %      The input matrix of the plant.
            %
            %   >> C (Matrix)
            %      The measurement matrix of the sensor.
            %
            %   >> Q (Positive semi-definite matrix)
            %      The state weighting matrix Q_k in the controller's underlying cost function.
            %      For the terminal state x_K, where K is the horizon
            %      length, the corresponding weighting matrix Q_K is chosen as the
            %      stabilizing solution of the DARE associated with A, B,
            %      Q, R. If this solution does not exist, Q is used instead.
            %
            %   >> R (Positive definite matrix)
            %      The input weighting matrix in the controller's underlying cost function.
            %
            %   >> caModeTransitionMatrix (Stochastic matrix, i.e. a square matrix with nonnegative entries whose rows sum to 1)
            %      The transition matrix of the mode theta_k of the augmented dynamics.
            %
            %   >> scDelayProb (Nonnegative vector)
            %      The vector describing the delay distribution of the
            %      SC-network.
            %
            %   >> sequenceLength (Positive integer)
            %      The length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %   >> maxMeasDelay (Nonnegative integer)
            %      The maximum delay a measurement may experience before it
            %      is discarded by the controller.
            %
            %   >> W (Square matrix)
            %      The covariance matrix of the plant noise.
            %
            %   >> V (Square matrix)
            %      The covariance matrix of the measurement noise.
            %
            %   >> horizonLength (Positive integer)
            %      The horizon length considered by the controller for optimization.
            %
            %   >> x0 (Vector)
            %      The controller's initial estimate of the plant state.
            %
            %   >> x0Cov (Positive definite matrix)
            %      The covariance matrix denoting the uncertainty of the controller's initial estimate.
            %
            %   >> useMexImplementation (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for the computation of the controller
            %      gains which is usually considerably faster than the Matlab
            %      implementation. 
            %      If left out, the default value true is used.
            %
            % Returns:
            %   << this (RecedingHorizonUdpLikeController)
            %      A new RecedingHorizonUdpLikeController instance.
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            this = this@SequenceBasedController(dimX, dimU, sequenceLength, false);
            
            % Q, R
            Validator.validateCostMatrices(Q, R, dimX, dimU);
            this.Q = Q;
            this.R = R;
                         
            numCaModes = sequenceLength + 1;
            Validator.validateTransitionMatrix(caModeTransitionMatrix, numCaModes);
            this.transitionMatrixCa = caModeTransitionMatrix;
              
            Validator.validateDiscreteProbabilityDistribution(scDelayProb);
            
            Validator.validateMeasurementMatrix(C, dimX);
            this.dimMeas = size(C, 1);
            assert(Checks.isNonNegativeScalar(maxMeasDelay) && mod(maxMeasDelay, 1) == 0, ...
                'RecedingHorizonUdpLikeController:InvalidMaxMeasDelay', ...
                ['** Input parameter <maxMeasDelay> (maximum measurement',...
                'delay (M)) must be a nonnegative integer **']);              
                         
            % check the noise covs
            Validator.validateSysNoiseCovarianceMatrix(W, dimX);
            Validator.validateMeasNoiseCovarianceMatrix(V, this.dimMeas);       

            % check the states
            assert(Checks.isVec(x0, dimX) && all(isfinite(x0)), ...
                'RecedingHorizonUdpLikeController:InvalidX0', ...
                '** Input parameter <x0> (initial estimate of plant state) must be a %d-dimensional vector **', dimX);
            assert(Checks.isCov(x0Cov, dimX), ...
                'RecedingHorizonUdpLikeController:InvalidX0Cov', ...
                '** Input parameter <x0Cov> (initial covariance) must be positive definite and %d-by%d-dimensional **', dimX, dimX);
            
            if nargin > 14
                this.useMexImplementation = useMexImplementation;
            end

            this.measBufferLength = maxMeasDelay + 1; % maxMeasDelay + 1 is buffer length
            this.dimAugmentedMeas = this.dimMeas * this.measBufferLength;
            this.dimEta = dimU * (sequenceLength * (sequenceLength - 1) / 2);
            % augmented state consists of x_k, x_k-1, ..., x_k-M and eta_k
            this.dimAugState = dimX * this.measBufferLength + this.dimEta;
             
            % and the augmented noise covs
            this.augW = blkdiag(W, zeros(maxMeasDelay*dimX+this.dimEta));
            this.augV = kron(eye(this.measBufferLength), V);           

            this.horizonLength = horizonLength;          
                        
            numModesSc = 3 ^ (this.measBufferLength);
            modesScDecimal = 0:3^(this.measBufferLength)-1; % decimal numbers        
            modesScBase3 = zeros(numModesSc, this.measBufferLength);
            for i = 0:maxMeasDelay
                modesScBase3(:, i+1) = mod(floor(modesScDecimal ./ (3 ^ i)), 3);
            end            
            this.computeAndSetScTransitionMatrix(scDelayProb, modesScBase3, modesScDecimal);
            this.initAugmentedMeasMatrices(modesScBase3, numModesSc, modesScDecimal, C);
            
            [F, G, H, ~] = Utility.createAugmentedPlantModel(sequenceLength, A, B);
            J = [speye(dimU) zeros(dimU, dimU * (sequenceLength - 1))];
            this.JRJ = J' * R * J;
            
            D = eye(dimX); 
           
            % A_tilde in the paper
            this.augA = repmat(blkdiag(A, diag(ones(dimX * (maxMeasDelay - 1), 1), -dimX), F), 1, 1, numCaModes);
            this.augA(dimX + 1:dimX * 2, 1:dimX, :) = repmat(D, 1, 1, numCaModes); % below A
            % H is empty for first and last mode
            % compute B * H(:, :,  i) for all modes (except first and last)
            % in one shot
            this.augA(1:dimX, end - this.dimEta + 1:end, 2:sequenceLength) = mtimesx(B, H(:, :, 2:sequenceLength));
            % B_tilde in the paper
            this.augB = repmat([zeros(this.dimAugState - size(G, 1), size(G, 2)); G], 1, 1, numCaModes);
            % augB is only different for first mode (multiply B by J)
            % B*J is zero for all modes but first
            %
            this.augB(1:dimX, :, 1) = B * J; 

            this.initAugmentedCostMatrices(numCaModes, H, A, B);
            this.initControllerGains();
            this.modeSc = 0;   % no measurements available initially
          
            this.setControllerPlantStateInternal(x0(:), x0Cov);
        end             
        
        %% setInitialControllerGains
        function setInitialControllerGains(this, initialM, initialK, initialL)
            % Set the controller gains M, K, and L to be
            % employed at the beginning of the time of operation.
            %
            % This function is intended to be used to set initial first
            % guess gains that serve as starting point for the first
            % optimization, or for testing purposes.
            %
            % Parameters:
            %   >> initialM (3D matrix)
            %      N initial gains M_0, M_1, ... M_{N-1} slice-wise arranged, where N is the optimization horizon.
            %
            %   >> initialK (3D matrix)
            %      N initial gains K_0, K_1, ... K_{N-1} slice-wise arranged, where N is the optimization horizon.
            %
            %   >> initialL (3D matrix)
            %      N initial gains L_0, L_1, ... L_{N-1} slice-wise arranged, where N is the optimization horizon.
            
            assert(Checks.isMat3D(initialM, this.dimAugState, this.dimAugState, this.horizonLength), ...
                'RecedingHorizonUdpLikeController:SetInitialControllerGains:InvalidM', ...
                '** <initialM> must be a 3D matrix of dimension %d-by-%d-by-%d **', ...
                this.dimAugState, this.dimAugState, this.horizonLength); 
            assert(Checks.isMat3D(initialK, this.dimAugState, this.dimAugmentedMeas, this.horizonLength), ...
                'RecedingHorizonUdpLikeController:SetInitialControllerGains:InvalidK', ...
                '** <initialK> must be a 3D matrix of dimension %d-by-%d-by-%d **', ...
                this.dimAugState, this.dimAugmentedMeas, this.horizonLength); 
             assert(Checks.isMat3D(initialL, this.dimPlantInput * this.sequenceLength, this.dimAugState, this.horizonLength), ...
                'RecedingHorizonUdpLikeController:SetInitialControllerGains:InvalidL', ...
                '** <initialL> must be a 3D matrix of dimension %d-by-%d-by-%d **', ...
                this.dimPlantInput * this.sequenceLength, this.dimAugState, this.horizonLength); 
            
            this.M = initialM;
            this.K = initialK;
            this.L = initialL;
        end
                      
        %% setControllerPlantState
        function setControllerPlantState(this, state)
            % Set the estimate of the plant state.
            %
            % This function is mainly used to set an initial estimate as
            % the controller's dynamics is used for updating the estimate of
            % the plant state.
            %
            % Parameters:
            %   >> state (Subclass of Distribution)
            %      The new system state.
            
            assert(Checks.isClass(state, 'Distribution') && state.getDim() == this.dimPlantState, ...
                'RecedingHorizonUdpLikeController:SetControllerPlantState:InvalidState', ...
                '** <state> must be a %d-dimensional Distribution **', this.dimPlantState);
            
            [x0, x0Cov] = state.getMeanAndCov();
            
            this.setControllerPlantStateInternal(x0, x0Cov);  
        end
        
        %% setEtaState
        function setEtaState(this, newEta)
            assert(Checks.isVec(newEta, this.dimEta) && all(isfinite(newEta)), ...
                'RecedingHorizonUdpLikeController:SetEtaState:InvalidEta', ...
                '** <newEta> must be a %d-dimensional vector **', ...
                this.dimEta);
            
            this.setEtaStateInternal(newEta(:));           
        end
        
        %% getControllerPlantState
        function plantState = getControllerPlantState(this)
            % Get the plant state as currently perceived by the controller, i.e., its internal estimate of the plant state.
            %
            % Returns:
            %   << plantState (Column vector)
            %      The controller's current estimate of the plant state.
            %
            plantState = this.augState(1:this.dimPlantState);
        end
        
        %% computeControlSequence
        function inputSequence = computeControlSequence(this, measurements, measurementDelays, modeMeas, modeDelays)
            if nargin == 1 || isempty(measurements)
                applicableMeas = [];
                applicableDelays = [];
            else
                this.checkMeasurementsAndDelays(measurements, measurementDelays);
                [applicableMeas, applicableDelays] = this.getApplicableMeasurements(measurements, measurementDelays);
            end
            if nargin < 4 || isempty(modeMeas)
                mostRecentMode = [];
            else
                this.checkModeMeasurementsAndDelays(modeMeas, modeDelays);
                % find the mode observation with smallest delay
                % get the ca mode + corresponding delay
                [modeDelay, idx] = min(modeDelays);
                mostRecentMode = modeMeas(idx);                 
            end           
            
            % we stack the available measurements
            % and we infer the SC mode 
            
            % convert current sc mode into base3 number
            measurementModes = mod(floor(this.modeSc ./ 3 .^ [0:this.measBufferLength-1]), 3);
            % now shift the whole vector and replace 1 with 2
            % because they have been processed in the last time step
            newMode = [0 measurementModes(1:end-1)];
            newMode(newMode == 1) = 2;
            % construct the augmented measurement, absence of measurement = 0
            augmentedMeas = zeros(this.dimAugmentedMeas, 1);
            if ~isempty(applicableMeas)
                for idx=1:numel(applicableDelays)
                    delay = applicableDelays(idx);
                    newMode(delay+1) = 1;
                    augmentedMeas(delay * this.dimMeas + 1:(delay + 1) * this.dimMeas) = applicableMeas(:, idx); 
                end
            end

            this.modeSc = sum(3 .^ [0:this.measBufferLength - 1] .* newMode); % new mode, as decimal number
            modeProbSc = [zeros(this.modeSc,1); 1; zeros(3 ^ this.measBufferLength - this.modeSc - 1, 1)];
            Se = this.computeExpectedMeasAvailabilityMatrices(modeProbSc);   

            % integrate the mode observation to update the ca mode, if
            % present
            if ~isempty(mostRecentMode) && modeDelay < this.lastModeObservationDelay
                %mostRecentModeProb = [zeros(mostRecentMode - 1, 1); 1; zeros(numCaModes - mostRecentMode, 1)]; 
                %this.modeCa = (this.transitionMatrixCa ^ modeDelay)' * mostRecentModeProb;
                % simply extract the corresponding row as mode
                % probabilities are unit vector
                probMat = this.transitionMatrixCa ^ modeDelay;
                this.modeCa = probMat(mostRecentMode, :)';                
                this.lastModeObservationDelay = modeDelay;
            end
            
            numCaModes = this.sequenceLength + 1;
            % predict the ca mode probs
            modeProbsCa = zeros(numCaModes, this.horizonLength + 1);
            modeProbsCa(:, 1) = this.modeCa;
            
            % initialization      
            for k = 1:this.horizonLength                
                % compute the ca mode probs
                modeProbsCa(:, k + 1) = this.transitionMatrixCa' * modeProbsCa(:, k);
            end
     
            if this.useMexImplementation
                [this.K, this.L, this.M, iternum] ...
                    = mex_RecedingHorizonUdpLikeController(this.augA, this.augB, this.augQ, this.transitionMatrixCa, ...
                        this.augW, this.augV, this.augmentedMeasMatrix, Se, this.JRJ, this.K, this.L, this.M, ...
                        this.augStateCov, this.augStateSecondMoment, this.terminalAugQ, modeProbsCa, this.maxNumIterations);
            else                
                oldCostsToGo = inf(1, this.horizonLength);
                oldCosts = inf;
                iternum = 1;
                while (iternum < this.maxNumIterations)
                    % forward pass: obtain the second moments for the given sequence of controller gains
                    [Xu, Xl] = this.predictXHorizon(modeProbsCa, Se);

                    Pu = repmat(this.terminalAugQ, 1, 1, numCaModes);
                    Pl = zeros(this.dimAugState, this.dimAugState, numCaModes);                
                    omega = zeros(1, numCaModes);

                    % backward pass to compute the gains
                    for k=this.horizonLength:-1:1                         
                        % compute the gain for the current k
                        [tempM, tempK, tempL, currPu_epsilon, currPl_epsilon] = this.computeGains(squeeze(Xu(:, :, k, :)), squeeze(Xl(:, :, k, :)), ...
                            Pu, Pl, modeProbsCa(:, k), Se(:, :, k));

                        [tempPu, tempPl, tempOmega] = this.computeCostate(currPu_epsilon, currPl_epsilon, omega, tempL, tempM, tempK, Se(:, :, k));
                        Xsum = squeeze(Xl(:, :, k, :) + Xu(:, :, k, :));
                        costsToGo = dot(tempOmega, modeProbsCa(:, k)) + trace(sum(mtimesx(tempPu, Xsum) + mtimesx(tempPl, squeeze(Xu(:, :, k, :))), 3));

                        if costsToGo < oldCostsToGo(k)                       
                            oldCostsToGo(k) = costsToGo;
                            this.M(:, :, k) = tempM;
                            this.L(:, :, k) = tempL;
                            this.K(:, :, k) = tempK;
                            Pu = tempPu;
                            Pl = tempPl;
                            omega = tempOmega;
                        else            
                            % recompute the costate with the old gain
                            [Pu, Pl, omega] = ...
                                this.computeCostate(currPu_epsilon, currPl_epsilon, omega, this.L(:, :, k), this.M(:, :, k), this.K(:, :, k), Se(:, :, k));
                        end
    %                     % compute the gain for the current k
    %                     [this.M{k}, this.K{k}, this.L{k}, currPu_epsilon, currPl_epsilon] = this.computeGains(squeeze(Xu(:, :, k, :)), squeeze(Xl(:, :, k, :)), ...
    %                         Pu, Pl, modeProbsCa(:, k), Se(:, :, k));
    % 
    %                     % update the costate using the computed gains
    %                     [Pu, Pl, omega] = ...
    %                         this.computeCostate(currPu_epsilon, currPl_epsilon, omega, this.L{k}, this.M{k}, this.K{k}, Se(:, :, k));
                    end
                    % the total costs are the costs-to-go/stage costs for the
                    % first time step
                    costs = oldCostsToGo(1);

                    if abs(oldCosts - costs) < RecedingHorizonUdpLikeController.convergenceDiff
                        %fprintf('Converged after %d iterations, current costs (bound): %f\n', iternum, costs);
                        break;
                    else
                        oldCosts = costs;
                        iternum = iternum + 1;
                    end                
                end
            end
            inputSequence = this.doControlSequenceComputation();
            % update the controller state (k+1)
            this.augState = this.M(:, :, 1) * this.augState + this.K(:, :, 1) * augmentedMeas;
            % also update the required moments
            this.updateControllerMoments(Se(:, :, 1));
            % also predict the ca mode
            this.modeCa = this.transitionMatrixCa' * this.modeCa;
            this.lastModeObservationDelay = this.lastModeObservationDelay + 1;
            
            this.lastNumIterations = iternum;
            
%             % shift the current gains
% %             celldisp(this.M)
%             this.M = circshift(this.M, -1, 2);
% %             celldisp(this.M)
%             this.L = circshift(this.L, -1, 2);
%             this.K = circshift(this.K, -1, 2);
%             
%             % and fill with a random matrix
%             this.M(:, :, end) = rand(this.dimAugState)-.5;
%             this.K(:, :, end) = rand(this.dimAugState,this.dimAugmentedMeas)-.5;
%             this.L(:, :, end) = rand(this.dimPlantInput * this.sequenceLength, this.dimAugState)-.5;
                        
            if nargin == 1
                this.lastNumUsedMeas = 0;
                this.lastNumDiscardedMeas = 0;
            else
                this.lastNumUsedMeas = size(applicableMeas, 2);
                this.lastNumDiscardedMeas = size(measurements, 2) - this.lastNumUsedMeas;
            end
        end
        
        %% getLastComputationMeasurementData
        function [numUsedMeas, numDiscardedMeas] = getLastComputationMeasurementData(this)
            % Get information on the last performed control sequence computation (due to a call of computeControlSequence).
            %
            % Returns:
            %   << numUsedMeas (Nonnegative integer)
            %      The number of measurements used by the last measurement update.
            %
            %   << numDiscardedMeas (Nonnegative integer)
            %      The number of measurements discarded during the last
            %      measurement update due to their delays.
            
            numUsedMeas = this.lastNumUsedMeas;
            numDiscardedMeas = this.lastNumDiscardedMeas;
        end
        
        %% reset
        function reset(this)
           
        end
    end
    
    methods (Access = private)
        %% initControllerGains
        function initControllerGains(this)
            % simple initialization procedure: pick initial gains randomly    
            this.M = rand(this.dimAugState, this.dimAugState, this.horizonLength) - 0.5;
            this.K = rand(this.dimAugState, this.dimAugmentedMeas, this.horizonLength) - 0.5;
            this.L = rand(this.dimPlantInput * this.sequenceLength, this.dimAugState, this.horizonLength) - 0.5;             
        end
        
        %% initAugmentedMeasMatrices
        function initAugmentedMeasMatrices(this, modesScBase3, numModesSc, modesScDecimal, C)
            if this.useMexImplementation                
                [this.S, this.augmentedMeasMatrix] = ...
                    mex_RecedingHorizonUdpLikeControllerMeasMatrices(modesScBase3, this.measBufferLength, C, this.dimEta);
            else
                % now the augmented meas matrix (sparse)
                this.augmentedMeasMatrix = [kron(speye(this.measBufferLength), C) zeros(this.dimAugmentedMeas, this.dimEta)];
                % now compute the S matrices, per mode
                this.S = zeros(this.measBufferLength * this.dimMeas, ...
                    this.measBufferLength * this.dimMeas, numModesSc);

                 for j=modesScDecimal
                    % only if measurement mode is 1 it is available for the
                    % controller
                    avail = modesScBase3(j+1, :) == 1;
                    this.S(:, :, j+1) = kron(diag(avail), speye(this.dimMeas));
                 end
            end
        end
        
        %% initAugmentedCostMatrices        
        function initAugmentedCostMatrices(this, numCaModes, H, A, B)
            this.computeAndSetTerminalAugQ(A, B);            
            this.computeAndSetAugQ(numCaModes, H);
        end
        
        %% computeAndSetTerminalAugQ
        function computeAndSetTerminalAugQ(this, A, B)
            [terminalQ, ~, ~, info] = idare(A, B, this.Q, this.R);
            if info.Report > 1
                % either report == 2: The solution is not finite
                % or report == 3: No solution found since the symplectic
                % matrix has eigenvalues on the unit circle
                terminalQ = this.Q;
                warning('RecedingHorizonUdpLikeController:ComputeAndSetTerminalAugQ:UnstabilizablePlant', ...
                    '** (A,B) seems not stabilizable, cannot use stabilizing solution of associated DARE as terminal weighting Q_N **');
            end            
            % terminalQ is the same for all modes
            this.terminalAugQ  = blkdiag(terminalQ, zeros(this.dimAugState - this.dimPlantState));
        end
        
        %% computeAndSetAugQ
        function computeAndSetAugQ(this, numCaModes, H)
            this.augQ = repmat(blkdiag(this.Q, zeros(this.dimAugState - this.dimPlantState)), 1, 1, numCaModes);
            idx = this.dimAugState - this.dimEta + 1;
            % add the mode dependent part
            % zero for first and last mode
            for i = 2:numCaModes - 1                 
                this.augQ(idx:this.dimAugState, idx:this.dimAugState, i) ...
                    = H(:, :, i)' * this.R * H(:, :, i); % H_tilde in the paper
            end
        end
        
        %% computeAndSetScTransitionMatrix
        function computeAndSetScTransitionMatrix(this, scDelayProb, modesScBase3, modesScDecimal)
            
            elementCount = numel(scDelayProb);
            if elementCount <= this.measBufferLength
                % fill up with zeros (row vector)
                trueScDelayProbs = [reshape(scDelayProb, 1, elementCount), zeros(1, this.measBufferLength + 1 - elementCount)];
            else                
                trueScDelayProbs = scDelayProb;
            end
            
            % ensure that the distribution is normalized to avoid NaNs
            idx = find(trueScDelayProbs <= 1e-15);
            normalizedProbs = trueScDelayProbs;
            if ~isempty(idx)
                normalizedProbs(idx) = 1e-15;
                normalizedProbs = normalizedProbs / sum(normalizedProbs); 
            end
            
            if this.useMexImplementation
                this.transitionMatrixSc_part = sparse(mex_RecedingHorizonUdpLikeControllerTransitionMatrix(modesScBase3, this.measBufferLength, normalizedProbs));
            else
                maxMeasDelay = this.measBufferLength - 1;
                % compute the individual transition matrices for i={0, M}
                measStateTransitionMatrices = repmat([0 0 0; 0 0 1; 0 0 1], 1, 1, this.measBufferLength);            
                sums = cumsum(normalizedProbs);
                a_i = normalizedProbs(2:maxMeasDelay+2) ./ (1 -sums(1:this.measBufferLength));
                measStateTransitionMatrices(1, 1, :) = 1 - a_i;
                measStateTransitionMatrices(1, 2, :) = a_i;            

                % only store the first third of the large matrix
                this.transitionMatrixSc_part = sparse(3 ^ maxMeasDelay, 3 ^ this.measBufferLength);
                for i=0:3 ^ maxMeasDelay-1
                    i_idx = modesScBase3(i+1, :); % old mode
                    for j=modesScDecimal
                        j_idx = modesScBase3(j+1, :); % new mode
                        prob= 0;
                        switch j_idx(1)
                            case 0
                                prob = 1-normalizedProbs(1);
                            case 1
                                prob = normalizedProbs(1);
                            otherwise
                                continue                   
                        end
                        for a=1:maxMeasDelay                        
                            prob = prob * measStateTransitionMatrices(i_idx(a)+1, j_idx(a+1)+1, a);
                        end
                        this.transitionMatrixSc_part(i+1,j+1) = prob;                    
                    end          
                end  
            end
%             Z = full([this.transitionMatrixSc_part; this.transitionMatrixSc_part; this.transitionMatrixSc_part])
% %             sum(Z, 2)
% %             %v = [1-trueScDelayProbs(1) trueScDelayProbs(1) 0; 0 0 1; 0 0 1];
%             v = repmat([1-trueScDelayProbs(1) trueScDelayProbs(1) 0]', 1, 3)
% %             %Z2 = kron(v', kron(measStateTransitionMatrices(:, :, 1), measStateTransitionMatrices(:, :, 2)))
%               Z2 = kron(v', kron(measStateTransitionMatrices(:, :, 1), measStateTransitionMatrices(:, :, 2)))
%               cell2mat(arrayfun(@(i) ismember(Z(:, i), Z2), 1:27, 'UniformOutput', false))
%             size(Z2)
%             sum(Z2, 2)
        end
        
        %% computeExpectedMeasAvailabilityMatrices
        function Se = computeExpectedMeasAvailabilityMatrices(this, currentModeSc)
            numCombinations = 3 ^ (this.measBufferLength);
       
            modeProbs = currentModeSc;
            part = numCombinations / 3;
            Se = zeros(this.measBufferLength * this.dimMeas, ...
                this.measBufferLength * this.dimMeas, this.horizonLength);
            % predict the mode probs
            for j=1:this.horizonLength % new mode
                Se(:, :, j) = sum(reshape(modeProbs, 1, 1, []) .* this.S, 3);
                
                modeProbs = this.transitionMatrixSc_part' ...
                    * (modeProbs(1:part) + modeProbs(part + 1:2 * part) ...
                    + modeProbs(2 * part+1:numCombinations));
            end
        end
        
        %% computeGains
        function [M, K, L, Pu_epsilon, Pl_epsilon] = computeGains(this, currXu, currXl, currPu, currPl, currModeProbs, currSe)            
            dimU = this.dimPlantInput;
 
            numCaModes = this.sequenceLength + 1;  
            SC = currSe * this.augmentedMeasMatrix;        
            noisePartV = reshape(currModeProbs, 1,1, numCaModes) .* this.augV;
            Xsums = currXu + currXl;     
            
            % compute the required matrices (they are large, due to kron)
            % initialize with first mode
            Pu_epsilon = cat(3, sum(reshape(this.transitionMatrixCa(1, :), 1, 1, []) .* currPu, 3), ...
                                zeros(this.dimAugState, this.dimAugState, numCaModes -1));
            Pl_epsilon = cat(3, sum(reshape(this.transitionMatrixCa(1, :), 1, 1, []) .* currPl, 3), ...
                                zeros(this.dimAugState, this.dimAugState, numCaModes -1));    
            
            currPu_epsilon = Pu_epsilon(:, :, 1);
            currPl_epsilon = Pl_epsilon(:, :, 1);
            
            Psum = currPu_epsilon + currPl_epsilon;

            Pl_epsilonB = currPl_epsilon * this.augB(:, :, 1);             
            
            SCXl = mtimesx(SC, currXl); % computes SC * curr(:, :, i) in one shot for all i
                      
            Psi = kron(currSe * (noisePartV(:, :, 1) ...
                + this.augmentedMeasMatrix * Xsums(:, :, 1) * this.augmentedMeasMatrix') * currSe', currPl_epsilon);
            Lambda = kron(currXl(:, :, 1), currPl_epsilon);                
            Phi = kron(currXl(:, :, 1), this.augB(:, :,1)' * Psum * this.augB(:, :,1) + this.JRJ);    
            Ypsilon = kron(SCXl(:, :, 1), Pl_epsilonB);                
            Sigma = kron(currXl(:, :, 1), Pl_epsilonB');                
            Pi = kron(SCXl(:, :, 1), currPl_epsilon);             
            
            for i=2:numCaModes                
                currPu_epsilon = sum(reshape(this.transitionMatrixCa(i, :), 1, 1, []) .* currPu, 3); 
                currPl_epsilon = sum(reshape(this.transitionMatrixCa(i, :), 1, 1, []) .* currPl, 3); 
                
                Psum = currPu_epsilon + currPl_epsilon;                  
                Pl_epsilonB = currPl_epsilon * this.augB(:, :, i);              
          
                Psi = Psi + kron(currSe * (noisePartV(:, :, i) ...
                    + this.augmentedMeasMatrix * Xsums(:, :, i) * this.augmentedMeasMatrix') * currSe', currPl_epsilon);
                
                Lambda = Lambda + kron(currXl(:, :, i), currPl_epsilon);

                Phi = Phi + kron(currXl(:, :, i), this.augB(:, :,i)' * Psum * this.augB(:, :,i));

                Ypsilon = Ypsilon + kron(SCXl(:, :, i), Pl_epsilonB);
                
                Sigma = Sigma + kron(currXl(:, :, i), Pl_epsilonB');
                
                Pi = Pi + kron(SCXl(:, :, i), currPl_epsilon);
                
                Pu_epsilon(:, :, i) = currPu_epsilon;
                Pl_epsilon(:, :, i) = currPl_epsilon;
            end
            
            % now come the vectors
            Pl_epsilonA = mtimesx(Pl_epsilon, this.augA);
            rho = sum(mtimesx(mtimesx(Pl_epsilonA, Xsums), SC, 'T'), 3);
            phi = sum(mtimesx(Pl_epsilonA, currXl), 3);
            gamma = sum(mtimesx(mtimesx(mtimesx(this.augB, 'T', Pl_epsilon + Pu_epsilon), this.augA), currXl), 3);
             
            % now solve the system of equations
            A = [Psi -Ypsilon Pi; -Ypsilon' Phi -Sigma; Pi' -Sigma' Lambda];
            A = (A + A') / 2; % ensure that A is symmetric (positive semidefiniteness can be lost due to numerical issues)
            gains = -pinv(A) * [-rho(:); gamma(:); -phi(:)];
            % maybe replace pinv by lsqminnorm
            % and extract the gain matrices
            startIdx = 1;
            endIdx = this.dimAugState * this.dimAugmentedMeas;
            K = reshape(gains(startIdx:endIdx), this.dimAugState, this.dimAugmentedMeas);
           
            startIdx = endIdx + 1;
            endIdx = endIdx + (dimU * this.sequenceLength * this.dimAugState);
            L = reshape(gains(startIdx:endIdx), dimU * this.sequenceLength, this.dimAugState);
           
            startIdx = endIdx + 1;
            M = reshape(gains(startIdx:end), this.dimAugState, this.dimAugState);
        end
        
        %% computeCostate
        function [Pu, Pl, omega] = computeCostate(this, currPu_epsilon, currPl_epsilon, currOmega, currL, currM, currK, Se)
            numCaModes = this.sequenceLength + 1;

            KS = currK * Se;
            KSC = KS * this.augmentedMeasMatrix;
            KSVSK = KS * this.augV * KS'; %E_tilde in the paper
            
            % J matrix is only nonzero for the first mode
            LJRJL = currL' * this.JRJ * currL;
            
            % compute this.Pl_epsilon(:, :, i) * KSVSK in one shot for all i
            tmp = mtimesx(currPl_epsilon, KSVSK) + mtimesx(currPl_epsilon + currPu_epsilon, this.augW);
            % evaluate the epsilon operator for all omegas
            currOmega_epsilon = dot(this.transitionMatrixCa, repmat(currOmega, numCaModes, 1), 2);
            % store omega for all modes as row vector
            omega = currOmega_epsilon' ...
                + arrayfun(@(i) trace(tmp(:, :, i)), 1:numCaModes);
            
            BL = mtimesx(this.augB, currL); % compute this.augB(:, :, i) * L in one shot for all i
            ABL = BL + this.augA;
            MBL = currM-BL;
            ABLMKSC = ABL - currM - KSC; % D_tilde in the paper
            
            Pu = this.augQ + mtimesx(mtimesx(ABL, 'T', currPu_epsilon), ABL) ...
                + mtimesx(mtimesx(ABLMKSC, 'T', currPl_epsilon), ABLMKSC);
            Pu(:, :, 1) =  Pu(:, :, 1) + LJRJL; % additional part for first mode
            
            Pl = mtimesx(mtimesx(BL, 'T', currPu_epsilon), BL) ...
                + mtimesx(mtimesx(MBL, 'T', currPl_epsilon), MBL);
            Pl(:, :, 1) = Pl(:, :, 1) + LJRJL; % additional part for first mode

            % ensure symmetry of costate
            Pu = (Pu + permute(Pu, [2 1 3])) / 2;
            Pl = (Pl + permute(Pl, [2 1 3])) / 2;
        end      
        
        %% updateControllerState
        function updateControllerState(this, augmentedMeas)
            
            % update controller state (k+1)
            this.augState = this.M(:, :, 1) * this.augState + this.K(:, :, 1) * augmentedMeas;
        end
        
        %% updateControllerMoments
        function updateControllerMoments(this, Se)
            numCaModes = this.sequenceLength + 1;
            
            % update the required moments
            newAugStateCov = zeros(this.dimAugState, this.dimAugState, numCaModes);
            newAugStateSecondMoment = zeros(this.dimAugState, this.dimAugState, numCaModes);
            KS = this.K(:, :, 1) * Se;
            KSC = KS * this.augmentedMeasMatrix;
            E_tilde = KS * this.augV * KS'; % E_tilde in the paper
            noise1 = reshape(this.modeCa, 1, 1, numCaModes) .* E_tilde;
            noise2 = noise1 + reshape(this.modeCa, 1, 1, numCaModes) .* this.augW;

            AKSC = this.augA - KSC;
            MKSC = KSC + this.M(:, :, 1);
            % computes this.augB(:, :, i) * L in one single shot for
            % all i
            AKSCMBL = AKSC + mtimesx(this.augB, this.L(:, :, 1)) - this.M(:, :, 1);
            
            AKSC_Xu_CSKA = mtimesx(mtimesx(AKSC, this.augStateCov), AKSC, 'T');
            AKSCMBL_Xl_LBMCSKA = mtimesx(mtimesx(AKSCMBL, this.augStateSecondMoment), AKSCMBL, 'T');
            MKSC_Xl_CSKM = mtimesx(mtimesx(MKSC, this.augStateSecondMoment), MKSC, 'T');
            KSC_Xu_CSK = mtimesx(mtimesx(KSC, this.augStateCov), KSC, 'T');
            
            for j=1:numCaModes % new mode
                newAugStateCov(:, :, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                    .* (AKSC_Xu_CSKA + AKSCMBL_Xl_LBMCSKA + noise2), 3);
                
                newAugStateSecondMoment(:, :, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                    .* (MKSC_Xl_CSKM + KSC_Xu_CSK + noise1), 3);
            end            

%             newAugStateSecondMoment = newAugStateCov + reshape(this.modeCa, 1, 1, numCaModes) .* (this.augState * this.augState');

            % ensure symmetry
            this.augStateCov = (newAugStateCov + permute(newAugStateCov, [2 1 3])) / 2;
            this.augStateSecondMoment = (newAugStateSecondMoment + permute(newAugStateSecondMoment, [2 1 3])) / 2;
        end
        
        %% predictXHorizon
        function [Xu, Xl] = predictXHorizon(this, modeProbsCa, Se)
            numCaModes = this.sequenceLength + 1;
            
            Xu = zeros(this.dimAugState, this.dimAugState, this.horizonLength + 1, numCaModes);
            Xl = zeros(this.dimAugState, this.dimAugState, this.horizonLength + 1, numCaModes);
            
            Xu(:, :, 1, :) = this.augStateCov; 
            Xl(:, :, 1, :) = this.augStateSecondMoment; 
                        
            KS = mtimesx(this.K, Se);
            KSC = mtimesx(KS, this.augmentedMeasMatrix);
            KSVSK = mtimesx(mtimesx(KS, this.augV), KS, 'T'); % E_tilde in the paper
            MKSC = KSC + this.M;
  
            for k=1:this.horizonLength
                AKSC = this.augA - KSC(:, :, k);
                % computes this.augB(:, :, i) * L in one single shot for
                % all i
                AKSCMBL = AKSC + mtimesx(this.augB, this.L(:, :, k)) - this.M(:, :, k);                
                E_tilde = reshape(modeProbsCa(:, k), 1, 1, numCaModes) .* KSVSK(:, :, k);
                noise2 = E_tilde + reshape(modeProbsCa(:, k), 1, 1, numCaModes) .* this.augW;
    
                AKSC_Xu_CSKA = mtimesx(mtimesx(AKSC, squeeze(Xu(:, :, k, :))), AKSC, 'T');
                AKSCMBL_Xl_LBMCSKA = mtimesx(mtimesx(AKSCMBL, squeeze(Xl(:, :, k, :))), AKSCMBL, 'T'); %D_tilde in the paper
                MKSC_Xl_CSKM = mtimesx(mtimesx(MKSC(:, :, k), squeeze(Xl(:, :, k, :))), MKSC(:, :, k), 'T');
                KSC_Xu_CSK = mtimesx(mtimesx(KSC(:, :, k), squeeze(Xu(:, :, k, :))), KSC(:, :, k), 'T');     
          
                for j=1:numCaModes % new mode                 
                    Xu(:, :, k+1, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                        .* (AKSC_Xu_CSKA + AKSCMBL_Xl_LBMCSKA + noise2), 3);
                    Xl(:, :, k+1, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                        .* (MKSC_Xl_CSKM + KSC_Xu_CSK + E_tilde), 3);
                end                
           
                % ensure symmetry
                Xu(:, :, k+1, :) = (squeeze(Xu(:, :, k+1, :)) + permute(squeeze(Xu(:, :, k+1, :)), [2 1 3])) / 2;
                Xl(:, :, k+1, :) = (squeeze(Xl(:, :, k+1, :)) + permute(squeeze(Xl(:, :, k+1, :)), [2 1 3])) / 2;
            end
        end
        
        %% checkMeasurementsAndDelays
        function checkMeasurementsAndDelays(this, measurements, measDelays)
            if ~isempty(measurements)
                assert(Checks.isFixedRowMat(measurements, this.dimMeas), ...
                    'RecedingHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeas', ...
                    '** Individual measurements must be %d-dimensional **', ...
                    this.dimMeas);

                numMeas = size(measurements, 2);
                assert(Checks.isNonNegativeVec(measDelays, numMeas) ...
                        && all(arrayfun(@(delay) mod(delay, 1) == 0, measDelays)), ...
                    'RecedingHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeasDelay', ...
                    '** Each measurement delay (%d in total) must be a nonnegative integer **', numMeas);
            end            
        end
        
        %% getApplicableMeasurements
        function [applicableMeas, applicableDelays] = getApplicableMeasurements(this, measurements, measDelays)
            numMeas = numel(measDelays);
            % the controller assumes only a single measurement per time step
            [~, uniqueIdx, ~] = unique(measDelays, 'stable');
            assert(numel(uniqueIdx) == numMeas, ...
                'RecedingHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsNotUnique', ...
                '** %s\nIgnoring %d of %d measurements. **', ...
                'Controller assumes only one measurement per time step', ...
                numMeas - numel(uniqueIdx), numMeas);
  
            % find the measurements with valid delays
            idx = find(measDelays(uniqueIdx) <= this.maxMeasurementDelay);
            applicableMeas = measurements(:, idx);
            applicableDelays = measDelays(idx);
            if numel(idx) ~= numel(uniqueIdx)
                warning('RecedingHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsDelayTooLarge', ...
                    '** %s\nIgnoring %d of %d measurements. **', ...
                    'Delay too large', ...
                    numel(uniqueIdx) - numel(idx), numel(uniqueIdx));
            end
        end
        
        %% checkModeMeasurementsAndDelays
        function checkModeMeasurementsAndDelays(this, trueModes, modeDelays)
            if ~isempty(trueModes)
                numCaModes = this.sequenceLength + 1;
                assert(Checks.isVec(trueModes) ...
                    && all(arrayfun(@(mode) Checks.isScalarIn(mode, 1, numCaModes) && mod(mode, 1) == 0, trueModes)), ...
                    'RecedingHorizonUdpLikeController:CheckModeMeasurementsAndDelays:InvalidModeObservations', ...
                    '** Mode observations must be a given as a vector of integers from {1, ..., %d} **', ...
                        numCaModes);
                
                assert(Checks.isVec(modeDelays, numel(trueModes)) ...
                    && all(arrayfun(@(x) Checks.isNonNegativeScalar(x)&& mod(x, 1) == 0, modeDelays)), ...
                    'RecedingHorizonUdpLikeController:CheckModeMeasurementsAndDelays:InvalidModeDelays', ...
                    '** Mode observation delays must be a given as a vector of %d nonnegative integers **', ...
                    numel(trueModes));
            end
        end        
      
        %% setControllerPlantStateInternal
        function setControllerPlantStateInternal(this, mean, cov)
            % no checks at all
            this.augState = [repmat(mean, this.measBufferLength, 1); zeros(this.dimEta, 1)];
            % x_o and "past states" are correlated, has to be reflected by cov
            covAug = blkdiag(kron(ones(this.measBufferLength), cov), zeros(this.dimEta));
            this.augStateCov = cat(3, zeros(this.dimAugState, this.dimAugState, this.sequenceLength), covAug);
            this.augStateSecondMoment = cat(3, zeros(this.dimAugState, this.dimAugState, this.sequenceLength), ...
                (this.augState * this.augState') + covAug);   

            this.modeCa = [zeros(this.sequenceLength, 1); 1]; 
        end
        
        %% setEtaStateInternal
        function setEtaStateInternal(this, etaState)
            this.augState(end-this.dimEta-1:end) = etaState;
            % the overall covariance remains unaffected, however update the
            % mode-conditioned second moments
            this.augStateSecondMoment = this.augStateCov + repmat(this.augState * this.augState', 1, 1, this.sequenceLength + 1);
        end
    end
    
    methods (Access = protected)
        %% doControlSequenceComputation
         function inputSequence = doControlSequenceComputation(this)
            inputSequence = this.L(:, :, 1) * this.augState;
         end
         
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, input, ~)
                        
            stageCosts = Utility.computeStageCosts(state, input, this.Q, this.R);
        end
         
         %% doCostsComputation
         function lQGCosts = doCostsComputation(this, stateTrajectory, appliedInputs)
            horizonlength = size(appliedInputs, 2);
            assert(size(stateTrajectory, 2) == horizonlength + 1, ...
                'RecedingHorizonUdpLikeController:DoCostsComputation', ...
                '** <stateTrajectory> is expected to have %d columns ', horizonlength + 1);
            
            lQGCosts = Utility.computeLQGCosts(horizonlength, stateTrajectory, appliedInputs, this.Q, this.R);
        end
    end
end

