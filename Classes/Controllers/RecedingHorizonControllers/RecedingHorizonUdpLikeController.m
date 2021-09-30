classdef RecedingHorizonUdpLikeController < SequenceBasedController & ModelParamsChangeable & CaDelayProbsChangeable & ScDelayProbsChangeable
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
    %    Copyright (C) 2018-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        defaultMaxNumIterations = 20;%30;%1000;
    end
         
    properties (Access = private, Constant)
        convergenceDiff = 1e-20;
        probBoundDelaysSc = 1e-50;
    end
    
    properties (SetAccess = immutable, GetAccess = private)
        dimMeas;
        dimAugmentedMeas;
        measBufferLength;
        Q;
        R;
        
        dimEta;
        
        dimAugState; % dimension of the controller's augmented state       
        JRJ; % only nonzero for the first mode
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
    
    properties (SetAccess = private, GetAccess = ?RecedingHorizonUdpLikeControllerTest)
        augA; % can be changed using changeModelParameters()
        augB;
        augW;
        
        S;
        augmentedMeasMatrix; % sparse!
        
        augQ;
        terminalAugQ; % the same for all modes
        
        % stores the most recent mode transition matrix T_k used for the
        % prediction of theta_k over the optimization horizon
        transitionMatrixCa;
        % history to store previous transition matrices T_{k-j} states (required for incorporation of mode observations)        
        % transitionMatrixCaHistory(:, :, 1) is thus T_{k-1}, used to propagate theta_{k-1} to theta_k
        % transitionMatrixCaHistory(:, :, 2) is thus T_{k-2}, used to propagate theta_{k-2} to theta_{k-1}
        % transitionMatrixCaHistory(:, :, 3) is thus T_{k-3}, used to propagate theta_{k-3} to theta_{k-2}
        % ...
        % transitionMatrixCaHistory(:, :, maxMeasurementDelay) is thus T_{k-maxMeasurementDelay}, used to propagate theta_{k-maxMeasurementDelay} to theta_{k-maxMeasurementDelay+1}
        transitionMatrixCaHistory;
        
        % stores the most recent transition matrix P_k used for the
        % prediction of measurement delay tau_k over the optimization horizon
        scDelayTransitionMatrix;
        % history to store previous transition matrices P_{k-j} states (required for incorporation of delayed measurements)
        % transitionMatrixScHistory(:, :, 1) is thus P_{k-1}, used to propagate tau_{k-1} to tau_k
        % transitionMatrixScHistory(:, :, 2) is thus P_{k-2}, used to propagate tau_{k-2} to tau_{k-1}
        % transitionMatrixScHistory(:, :, 3) is thus P_{k-3}, used to propagate tau_{k-3} to tau_{k-2}
        % ...
        % transitionMatrixScHistory(:, :, maxMeasurementDelay) is thus T_{k-maxMeasurementDelay}, used to propagate tau_{k-maxMeasurementDelay} to tau_{k-maxMeasurementDelay+1}
        transitionMatrixScHistory;

        % we store the previous estimate of the measurement delays tau_k to tau_k-M (unit vector if measurement is available)
        measDelayProbs;
        % measDelayProbs(:, 1) = current delay probs (time k)
        % measDelayProbs(:, 2) = previous delay probs(time k-1)
        % ...
        % measDelayProbs(:, M+1) = delay probs at time k-M, where M is maxMeasDelay
        measAvailabilityStates; % as binary number, ordered according to de2bi
        measAvailabilityIdx; % indices of meas delay probs corresponding to availability states
        
        lastNumUsedMeas = 0;
        lastNumDiscardedMeas = 0;        
    
        modeCa; % this one is only partially observed (stored as probability (column) vec)
        lastModeObservationDelay = inf; % store the time delay to the last mode observation
        
        L;
        K;
        augState; % the controller's augmented plant state
        augStateSecondMoment; % internal estimate of second moment per ca mode
        augStateCov; % internal estimate of error covariance per ca mode
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
            this.dimEta = dimU * (sequenceLength * (sequenceLength - 1) / 2);
            % augmented state consists of x_k, x_k-1, ..., x_k-M and eta_k
            this.dimAugState = dimX * this.measBufferLength + this.dimEta;
            this.dimAugmentedMeas = this.dimMeas * this.measBufferLength;
            
            % and the augmented noise covs
            this.augW = blkdiag(W, zeros(maxMeasDelay*dimX+this.dimEta));
            this.augV = kron(eye(this.measBufferLength), V);

            this.horizonLength = horizonLength;           
                        
            [F, G, H, ~] = Utility.createActuatorMatrices(sequenceLength, dimU);
            J = [speye(dimU) zeros(dimU, dimU * (sequenceLength - 1))];
            this.JRJ = J' * R * J;
            
            D = eye(dimX); 
            if maxMeasDelay > 0
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
            else
                % corner case: max measurement delay = 0
                [~,~,~,~,this.augA, this.augB] = Utility.createAugmentedPlantModel(sequenceLength, A, B);
            end

            this.initAugmentedCostMatrices(numCaModes, H, A, B);

            this.setControllerPlantStateInternal(x0(:), x0Cov);
            
            this.transitionMatrixCaHistory = repmat(this.transitionMatrixCa, 1, 1, this.maxMeasurementDelay);
            
            % measurement availability is the second mode of the system
            % gamma_k+j|k-i = 1: measurement y_{k-i} can be processed at time k+j (i.e., it arrives at time k+j, thus has a delay of j+i time steps)
            % gamma_k_k+j|k-i = 0: measurement y_{k-i} cannot be processed at time k+j (i.e., it has already been processed, and therefore arrived at a previous time, or is simply not yet available)
            % -> delay of y_{k-i} is not i+j time steps
            % so two possible states per relevant measurement, in total
            % 2^M+1, with M the maxMeasurement delay to be maintained at
            % time k: gamma_k|k, gamma_k|k-1, ..., gamma_k|k-M
            this.measAvailabilityStates = de2bi(0:2^this.measBufferLength-1);
            
            % fix the number of elements of the delay probs: M+2
            normalizedScDelayProbs = Utility.normalizeProbabilities(Utility.truncateDiscreteProbabilityDistribution(scDelayProb, maxMeasDelay + 2), ...
                RecedingHorizonUdpLikeController.probBoundDelaysSc);
            this.scDelayTransitionMatrix = repmat(reshape(normalizedScDelayProbs, 1, []), numel(normalizedScDelayProbs), 1);
            this.transitionMatrixScHistory = repmat(this.scDelayTransitionMatrix, 1, 1, this.maxMeasurementDelay);
            
            % extract the corresponding indices
            allIdx = 1:2^this.measBufferLength;
            this.measAvailabilityIdx = cell(this.measBufferLength, 2^this.measBufferLength);
            k = 1;
            for j=1:this.measBufferLength
                zeroIdx = [];
                for i=1:2^(k-1)        
                    zeroIdx = [zeroIdx, i:2^k:2^this.measBufferLength];
                end
                k = k + 1;
                [this.measAvailabilityIdx{j,zeroIdx}] = deal([1:j-1, j+1:numel(normalizedScDelayProbs)]);
                [this.measAvailabilityIdx{j, setdiff(allIdx, zeroIdx)}] = deal(j);
            end           
            % initial guess: we don't know anything about the measurement delays
            this.measDelayProbs = zeros(numel(normalizedScDelayProbs), this.measBufferLength) + 1/numel(normalizedScDelayProbs);

            this.initAugmentedMeasMatrices(C);
            this.initControllerGains();
        end
        
        %% setInitialControllerGains
        function setInitialControllerGains(this, initialK, initialL)
            % Set the controller gains K and L to be employed at the beginning of the time of operation.
            %
            % This function is intended to be used to set initial first
            % guess gains that serve as starting point for the first
            % optimization, or for testing purposes.
            %
            % Parameters:
            %   >> initialK (3D matrix)
            %      N initial gains K_0, K_1, ... K_{N-1} slice-wise arranged, where N is the optimization horizon.
            %
            %   >> initialL (3D matrix)
            %      N initial gains L_0, L_1, ... L_{N-1} slice-wise arranged, where N is the optimization horizon.
            
            assert(Checks.isMat3D(initialK, this.dimAugState, this.dimAugmentedMeas, this.horizonLength), ...
                'RecedingHorizonUdpLikeController:SetInitialControllerGains:InvalidK', ...
                '** <initialK> must be a 3D matrix of dimension %d-by-%d-by-%d **', ...
                this.dimAugState, this.dimAugmentedMeas, this.horizonLength); 
             assert(Checks.isMat3D(initialL, this.dimPlantInput * this.sequenceLength, this.dimAugState, this.horizonLength), ...
                'RecedingHorizonUdpLikeController:SetInitialControllerGains:InvalidL', ...
                '** <initialL> must be a 3D matrix of dimension %d-by-%d-by-%d **', ...
                this.dimPlantInput * this.sequenceLength, this.dimAugState, this.horizonLength); 
            
            this.K = initialK;
            this.L = initialL;
        end
        
        %% changeCaDelayProbs
        function changeCaDelayProbs(this, newCaDelayProbs)
            % Change the distribution of the control packet delays to be
            % assumed by the controller during optimization.
            %
            % Parameters:
            %  >> newCaDelayProbs (Nonnegative vector)
            %     Vector specifiying the new delay distribution.
            %
            newMat = Utility.calculateDelayTransitionMatrix(...
                Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.sequenceLength + 1));
            if ~isequal(newMat, this.transitionMatrixCa)
                this.transitionMatrixCa = newMat;               
            end            
        end
        
        %% changeScDelayProbs
        function changeScDelayProbs(this, newScDelayProbs)
            % Change the distribution of the measurement packet delays to be
            % assumed by the controller during optimization.
            %
            % Parameters:
            %  >> newScDelayProbs (Nonnegative vector)
            %     Vector specifiying the new delay distribution.
            %
            
            % we must truncate the distribution so that the number of
            % element is conformable, then normalize             
            normalizedScDelayProbs = Utility.normalizeProbabilities(...
                Utility.truncateDiscreteProbabilityDistribution(newScDelayProbs, size(this.scDelayTransitionMatrix, 1)), ...
                RecedingHorizonUdpLikeController.probBoundDelaysSc);
            newMat = repmat(reshape(normalizedScDelayProbs, 1, []), numel(normalizedScDelayProbs), 1);            
            if ~isequal(newMat, this.scDelayTransitionMatrix)
                this.scDelayTransitionMatrix = newMat;
            end            
        end
        
        %% changeModelParameters
        function changeModelParameters(this, newA, newB, newW)
            Validator.validateSystemMatrix(newA, this.dimPlantState);
            Validator.validateInputMatrix(newB, this.dimPlantState, this.dimPlantInput);
            Validator.validateSysNoiseCovarianceMatrix(newW, this.dimPlantState);
            
            % and the augmented plant noise cov changes
            this.augW(1:this.dimPlantState, 1:this.dimPlantState) = newW;           
            
            [~, ~, H, J] = Utility.createActuatorMatrices(this.sequenceLength, this.dimPlantInput);
            
            % A_tilde in the paper: update the A portion of A_tilde
            this.augA(1:this.dimPlantState, 1:this.dimPlantState, 1) = newA;
            for i=2:this.sequenceLength
                this.augA(1:this.dimPlantState, 1:this.dimPlantState, i) = newA;
                % update B*H for all but first and last mode
                % H is empty for first and last mode
                this.augA(1:this.dimPlantState, end-this.dimEta +1:end, i) = newB * H(:, :, i);
            end
            this.augA(1:this.dimPlantState, 1:this.dimPlantState, this.sequenceLength + 1) = newA;            
            % augB is only different for first mode (multiply B by J)
            % B*J is zero for all modes but first
            %
            this.augB(1:this.dimPlantState, :, 1) = newB * J(:, :, 1); 
            % terminal state weighting Q_K has to be adapted,
            % dynamics-dependent
            this.computeAndSetTerminalAugQ(newA, newB);
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
            mostRecentMode = [];
            if nargin >= 4 && ~isempty(modeMeas)
                this.checkModeMeasurementsAndDelays(modeMeas, modeDelays);
                % find the mode observation with smallest delay
                % get the ca mode + corresponding delay
                [modeDelay, idx] = min(modeDelays);
                if modeDelay <= this.maxMeasurementDelay
                    mostRecentMode = modeMeas(idx);
                end
            end
            
            % we stack the available measurements 
            % construct the augmented measurement, absence of measurement = 0 w.l.o.g.
            augmentedMeas = zeros(this.dimAugmentedMeas, 1);
            newMeasMode = zeros(1, this.measBufferLength); % row vector (gamma_k)
            for j=this.maxMeasurementDelay:-1:0
                % check if measurement y_{k-j} is available
                [avail, idx]= ismember(j, applicableDelays);
                if avail
                    % measurement is available for processing, it arrived
                    % at time k, so delay is j
                    augmentedMeas(j * this.dimMeas + 1:(j + 1) * this.dimMeas) = applicableMeas(:, idx);
                    % update the corresponding delay prob (1 at position j+1, zero elsewhere)
                    this.measDelayProbs(1:j, j+1) = 0;
                    this.measDelayProbs(j+1 , j+1) = 1;
                    this.measDelayProbs(j+2:end, j +1) = 0;
                    newMeasMode(j + 1) = 1; % indicate availability
                elseif ~any(this.measDelayProbs(1:j, j+1) == 1)
                    % measurement has not arrived, so (1) its delay must be > j
                    % or (2) it has already arrived in a previous time step
                    % in the latter case (2), the delay probs are already a
                    % unit vector with one at position 1<= i <=j+1
                    % only in the former case (1) we can update the delay probs
                    this.measDelayProbs(1:j+1, j+1) = 0;
                    % and normalize (implying that we divide by Pr[delay > j])
                    this.measDelayProbs(:, j+1) = this.measDelayProbs(:, j+1) ./ sum(this.measDelayProbs(:, j+1));
                end
                % this also affects the next state, so we can update our delay estimates
                if j~=0
                   this.measDelayProbs(:, j) = this.transitionMatrixScHistory(:, :, j)' * this.measDelayProbs(:, j+1);
                end
            end           
            % get S_exp=E[S(gamma_k+n)|I_k] for the whole horizon based on the current available information I_k
            Sexp = this.computeExpectedMeasAvailabilityMatrices(newMeasMode);
            % integrate the mode observation to update the ca mode theta_k, if
            % present
            if ~isempty(mostRecentMode) && modeDelay < this.lastModeObservationDelay                
                % simply extract the corresponding row as mode
                % probabilities are unit vector
                probMat = eye(this.sequenceLength + 1);
                for j=modeDelay:-1:1
                    probMat = probMat * this.transitionMatrixCaHistory(:, :, j); 
                end                
                this.modeCa = probMat(mostRecentMode, :)';                
                this.lastModeObservationDelay = modeDelay;
            end
            
            numCaModes = this.sequenceLength + 1;
            
            if this.useMexImplementation
                [this.K, this.L, iternum] ...
                    = mex_RecedingHorizonUdpLikeController(this.augA, this.augB, this.augQ, this.transitionMatrixCa, ...
                        this.augW, this.augV, this.augmentedMeasMatrix, Sexp, this.JRJ, this.K, this.L, ...
                        this.augStateCov, this.augStateSecondMoment, this.terminalAugQ, this.modeCa, this.maxNumIterations);
                Aexp =  sum(reshape(this.modeCa, 1, 1, numCaModes) .* this.augA, 3);
                Bexp =  sum(reshape(this.modeCa, 1, 1, numCaModes) .* this.augB, 3);
            else
                % predict the ca mode probs over the horizon from theta_k to theta_{k+K}
                modeProbsCa = zeros(numCaModes, this.horizonLength + 1);
                modeProbsCa(:, 1) = this.modeCa; %theta_k based on current information set I_k
                Aexp = zeros(this.dimAugState, this.dimAugState, this.horizonLength);
                Bexp = zeros(this.dimAugState, this.dimPlantInput * this.sequenceLength, this.horizonLength);
                % initialization
                for j = 1:this.horizonLength
                    % compute the expected dynamic matrices A_hat = E[augA(theta_k+j)|I_k] and
                    % B_hat = E[augB(theta_k+j)|I_k]
                    Aexp(:, :, j) = sum(reshape(modeProbsCa(:, j), 1, 1, numCaModes) .* this.augA, 3); % Ahat_{k+j}
                    Bexp(:, :, j) = sum(reshape(modeProbsCa(:, j), 1, 1, numCaModes) .* this.augB, 3); % Bhat_{k+j}
                    % predict the ca mode probs
                    modeProbsCa(:, j + 1) = this.transitionMatrixCa' * modeProbsCa(:, j); %theta_{k+j}
                end

                oldCostsToGo = inf(1, this.horizonLength);
                oldCosts = inf;
                iternum = 0;
                while (iternum < this.maxNumIterations)
                    iternum = iternum + 1;
                    % forward pass: obtain the second moments for the given sequence of controller gains
                    [Xu, Xl] = this.predictXHorizon(modeProbsCa, Sexp, Aexp, Bexp);

                    Pu = repmat(this.terminalAugQ, 1, 1, numCaModes);
                    Pl = zeros(this.dimAugState, this.dimAugState, numCaModes);                
                    omega = zeros(1, numCaModes);

                    % backward pass to compute the gains
                    for k=this.horizonLength:-1:1                         
                        % compute the gains (K_k,L_k) for the current k
                        [tempK, tempL, currPu_epsilon, currPl_epsilon] = this.computeGains(squeeze(Xu(:, :, k, :)), squeeze(Xl(:, :, k, :)), ...
                            Pu, Pl, modeProbsCa(:, k), Sexp(:, :, k), Aexp(:, :, k), Bexp(:, :, k));
                        
                        [tempPu, tempPl, tempOmega] = this.computeCostate(currPu_epsilon, currPl_epsilon, omega, tempL, tempK, Sexp(:, :, k), Aexp(:, :, k), Bexp(:, :, k));
                        Xsum = squeeze(Xl(:, :, k, :) + Xu(:, :, k, :));
                        costsToGo = dot(tempOmega, modeProbsCa(:, k)) + trace(sum(mtimesx(tempPu, Xsum) + mtimesx(tempPl, squeeze(Xu(:, :, k, :))), 3));
                                                
                        if costsToGo < oldCostsToGo(k)
                            oldCostsToGo(k) = costsToGo;
                            this.L(:, :, k) = tempL;
                            this.K(:, :, k) = tempK;
                            Pu = tempPu;
                            Pl = tempPl;
                            omega = tempOmega;
                        else            
                            % recompute the costate with the old gain
                            [Pu, Pl, omega] = ...
                                this.computeCostate(currPu_epsilon, currPl_epsilon, omega, this.L(:, :, k), this.K(:, :, k), Sexp(:, :, k), Aexp(:, :, k), Bexp(:, :, k));
                        end
                    end
                    % the total costs are the costs-to-go/stage costs for the
                    % first time step
                    costs = oldCostsToGo(1);
                    if abs(oldCosts - costs) < RecedingHorizonUdpLikeController.convergenceDiff
                        %fprintf('Converged after %d iterations, current costs (bound): %f\n', iternum, costs);                        
                        break;
                    else
                        oldCosts = costs;                        
                    end                
                end
            end
            
            inputSequence = this.doControlSequenceComputation();            
            innovation = augmentedMeas - Sexp(:, :, 1) * this.augmentedMeasMatrix * this.augState;
                        
            % update the controller state (k+1)            
            this.augState = Aexp(:, :, 1) * this.augState + Bexp(:, :, 1) * inputSequence + this.K(:, :, 1) * innovation;            
            
            % also update the required moments
            this.updateControllerMoments(Sexp(:, :, 1), Aexp(:, :, 1), Bexp(:, :, 1));
            % also predict the ca mode
            this.modeCa = this.transitionMatrixCa' * this.modeCa;
            this.lastModeObservationDelay = this.lastModeObservationDelay + 1;
            this.transitionMatrixCaHistory  = cat(3, this.transitionMatrixCa, ...
                this.transitionMatrixCaHistory(:, :, 1:end-1));

            % we also have to update the measurement delay probs
            % predict the meas delay probs (k+1)            
            if this.measBufferLength == 1
                % corner case
                this.measDelayProbs = this.scDelayTransitionMatrix' * this.measDelayProbs;
            else
                this.measDelayProbs = circshift(this.measDelayProbs, 1, 2);
                this.measDelayProbs(:, 1) = this.scDelayTransitionMatrix' * this.measDelayProbs(:, 2);
            end

            this.transitionMatrixScHistory = cat(3, this.scDelayTransitionMatrix, ...
                this.transitionMatrixScHistory(:, :, 1:end-1));

            this.lastNumIterations = iternum;
            
            % shift the current gains
%             this.M = circshift(this.M, -1, 3);
%             this.L = circshift(this.L, -1, 3);
%             this.K = circshift(this.K, -1, 3);
            
            % and fill with a random matrix
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
            this.K = rand(this.dimAugState, this.dimAugmentedMeas, this.horizonLength);
            this.L = rand(this.dimPlantInput * this.sequenceLength, this.dimAugState, this.horizonLength);            
        end
        
        
        %% initAugmentedMeasMatrices
        function initAugmentedMeasMatrices(this, C)
            if this.useMexImplementation                
                [this.S, this.augmentedMeasMatrix] = ...
                    mex_RecedingHorizonUdpLikeControllerMeasMatrices(this.measAvailabilityStates, C, this.dimEta);                
            else
                % now the augmented meas matrix (sparse)
                this.augmentedMeasMatrix = [kron(speye(this.measBufferLength), C) zeros(this.dimAugmentedMeas, this.dimEta)];
                
                % now compute the S matrices, per mode
                this.S = zeros(this.dimAugmentedMeas,this.dimAugmentedMeas, 2^this.measBufferLength);
         
                for j=1:size(this.measAvailabilityStates, 1)
                    % only if measurement mode is 1 it is available for the
                    % controller
                    avail = this.measAvailabilityStates(j, :) == 1;
                    this.S(:, :, j) = kron(diag(avail), speye(this.dimMeas));                    
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
        
        %% computeExpectedMeasAvailabilityMatrices
        function Sexp = computeExpectedMeasAvailabilityMatrices(this, currentMeasAvailability)            
            % at the current time k, we know the availability of the measurements (gamma_k)
            % simply check which measurements are available
            Sexp = cat(3, this.S(:, :, bi2de(currentMeasAvailability) + 1), ...
                zeros(this.dimAugmentedMeas, this.dimAugmentedMeas, this.horizonLength -1));
            
            predMeasDelayProbs = this.measDelayProbs; % at time k, last column is delay probs for y_k-M
            % we use the whole matrix of delay probs here and do not predict with the second-last column 
            % corresponding to y_k-M+1 (that column is needed to compute availability probs gamma_k+1)
            % tau_k-M+2 is needed to compute probs gamma_k+2
            % ...
            % tau_k is needed to compute probs gamma_k+M
            % and so forth
            % tau_k-M+n is needed to compute probs gamma_k+n
            % reason why we use the whole matrix is that we might have updated several delay probs 
            % due to received measurements or knowledge that a measurement has not arrived yet
            for n=1:this.horizonLength-1
                % predict the meas delay probs (k+1)
                % here we make a small error if delay transition matrix
                % changed in the last M time steps
                % we always predict the delay probabilities with the most recent transition matrix
                % predMeasDelayProbs(:, 1) contains the delay probs for time k+n, i.e. at stage n of the horizon
                predMeasDelayProbs(:, end) = this.scDelayTransitionMatrix' * predMeasDelayProbs(:, 1); % delay probs for stage n+1
                predMeasDelayProbs = circshift(predMeasDelayProbs, 1, 2); % shift columns to the right                
                % now we can predict the measurement availability at stage n+1 (time k+n+1)
                availProbs = this.computeMeasAvailabilityProbs(predMeasDelayProbs(:, end)); 
                % this gives us the expected measurement availability matrix
                Sexp(:, :, n+1) = sum(reshape(availProbs, 1, 1, []) .* this.S, 3);                
            end  
        end
        
        %% computeGains
        function [K, L, Pu_epsilon, Pl_epsilon] = computeGains(this, currXu, currXl, currPu, currPl, ...
                currModeProbs, currSexp, currAexp, currBexp)
                                    
            numCaModes = this.sequenceLength + 1;  
            CS = (currSexp * this.augmentedMeasMatrix)'; %SC'
            noisePartV = reshape(currModeProbs, 1,1, numCaModes) .* this.augV;            
            
            % compute the required matrices (they are large, due to kron)
            % initialize with first mode
            Pu_epsilon = cat(3, sum(reshape(this.transitionMatrixCa(1, :), 1, 1, []) .* currPu, 3), ...
                                zeros(this.dimAugState, this.dimAugState, numCaModes -1));
            Pl_epsilon = cat(3, sum(reshape(this.transitionMatrixCa(1, :), 1, 1, []) .* currPl, 3), ...
                                zeros(this.dimAugState, this.dimAugState, numCaModes -1));   
            
            currPu_epsilon = Pu_epsilon(:, :, 1);
            currPl_epsilon = Pl_epsilon(:, :, 1);
            
            % quadratic in K
            Psi = kron(currSexp * (noisePartV(:, :, 1) ...
                + this.augmentedMeasMatrix * currXu(:, :, 1) * this.augmentedMeasMatrix') * currSexp', currPl_epsilon);
            % quadratic in L
            part = this.JRJ + this.augB(:, :, 1)' * (currPu_epsilon)* this.augB(:, :, 1) ...
                + (this.augB(:, :, 1)-currBexp)' * currPl_epsilon * (this.augB(:, :, 1)-currBexp);
            Phi = kron(currXl(:, :, 1), part);            
                      
            % linear in vec(k)
            % negative!
            rho = currPl_epsilon * this.augA(:, :, 1) * currXu(:, :, 1) * CS;
            
            % linear in vec (L)
            gamma = (this.augB(:, :, 1)' * currPu_epsilon * this.augA(:, :, 1) ...
                + (this.augB(:, :, 1) - currBexp)' * currPl_epsilon * (this.augA(:, :, 1) - currAexp)) ...
                * currXl(:, :, 1);
            
            for i=2:numCaModes                                
                currPu_epsilon = sum(reshape(this.transitionMatrixCa(i, :), 1, 1, []) .* currPu, 3); 
                currPl_epsilon = sum(reshape(this.transitionMatrixCa(i, :), 1, 1, []) .* currPl, 3); 
                                
                Pu_epsilon(:, :, i) = currPu_epsilon;
                Pl_epsilon(:, :, i) = currPl_epsilon;                         
          
                Psi = Psi + kron(currSexp * (noisePartV(:, :, i) ...
                    + this.augmentedMeasMatrix * currXu(:, :, i) * this.augmentedMeasMatrix') * currSexp', currPl_epsilon);
                
                part = this.augB(:, :, i)' * (currPu_epsilon)* this.augB(:, :, i) ...
                    + (this.augB(:, :, i)-currBexp)' * currPl_epsilon * (this.augB(:, :, i)-currBexp);
                Phi = Phi + kron(currXl(:, :, i), part);               
                
                rho = rho + currPl_epsilon * this.augA(:, :, i) * currXu(:, :, i) * CS;
                                
                gamma = gamma + (this.augB(:, :, i)' * currPu_epsilon * this.augA(:, :, i) ...
                    + (this.augB(:, :, i) - currBexp)' * currPl_epsilon * (this.augA(:, :, i) - currAexp)) ...
                    * currXl(:, :, i);
            end            
            % K, L
            K = reshape(lsqminnorm(Psi, rho(:)), this.dimAugState, this.dimAugmentedMeas);
            L = reshape(lsqminnorm(Phi, -gamma(:)), this.dimPlantInput * this.sequenceLength, this.dimAugState);
        end

        %% computeCostate
        function [Pu, Pl, omega] = computeCostate(this, currPu_epsilon, currPl_epsilon, currOmega, currL, currK, currSexp, currAexp, currBexp)
            numCaModes = this.sequenceLength + 1;

            KS = currK * currSexp;
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
            AexpBexpL = currAexp + currBexp * currL;            
            ABL = BL + this.augA;
            diffPart = ABL - AexpBexpL; % for all modes
            O_tilde = AexpBexpL-BL - KSC; %Aexp+Bexp*L-K*S*C-B(i)*L
                  
            Pu = this.augQ + mtimesx(mtimesx(ABL, 'T', currPu_epsilon), ABL) ...
                + mtimesx(mtimesx(diffPart, 'T', currPl_epsilon), diffPart);
            Pu(:, :, 1) =  Pu(:, :, 1) + LJRJL; % additional part for first mode
            
            Pl = mtimesx(mtimesx(BL, 'T', currPu_epsilon), BL) ...
                + mtimesx(mtimesx(O_tilde, 'T', currPl_epsilon), O_tilde);
            Pl(:, :, 1) = Pl(:, :, 1) + LJRJL; % additional part for first mode

            % ensure symmetry of costate
            Pu = (Pu + permute(Pu, [2 1 3])) / 2;
            Pl = (Pl + permute(Pl, [2 1 3])) / 2;
        end
                
        %% updateControllerMoments
        function updateControllerMoments(this, Sexp, Aexp, Bexp)
            numCaModes = this.sequenceLength + 1;
            
            % update the required moments
            newAugStateCov = zeros(this.dimAugState, this.dimAugState, numCaModes);
            newAugStateSecondMoment = zeros(this.dimAugState, this.dimAugState, numCaModes);
            KS = this.K(:, :, 1) * Sexp;
            KSC = KS * this.augmentedMeasMatrix;
            E_tilde = KS * this.augV * KS'; % E_tilde in the paper
            noise1 = reshape(this.modeCa, 1, 1, numCaModes) .* E_tilde;
            noise2 = noise1 + reshape(this.modeCa, 1, 1, numCaModes) .* this.augW;

            AKSC = this.augA - KSC;
            AexpBexpL = Aexp + Bexp * this.L(:, :, 1);
            diffPart = this.augA - Aexp + mtimesx(this.augB - Bexp, this.L(:, :, 1)); % for all modes
            AABBL_Xl_LBBAA = mtimesx(mtimesx(diffPart, this.augStateSecondMoment), diffPart, 'T');
            AKSC_Xu_CSKA = mtimesx(mtimesx(AKSC, this.augStateCov), AKSC, 'T');
                        
            KSC_Xu_CSK = mtimesx(mtimesx(KSC, this.augStateCov), KSC, 'T');     
            ABL_Xl_BLA = mtimesx(mtimesx(AexpBexpL, this.augStateSecondMoment), AexpBexpL, 'T');
            
            for j=1:numCaModes % new mode
                newAugStateCov(:, :, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                    .* (AKSC_Xu_CSKA + AABBL_Xl_LBBAA + noise2), 3);
                
                newAugStateSecondMoment(:, :, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                    .* (ABL_Xl_BLA + KSC_Xu_CSK + noise1), 3);
            end

            % ensure symmetry
            this.augStateCov = (newAugStateCov + permute(newAugStateCov, [2 1 3])) / 2;
            this.augStateSecondMoment = (newAugStateSecondMoment + permute(newAugStateSecondMoment, [2 1 3])) / 2;
        end
        
        %% predictXHorizon
        function [Xu, Xl] = predictXHorizon(this, modeProbsCa, Sexp, Aexp, Bexp)
            numCaModes = this.sequenceLength + 1;
            
            Xu = zeros(this.dimAugState, this.dimAugState, this.horizonLength + 1, numCaModes);
            Xl = zeros(this.dimAugState, this.dimAugState, this.horizonLength + 1, numCaModes);
            
            Xu(:, :, 1, :) = this.augStateCov; 
            Xl(:, :, 1, :) = this.augStateSecondMoment;                         
            
            KS = mtimesx(this.K, Sexp);
            KSC = mtimesx(KS, this.augmentedMeasMatrix);
            KSVSK = mtimesx(mtimesx(KS, this.augV), KS, 'T'); % E_tilde in the paper
            AexpBexpL = Aexp + mtimesx(Bexp, this.L); % for all k
            for k=1:this.horizonLength
                AKSC = this.augA - KSC(:, :, k);
                diffPart = this.augA - Aexp(:, :, k) + mtimesx(this.augB - Bexp(:, :, k), this.L(:, :, k));
                AABBL_Xl_LBBAA = mtimesx(mtimesx(diffPart, squeeze(Xl(:, :, k, :))), diffPart, 'T');
                AKSC_Xu_CSKA = mtimesx(mtimesx(AKSC, squeeze(Xu(:, :, k, :))), AKSC, 'T');
                
                E_tilde = reshape(modeProbsCa(:, k), 1, 1, numCaModes) .* KSVSK(:, :, k);
                noise2 = E_tilde + reshape(modeProbsCa(:, k), 1, 1, numCaModes) .* this.augW;                
                                                
                KSC_Xu_CSK = mtimesx(mtimesx(KSC(:, :, k), squeeze(Xu(:, :, k, :))), KSC(:, :, k), 'T');     
                ABL_Xl_BLA = mtimesx(mtimesx(AexpBexpL(:, :, k), squeeze(Xl(:, :, k, :))), AexpBexpL(:, :, k), 'T');
                
                for j=1:numCaModes % new mode                 
                    Xu(:, :, k+1, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                        .* (AKSC_Xu_CSKA + AABBL_Xl_LBBAA + noise2), 3);
                    Xl(:, :, k+1, j) = sum(reshape(this.transitionMatrixCa(:, j), 1, 1, numCaModes) ...
                        .* (ABL_Xl_BLA + KSC_Xu_CSK + E_tilde), 3);
                end               
            end           
            % ensure symmetry
            Xu = (Xu + permute(Xu, [2 1 3 4])) / 2;
            Xl = (Xl + permute(Xl, [2 1 3 4])) / 2;
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
            this.modeCa = [zeros(this.sequenceLength, 1); 1]; 
            % no checks at all
            this.augState = [repmat(mean, this.measBufferLength, 1); zeros(this.dimEta, 1)];
            % x_o and "past states" are correlated, has to be reflected by cov
            covAug = blkdiag(kron(ones(this.measBufferLength), cov), zeros(this.dimEta));
            this.augStateCov = cat(3, zeros(this.dimAugState, this.dimAugState, this.sequenceLength), covAug);
            this.augStateSecondMoment = cat(3, zeros(this.dimAugState, this.dimAugState, this.sequenceLength), ...
                this.augState * this.augState' + covAug);  
        end
        
        %% setEtaStateInternal
        function setEtaStateInternal(this, etaState)
            this.augState(end-this.dimEta-1:end) = etaState;
            % the overall covariance remains unaffected, however update the
            % mode-conditioned second moments
            this.augStateSecondMoment = this.augStateCov + repmat(this.augState * this.augState', 1, 1, this.sequenceLength + 1);
        end
        
        %% computeMeasAvailabilityProbs
        function probs = computeMeasAvailabilityProbs(this, delayProbs)
            if this.useMexImplementation
                probs = mex_RecedingHorizonUdpLikeControllerMeasAvailProbs(this.scDelayTransitionMatrix, this.measAvailabilityIdx, delayProbs);
            else               
                switch this.measBufferLength
                    case 1 % no measurement delay allowed
                        probs = zeros(1, 2);
                        probs(1) = sum(delayProbs(this.measAvailabilityIdx{1, 1}));
                        probs(2) = sum(delayProbs(this.measAvailabilityIdx{1, 2}));
                    case 2
                        % delay probs is a column vector
                        probs = zeros(1, 4);
                        for p=1:4                            
                            localIdx = this.measAvailabilityIdx{2, p};
                            for l=this.measAvailabilityIdx{1, p}        
                                probs(p) = probs(p) + sum(this.scDelayTransitionMatrix(localIdx, l) .* delayProbs(localIdx));
                            end
                        end
                    case 3
                        % delay probs is a column vector
                        probs = zeros(1, 8);
                        for p=1:8
                            localIdx2 = this.measAvailabilityIdx{3, p};
                            for l=this.measAvailabilityIdx{1, p}
                                for m=this.measAvailabilityIdx{2, p}
                                    probs(p) = probs(p) + this.scDelayTransitionMatrix(m, l) * sum(this.scDelayTransitionMatrix(localIdx2, m) .* delayProbs(localIdx2));    
                                end                                 
                            end
                        end
                    case 4
                        % delay probs is a column vector
                        probs = zeros(1, 16); % 2^4
                        for p=1:16
                            localIdx3 = this.measAvailabilityIdx{4, p};
                            for l=this.measAvailabilityIdx{1, p}
                                for m=this.measAvailabilityIdx{2, p}
                                    for n=this.measAvailabilityIdx{3, p}
                                        probs(p) = probs(p) +  this.scDelayTransitionMatrix(m, l) * this.scDelayTransitionMatrix(n, m) ...
                                            * sum(this.scDelayTransitionMatrix(localIdx3, n) .* delayProbs(localIdx3));    
                                    end
                                end
                            end
                        end
                    case 5
                        % delay probs is a column vector
                        probs = zeros(1, 32); % 2^5
                        for p=1:32
                            localIdx4 = this.measAvailabilityIdx{5, p};
                            for l=this.measAvailabilityIdx{1, p}
                                for m=this.measAvailabilityIdx{2, p}
                                    for n=this.measAvailabilityIdx{3, p}
                                        for r=this.measAvailabilityIdx{4, p}
                                            probs(p) = probs(p) + this.scDelayTransitionMatrix(m, l) * this.scDelayTransitionMatrix(n, m) ...
                                                * this.scDelayTransitionMatrix(r, n) ...    
                                                * sum(this.scDelayTransitionMatrix(localIdx4, r) .* delayProbs(localIdx4));    
                                        end
                                    end
                                end
                            end
                        end
                    otherwise
                        % recursion is VERY slow in Matlab
                        probs = this.computeMeasAvailabilityProbsRecursive(delayProbs);
                end
            end                  
        end
        
        %% computeMeasAvailabilityProbsRecursive
        function probs = computeMeasAvailabilityProbsRecursive(this, delayProbs)
            % delay probs is a column vector
            probs = zeros(1, 2^this.measBufferLength);
            for k=1:2^this.measBufferLength
                currIdx = this.measAvailabilityIdx{1, k};
                for l=currIdx        
                    probs(k) = probs(k) + computeRecursive(2, l, delayProbs);
                end
            end
            
            function e = computeRecursive(depth, j, delayProbs)                
                localIdx = this.measAvailabilityIdx{depth, k};                
                if depth== this.measBufferLength  
                    e = sum(this.scDelayTransitionMatrix(localIdx, j) .* delayProbs(localIdx));
                else
                    e = 0;
                    for i=localIdx
                        e = e + this.scDelayTransitionMatrix(i, j) * computeRecursive(depth + 1, i, delayProbs);     
                    end                    
                end
            end     
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

