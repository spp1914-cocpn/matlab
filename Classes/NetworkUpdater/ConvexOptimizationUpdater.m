classdef ConvexOptimizationUpdater < Updater
    % Class for estimating the transition probability matrix (TPM) of a
    % network described based on a Markov chain. The filter utilizes the
    % probabilities and likelihoods for a specific mode of the Markov jump
    % system as well as covariances and means calculated from the system, 
    % measurement, and input matrix of the system to be controlled.
    %
    % Literature:
    % Wang, Gang. "ML estimation of transition probabilities in jump 
    % Markov systems via convex optimization." IEEE Transactions on 
    % Aerospace and Electronic Systems 46.3 (2010): 1492-1502.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2019  Joanna Mueller <joanna.mueller@student.kit.edu>
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
    properties (Access = private)
        % Integer describing how much time has passed since the filter was
        % initilized.
        timeStep = 1;
        % Integer describing how many time steps are incorporated in the
        % estimation.
        horizon;
        % Array of decimals describing the probabilities of the modes of
        % the Markov jump system at a specific time step.
        modeProbabilities = [];
        % Number of modes describing the Markov jump system.
        numModes;
        % The newest values of updateHistory, likelihoodHistory, and 
        % probabilityHistory (of time step = k) are the first element. 
        % The following values are from time step k + 1 - index.
        % UpdateHistory contains an array of integers which describe a
        % calculation step of TPM update of the last several steps.
        updateHistory = [];
        % likelihoodHistory contains an array of integers which describe
        % the mode likelihoods of the last several steps.
        likelihoodHistory = [];
        % probabilityHistory contains an array of integers which describe
        % the mode probabilities of the last several steps.
        probabilityHistory = [];
        % Vector of Integer describing the measurement mean.
        currentMeasurementMean = [];
        % Matrix describing the measurement covariance.
        currentMeasurementCovariance = [];
        % Matrix describing the relation of state vector and measurements
        % partly.
        measurementMatrix = [];     
        % Matrix describing the relation of the current state and the next
        % state partly.
        systemMatrix = [];
        % Matrix describing the influence of input of the systems for the
        % next state.
        inputMatrix = [];
        % Matrix describing the covariances of the system state.
        stateCovariance = [];
        % Matrix describing the covariances of measurements.
        measurementCovariance = [];
        % Parameter for minimization of mode likelihoods and calculation of
        % TPM update. Integer describing the maximum distance between the
        % start vector and the result vector of the minimization.
        diff;
        % Optimality tolerance of the optimization algorithm.
        opt;
    end
    
    methods (Access = public)
        function this = ConvexOptimizationUpdater( ...
                horizon, numModes, measurementMatrix, systemMatrix, ...
                inputMatrix, stateCovariance, measurementCovariance)
            % Constructor of ConvexOpimitzationUpdater.
            %
            % parameters:
            %   >> horizon: number of time steps to be considered in the
            %   calculation
            %   >> numModes: number of modes of the Markov Jump System as
            %   well as maximum delay time considered + 1
            %   >> measurementMatrix: two-dimensional matrix modelling the
            %   relation of state and measurements
            %   >> systemMatrix: two-dimensional matrix modelling the
            %   relation of the current state and the next state
            %   >> inputMatrix: three-dimensional matrix describing the
            %   influence of control inputs in the next state.
            %   >> stateCovariance: two-dimensional matrix
            %   >> measurementCovariance: two-dimensional matrix
            % Returns:
            %   << this: object of type ConvexOptimizationUpdater
            this.horizon = horizon;
            this.numModes = numModes;
            % Holds the last few values of the sum. The history arrays
            % include the current value as well.
            this.updateHistory = zeros((horizon + 1));
            this.likelihoodHistory = zeros(numModes, (horizon + 1));
            this.probabilityHistory = zeros(numModes, (horizon + 1));
            
            this.measurementMatrix = measurementMatrix;
            this.systemMatrix = systemMatrix;
            this.inputMatrix = inputMatrix;
            this.stateCovariance = stateCovariance;
            this.measurementCovariance = measurementCovariance;
            
            this.opt = 1e-8;
            this.diff = 0.05;
        end
        
        function initializeUpdateStep(this, stateMean, stateCovariance, ...
                modeProbabilities, inputs)
            % Sets and calculates values needed for next update
            %
            % Parameters:
            %   >> stateMean: vector of doubles
            %   >> stateCovariance: matrix of doubles
            %   >> modeProbabilities: mode weights
            %   >> input: vector of doubles
            this.modeProbabilities = ...
                ConvexOptimizationUpdater.boundProbabilities(modeProbabilities);

            this.calculateMean(inputs, stateMean);
            this.calculateCovariance(stateCovariance);
        end
        
        function newTransitionMatrix = updateTransitionProbabilityMatrix(...
                this, measurement, previousTransitionMatrix)
            % Updates transition probability matrix based 
            % on previous probability matrix. 
            %
            % Parameters:
            %   >> measurements: measurements which arrived at the filter
            %   in the current time step.
            %   If there are not enough measurements, use alternative 
            %   calculation
            %   >> previousTransitionMatrix: current two-dimensional 
            %   transition probability matrix
            % Returns: 
            %   << newTransitionMatrix: two-dimentional transition
            %   probability matrix for the next time step            
            newTransitionMatrix = previousTransitionMatrix;
            
            % modeProbabilities = "a priori mode probabilities" = mode
            % weights
            modeLikelihoods = this.getModeLikelihoods(measurement, ...
                this.numModes);
            
            % Calculate the likelihoods of the transition matrix
            % Update the histories.
            this.updateHistories(modeLikelihoods, ...
                previousTransitionMatrix);
            
            if this.timeStep > this.horizon + 1
                x0 = previousTransitionMatrix(this.numModes, :) ...
                    + rand([1, this.numModes]) * 1e-12;
                x0 = x0 / sum(x0);
                options = optimoptions('fmincon', 'Display', 'off', ...
                    'SpecifyObjectiveGradient', false, ...
                    'Algorithm', 'interior-point', ...
                    'OptimalityTolerance', this.opt, ...
                    'MaxIterations', 20, ...
                    'DiffMaxChange', this.diff);
                newTransitionMatrix(this.numModes, :) = ...
                    fmincon( @(x)this.costFunction(x) , ...
                    x0, [], ... % x0, A
                    [], ones( 1, this.numModes ), 1, zeros( 1, this.numModes ), ... % b, Aeq, beq, lb
                    ones( 1, this.numModes ), [], options); % ub

                newTransitionMatrix = RKLUpdater.project(newTransitionMatrix);
            end
        end
        
        function incrementk(this)
            % Function to inform filter that one time step has passed. The
            % time step is used to state if enough measurements were
            % gathered to reasonably calculate mode probabilities and in
            % result the TPM update.
            this.timeStep = this.timeStep + 1;
        end
        
        function updateHistories(this, modeLikelihoods, ...
                previousTransitionMatrix)
            % Only needs to be calculated with respect to last row of the
            % transition probability matrix because the rest of the new 
            % transition matrix will be calculated from that.
            % The newest values (of timeStep = k) is the first element
            this.updateHistory = circshift(this.updateHistory, 1);
            this.probabilityHistory = circshift(this.probabilityHistory, 1, 2);
            this.likelihoodHistory = circshift(this.likelihoodHistory, 1, 2);

            this.updateHistory(1) = dot( this.modeProbabilities(1:(this.numModes-1)), ...
                (previousTransitionMatrix(1:(this.numModes-1), :) * ...
                modeLikelihoods.') );
            this.probabilityHistory(:, 1) = this.modeProbabilities;
            this.likelihoodHistory(:, 1) = modeLikelihoods;
        end
    end
    
    methods (Access = private)  
        %% calculateLikelihoods
        function likelihoods = getModeLikelihoods(this, measurement, ...
                numModes)
            % Calculates the mode likelihoods from the measurement which
            % arrived in the last time step.
            %
            % Parameters:
            %   >> measurements: measurements which arrived at the filter
            %   in the current time step.
            %   If there are not enough measurements, use alternative 
            %   calculation
            %   >> numModes: number of modes of the Markov jump system. 
            % Returns: 
            %   << likelihoods: one-dimentional vector with a likelihood
            %   per mode.
            logLikelihoods = zeros(1, numModes);
            dist = Gaussian();
            for i=1:numModes
                % measMean und MeasCov von RKL mit i = j
                % assume Gaussian likelihood: set the moments
                dist.set(this.currentMeasurementMean(:, i), ...
                    this.currentMeasurementCovariance(:, :, i));
                logLikelihoods(i) = dist.logPdf(measurement);
            end
            % for numerical stability
            maxLogValue = max(logLikelihoods);
            logLikelihoods = logLikelihoods - maxLogValue;
            
            likelihoods = exp(logLikelihoods);
        end
        
        %% OptimizationFunction
        function costFunction = costFunction(this, nextTransitionMatrix)
            % Cost function to calculate the update for the next transition
            % probability matrix.
            % 
            % Parameters:
            %   >> nextTransitionMatrix (1xthis.numModes array): last row
            %   of the new transition matrix
            % Returns:
            %   << fun: cost function to calculate update of the transition
            %   probability matrix
            costFunction = 1;
            for n = 1:this.horizon
                costFunction = dot( costFunction, ...
                       ( ...
                            ( ... 
                                this.probabilityHistory(this.numModes, n) * ...
                                dot(nextTransitionMatrix, this.likelihoodHistory(:, n)) ...
                            ) ...
                            + this.updateHistory(n) ...
                        ) ...
                      );
            end
            costFunction = -costFunction;
        end
        
        function calculateMean(this, inputs, stateEstimate)
            % Calculation of measurement mean for the update of the
            % transition probability matrix.
            % 
            % Parameters:
            %   >> inputs: possible control inputs of the current time
            %   step
            %   >> stateEstimate: estimate of the current one-dimensional
            %   state vector
            lengthInputs = this.numModes;
            
            this.currentMeasurementMean = zeros(size(this.measurementMatrix, 1), ...
              lengthInputs);
      
            if isempty(inputs)
                inputs = zeros(1, lengthInputs);
            end
            
            for i = 1:lengthInputs
                A = this.systemMatrix * stateEstimate(:, i) ...
                    + this.inputMatrix * inputs(i);
                mean = this.measurementMatrix * A;
                this.currentMeasurementMean(:, i) = mean;
            end
        end

        function calculateCovariance(this, ...
              covarianceStateEstimates)
            % Calculation of measurement covariance for computing the 
            % likelihood of the measurement at a specific time step
            %
            % Parameters:
            %   >> covarianceStateEstimates: two-dimensional matrix 
            this.currentMeasurementCovariance = zeros(size(this.measurementMatrix, 1),...
              size(this.measurementMatrix, 1), this.numModes);
            
            for i = 1:this.numModes         
                A = this.systemMatrix;

                covariance = this.measurementMatrix * ( ...
                    A * ...
                    covarianceStateEstimates(:,:,i) * A.' + ...
                    this.stateCovariance ) * this.measurementMatrix.' + ...
                    this.measurementCovariance;
                this.currentMeasurementCovariance(:,:,i) = covariance;
            end
        end
    end
end