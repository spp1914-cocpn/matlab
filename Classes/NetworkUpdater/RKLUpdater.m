classdef RKLUpdater < Updater
   % Implementation of a recursive Kullback-Leibler (RKL) estimator based 
   % on Orguner and Demirekler
   %
   % Literature: 
   %    Orguner, Umut & Demirekler, Mubeccel. (2006). An online sequential 
   %    algorithm for the estimation of transition probabilities for 
   %    jump Markov linear systems. Automatica. 42. 1735-1744.
   
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
        % Array of decimals describing the probabilities of the modes of
        % the Markov jump system at a specific time step.
        currentModeWeights = [];
        
        % scale factor, needs to fullfill the following contraints:
        %   infinite sum of every step size needs to be infinite
        %   infinite sum of every squared step size needs to be smaller
        %   than infinity
        %   => stepSize needs to decrease with every recursion step
        %       this decreasing needs to be slow for the algorithm to
        %       converge. 
        stepSize;
        
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
        % state estimate of the current time step to calculate measurement
        % mean and covariance.
        stateEstimate = [];
        % Number of modes describing the Markov jump system.
        numModes = [];
   end

   methods (Access = public)
        function this = RKLUpdater(initialModeWeights, ...
                stepSize, measurementMatrix, systemMatrix, ...
                inputMatrix, stateCovariance, measurementCovariance)
            % Constructor of RKLUpdater.
            %
            % parameters:
            %   >> initialModeWeights: intial probabilities for the modes
            %   of the Markov jump system
            %   >> stepSize: update weight to determine the speed of the
            %   gradient descend.
            %   >> measurementMatrix: two-dimensional matrix modelling the
            %   relation of state and measurements
            %   >> systemMatrix: two-dimensional matrix modelling the
            %   relation of the current state and the next state
            %   >> inputMatrix: three-dimensional matrix describing the
            %   influence of control inputs in the next state.
            %   >> stateCovariance: two-dimensional matrix
            %   >> measurementCovariance: two-dimensional matrix
            % Returns:
            %   << this: object of type RKLUpdater
            this.currentModeWeights = initialModeWeights;
            this.numModes = length(initialModeWeights);
            
            % magic number taken over from examples of the authors.
            this.stepSize = stepSize;% 0.02; 
            
            this.measurementMatrix = measurementMatrix;
            this.systemMatrix = systemMatrix;
            this.inputMatrix = inputMatrix;
            this.stateCovariance = stateCovariance;
            this.measurementCovariance = measurementCovariance;
        end
       
        function newTransitionMatrix = updateTransitionProbabilityMatrix(...
            this, measurement, previousTransitionMatrix)
            % Updates transition probability matrix based on previous 
            % probability matrix and the current measurements. 
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
            estimations = previousTransitionMatrix;
            dimProbabilityMatrix = size( estimations, 1 );
           
            loglikelihoods = zeros(dimProbabilityMatrix);
            likelihoods = zeros(dimProbabilityMatrix);
            for i = 1:dimProbabilityMatrix
                for j = 1:dimProbabilityMatrix
                    % MeasurementMean is only calculated for i =
                    % dimProbabilityMatrix
                    covariance = this.currentMeasurementCovariance(:, :, i, j);
                    gaussian = Gaussian(this.currentMeasurementMean(:, i, j), ...
                        covariance );

                    loglikelihoods(i, j) = gaussian.logPdf( measurement );
                end
                
                maxLogValue = max( loglikelihoods(i, :) );
                likelihoods(i, :) = exp( loglikelihoods(i, :) - maxLogValue );
                likelihoods(i, :) = Updater.boundProbabilities(likelihoods(i, :));    
            end
            
            denominator = sum(likelihoods .* estimations, 'all');
            
            newTransitionMatrix = estimations + this.stepSize * likelihoods / denominator;
            
            % Projects TPM estimation to ensure satisfying the constraints
            newTransitionMatrix = RKLUpdater.project(newTransitionMatrix);
        end
        
        function initializeUpdateStep(this, stateMean, stateCovariance, modeWeights, ...
                inputs)    
            % Sets and calculates values needed for next update
            %
            % Parameters:
            %   >> stateMean: vector of doubles
            %   >> stateCovariance: matrix of doubles
            %   >> modeWeights: vector of doubles
            %   >> input: vector of doubles
            this.currentModeWeights = Updater.boundProbabilities(...
                modeWeights); 
            
            this.stateEstimate = stateMean;
            this.calculateMean(inputs);
            this.calculateCovariance(stateCovariance);
        end
        
        function calculateMean(this, inputs)
            % Calculation of measurement mean for the update of the
            % transition probability matrix.
            % 
            % Parameters:
            %   >> inputs: possible control inputs of the current time
            %   step
            lengthInputs = this.numModes;
            
            this.currentMeasurementMean = zeros(size(this.measurementMatrix, 1), ...
              lengthInputs, lengthInputs);
      
            if isempty(inputs)
                inputs = zeros(1, lengthInputs);
            end
            
            for i = 1:lengthInputs
                for j = 1:lengthInputs
                    A = this.systemMatrix * this.stateEstimate(:, i) ...
                        + this.inputMatrix * inputs(j);
                    mean = this.measurementMatrix * A;
                    this.currentMeasurementMean(:, i, j) = mean;
                end
            end
        end

        function calculateCovariance(this, ...
              covarianceStateEstimates)
            % Calculation of measurement covariance for computing the 
            % likelihood of the measurement at a specific time step
            %
            % Parameters:
            %   >> covarianceStateEstimates: two-dimensional matrix 
            lengthInputs = this.numModes;
            
            this.currentMeasurementCovariance = zeros(size(this.measurementMatrix, 1),...
              size(this.measurementMatrix, 1), lengthInputs, lengthInputs);
            
            for i = 1:lengthInputs
                for j = 1:lengthInputs            
                    A = this.systemMatrix;

                    covariance = this.measurementMatrix * ( ...
                        A * ...
                        covarianceStateEstimates(:,:,i) * A.' + ...
                        this.stateCovariance ) * this.measurementMatrix.' + ...
                        this.measurementCovariance;
                    this.currentMeasurementCovariance(:,:,i,j) = covariance;
                end
            end
        end
   end
   
   methods (Static)   
        function projection = project(matrix)
            % Projects input transition probability matrix (TPM) to ensure 
            % satisfying the contraints.
            % Returns matrix itself if constrains are already satisfied.
            %
            % Parameters:
            %   >> matrix (NxN matrix):
            %      Updated TPM which needs to be postprocessed to make it
            %      valid
            %
            % Returns:
            %   << projection (NxN matrix, projection(i, j) in [0,1] 
            %       for all 1 <= i, j <= dimension):
            %      projection of matrix to ensure that contrains of TPM are
            %      satisfied.
            projection = matrix;
            dimension = size( matrix, 1 );
            
            if ~RKLUpdater.isValid( projection( dimension, : ) )
                %[projection( dimension, : ), wasnan] = RKLUpdater.project_vector(...
                %    projection( dimension, : ), dimension, distribution );
                projection( dimension, : ) = RKLUpdater.projsplx(...
                    projection( dimension, : ) );

            end
                
            projection( dimension, : ) = ...
                Updater.boundProbabilities(projection( dimension, : ));
            
            projection( dimension, : ) = projection( dimension, : ) / ...
                sum( projection(dimension, :) );
            
            if dimension > 1
                projection( (dimension - 1), : ) = ...
                    projection( dimension, : );

                for i = 2:(dimension - 1)
                    index = dimension - i;
                    projection( index, : ) = ...
                        projection( ( index + 1 ), : );
                    projection( index, ( index + 1 ) ) = ...
                        projection( index, ( index + 1 ) ) ...
                        + projection( index, ( index + 2 ) );
                    projection( index, ( index + 2 ) ) = 0;
                end
            end
        end
        
        function x = projsplx(y)
            % project an n-dim vector y to the simplex Dn
            % Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
            % (c) Xiaojing Ye
            % xyex19@gmail.com
            %
            % Algorithm is explained as in the linked document
            % http://arxiv.org/abs/1101.6081
            % or
            % http://ufdc.ufl.edu/IR00000353/
            %
            % Jan. 14, 2011.
            m = length(y); bget = false;
            s = sort(y,'descend'); tmpsum = 0;
            for ii = 1:m-1
                tmpsum = tmpsum + s(ii);
                tmax = (tmpsum - 1)/ii;
                if tmax >= s(ii+1)
                    bget = true;
                    break;
                end
            end

            if ~bget, tmax = (tmpsum + s(m) -1)/m; 
            end
         x = max(y-tmax,0);
        end
        
        function [fun, der] = costFunctionAndDerivative(x,p)
            % Cost function and corresponding derivative to calculate the 
            % update for the next transition probability matrix.
            % 
            % Parameters:
            %   >> x: start value for optimization
            %   >> p: one-dimensional vector which is the last row of the
            %   new transition matrix  
            % Returns:
            %   << fun: cost function to calculate update of the transition
            %   probability matrix
            fun = norm(x-p);
            
            if nargout > 1
                der = (1 / 2) * x - x;
            end
        end
        
        function valid = isValid(modeTransitionVector)
            % Tests if every value of vector is in [0,1] and the sum of all
            % values equals one.
            %
            % Parameters:
            %   >> vector: one-dimensional vector
            %
            % Returns:
            %   << valid: either 1 or 0       
            valid = 1;
            if ~isequal(modeTransitionVector <= 1 & modeTransitionVector >= 0, ...
                        ones(size(modeTransitionVector))) ...
                    || ~isequal(round(sum(modeTransitionVector) * 1e8), ...
                        ones(size(modeTransitionVector, 1), 1) * 1e8) ...
                        
                    valid = 0;
            end
        end
   end
end