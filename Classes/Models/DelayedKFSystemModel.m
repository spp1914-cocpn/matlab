classdef DelayedKFSystemModel < LinearSystemModel
    % Implementation of the augmented system model
    % to cope with delayed/lost measurements and uncertain/delayed inputs which is used by the Kalman
    % filter developed by Moayedia et al. for networked control systems.
    
    % Literature: 
    %   Maryam Moayedia, Yung Kuan Foo and Yeng Chai Soha
    %   Filtering for networked control systems with single/multiple measurement packets
    %   subject to multiple-step measurement delays and multiple packet dropouts
    %   International Journal of Systems Science (2011).
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2018  Florian Rosenthal <florian.rosenthal@kit.edu>
    %
    %                        Institute for Anthropomatics and Robotics
    %                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
    %                        Karlsruhe Institute of Technology (KIT), Germany
    %
    %                        http://isas.uka.de
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
    
    properties (SetAccess = immutable, GetAccess = private)
        weights;
        numModes;
        noiseMean;
        noiseCov;
        plantInputMatrix;
    end
      
    properties (Dependent, Access = private)
        dimInput;
    end
    
   methods
       function dimU = get.dimInput(this)
           dimU = size(this.plantInputMatrix, 2);
       end
   end
   
    methods (Access = public)
        %% DelayedKFSystemModel
        function this = DelayedKFSystemModel(systemMatrix, inputMatrix, systemNoise, numModes, ...
                maxMeasDelay, delayWeights)
            % Class constructor.
            %
            % Parameters:
            %   >> systemMatrix (Square matrix)
            %      The system matrix of the underlying linear plant.
            %
            %   >> inputMatrix (Matrix)
            %      The plant input matrix.
            %
            %   >> systemNoise (Subclass of Distribution)
            %      The noise affecting the plant.
            %
            %   >> numModes (Positive integer)
            %      The number of modes, i.e., the number of possibly active control inputs (including a default input) per time step.   
            %
            %   >> maxMeasDelay (Nonnegative integer)
            %      The maximum allowed measurement delay this instance can cope with.
            %
            %   >> delayWeights (Nonnegative Vector)
            %      A vector where the i-th elements denotes the probability
            %      of the i-th possible input to be applied.
            %          
            % Returns:
            %   << this (DelayedKFSystemModel)
            %      A new DelayedKFSystemModel instance.
            
                     
            Validator.validateSystemMatrix(systemMatrix);
            dimState = size(systemMatrix, 1);
            
            Validator.validateInputMatrix(inputMatrix, dimState);
          
            assert(Checks.isNonNegativeScalar(maxMeasDelay) && mod(maxMeasDelay, 1) == 0, ...
                'DelayedKFSystemModel:InvalidMaxMeasDelay', ...
                '** Maximum measurement delay <maxMeasDelay> must be a nonnegative integer **');
            
            assert(Checks.isClass(systemNoise, 'Distribution'), ...
                'DelayedKFSystemModel:InvalidSystemNoise', ...
                '** System noise <systemNoise> must be a subclass of Distribution **');
            
            assert(Checks.isPosScalar(numModes) && mod(numModes, 1) == 0, ...
                'DelayedKFSystemModel:InvalidNumModes', ...
                '** Number of modes (i.e., the number of possible inputs) must be a positive integer **');
                        
            Validator.validateDiscreteProbabilityDistribution(delayWeights, numModes);
        
            zeroMatrix = zeros(dimState, dimState * maxMeasDelay);
            % augmented system noise matrix
            augmentedNoiseMatrix = [eye(dimState); zeroMatrix'];            
            
            this = this@LinearSystemModel([systemMatrix, zeroMatrix; ...
                Utils.blockDiag(eye(dimState), maxMeasDelay), zeroMatrix'], ...
               ... % augmented system noise matrix
                augmentedNoiseMatrix);
            
            this.plantInputMatrix = inputMatrix;
            % moment matching to get Gaussian
            [this.noiseMean, this.noiseCov] = systemNoise.getMeanAndCov();
            this.resetNoise();
                       
            this.numModes = numModes;
            % store weights as row vector, i.e., column-wise arranged
            this.weights = reshape(delayWeights, 1, numModes);
        end
        
        %% reset
        function reset(this)
            % Resets the model by removing the additional system noise due
            % to uncertain inputs.
            this.resetNoise();
        end
        
        %% setSystemMatrix
        function setSystemMatrix(this, sysMatrix)
            % Set the system matrix.
            %
            % By default, the system matrix is an empty matrix
            % (i.e., an identity matrix of appropriate dimensions).
            %
            % Parameters:
            %   >> sysMatrix (Square matrix or empty matrix)
            %      The new system matrix.
            %      An empty matrix means the identity matrix of appropriate dimensions.
            
            dimState = size(this.plantInputMatrix, 1);
            newSysMatrix = this.sysMatrix;
            if isempty(sysMatrix)
               newSysMatrix(1:dimState,1:dimState) = eye(dimState);
            else
                Validator.validateSystemMatrix(sysMatrix, dimState);
                newSysMatrix(1:dimState,1:dimState) = sysMatrix;
            end

            setSystemMatrix@LinearSystemModel(this, newSysMatrix);
        end
        
       function setSystemNoiseMatrix(this, sysNoiseMatrix)
            % Set the system noise matrix.
            %
            % This is not the system noise covariance matrix!
            %
            % By default, the noise matrix is an empty matrix
            % (i.e., an identity matrix of appropriate dimensions).
            %
            % Parameters:
            %   >> sysNoiseMatrix (Matrix or empty matrix)
            %      The new system noise matrix.
            %      An empty matrix means the identity matrix of appropriate dimensions.
            
            dimState = size(this.plantInputMatrix, 1);
            newSysNoiseMatrix = this.sysNoiseMatrix;
            if isempty(sysNoiseMatrix)
                newSysNoiseMatrix(1:dimState, :) = eye(dimState);
            else
                assert(Checks.isSquareMat(sysNoiseMatrix, dimState) && all(isfinite(sysNoiseMatrix(:))), ...
                    'DelayedKFSystemModel:SetSystemNoiseMatrix:InvalidDimensions', ...
                    '** System noise matrix <sysNoiseMatrix> must be square (%d-by-%d) and real-valued **', ...
                    dimState, dimState);
                
                newSysNoiseMatrix(1:dimState, :) = sysNoiseMatrix;
            end
            setSystemNoiseMatrix@LinearSystemModel(this, newSysNoiseMatrix);
        end
        
        %% setSystemInput
        function setSystemInput(this, possibleInputs)
            % Set the system inputs that might be currently applied to the
            % plant, adding additional uncertainty to the system.
            %
            % By default, the system input is an empty matrix.
            %
            % Parameters:
            %   >> sysInput (Matrix)
            %      A matrix with column-wise arranged inputs that might be
            %      currently applied at the plant. 
            %      An empty matrix means no input vector, and, hence, no
            %      addtional uncertainty.
            if isempty(possibleInputs)
                if ~isempty(this.noiseMean) % should only be empty if setSystemInput is called in ctor of base class
                    this.resetNoise();
                end
                return;
            end
            assert(Checks.isMat(possibleInputs, this.dimInput, this.numModes) && all(isfinite(possibleInputs(:))), ...
                'DelayedKFSystemModel:InvalidInput', ...
                '** <possibleInputs> must be given as a real-valued %d-by-%d matrix (i.e., the possibly active inputs must be column-wise arranged) **',  ...
                this.dimInput, this.numModes);
            
            % last element is considered the default input (e.g., zero) applied by actuator in
            % case of packet loss
            % compute expected, but uncertain input
            [inputMean, inputCov] = Utils.getMeanAndCov(possibleInputs, this.weights);
            virtualInput = this.plantInputMatrix * inputMean;
            virtualInputCov = this.plantInputMatrix * inputCov * this.plantInputMatrix';
            this.noise.set(this.noiseMean + virtualInput, this.noiseCov + virtualInputCov);
        end
    end
     
    methods (Access = private)
        function resetNoise(this)
            if isempty(this.noise)
                this.setNoise(Gaussian(this.noiseMean, this.noiseCov));
            else
                this.noise.set(this.noiseMean, this.noiseCov);
            end
        end
    end
end

