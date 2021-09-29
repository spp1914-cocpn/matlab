classdef LinearPlant < LinearSystemModel & matlab.mixin.Copyable
    % Implementation (based on Maxim Dolgov's original one) of a linear (potentially time-varying) plant with white and zero mean Gaussian
    % process noise, i.e. x_{k+1} = A_{k}x_{k} + B_{k}u_{k} + G_{k}w_{k}.
    % By default, the system noise matrix G_{k} is the identity matrix.
    % It can be changed, for instance, if the noise affects only a subset of the state variables.    
    %
    % This implementation is based on the original one by Jörg Fischer.
    %
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        % input matrix (B); (matrix of dimension <dimState> x <dimInput>)
        inputMatrix = [];
        % state constraints: state must be within [lb, ub]
        stateConstraints = []; % column wise arranged: [lb ub]
    end

    properties (SetAccess = immutable, GetAccess = public)
        % control input dimension (<dimInput>); (positive integer)
        dimInput= [];
        % system state dimension; (positive integer)
        dimState = [];
    end
    
    properties (Access = private)
        
        % applied control input; (col vector of dimension <dimInput>)
        input;
    end
    
    methods (Access = public)
        %% LinearPlant
        function this = LinearPlant(systemMatrix, inputMatrix, noiseCovMatrix, sysNoiseMatrix)
            % Class constructor.
            %
            % Parameters:
            %   >> systemMatrix (Square matrix (<dimState> x <dimState>))
            %      The system matrix.
            %   >> inputMatrix (Matrix (<dimState> x <dimInput>))
            %      The input matrix.
            %   >> noiseCovMatrix (Positive definite matrix or vector)
            %      Covariance matrix (<dimNoise> x <dimNoise>) of the noise distribution.
            %      If a vector (<dimNoise>) is passed, 
            %      its values are interpreted as the variances of a diagonal covariance matrix.
            %   >> sysNoiseMatrix (Matrix, optional)
            %      The system noise matrix G.
            %      If left out, or the empty matrix is passed, the identity
            %      matrix of dimension <dimState>  is assumed, so that <dimNoise> = <dimState>.            
            %
            % Returns:
            %   << this (LinearPlant)
            %      A new LinearPlant instance.
            
            if nargin == 3
                sysNoiseMatrix = [];
            end
            this@LinearSystemModel(systemMatrix, sysNoiseMatrix);
            this.dimState = size(this.sysMatrix, 1);
            
            Validator.validateInputMatrix(inputMatrix, this.dimState);
            this.inputMatrix = inputMatrix;
            this.dimInput = size(inputMatrix, 2);
            
            this.setNoise(noiseCovMatrix);
        end
        
        %% setNoise
        function setNoise(this, covariance)
            % Set the system noise to a zero mean Gaussian with the
            % specified covariance matrix.
            %
            %   >> covariance (Positive definite matrix or vector)
            %      Covariance matrix (<dimNoise> x <dimNoise>) of the noise distribution.
            %      If a vector (<dimNoise>) is passed, 
            %      its values are interpreted as the variances of a diagonal covariance matrix.
            %      Here, <dimNoise> indicates the dimension of the noise,
            %      which might be unequal to <dimState> in case a system
            %      noise matrix G is set.
            
            if isempty(this.sysNoiseMatrix)
                assert(Checks.isSquareMat(covariance, this.dimState) || Checks.isVec(covariance, this.dimState), ...
                    'LinearPlant:SetNoise:InvalidCovariance', ...
                    '** <covariance> must be a %d-by-%d matrix as no system noise matrix G is set', ...
                        this.dimState, this.dimState);
                if isempty(this.noise)
                    % no system noise matrix G -> identity matrix, no
                    % subspace noise
                    setNoise@LinearSystemModel(this, Gaussian(zeros(this.dimState, 1), covariance));
                else
                   this.noise.set(zeros(this.dimState, 1), covariance);
                end
            else
                % we check if there is a match to the sys noise matrix
                    dimNoise = size(this.sysNoiseMatrix, 2);
                    assert(Checks.isSquareMat(covariance, dimNoise) || Checks.isVec(covariance, dimNoise), ...
                        'LinearPlant:SetNoise:InvalidCovariance', ...
                        '** <covariance> must be a %d-by-%d matrix as system noise matrix is %d-by-%d', ...
                        dimNoise, dimNoise, this.dimState, dimNoise);
                if isempty(this.noise)                    
                    setNoise@LinearSystemModel(this, Gaussian(zeros(dimNoise, 1), covariance));
                else
                    this.noise.set(zeros(dimNoise, 1), covariance);
                end
            end
        end
        
        %% setSystemInput
        function setSystemInput(this, sysInput)
            % Set the system input vector.
            %
            % By default, the system input is an empty matrix.
            %
            % Parameters:
            %   >> sysInput (Column vector or empty matrix)
            %      The new system input vector.
            %      An empty matrix means no input vector.
            
            newInput = [];
            if ~isempty(sysInput)
                assert(Checks.isColVec(sysInput, this.dimInput) && all(isfinite(sysInput)), ...
                    'LinearPlant:InvalidInput', ...
                    '** Control input to apply must be a real-valued %d-dimensional column vector  **', ...
                    this.dimInput);
                
                newInput = this.inputMatrix * sysInput;
            end
            setSystemInput@LinearSystemModel(this, newInput);
            this.input = sysInput;
        end
        
        %% getSystemInput
        function input = getSystemInput(this)
            % Get the system input vector.
            %
            % Returns:
            %   << input (Column vector of dimension <dimU> or empty matrix)
            %      The currently applied input vector.
            %      An empty matrix means no input vector.
            input = this.input;
        end
        
        %% setSystemInputMatrix
        function setSystemInputMatrix(this, inputMatrix)
            % Set the system input matrix (B).
            %
            % Parameters:
            %   >> inputMatrix (Matrix with <dimX> rows)
            %      The new system input matrix (B).
            
            Validator.validateInputMatrix(inputMatrix, this.dimState, this.dimInput);
            this.inputMatrix = inputMatrix;
            
            % side effect in case system input is present
            this.setSystemInput(this.input);
         end
        
        %% setSystemMatrix
        function setSystemMatrix(this, sysMatrix)
            % Set the system matrix.
            %
            %
            % Parameters:
            %   >> sysMatrix (Square matrix)
            %      The new system matrix.
            
            Validator.validateSystemMatrix(sysMatrix, this.dimState);
            setSystemMatrix@LinearSystemModel(this, sysMatrix);            
        end
        
        %% setStateConstraints
        function setStateConstraints(this, lowerBound, upperBound)
             assert(Checks.isVec(lowerBound, this.dimState) && all(~isnan(lowerBound)), ...
                'LinearPlant:SetStateConstraints:InvalidLowerBound', ...
                '** <state> must be a %d-dimensional vector (-inf/inf allowed) **', this.dimState);
            assert(Checks.isVec(upperBound, this.dimState) && all(~isnan(lowerBound)), ...
                'LinearPlant:SetStateConstraints:InvalidUpperBound', ...
                '** <state> must be a %d-dimensional vector (-inf/inf allowed) **', this.dimState);
            
            this.stateConstraints = [lowerBound(:) upperBound(:)];
        end
        
        %% isValidState
        function isValid = isValidState(this, state)
            % Function to check whether a given plant state is valid (e.g.,
            % does not violate constraints).
            %
            % Parameters:
            %   >> state (Column vector)
            %      The system state to check.
            %
            % Returns:
            %   << isValid (Flag, i.e., boolean)
            %      Flag to indicate whether the given state is valid (e.g.,
            %      admissible with regards to contraints).
            
            assert(Checks.isColVec(state, this.dimState) && all(isfinite(state)), ...
                'LinearPlant:IsValidState:InvalidSystemState', ...
                '** <state> must be a real-valued %d-dimensional column vector  **', this.dimState);
            
            isValid = true;
            if ~isempty(this.stateConstraints)
                isValid = all(this.stateConstraints(:, 1) <= state) && all(state <= this.stateConstraints(:, 2));
            end
        end
    end
    
    methods (Access = protected)
        %% copyElement
        function copyObj = copyElement(this)
            % noise cannot be empty
            [~,noiseCov] = this.noise.getMeanAndCov();
            copyObj = LinearPlant(this.sysMatrix, this.inputMatrix, noiseCov, this.sysNoiseMatrix);
            copyObj.setSystemInput(this.input); % this.input can be empty
            copyObj.stateConstraints = this.stateConstraints; % they can be copied
        end
    end
end
