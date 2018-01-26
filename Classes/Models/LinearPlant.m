classdef LinearPlant < LinearSystemModel
    % Implemententation of a linear (potentially time-varying) plant with white and zero mean Gaussian
    % process noise, i.e. x_{k+1} = A_{k}x_{k} + B_{k}u_{k} + G_{k}w_{k}.
    % By default, the system noise matrix G_{k} is the identity matrix.
    % It can be changed, for instance, if the noise affects only a subset of the state variables.
    % Then, G_{k} should be such that G_{k}*G_{k}' equals the desired noise
    % covariance, and setNoise(..) should be called with the identity
    % matrix.
    %
    %
    % AUTHOR:       Maxim Dolgov, 12.04.2013
    % LAST UPDATE:  Maxim Dolgov, 12.04.2013
    %               Florian Rosenthal, 05.10.2017
    
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
        % system state dimension; (positive integer)
        dimState = [];
        % input matrix (B); (matrix of dimension <dimX> x <dimU>)
        inputMatrix = [];
    end

    properties (Access = private)
        % control input dimension (<dimU>); (positive integer)
        dimInput= [];
        % applied control input; (col vector of dimension <dimU>)
        input;
    end
    
    methods (Access = public)
        %% LinearPlant
        function this = LinearPlant(systemMatrix, inputMatrix, noiseCovMatrix)
            % Class constructor.
            %
            % Parameters:
            %   >> systemMatrix (Square matrix (<dimState> x <dimState>))
            %      The system matrix.
            %   >> inputMatrix (Matrix (<dimState> x <dimInput>))
            %      The input matrix.
            %   >> noiseCovMatrix (Positive definite matrix or vector)
            %      Covariance matrix (<dimState> x <dimState>) of the noise distribution.
            %      If a vector (<dimState>) is passed, 
            %      its values are interpreted as the variances of a diagonal covariance matrix.
            %
            % Returns:
            %   << this (LinearPlant)
            %      A new LinearPlant instance.
            
            this@LinearSystemModel(systemMatrix);
            this.dimState = size(this.sysMatrix, 1);
            this.setSystemInputMatrix(inputMatrix); 
            this.setNoise(noiseCovMatrix);
        end
        
        %% setNoise
        function setNoise(this, covariance)
            % Set the system noise to a zero mean Gaussian with the
            % specified covariance matrix.
            %
            %   >> covariance (Positive definite matrix or vector)
            %      Covariance matrix (<dimState> x <dimState>) of the noise distribution.
            %      If a vector (<dimState>) is passed, 
            %      its values are interpreted as the variances of a diagonal covariance matrix.
            
            setNoise@LinearSystemModel(this, Gaussian(zeros(this.dimState, 1), covariance));      
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
                if ~Checks.isColVec(sysInput, this.dimInput) || ~all(isfinite(sysInput))
                    error('LinearPlant:InvalidInput', ...
                        '** Control input to apply must be a real-valued %d-dimensional column vector  **', this.dimInput);
                end
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
                       
            Validator.validateInputMatrix(inputMatrix, this.dimState);
            this.inputMatrix = inputMatrix;
            this.dimInput = size(inputMatrix, 2);
         end
        
        %% setSystemMatrix
        function setSystemMatrix(this, sysMatrix)
            % Set the system matrix.
            %
            %
            % Parameters:
            %   >> sysMatrix (Square matrix)
            %      The new system matrix.
            
            Validator.validateSystemMatrix(sysMatrix);
            setSystemMatrix@LinearSystemModel(this, sysMatrix);
            this.dimState = size(this.sysMatrix, 1);
        end
   end % methods public
end % classdef