classdef Validator
    % Utility class containing convenience functions for parameter
    % validation.
        
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
    
    methods (Access = private)
        function this = Validator()
        end
    end
    
    methods (Static, Access = public)
        
        %% validateSystemMatrix
        function validateSystemMatrix(A, dimX)
            % Convenience function to validate the system matrix A.
            % In particular, this function errors if the matrix is not
            % square (and of appropriate dimension) or contains NaN or inf values.
            % 
            %
            % Parameters:
            %  >> A (Square Matrix)
            %     The system matrix of the plant.
            %
            %   >> dimX (Positive integer, Optional)
            %      The dimension of the plant's state.
            %      If a value is specified here, the dimensions of the given matrix are also checked. 
            %
            if nargin == 1
                assert(Checks.isSquareMat(A) && all(isfinite(A(:))), ...
                    'Validator:ValidateSystemMatrix:InvalidMatrix', ...
                    '** System matrix <A> must be square and real-valued **');
            elseif nargin == 2
                assert(Checks.isSquareMat(A, dimX) && all(isfinite(A(:))), ...
                    'Validator:ValidateSystemMatrix:InvalidDimensions', ...
                    '** System matrix <A> must be square (%d-by-%d) and real-valued **', dimX, dimX);
            end
        end
        
        %% validateInputMatrix
        function validateInputMatrix(B, dimX)
            % Convenience function to validate the input matrix B.
            % In particular, this function errors if the matrix is not
            % of appropriate dimensions (i.e., does not possess <dimX> rows) or contains NaN or inf values.
            % 
            %
            % Parameters:
            %   >> B (Matrix)
            %      The input matrix of the plant.
            %
            %   >> dimX (Positive integer)
            %      The dimension of the plant's state.
            %
            assert(Checks.isFixedRowMat(B, dimX) && all(isfinite(B(:))), ...
                'Validator:ValidateInputMatrix:InvalidInputMatrix', ...
                '** Input matrix <B> must be a real-valued matrix with %d rows **', dimX); 
        end
        
        %% validateMeasurementMatrix
        function validateMeasurementMatrix(C, dimX, dimY)
            % Convenience function to validate the measurement matrix C.
            % In particular, this function errors if the matrix is not
            % of appropriate dimensions (i.e., does not possess <dimX> cols, 
            % or, if dimY is specified, is not <dimY>-by-<dimX>) or contains NaN or inf values.
            % 
            % Parameters:
            %   >> C (Matrix)
            %      The measurement matrix of the linear sensor.
            %
            %   >> dimX (Positive integer)
            %      The dimension of the plant's state.
            % 
            %   >> dimY (Positive integer, Optional)
            %      The dimension of the measurements.
            %
            if nargin == 2
                assert(Checks.isFixedColMat(C, dimX) && all(isfinite(C(:))), ...
                    'Validator:ValidateMeasurementMatrix:InvalidMeasMatrix', ...
                    '** Measurement matrix <C> must be a real-valued matrix with %d cols **', dimX); 
            elseif nargin == 3
                assert(Checks.isMat(C, dimY, dimX) && all(isfinite(C(:))), ...
                    'Validator:ValidateMeasurementMatrix:InvalidMeasMatrixDims', ...
                    '** Measurement matrix <C> is expected to be real-valued %d-by-%d matrix**', dimY, dimX);
            end
        end
        
        %% validateMeasNoiseCovarianceMatrix
        function V_sqrt = validateMeasNoiseCovarianceMatrix(V, dimY)
            % Convenience function to validate the covariance matrix V of the measurement noise.
            % In particular, this function errors if the matrix is not
            % positive definite, or of appropriate dimensions (i.e., if dimY is specified, is not <dimY>-by-<dimY>).
            % 
            % Parameters:
            %   >> V (Matrix)
            %      The covariance matrix of the measurement noise.
            %
            % 
            %   >> dimY (Positive integer, Optional)
            %      The dimension of the measurements.
            %
            % Returns:
            %   >> V_sqrt (Square matrix, optional) 
            %      The lower Cholesky factor of V.
            %
            if nargin == 2
                [isCov, tmp] = Checks.isCov(V, dimY);
                assert(isCov, 'Validator:ValidateMeasNoiseCovarianceMatrix:InvalidCovDim', ...
                    '**  <V> (measurement noise covariance) must be a real-valued, positive definite %d-by-%d matrix **', dimY, dimY);                
            else
                [isCov, tmp] = Checks.isCov(V);
                assert(isCov, 'Validator:ValidateMeasNoiseCovarianceMatrix:InvalidCov', ...
                    '**  <V> (measurement noise covariance) must be a real-valued, positive definite matrix **');                
            end
             
            if nargout == 1
                V_sqrt = tmp;
            end
        end
        
        %% validateSysNoiseCovarianceMatrix
        function W_sqrt = validateSysNoiseCovarianceMatrix(W, dimW)
            % Convenience function to validate the covariance matrix W of the plant noise.
            % In particular, this function errors if the matrix is not
            % positive definite, or of appropriate dimensions (i.e., if dimW is specified, is not <dimW>-by-<dimW>).
            % 
            % Parameters:
            %   >> W (Matrix)
            %      The covariance matrix of the plant noise.
            %
            % 
            %   >> dimW (Positive integer, Optional)
            %      The dimension of the noise, which is usually the dimension of the plant state.
            %
            % Returns:
            %   >> W_sqrt (Square matrix, optional) 
            %      The lower Cholesky factor of W.
            %
            if nargin == 2
               [isCov, tmp] = Checks.isCov(W, dimW);
               assert(isCov, 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCovDim', ...
                    '**  <W> (process noise covariance) must be a real-valued, positive definite %d-by-%d matrix **', dimW, dimW);
            else
                [isCov, tmp] = Checks.isCov(W);
                assert(isCov, 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCov', ...
                    '**  <W> (process noise covariance) must be a real-valued, positive definite matrix **');
            end
            
            if nargout == 1
                W_sqrt = tmp;
            end
        end
        
        %% validateCostMatrices
        function validateCostMatrices(Q, R, dimX, dimU)
            % Convenience function to validate the weighting matrices Q and
            % R of a quadratic cost function.
            % In particular, this function errors if the state weighting
            % matrix Q is of wrong dimensions and not positive
            % semi-definite, or the input weighting matrix R is not
            % positive definite.
            % 
            %
            % Parameters:
            %   >> Q (Square matrix, dimX-by-dimX)
            %      The time-invariant state or plant output weighting matrix.
            %
            %   >> R (Square Matrix, dimU-by-dimU)
            %      The time-invariant input weighting matrix.
            %
            %   >> dimX (Positive integer)
            %      The dimension of the plant's state.
            %
            %   >> dimU (Positive integer)
            %      The dimension of the inputs applied to the plant.
            
            % state weighting matrix, Q, must be positive-semidefinite
            if ~Checks.isSquareMat(Q, dimX) || any(~isfinite(Q(:)))
                error('Validator:ValidateCostMatrices:InvalidQMatrix', ...
                    '** Cost matrix <Q> must be a real-valued %d-by-%d matrix **', dimX, dimX);
            elseif ~issymmetric(Q) || any(eig(Q) < 0)
                error('Validator:ValidateCostMatrices:InvalidQMatrixPSD', ...
                    '** Cost matrix <Q> must be positive semi-definite, i.e., symmetric with nonnegative eigenvalues **');
            end
            
            % input weighting matrix, R, must be positive-definite
            assert(Checks.isCov(R, dimU) && all(isfinite(R(:))), ...
                'Validator:ValidateCostMatrices:InvalidRMatrix', ...
                '** Cost matrix <R> must be a real-valued, symmetric and invertible %d-by-%d matrix **', dimU, dimU);
        end
        
        %% validateDiscreteProbabilityDistribution           
        function validateDiscreteProbabilityDistribution(probDist, numElements)
            % Convenience function to validate a vector defining a discrete probability distribution.
            % In particular, this function errors if entries of the vector
            % are not nonnegative or do not sum up to 1.
            % 
            % Parameters:
            %   >> probDist (Vector)
            %      The vector to describe a probability distribution.
            %
            % 
            %   >> numElements (Positive integer, Optional)
            %      The expected number of elements of the given vector.
            %
            if nargin == 1
                assert(Checks.isNonNegativeVec(probDist) && all(isfinite(probDist)) && round(sum(probDist * 1e8)) == 1e8, ...
                    'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs', ...
                    '** <probDist> must be a nonnegative vector whose elements sum up to 1 **');
            else
                assert(Checks.isNonNegativeVec(probDist, numElements) && all(isfinite(probDist)) && round(sum(probDist * 1e8)) == 1e8, ...
                    'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbsNum', ...
                    '** <probDist> must be a nonnegative vector with %d elements summing up to 1 **', numElements);
            end
        end
        
        %% validateTransitionMatrix
        function validateTransitionMatrix(modeTransitionMatrix, numModes)
            % Convenience function to validate whether a given matrix is valid transition matrix of a Markov chain.
            % In particular, this function errors if the given matrix is not square, any entry is not in
            % [0, 1], or any row does not sum up to 1.
            % 
            % Parameters:
            %   >> modeTransitionMatrix (Square matrix)
            %      The mode transition matrix to check.
            %
            % 
            %   >> numModes (Positive integer, Optional)
            %      The number of modes of the corresponding Markov chain.
            %
            if nargin == 1
                assert(Checks.isSquareMat(modeTransitionMatrix) ...
                    && isequal(modeTransitionMatrix <= 1 & modeTransitionMatrix >= 0, ones(size(modeTransitionMatrix))) ...
                    && isequal(round(sum(modeTransitionMatrix, 2) * 1e8), ones(size(modeTransitionMatrix, 1), 1) * 1e8), ...
                    'Validator:ValidateTransitionMatrix:InvalidTransitionMatrix', ...
                    ['** <modeTransitionMatrix> must be square and stochastic, i.e. each entry must be in [0,1] ' ...
                    'and all rows must sum up to 1 **']);
            else
                assert(Checks.isSquareMat(modeTransitionMatrix, numModes) ...
                    && isequal(modeTransitionMatrix <= 1 & modeTransitionMatrix >= 0, ones(numModes)) ...
                    && isequal(round(sum(modeTransitionMatrix, 2) * 1e8), ones(numModes, 1) * 1e8), ...
                    'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim', ...
                    ['** <modeTransitionMatrix> must be %d-dimensional and stochastic, i.e. each entry must be in [0,1] ' ...
                    'and all rows must sum up to 1 **'], numModes); 
            end           
        end
        
        %% validateHorizonLength
        function validateHorizonLength(horizonLength)
            assert(Checks.isPosScalar(horizonLength) && mod(horizonLength, 1) == 0, ...
                'Validator:ValidateHorizonLength:InvalidHorizonLength', ...
                '** <horizonLength> must be a positive integer **');
        end
        
        %% validateSequenceLength
        function validateSequenceLength(sequenceLength)
            assert(Checks.isPosScalar(sequenceLength) && mod(sequenceLength, 1) == 0, ...
                'Validator:ValidateSequenceLength:InvalidSequenceLength', ...
                '** Input parameter <sequenceLength> (control sequence length) must be a positive integer **');
        end
    end
end

