classdef Utility < handle
    % This class provides various utility functions.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
    %                             Fabio Broghammer <fabio.broghammer@student.kit.edu>
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
        function this = Utility()
            
        end
    end
    
    methods (Static, Access = public)    
        %% createBlockdiagonalMatrix
        function diagMatrix = createBlockdiagonalMatrix(matrix)
            % Creates a block diagonal matrix out of a 3D-Matrix.
            %
            % Parameters:
            %   >> matrix (3D-Matrix, m-by-n-by-k)
            %      The 3D-Matrix.
            %
            % Returns:
            %   << diagMatrix (Matrix of shape m*k-by-n*k)
            %      A block diagonal matrix with <matrix> (M) on the diagonal:
            %      |M(:, :, 1)      0           0       ...     0     |
            %      |    0       M(:,:, 2)               ...     0     |
            %      |    0           0       M(:, :, 3)  ...     0     |
            %      |    .           .           .        .      0     |
            %      |    0           0           0       ... M(:, :, k)|

            [rows, columns, numberOfMatrices] = size(matrix);            

            diagMatrix = zeros(rows * numberOfMatrices, columns * numberOfMatrices);

            for i = 1:numberOfMatrices
                diagMatrix((i - 1) * rows + 1: i * rows, (i - 1) * columns + 1:i * columns) = matrix(:, :, i);
            end
        end
        
        %% computeStageCosts
        function stageCosts = computeStageCosts(state, appliedInput, Q, R)
            % Compute the stage-costs of a quadratic cost function,
            % i.e., compute x'*Q*x + u'*R*u for a given pair of state x and input u.
            %
            % Parameters:
            %
            %   >> state (Column vector)
            %      The state vector at the current stage.
            %  
            %   >> appliedInput (Column vector)
            %      The input at the current stage.
            %
            %   >> Q (positive-semidefinite matrix)
            %      The state weighting matrix, expected to be at least
            %      positive-semidefinite and of appropriate dimension.
            %
            %   >> R (positive-definite matrix)
            %      The input weighting matrix, expected to be 
            %      positive-definite and of appropriate dimension.
            %
            % Returns:
            %   << stageCosts (Nonnegative scalar)
            %      The stage-costs, coresponding to the given parameters.
            
            stageCosts = state' * Q * state + appliedInput' * R * appliedInput;
        end
        
        %% computeLQGCosts
        function costs = computeLQGCosts(horizonLength, stateTrajectories, appliedInputs, Q, R)
            % Evaluate a quadratic cost function for a finite horizon N,
            % i.e., compute x_0'*Q*x_0 + u_0'*R*u_0 + ... + x_{N-1}'*Q*x_{N-1} + u_{N-1}'*R*u_{N-1} 
            % + x_N'*Q*x_N, for given pairs of state and input
            % trajectories.
            %
            % Parameters:
            %   >> horizon length (Positive integer)
            %      The length N of the horizon to be considered.
            %
            %   >> stateTrajectories (3D matrix)
            %      A 3D matrix with the different state trajectories     
            %      slice-wise arranged. Each state trajectory is expected
            %      to be column-wise arranged and to possess at least N+1 elements.
            %  
            %   >> appliedInputs (3D matrix)
            %      A 3D matrix with the different input trajectories     
            %      slice-wise arranged. Each input trajectory is expected
            %      to be column-wise arranged and to possess at least N elements.
            %
            %   >> Q (positive-semidefinite matrix)
            %      The state weighting matrix, expected to be at least
            %      positive-semidefinite and of appropriate dimension.
            %
            %   >> R (positive-definite matrix)
            %      The input weighting matrix, expected to be 
            %      positive-definite and of appropriate dimension.
            %
            % Returns:
            %   << costs (Row vector)
            %      The value of the cost function, evaluated for each pair
            %      of states and inputs, column-wise arranged.
                        
            % stateTrajectories: dimState x horizonLength+1 x n  (true states or
            % estimates by controller)
            % applied inputs: dimInput x horizonLength x n
            % n: number of controller algorithms to compare
            % Q: psd weighting matrix (dimState x dimState)
            % R: pd weighting matrix (dimInput x dimInput)

            % compute x'*Q*x + u'*R*u for every time step and sum up
            s1 = size(stateTrajectories);
            if numel(s1) == 2
                numControllers = 1;
            else
                numControllers = s1(3);
            end
            if ndims(appliedInputs) == ndims(stateTrajectories) && size(appliedInputs, 2) == s1(2) - 1 ...
                    && s1(2) == horizonLength + 1
                costs = zeros(1, numControllers);
                Rsqrt = chol(R); % upper Cholesky factor of R, R = Rsqrt' * Rsqrt
                for i=1:numControllers
                    costs(i) = sum(sum(stateTrajectories(:,:, i) .* (Q * stateTrajectories(:,:, i)))) ...
                        + sum(sum((Rsqrt * appliedInputs(:,:, i)).^ 2));
                end              
            else
                costs = [];
            end
        end
        
        %% createActuatorMatrices
        function [F, G, H, J] = createActuatorMatrices(controlSequenceLength, dimPlantInput)
            % Create the matrices governing the dynamics of eta_k in the augmented state space describing the
            % buffering procedure employed at the actuator.
            %
            % Literature: 
            %   Jörg Fischer, Achim Hekler and Uwe D. Hanebeck,
            %   State Estimation in Networked Control Systems,
            %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
            %   Singapore, July 2012.
            %
            %   Jörg Fischer, Achim Hekler, Maxim Dolgov and Uwe D. Hanebeck,
            %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject
            %   to Random Transmission Delays and Packet Losses,
            %   Proceedings of the 2013 American Control Conference (ACC 2013), 
            %   Washington D. C., USA, June 2013.
            %
            % Parameters:
            %   >> controlSequenceLength (Positive integer)
            %      The length of a control sequence (without default input)
            %      generated by the controller.
            %
            %   >> dimPlantInput (Positive integer)
            %      The dimension of a single control input.
            %
            % Returns:
            %   << F (Square matrix)
            %      System matrix of the dynamics of the network-actuator
            %      subsystem, which governs the control input selection at the actuator.
            %
            %   << G (Matrix)
            %      Input matrix of the network-actuator subsystem, that maps the entries of the current control sequence
            %      which are meant to be applied to the plant at future time steps to the
            %      vector containing possible future control inputs.
            %
            %   << H (3D matrix, one for each mode of the MJLS)
            %      H selects a control input from the vector of possible control inputs according to the state of the Markov
            %      chain which governs the packet arrival process at the actuator.
            %
            %   << J (3D matrix, one for each mode of the MJLS)
            %      J selects the first control input from the most recent control input sequence 
            %      according to the state of the Markov chain which governs the packet arrival at the actuator.
            
            if controlSequenceLength == 1
                F = [];
                G = zeros(0, dimPlantInput); % so that caller has not to perform a check for length of sequence
                H = zeros(dimPlantInput, 0, 2); % so that caller has not to perform a check
                J = cat(3, eye(dimPlantInput), zeros(dimPlantInput));                
            else
                dimInputSequenceState = dimPlantInput * (controlSequenceLength * (controlSequenceLength - 1) / 2); % the 'input sequence part' of the controller state;
                % create F
                F = zeros(dimInputSequenceState);
                idx = [1:controlSequenceLength-2]; %#ok
                dims = idx * dimPlantInput;
                sums = cumsum(dims);
                for i=idx 
                    rowIndex = dimInputSequenceState - sums(i) + 1; 
                    colIndex = (controlSequenceLength - i - 1) * dimPlantInput + (sums(end) - sums(i)) + 1;
                    F(rowIndex:rowIndex + dims(i) - 1, colIndex:colIndex + dims(i) - 1) = speye(dims(i));
                end
    
                % create G
                dimEye = dimPlantInput * (controlSequenceLength - 1);
                G = zeros(dimEye * controlSequenceLength / 2, dimPlantInput * controlSequenceLength);
                G(1:dimEye, (dimPlantInput + 1):end) = speye(dimEye);

                % create H matrices
                numCols = dimInputSequenceState;
                numRows = dimPlantInput;
                H = zeros(numRows, numCols, controlSequenceLength + 1);
                sums = dimPlantInput * cumsum(1:controlSequenceLength -1);
                % for mode 1 (no delay): matrix remains zero
                % for mode N + 1 (packetDelay of N): matrix remains zero
                for i = (controlSequenceLength - 1):-1:1
                    colIdx = numCols - sums(controlSequenceLength - i) + 1;
                    H(:, colIdx:colIdx + dimPlantInput -1, i + 1) = speye(dimPlantInput);
                end

                % create J matrices, zeros for all but first mode
                J = zeros(dimPlantInput, dimPlantInput * controlSequenceLength, controlSequenceLength + 1);
                J(:, 1:dimPlantInput, 1) = speye(dimPlantInput);
            end
        end
        
        %% createAugmentedPlantModel
        function [F, G, H, J, augA, augB] ...
            = createAugmentedPlantModel(controlSequenceLength, A, B)
            % Create the matrices of the augmented state space which describes the
            % linear NCS (including the network-actuator subsystem), where
            % a sequence-base controller and zero default input are
            % employed, in terms of an MJLS.
            %
            % Literature: 
            %   Jörg Fischer, Achim Hekler and Uwe D. Hanebeck,
            %   State Estimation in Networked Control Systems,
            %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
            %   Singapore, July 2012.
            %
            %   Jörg Fischer, Achim Hekler, Maxim Dolgov and Uwe D. Hanebeck,
            %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject
            %   to Random Transmission Delays and Packet Losses,
            %   Proceedings of the 2013 American Control Conference (ACC 2013), 
            %   Washington D. C., USA, June 2013.
            %
            % Parameters:
            %   >> controlSequenceLength (Positive integer)
            %      The length of a control sequence (without default input)
            %      generated by the controller.
            %
            %   >> A (Square matrix, dimPlantState-by-dimPlantState)
            %      The time-invariant system matrix of the plant dynamics.
            %
            %   >> B (Matrix, dimPlantState-by-dimPlantInput)
            %      The time-invariant input matrix of the plant dynamics.
            %
            % Returns:
            %   << F (Square matrix)
            %      System matrix of the dynamics of the network-actuator
            %      subsystem, which governs the control input selection at the actuator.
            %
            %   << G (Matrix)
            %      Input matrix of the network-actuator subsystem, that maps the entries of the current control sequence
            %      which are meant to be applied to the plant at future time steps to the
            %      vector containing possible future control inputs.
            %
            %   << H (3D matrix, one for each mode of the MJLS)
            %      H selects a control input from the vector of possible control inputs according to the state of the Markov
            %      chain which governs the packet arrival process at the actuator.
            %
            %   << J (3D matrix, one for each mode of the MJLS)
            %      J selects the first control input from the most recent control input sequence 
            %      according to the state of the Markov chain which governs the packet arrival at the actuator.
            %
            %   << augA (3D matrix, one for each mode of the MJLS (Optional output parameter))
            %      System matrix of the augmented open-loop dynamics, one per
            %      mode.
            %
            %   << augB (3D matrix, one for each mode of the MJLS (Optional output parameter))
            %      Input matrix of the augmented open-loop dynamics, one per
            %      mode.

            dimPlantState = size(A, 1);
            dimPlantInput = size(B, 2);
            
            [F, G, H, J] = Utility.createActuatorMatrices(controlSequenceLength, dimPlantInput);
            
            if nargout == 6
                if controlSequenceLength == 1
                    augA = repmat(A, 1, 1, 2);
                    augB = cat(3, B, zeros(dimPlantState, dimPlantInput));
                else
                    augA = repmat(blkdiag(A, F), 1, 1, controlSequenceLength + 1);
                    augB = repmat([zeros(dimPlantState, size(G,2)); G], 1, 1, controlSequenceLength + 1);
                    % upper part of augB is only non-zero for the first mode
                    augB(1:dimPlantState, 1:dimPlantInput, 1) = B;
                    % augA mode-dependent, but H is zero for first and last mode
                    augA(1:dimPlantState, dimPlantState + 1:end, 2:controlSequenceLength) = mtimesx(B, H(:, :, 2:controlSequenceLength));
                end
            end            
        end
        
        %% createAugmentedCostModel
        function [augQ, augR] = createAugmentedCostModel(controlSequenceLength, Q, R, Z)
            % Create augmented cost matrices to express the quadratic cost
            % function in terms of the augmented state space.
            %
            % Literature: 
            %   Jörg Fischer,
            %   Optimal sequence-based control of networked linear systems,
            %   Karlsruhe series on intelligent sensor-actuator-systems, Volume 15,
            %   KIT Scientific Publishing, 2015.
            %
            %   Jörg Fischer, Achim Hekler, Maxim Dolgov and Uwe D. Hanebeck,
            %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject
            %   to Random Transmission Delays and Packet Losses,
            %   Proceedings of the 2013 American Control Conference (ACC 2013), 
            %   Washington D. C., USA, June 2013.
            %
            % Parameters:
            %   >> controlSequenceLength (Nonnegative integer)
            %      The length of a control sequence (without default input)
            %      generated by the controller.
            %
            %   >> Q (Square matrix, dimPlantState-by-dimPlantState)
            %      The time-invariant state or plant output weighting matrix.
            %
            %   >> R (Square Matrix, dimPlantInput-by-dimPlantInput)
            %      The time-invariant input weighting matrix.
            %
            %   >> Z (Matrix, n-by-dimPlantState, optional)
            %      The time-invariant plant output (performance) matrix, 
            %      i.e., z_k = Z*x_k
            %      By default, Z = I (identity matrix).
            %
            % Returns:
            %   << augQ (3D matrix, one for each mode of the MJLS)
            %      Weighting matrix of the augmented state, one per mode.
            %
            %   << augR (3D matrix, one for each mode of the MJLS)
            %      Weighting matrix of the augmente input (control sequence), one per mode.
            %
            dimPlantInput = size(R, 2);
            if nargin == 4 
                stateWeighting = Z' * Q * Z; % performance output for tracking
            else
                stateWeighting = Q;
            end
            if controlSequenceLength == 1
                augQ = repmat(stateWeighting, 1, 1, 2);
                augR = cat(3, R, zeros(dimPlantInput));
            else
                dimState = size(stateWeighting, 1);
                dimAugmentedState = dimState ...
                    + dimPlantInput * (controlSequenceLength * (controlSequenceLength - 1) / 2);
                 dimEta = dimAugmentedState - dimState;
                
                augR = zeros(dimPlantInput * controlSequenceLength, ...
                    dimPlantInput * controlSequenceLength, controlSequenceLength + 1);
                % augR is only non-zero for first mode
                augR(1:dimPlantInput, 1:dimPlantInput, 1) = R;
                                
                augQ = repmat(blkdiag(stateWeighting, zeros(dimEta)), ...
                    1, 1, controlSequenceLength + 1);
                % cost matrices per mode
                % for first and last mode, no work work required
                sums = dimPlantInput * cumsum(1:controlSequenceLength -1);
                for i = (controlSequenceLength - 1):-1:1
                    startIdx = dimState + dimEta - sums(controlSequenceLength - i) + 1;
                    endIdx = startIdx + dimPlantInput - 1;
                    augQ(startIdx:endIdx, startIdx:endIdx, i + 1) = R;
                end
            end
        end
        
        %% performModelAugmentation
        function [dimAugmentedState, F, G, augA, augB, augQ, augR] ...
            = performModelAugmentation(sequenceLength, dimX, dimU, A, B, Q, R, Z)
            % Convenience function to compute the matrices of 1) the augmented state space which describes the
            % linear NCS (including the network-actuator subsystem) in
            % terms of an MJLS, and 2) the correspondingly augmented cost matrices. 
            % Note that no parameter checks are carried out.
            %
            % Literature: 
            %   Jörg Fischer,
            %   Optimal sequence-based control of networked linear systems,
            %   Karlsruhe series on intelligent sensor-actuator-systems, Volume 15,
            %   KIT Scientific Publishing, 2015.
            %
            %   Jörg Fischer, Achim Hekler, Maxim Dolgov and Uwe D. Hanebeck,
            %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject
            %   to Random Transmission Delays and Packet Losses,
            %   Proceedings of the 2013 American Control Conference (ACC 2013), 
            %   Washington D. C., USA, June 2013.
            %
            % Parameters:
            %   >> sequenceLength (Positive integer)
            %      The length of a control sequence (without default input)
            %      generated by the controller.
            %
            %   >> dimX (Positive integer)
            %      The dimension of the plant's state.
            %
            %   >> dimU (Positive integer)
            %      The dimension of the inputs applied to the plant.
            %
            %   >> A (Square matrix, dimX-by-dimX)
            %      The time-invariant system matrix of the plant dynamics.
            %
            %   >> B (Matrix, dimX-by-dimU)
            %      The time-invariant input matrix of the plant dynamics.
            %
            %   >> Q (Square matrix, dimX-by-dimX)
            %      The time-invariant state or plant output weighting matrix.
            %
            %   >> R (Square Matrix, dimU-by-dimU)
            %      The time-invariant input weighting matrix.
            %
            %   >> Z (Matrix, n-by-dimX, optional)
            %      The time-invariant plant output (performance) matrix, 
            %      i.e., z_k = Z*x_k
            %      By default, Z = I (identity matrix).
            %
            % Returns:
            %   << dimAugmentedState (Positive integer)
            %      The dimension of the augmented state which consists of the plant
            %      state and the applicable control inputs.
            %
            %   << F (Square matrix)
            %      System matrix of the dynamics of the network-actuator
            %      subsystem, which governs the control input selection at the actuator.
            %
            %   << G (Matrix)
            %      Input matrix of the network-actuator subsystem, that maps the entries of the current control sequence
            %      which are meant to be applied to the plant at future time steps to the
            %      vector containing possible future control inputs.
            %
            %   << augA (3D matrix, one for each mode of the MJLS)
            %      System matrix of the augmented open-loop dynamics, one per
            %      mode.
            %
            %   << augB (3D matrix, one for each mode of the MJLS)
            %      Input matrix of the augmented open-loop dynamics, one per
            %      mode.
            %
            %   << augQ (3D matrix, one for each mode of the MJLS)
            %      Weighting matrix of the augmented state, one per mode.
            %
            %   << augR (3D matrix, one for each mode of the MJLS)
            %      Weighting matrix of the augmente input (control sequence), one per mode.
            %
            if sequenceLength == 1
                dimAugmentedState = dimX;
            else
                dimAugmentedState = dimX + dimU * (sequenceLength * (sequenceLength - 1) / 2);
            end 
            [F, G, ~, ~, augA, augB] = Utility.createAugmentedPlantModel(sequenceLength, A, B);
            if nargin == 8
                [augQ, augR] = Utility.createAugmentedCostModel(sequenceLength, Q, R, Z);
            else
                [augQ, augR] = Utility.createAugmentedCostModel(sequenceLength, Q, R);
            end
        end
        
        %% createAugmentedLinearIntegralConstraints
        function [augStateWeightings, augInputWeightings] = ...
            createAugmentedLinearIntegralConstraints(sequenceLength, stateWeightings, inputWeightings)
            % Create augmented linear integral constraints to express the
            % original linear integral state and input constraints in terms of the augmented state space.
            %
            % Literature: 
            %   Maxim Dolgov, Jörg Fischer, and Uwe D. Hanebeck,
            %   Sequence-based LQG Control over Stochastic Networks
            %   with Linear Integral Constraints,
            %   Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014), 
            %   Los Angeles, California, USA, December 2014.
            %
            % Parameters:
            %   >> sequenceLength (Nonnegative integer)
            %      The length of a control sequence (without default input)
            %      generated by the controller.
            %
            %   >> stateWeightings (3D matrix, dimPlantState-by-horizonLength-by-numConstraints)
            %      The state weightings (without those for the terminal state) of the individual constraints, 
            %      slice-wise arranged.
            %
            %   >> inputWeightings (3D matrix, dimPlantInput-by-horizonLength-by-numConstraints)
            %      The input weightings of the individual constraints, slice-wise arranged.
            %
            % Returns:
            %   << augStateWeightings (4D matrix, dimAugmentedState-by-horizonLength-by-numConstraints-by-numModes)
            %      State weightings expressed in terms of the augmented state, one per mode.
            %
            %   << augInputWeightings (4D matrix, dimAugmentedInput-by-horizonLength-by-numConstraints-by-numModes)
            %      Input weightings expressed in terms of the augmented input (control sequence), one per mode.
            %
                 
            % weightings are 3d matrices whith number of slices identical
            if sequenceLength == 1
                % we have two modes
                augStateWeightings = repmat(stateWeightings, 1, 1, 1, 2);
                augInputWeightings = cat(4, inputWeightings, zeros(size(inputWeightings)));
            else
                % translate the weightings, mode-dependent
                dimPlantInput = size(inputWeightings, 1);
                dimSequence = sequenceLength * dimPlantInput;
                dimPlantState = size(stateWeightings, 1);
                dimAugState = dimPlantState + dimPlantInput * (sequenceLength * (sequenceLength - 1) / 2);
                dimEta = dimAugState - dimPlantState;
                numModes = sequenceLength + 1;
             
                augInputWeightings = zeros(dimSequence, size(inputWeightings, 2), size(inputWeightings, 3), numModes);
                augInputWeightings(1:dimPlantInput, :, :, 1) = inputWeightings;
                
                augStateWeightings = repmat([stateWeightings; zeros(dimEta, size(stateWeightings, 2), size(stateWeightings, 3))],...
                                        1, 1, 1, numModes);
                % for first and last mode, no work work required
                sums = dimPlantInput * cumsum(1:sequenceLength -1);
                for i = (sequenceLength - 1):-1:1
                    startIdx = dimAugState - sums(sequenceLength - i) + 1;
                    endIdx = startIdx + dimPlantInput - 1;
                    augStateWeightings(startIdx:endIdx, :, :, i + 1) = inputWeightings;
                end
            end
        end
        
        %% normalizeProbabilities
        function normalizedProbs = normalizeProbabilities(probs, lowerBound)
            % Convenience function to normalize a stochastic vector or column-stochastic matrix.            
            % Here, normlization means that a lower bound on the
            % probabilities is imposed, thus avoiding zero entries.
            %
            % Parameters:
            %   >> probs (Vector or matrix)
            %      Either a stochastic vector, i.e., a vector with nonnegative
            %      entries where each row sums to 1, or a not necessarily
            %      square matrix with nonnegative entries whose columns
            %      sum to 1.
            %      
            %   >> lowerBound (Nonnegative scalar, optional)
            %      The lower bound for the entries of the given vector or matrix. 
            %      If 0 is passed here, this function has no effect.
            %      If no value is passed, the dafult value 1e-12 is used.
            %
            % Returns:
            %   << normalizedProbs (Vector with nonnegative entries)
            %      The given stochastic vector or column-stochastic matrix after normalization,
            %      such that no entries is less than the given lower bound.            
            %
            if nargin < 2
                lowerBound = 1e-12;
            end           
            
            idx = find(probs <= lowerBound);
            normalizedProbs = probs;
            if ~isempty(idx)
                normalizedProbs(idx) = lowerBound;                
            end
            normalizedProbs = normalizedProbs ./ sum(normalizedProbs);
        end
        
        %% normalizeTransitionMatrix
        function normalizedTransMat = normalizeTransitionMatrix(transMat, lowerBound)
            % Convenience function to normalize a row-stochastic matrix,
            % such as as transition matrix of a Markov chain.
            % Here, normlization means that a lower bound on the
            % probabilities is imposed, thus avoiding zero entries.
            %
            % Parameters:
            %   >> transMat (Matrix)
            %      A row-stochastic matrix, i.e., a matrix with nonnegative
            %      entries where each row sums to 1, but not necessarily
            %      square.
            %      
            %   >> lowerBound (Nonnegative scalar, optional)
            %      The lower bound for the entries of the given matrix. 
            %      If 0 is passed here, this function has no effect.
            %      If no value is passed, the dafult value 1e-12 is used.
            %
            % Returns:
            %   << normalizedTransMat (Vector with nonnegative entries)
            %      The given row-stochastic matrix after normalization,
            %      such that no entries is less than the given lower bound.            
            %
            if nargin < 2
                lowerBound = 1e-12;
            end           
            
            idx = find(transMat <= lowerBound);
            normalizedTransMat = transMat;
            if ~isempty(idx)
                normalizedTransMat(idx) = lowerBound;                
            end
            normalizedTransMat = normalizedTransMat ./ sum(normalizedTransMat, 2);
        end
        
        %% calculateDelayTransitionMatrix
        function transitionMatrix = calculateDelayTransitionMatrix(delayProbs)
            % Compute the (possibly time-varying) transition matrix of the Markov chain that occurs in the combined stochastic model 
            % of communication network and actuator. 
            % Note that no parameter checks are carried out.
            %
            % Literature: 
            %   Jörg Fischer, Achim Hekler, and Uwe D. Hanebeck,
            %   State Estimation in Networked Control Systems,
            %   Proceedings of the 15th International Conference on Information Fusion (Fusion 2012),
            %   Singapore, July 2012.
            %
            %   Jörg Fischer, Achim Hekler, Maxim Dolgov and Uwe D. Hanebeck,
            %   Optimal Sequence-Based LQG Control over TCP-like Networks Subject
            %   to Random Transmission Delays and Packet Losses,
            %   Proceedings of the 2013 American Control Conference (ACC 2013), 
            %   Washington D. C., USA, June 2013.
            %
            %   Florian Rosenthal and Uwe D. Hanebeck,
            %   Stability Analysis of Polytopic Markov Jump Linear Systems 
            %   with Applications to Sequence-Based Control over Networks,
            %   Proceedings of the 21st IFAC World Congress (IFAC 2020),
            %   Berlin, Germany, July 2020.
            %
            % Parameters:
            %   >> delayProbs (Nonnegative vector or matrix)
            %      -If a vector is passed:
            %       A vector with nonnegative entries that sum up to 1
            %       that describes the probability distribution specifying
            %       the (time-invariant) delays of the control sequences.
            %       The i-th entry denotes the probability that a delay of
            %       (i-1) time steps occurs, and the last entry specifies
            %       the loss probability (i.e., infinite delay).
            %      -If a matrix is passed:
            %       A matrix where each column is a vector with nonnegative 
            %       entries that sum up to 1 to describe the probability
            %       distributions of the (time-varying) delays of the
            %       control sequences.
            %       In particular, the first column specifies the delay
            %       distribution of the control sequence from the current
            %       time step, the second column the delay distribution of
            %       the control sequence from the previous time step, and
            %       so forth. In total, N-1 distributions (i.e. columns) must
            %       be given, where N-1 is the sequence length.
            %       The i-th entry of each column denotes the probability that a delay of
            %       (i-1) time steps occurs, and the last entry specifies
            %       the loss probability (i.e., infinite delay).
            %
            % Returns:
            %   << transitionMatrix (Square matrix)
            %      The transition matrix of the resulting Markov chain,
            %      which is of dimension N-by-N where N is the number of
            %      elements of the given delay distribution.
            %      Note that this matrix is always lower Hessenberg.
                       
            normalizedDelayProbs = Utility.normalizeProbabilities(delayProbs);
            sums = cumsum(normalizedDelayProbs);
            if isvector(normalizedDelayProbs)
                % probability vector given -> time-invariant transition matrix
                transitionMatrix = diag(1 - sums(1:end -1), 1);
                % calculate Transition matrix
                for j=1:numel(normalizedDelayProbs)
                    transitionMatrix(j:end, j) = normalizedDelayProbs(j);
                end
            else
                % matrix given, with column-wise arranged delay
                % distributions: the one for packet from the current time step in the first
                % colum, and so forth
                [numModes, seqLength] = size(normalizedDelayProbs);
                % number of rows =  modes 
                % number of columns should be numModes-1 = sequence length
                
                % compute q_tilde(m)
                q_tilde = diag(normalizedDelayProbs) ./ [1; 1 - diag(sums, 1)];
                prods = cumprod(1-q_tilde);
                transitionMatrix = diag(prods, 1);
                transitionMatrix(1:seqLength, 1) = q_tilde(1);
                for j = 2:seqLength
                    transitionMatrix(j:seqLength, j) = q_tilde(j) * prods(j-1);
                end
                transitionMatrix(numModes, :) = transitionMatrix(seqLength, :);
            end
        end        
        
        %% calculateTransitionMatrixCorrelatedDelays
        function transitionMatrix = calculateTransitionMatrixCorrelatedDelays(delayTransitionMatrix, controlSeqLength, useMex, useStatDist)
            % Compute the transition matrix of the Markov chain that occurs in the combined stochastic model of communication network and actuator
            % for the case of correlated packet delays and losses, modeled
            % in the form of a Markov chain tau_k.
            % The computation is based on a lumping of the aggregated Markov
            % chain tau_k, tau_{k-1}, tau_{k-(N-1)}, which is easily
            % created given the transition matrix of tau_k and where N is
            % the desired control sequence length.
            % Note that, however, this aggregated chain is generally not
            % exactly lumpable with respect to the clusters theta_k (which are the modes of
            % the augmented dynamics).
            % Note also that no parameter checks are carried out.
            %
            % Literature: 
            %   John G. Kemeny and J. Laurie Snell,
            %   Finite Markov Chains,
            %   Sections 6.3-6.4,
            %   Van Nostrand, Princeton, NJ, USA, 1960.
            %
            % Parameters:
            %   >>  delayTransitionMatrix (Stochastic matrix)
            %       The transition matrix of the Markov chain governing the
            %       packet delays and losses (tau_k) in the communication
            %       between controller and actuator. The dimension of this
            %       matrix must be at least (N+1)-by-(N+1), where N is the
            %       passed length of the control sequences.
            %
            %   >> controlSeqLength (Positive integer)
            %      The length of a control sequence (without default input)
            %      generated by the controller.
            %
            %   >> useMex (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether the C++ (mex) implementation
            %      shall be used for the computation of the transition matrix. 
            %      If left out, the default value false is used.
            %
            %   >> useStatDist (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether stationary distribution of the
            %      augmented delay Markov chain shall be used for
            %      constructing the distribution matrix U, which ensures
            %      that the resulting lumped chain has the same convergence
            %      properties.
            %      If false is passed here, the distribution matrix is a
            %      constructed such that each entry is 1/d_l, with d_l the
            %      number of elements in the cluster corresponding to that
            %      entry, resulting in a transition matrix for theta that
            %      has the same properties as if the delays were a white
            %      process.
            %      If left out, the default value true is used.
            %
            % Returns:
            %   << transitionMatrix (Square matrix)
            %      The transition matrix of the resulting Markov chain,
            %      which is of dimension (N+1)-by-(N+1) where N is the given sequence length.
            %      Note that this matrix is always lower Hessenberg.
            
            if nargin < 3
                useMex = false;
                useStatDist = true;
            elseif nargin < 4
                useStatDist = true;  
            end
            
            if ~useMex
                numDelays = size(delayTransitionMatrix, 1); % last state also comprises packet losses
                numCaModes = controlSeqLength + 1;
                assert(Checks.isSquareMat(delayTransitionMatrix, numDelays) && numDelays >= numCaModes, ...
                    'Utility:CalculateTransitionMatrixCorrelatedDelays:InvalidTransitionMatrix', ...
                    '** Transition matrix of the delay Markov chain must be at least %d-by-%d **', numCaModes, numCaModes);

                numAggregations = controlSeqLength;
                numAggregatedStates = numDelays^numAggregations;
                % the aggregated chain is not regular
                aggregatedStatesDecimal = 0:numAggregatedStates-1; % decimal numbers
                statesBaseN = zeros(numAggregatedStates, numAggregations);
                for l = 0:numAggregations-1
                    statesBaseN(:, l+1) = mod(floor(aggregatedStatesDecimal ./ (numDelays ^ l)), numDelays);
                end
                % use a sparse matrix for the transitions of the aggregated
                % chain
                % matrix has repetitive structure (regarding the rows)
                numRows = numAggregatedStates / numDelays;
                numNonzeroEntries = numDelays ^ numAggregations; % maximum number of nonzero transition probs of the aggregated chain
                % first column indicates rows
                % second column indicates cols
                % third column contains nonzero values
                sparseEntries(:, 3) = repmat(reshape(delayTransitionMatrix', numDelays * numDelays, 1), ...
                    numNonzeroEntries / (numDelays * numDelays), 1);
                sparseEntries(:, 2) = 1:numNonzeroEntries;
                sparseEntries(:, 1) = kron(1:numRows, ones(1, numDelays));

                aggT = sparse(sparseEntries(:, 1), sparseEntries(:, 2), sparseEntries(:, 3), numRows, numAggregatedStates);
                
                % compute the indices of the clusters
                % states of aggregated chain are lumped into clusters-> mode of
                % augmented system
                numClusters = numCaModes; % equals the number of modes
                clusterIdx = cell(1, numClusters);
                % first cluster is easy: tau_k = 0
                clusterIdx{1} = find(~(statesBaseN(:, 1)));
                % i-th cluster: tau_k > 0, tau_{k-1} > 1, ... tau_{k-(i-1)} > i-1, tau_{k-i} <= i
                lastIdx = 1:numAggregatedStates;
                usedIdx = clusterIdx{1};
                for j=2:numClusters-1
                    % create the next cluster by concatenating function handles
                    fun = @() statesBaseN(:, 1) > 0; % tau_k > 0
                    for k=2:j-1
                        fun = @() fun() & statesBaseN(:, k) > k-1;
                    end
                    % finally the last one
                    fun = @() fun() & statesBaseN(:, j) <= j-1;

                    clusterIdx{j} = find(fun());
                    usedIdx = union(usedIdx, clusterIdx{j});
                end
                %last cluster: all indices not yet assigned
                clusterIdx{end} = find(~ismember(lastIdx, usedIdx));
                % construct matrices U, V required for clustering into modes
                                
                U = zeros(numClusters, numAggregatedStates);
                V = zeros(numAggregatedStates, numClusters);
                if useStatDist
                    % use the fact that stationary distribution can be computed
                    % by means of stationary distribution of delayTransitionMatrix

                    % (theta_k) (MC lumping)                
                    %[eigVec, ~] = eigs(repmat(aggT, numDelays, 1)', 1); % compute Perron eigenvalue (1) and corresponding left eigenvector
                    %statDist = eigVec ./ sum(eigVec);
                    % generalization of Theorem 6.5.2 in Kemeny and
                    % Small's book to obtain stationary distribution of
                    % aggT based stationary distribution of delayTransitionMatrix
                    % stat (i,j,l) = statP(l) * t(l,j)*t(j,i) and so forth
                    % where statP is stationary distribution of delay MC
                    % and t are the entries of delayTransitionMatrix
                    tmp = delayTransitionMatrix;
                    for n=1:numAggregations-2
                        tmp = delayTransitionMatrix .* permute(tmp, [3 1 2]);                        
                        tmp = reshape(permute(tmp, [1 3 2]), numDelays, []);                        
                    end
                    statDist = reshape(transpose(Utility.computeStationaryDistribution(delayTransitionMatrix) .* tmp), numAggregatedStates, 1);
                
                    %                     a = Utility.computeStationaryDistribution(delayTransitionMatrix);
                    %                     for i=1:numDelays
                    %                         for j=1:numDelays                            
                    %                             for l=1:numDelays
                    %                                 %res(l,j,i) = [a(i)*delayTransitionMatrix(i,j) * delayTransitionMatrix(j,l)];
                    %                                 for m=1:numDelays
                    %                                 %res2 = [res2; delayTransitionMatrix(i, :) .* delayTransitionMatrix(:,l)];
                    %                                     res(i, j, l, m) = [a(m)*delayTransitionMatrix(m,l) * delayTransitionMatrix(l,j) * delayTransitionMatrix(j,i)];
                    %                                 end
                    %                             end
                    %                         end
                    %                     end

                    for j=1:numClusters                        
                        U(j, clusterIdx{j}) = statDist(clusterIdx{j}) ./ sum(statDist(clusterIdx{j}));
                        V(clusterIdx{j}, j) = 1;
                    end
                else
                    % yields a transition matrix for theta that has the familiar structure: 
                    % last 2 rows equal, also entries below main diagonal the same per column 
                    % (including diagonal entrystarting from column 2)
                    % also, requirement about delay chain and its properties (irreducible etc) not needed)
                    for j=1:numClusters
                        U(j, clusterIdx{j}) = 1/ numel(clusterIdx{j});
                        V(clusterIdx{j}, j) = 1; 
                    end
                end       
                transitionMatrix = U * repmat(aggT, numDelays, 1) * V;
            else
                transitionMatrix = mex_CalculateTransitionMatrixCorrDelays(delayTransitionMatrix, controlSeqLength, useStatDist);
            end
        end
        
        %% computeStationaryDistribution
        function stationaryDist = computeStationaryDistribution(transitionMatrix, useEig)
            % Convenience function to compute the stationary distribution of a
            % time-homogeneous Markov chain with given transition matrix.
            % Note that no parameter checks are carried out.
            %
            % Parameters:
            %   >> transitionMatrix (Square matrix)
            %      The transition matrix of the Markov chain.
            %
            %   >> useEig (Logical, optional)
            %      If true is passed (default), then an eigenvalue decompostion is used to compute the stationary distribution.
            %      More precisely, the stationary distribution p is given by the
            %      equation p*T=p, where T is the transition matrix.
            %      That is, p is a (normalized) left eigenvector corresponding to the eigenvalue 1.
            %      If false is passed, the above fixed-point iteration is
            %      solved until convergence is detected, or a maximum number of
            %      iterations (100000) was reached.
            %
            % Returns:
            %   << stationaryDist (Vector with nonnegative entries)
            %      The stationary distribution of the Markov chain.
            %
            
            if nargin == 1
                useEig = true;
            end
            if useEig
                % get left eigenvectors of transition matrix which are the usual
                % eigenvectors of transposed transition matrix -> stationary
                % distribution of the Markov chain
                [eigenVecs, eigenValues] = eig(transitionMatrix', 'vector');
                [maxVal,idx] = max(eigenValues);
                if round(maxVal * 1e8) ~= 1e8
                    error('Utility:ComputeStationaryDistribution:InvalidTransitionMatrix', ...
                        ['** Transition matrix of the Markov chain ' ...
                        ' is supposed to have (max) eigenvalue 1 **']);
                end
                % compute stationary distribution, which are just the normalized
                % eigenvectors of eigenvalue 1
                stationaryDist = eigenVecs(:, idx) ./ sum(eigenVecs(:, idx));
            else
                % compute stationary distribution by iteration until convergence
                maxIterNum = 100000;
                convergenceDiff = 1e-8;
                % start with uniform distribution
                stationaryDist = ones(1, size(transitionMatrix, 1)) / size(transitionMatrix, 1);
                old = inf(1, size(transitionMatrix, 1));
                i = 1;
                while i <= maxIterNum && norm (old - stationaryDist) > convergenceDiff
                    old = stationaryDist;
                    stationaryDist = old * transitionMatrix;
                    i = i + 1;
                end
            end
        end
        
        %% truncateDiscreteProbabilityDistribution
        function distribution = truncateDiscreteProbabilityDistribution(probs, numElements)
            % Convenience function to truncate or extend a given discrete
            % probability distribution to a given number of elements.
            % Note that no parameter checks are carried out.
            %
            % Parameters:
            %   >> probs (Nonnegative vector)
            %      A vector with nonnegative entries that sum up to 1
            %      that describes a discrete probability distribution.
            %
            %   >> numElements (Positive integer)
            %      The desired number of elements in the resulting
            %      probability vector.
            %
            % Returns:
            %   << distribution (Nonnegative column vector)
            %      A column vector with nonnegative entries that sum up to 1
            %      that describes a discrete probability distribution with the desired number of elements.
            
            numProbs = length(probs);
            
            if numProbs <= numElements
                % fill up with zeros
                distribution = [probs(:); zeros(numElements - numProbs, 1)];
            else
                % cut the distribution and fill up with an entry so that
                % sum is 1 again
                distribution = [reshape(probs(1:1:numElements -1), numElements -1, 1); 1 - sum(probs(1:1:numElements -1))];            
            end
        end
        
        %% calculateRMSE
        function rmse = calculateRMSE(trueStates, estimates)
            % Calculate the Root mean squared error (RMSE) of a collection of
            % estimates resulting from a Monte Carlo simulation.
            %
            % Parameters:
            %   >> trueStates (3D-matrix or 2D-matrix)
            %      The true states at each time step (column-wise) and per run (page-wise).
            %      In case a 2D-matrix is passed, then only 1 run is assumed.
            %
            %   >> estimates (3D-matrix or 2D-matrix)
            %      The estimates produced by the filter, must be a matrix of the
            %      same size as the true states.
            %
            % Returns:
            %   >> rmse (Row vector) 
            %      The computed RMSE, where the RMSE at time step i is stored in the
            %      i-th component of the row vector.
    
            s1 = size(trueStates); % first dim: state, second dim: numTimeSteps, third dim: numRuns
            s2 = size(estimates);

            if isequal(s1, s2)
                numTimeSteps = s1(2);
                if numel(s1) == 3
                    numRuns = s1(3);
                else
                    numRuns = 1;
                end
                rmse(numTimeSteps) = nan;
                squaredErrors = (estimates - trueStates) .^ 2;
                for k=1:numTimeSteps
                    rmse(k) = sqrt(sum(sum(squaredErrors(:, k, :))) / numRuns);
                end
            else
                rmse = [];
            end
        end
    end    
end

