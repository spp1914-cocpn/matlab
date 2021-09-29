classdef RecedingHorizonUdpLikeControllerTest < matlab.unittest.TestCase
    % Test cases for RecedingHorizonUdpLikeController.
    
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
    
    properties (Constant)
        absTol = 5 * 1e-3;
        numIterations = 100;
    end
    
     properties (Access = private)
        A;
        B;
        Q;
        R;
        C;
        J;
        sequenceLength;
        dimX;
        dimU;
        dimY;
        
        modeTransitionMatrix;
                        
        W;
        V;

        x0;
        x0Cov;
        
        caDelayProbs;
        scDelayProbs;
        scDelayTransMat;
        maxMeasDelay;
        truncatedScDelayProbs;
        truncatedScDelayTransMat;
        
        horizonLength;
        
        augA;
        augB;
        augC;
        augW;
        augV;        
        augQ;
        augX0;
        augX0Cov;
        
        xUpperbar;
        xUnderbar;
        
        controllerUnderTest;
        
        S; 

        initialL;
        initialK;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            % use (noise-free) stirred tank example (Example 6.15, p. 500-501) from
            %
            % Huibert Kwakernaak, and Raphael Sivan, 
            % Linear Optimal Control Systems,
            % Wiley-Interscience, New York, 1972.
            %
            
            this.dimX = 2;
            this.dimU = 2;
            this.dimY = 2;                             
            this.maxMeasDelay = 1;            
            this.horizonLength = 3;
            
            this.sequenceLength = 2; % so we have three modes
            this.caDelayProbs = ones(1, 5) / 5;
            this.scDelayProbs = ones(1, 6) / 6;
            this.scDelayTransMat = repmat(reshape(this.scDelayProbs, 1, []), numel(this.scDelayProbs), 1);
            
            this.truncatedScDelayProbs = Utility.truncateDiscreteProbabilityDistribution(this.scDelayProbs, this.maxMeasDelay + 2);
            this.truncatedScDelayTransMat = repmat(reshape(this.truncatedScDelayProbs, 1, []), numel(this.truncatedScDelayProbs), 1);
            
            % from the delay probs and the 3 modes we get the following
            % transition matrix
            this.modeTransitionMatrix = [1/5 4/5 0; 1/5 1/5 3/5; 1/5 1/5 3/5];
            
            this.A = diag([0.9512, 0.9048]);           
            this.B = [4.877 4.877; -1.1895 3.569];
            this.C = 2 * ones(this.dimY, this.dimX);
            
            this.Q = diag([0.01, 1]) * diag([50 0.02]) * diag([0.01, 1]); % R_3 in the book
            this.R = diag([1/3, 3]); % R_2 in the book
         
            this.W = eye(this.dimX) * 0.01;
            this.V = eye(this.dimY) * 0.001;
            
            this.x0 = ones(this.dimX,1);
            this.x0Cov = eye(this.dimX) * 0.2;
            
            this.initAdditionalProperties();

            this.controllerUnderTest = RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, ...
                this.horizonLength, this.x0, this.x0Cov);
            
            this.controllerUnderTest.maxNumIterations = RecedingHorizonUdpLikeControllerTest.numIterations;
        end
    end
    
    methods (Access = private)
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            numModes = this.sequenceLength + 1;
            dimEta = 2;
            dimAugState = 2 * this.dimX + dimEta;
            numSMatrices = 4; %2^M, where M is max meas delay + 1 
            dimAugmentedMeas = 2 * this.dimY;         
            
            this.augC = [kron(eye(2), this.C),...
                zeros(dimAugmentedMeas, dimEta)];

            this.S = zeros(dimAugmentedMeas, dimAugmentedMeas, numSMatrices);
            % ordered according to de2bi([0 1 2 3])
            % 0 0
            % 1 0
            % 0 1
            % 1 1
            % here 1 indicates that measurement is available for processing
            this.S(:, :, 2) = blkdiag(eye(this.dimY), zeros(this.dimY));
            this.S(:, :, 3) = blkdiag(zeros(this.dimY), eye(this.dimY));
            this.S(:, :, 4) = blkdiag(eye(this.dimY), eye(this.dimY));

            [F, G, H, this.J] = Utility.createAugmentedPlantModel(this.sequenceLength, this.A, this.B);           
            D = eye(this.dimX);            
            
            this.augA = zeros(dimAugState, dimAugState, numModes);
            Abar = zeros(dimAugState - dimEta);
            Abar(:, 1:this.dimX) = [this.A; D];
            for j=1:numModes
                this.augA(:, :, j) = blkdiag(Abar, F);
                this.augA(1:this.dimX, end-dimEta+1:end, j) = this.B * H(:, :, j);
            end        
            this.augB = repmat([zeros(dimAugState - size(G, 1), size(G, 2)); G], 1, 1, numModes);
            this.augB(1:this.dimX, :, :) = mtimesx(this.B, this.J);
               
            this.augQ = repmat(blkdiag(this.Q, zeros(dimAugState - this.dimX)), 1, 1, numModes);
            idx = dimAugState - dimEta + 1;
            % add the mode dependent part
            for i = 1:numModes                
                this.augQ(idx:dimAugState, idx:dimAugState, i) ...
                    = H(:, :, i)' * this.R * H(:, :, i); % H_tilde in the paper
            end
            
            this.augX0 = [this.x0(:); this.x0(:); zeros(dimEta, 1)];
            this.augX0Cov = blkdiag(kron(ones(2,2), this.x0Cov), zeros(dimEta));
            
            this.augW = blkdiag(this.W, zeros(this.dimX), zeros(dimEta));
            this.augV = blkdiag(this.V, this.V);
            
            this.xUpperbar = zeros(dimAugState, dimAugState, numModes);
            this.xUpperbar(:, :, 2) = this.augX0Cov;
            this.xUnderbar = zeros(dimAugState, dimAugState, numModes);
            this.xUnderbar(:, :, 2) = this.augX0 * this.augX0' + this.augX0Cov;
            
            % pick some initial gains        
            this.initialK = ones(dimAugState, 2 * this.dimY, this.horizonLength);
            this.initialL = ones(this.dimU * this.sequenceLength, dimAugState, this.horizonLength);
        end
        
        %% computeGains
        function [K, L, Aexp, Bexp] = computeGains(this, caModeProbs, S_tilde)
            dimEta = 2;
            dimAugState = 2 * this.dimX + dimEta;
            numCaModes = this.sequenceLength + 1;
            
            K = this.initialK;
            L = this.initialL;
            
            % first compute Aexp and Bexp for the whole horizon
            Aexp = zeros(dimAugState, dimAugState, this.horizonLength);
            Bexp = zeros(dimAugState, this.dimU * this.sequenceLength, this.horizonLength);
            for k=1:this.horizonLength
                for i=1:numCaModes
                    Aexp(:, :, k) = Aexp(:, :, k) + caModeProbs(i, k) * this.augA(:, :, i);
                    Bexp(:, :, k) = Bexp(:, :, k) + caModeProbs(i, k) * this.augB(:, :, i);
                end
            end
   
            for iter=1:RecedingHorizonUdpLikeControllerTest.numIterations            
                XUpper = zeros(dimAugState, dimAugState, numCaModes, this.horizonLength + 1);
                XUpper(:, :, :, 1) = this.xUpperbar;
                XUnder = zeros(dimAugState, dimAugState, numCaModes, this.horizonLength + 1);
                XUnder(:, :, :, 1) = this.xUnderbar;            

                % predict XUpper and Xunder over the complete horizon            
                for k=1:this.horizonLength                   
                    for j=1:numCaModes
                        X_temp = zeros(dimAugState, dimAugState);
                        X_temp2 = zeros(dimAugState, dimAugState);
                        for i=1:numCaModes
                            X_temp = X_temp + this.modeTransitionMatrix(i, j) * ((this.augA(:, :, i) - K(:, :, k)* S_tilde(:, :, k)*this.augC) * XUpper(:, :, i, k) * (this.augA(:, :, i) - K(:, :, k)* S_tilde(:, :, k)*this.augC)' ...
                                + (this.augA(:, :, i) -Aexp(:, :, k) + this.augB(:, :, i) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k)) * XUnder(:, :, i, k) * (this.augA(:, :, i) -Aexp(:, :, k) + this.augB(:, :, i) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k))' ...
                                + caModeProbs(i, k) * (this.augW + (K(:, :, k)* S_tilde(:, :, k)) * this.augV * (K(:, :, k)* S_tilde(:, :, k))'));

                            X_temp2 = X_temp2 + this.modeTransitionMatrix(i, j) * ((K(:, :, k)* S_tilde(:, :, k)*this.augC) * XUpper(:, :, i, k) * (K(:, :, k)* S_tilde(:, :, k)*this.augC)' ...
                                + (Aexp(:, :, k) + Bexp(:, :, k)* L(:, :, k)) * XUnder(:, :, i, k) * (Aexp(:, :, k) + Bexp(:, :, k)* L(:, :, k))' ...
                                + caModeProbs(i, k) * K(:, :, k)* S_tilde(:, :, k) * this.augV * (K(:, :, k)* S_tilde(:, :, k))');
                        end
                        XUpper(:, :, j, k+1) = X_temp; 
                        XUnder(:, :, j, k+1) = X_temp2;
                    end
                end
            
                % backward pass for the gains and the costate
                terminalQ = idare(this.A, this.B, this.Q, this.R);
                PUpper = repmat(blkdiag(terminalQ, zeros(dimAugState - this.dimX)), 1, 1, numCaModes); % terminal Q
                PUnder = zeros(dimAugState, dimAugState, numCaModes);
                omega = zeros(1, numCaModes);
                oldCosts = inf;
                for k=this.horizonLength:-1:1
                    omega_sum1 = this.modeTransitionMatrix(1, 1) * omega(1) + this.modeTransitionMatrix(1, 2) * omega(2) ...
                        + this.modeTransitionMatrix(1, 3) * omega(3);
                    omega_sum2 = this.modeTransitionMatrix(2, 1) * omega(1) + this.modeTransitionMatrix(2, 2) * omega(2) ...
                        + this.modeTransitionMatrix(2, 3) * omega(3);
                    omega_sum3 = this.modeTransitionMatrix(3, 1) * omega(1) + this.modeTransitionMatrix(3, 2) * omega(2) ...
                        + this.modeTransitionMatrix(3, 3) * omega(3);
            
                    P_sum11 = this.modeTransitionMatrix(1, 1) * PUnder(:, :, 1) + this.modeTransitionMatrix(1, 2) * PUnder(:, :, 2) ...
                    + this.modeTransitionMatrix(1, 3) * PUnder(:, :, 3);
                    P_sum12 = this.modeTransitionMatrix(1, 1) * PUpper(:, :, 1) + this.modeTransitionMatrix(1, 2) * PUpper(:, :, 2) ...
                        + this.modeTransitionMatrix(1, 3) * PUpper(:, :, 3);
                    P_sum21 = this.modeTransitionMatrix(2, 1) * PUnder(:, :, 1) + this.modeTransitionMatrix(2, 2) * PUnder(:, :, 2) ...
                    + this.modeTransitionMatrix(2, 3) * PUnder(:, :, 3);
                    P_sum22 = this.modeTransitionMatrix(2, 1) * PUpper(:, :, 1) + this.modeTransitionMatrix(2, 2) * PUpper(:, :, 2)...
                        + this.modeTransitionMatrix(2, 3) * PUpper(:, :, 3);
                    P_sum31 = this.modeTransitionMatrix(3, 1) * PUnder(:, :, 1) + this.modeTransitionMatrix(3, 2) * PUnder(:, :, 2) ...
                        + this.modeTransitionMatrix(3, 3) * PUnder(:, :, 3);
                    P_sum32 = this.modeTransitionMatrix(3, 1) * PUpper(:, :, 1) + this.modeTransitionMatrix(3, 2) * PUpper(:, :, 2)...
                        + this.modeTransitionMatrix(3, 3) * PUpper(:, :, 3);                 
                    
                    % psi matrix
                    Psi = kron(S_tilde(:, :, k) * (caModeProbs(1, k) * this.augV + this.augC * XUpper(:, :, 1, k) * this.augC') * S_tilde(:, :, k)', P_sum11) ...
                        + kron(S_tilde(:, :, k) * (caModeProbs(2, k) * this.augV + this.augC * XUpper(:, :, 2, k) * this.augC') * S_tilde(:, :, k)', P_sum21) ...
                        + kron(S_tilde(:, :, k) * (caModeProbs(3, k) * this.augV + this.augC * XUpper(:, :, 3, k) * this.augC') * S_tilde(:, :, k)', P_sum31);
                    
                    % phi matrix
                    Phi = kron(XUnder(:, :, 1, k), ...
                            this.augB(:, :, 1)' * P_sum12 * this.augB(:, :, 1) + this.J(:, :, 1)' * this.R * this.J(:, :, 1) ...
                            + (this.augB(:, :, 1)-Bexp(:, :, k))' * P_sum11 * (this.augB(:, :, 1)-Bexp(:, :, k)));
                    Phi = Phi + kron(XUnder(:, :, 2, k), ...
                            this.augB(:, :, 2)' * P_sum22 * this.augB(:, :, 2) + this.J(:, :, 2)' * this.R * this.J(:, :, 2) ...
                            + (this.augB(:, :, 2)-Bexp(:, :, k))' * P_sum21 * (this.augB(:, :, 2)-Bexp(:, :, k)));
                    Phi = Phi + kron(XUnder(:, :, 3, k), ...
                            this.augB(:, :, 3)' * P_sum32 * this.augB(:, :, 3) + this.J(:, :, 3)' * this.R * this.J(:, :, 3) ...
                            + (this.augB(:, :, 3)-Bexp(:, :, k))' * P_sum31 * (this.augB(:, :, 3)-Bexp(:, :, k)));                   
                    
                    % rho vector
                    rho = P_sum11 * this.augA(:, :, 1) * XUpper(:, :, 1, k) * this.augC' * S_tilde(:, :, k)' ...
                        + P_sum21 * this.augA(:, :, 2) * XUpper(:, :, 2, k) * this.augC' * S_tilde(:, :, k)' ...
                        + P_sum31 * this.augA(:, :, 3) * XUpper(:, :, 3, k) * this.augC' * S_tilde(:, :, k)';
                    % gamma vector
                    gamma = this.augB(:, :, 1)' * P_sum12 * this.augA(:, :, 1) * XUnder(:, :, 1, k) ...
                        + (this.augB(:, :, 1)-Bexp(:, :, k))' * P_sum11 * (this.augA(:, :, 1) - Aexp(:, :, k)) * XUnder(:, :, 1, k) ...
                        + this.augB(:, :, 2)' * P_sum22 * this.augA(:, :, 2) * XUnder(:, :, 2, k) ...
                        + (this.augB(:, :, 2)-Bexp(:, :, k))' * P_sum21 * (this.augA(:, :, 2) - Aexp(:, :, k)) * XUnder(:, :, 2, k) ...
                        + this.augB(:, :, 3)' * P_sum32 * this.augA(:, :, 3) * XUnder(:, :, 3, k) ...
                        + (this.augB(:, :, 3)-Bexp(:, :, k))' * P_sum31 * (this.augA(:, :, 3) - Aexp(:, :, k)) * XUnder(:, :, 3, k);           
                    
                    % construct system of linear equations to solve (Amat *x = b)                    
                    xk = lsqminnorm(Psi, rho(:));                    
                    K(:, :, k) = reshape(xk, dimAugState, []);                    
                  
                    xl = lsqminnorm(Phi, -gamma(:));                    
                    L(:, :, k) = reshape(xl, 2* this.dimU, []);                    
                    
                    % now we can update the costate with the computed gains
                    PUpper(:, :, 1) = L(:, :, k)' * this.J(:, :, 1)' * this.R * this.J(:, :, 1) * L(:, :, k) + this.augQ(:, :, 1) ...
                        + (this.augA(:, :, 1) -Aexp(:, :, k) + this.augB(:, :, 1) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k))' * P_sum11 * (this.augA(:, :, 1) -Aexp(:, :, k) + this.augB(:, :, 1) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k)) ...
                        + (this.augA(:, :, 1) + this.augB(:, :, 1) * L(:, :, k))' * P_sum12 * (this.augA(:, :, 1) + this.augB(:, :, 1) * L(:, :, k));                    
                    PUnder(:, :, 1) = (Aexp(:, :, k) + Bexp(:, :, k) * L(:, :, k) - K(:, :, k) * S_tilde(:, :, k) * this.augC - this.augB(:, :, 1) * L(:, :, k))' * P_sum11 * (Aexp(:, :, k) + Bexp(:, :, k) * L(:, :, k) - K(:, :, k) * S_tilde(:, :, k) * this.augC - this.augB(:, :, 1) * L(:, :, k)) ...
                        + (this.augB(:, :, 1) * L(:, :, k))' * P_sum12 * (this.augB(:, :, 1) * L(:, :, k)) ...
                        + L(:, :, k)' * this.J(:, :, 1)' * this.R * this.J(:, :, 1) * L(:, :, k);
                    omega(1) = omega_sum1 + trace((P_sum11 + P_sum12) * this.augW) ...
                        + trace(P_sum11 *  K(:, :, k)* S_tilde(:, :, k) * this.augV * (K(:, :, k)* S_tilde(:, :, k))');
                    % mode 2
                    PUpper(:, :, 2) = L(:, :, k)' * this.J(:, :, 2)' * this.R * this.J(:, :, 2) * L(:, :, k) + this.augQ(:, :, 2) ...
                        + (this.augA(:, :, 2) -Aexp(:, :, k) + this.augB(:, :, 2) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k))' * P_sum21 * (this.augA(:, :, 2) -Aexp(:, :, k) + this.augB(:, :, 2) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k)) ...
                        + (this.augA(:, :, 2) + this.augB(:, :, 2) * L(:, :, k))' * P_sum22 * (this.augA(:, :, 2) + this.augB(:, :, 2) * L(:, :, k));
                    PUnder(:, :, 2) = (Aexp(:, :, k) + Bexp(:, :, k) * L(:, :, k) - K(:, :, k) * S_tilde(:, :, k) * this.augC - this.augB(:, :, 2) * L(:, :, k))' * P_sum21 * (Aexp(:, :, k) + Bexp(:, :, k) * L(:, :, k) - K(:, :, k) * S_tilde(:, :, k) * this.augC - this.augB(:, :, 2) * L(:, :, k)) ...
                        + (this.augB(:, :, 2) * L(:, :, k))' * P_sum22 * (this.augB(:, :, 2) * L(:, :, k)) ...
                        + L(:, :, k)' * this.J(:, :, 2)' * this.R * this.J(:, :, 2) * L(:, :, k);
                    omega(2) = omega_sum2 + trace((P_sum21 + P_sum22) * this.augW) ...
                        + trace(P_sum21 *  K(:, :, k)* S_tilde(:, :, k) * this.augV * (K(:, :, k)* S_tilde(:, :, k))');
                    % mode 3
                    PUpper(:, :, 3) = L(:, :, k)' * this.J(:, :, 3)' * this.R * this.J(:, :, 3) * L(:, :, k) + this.augQ(:, :, 3) ...
                        + (this.augA(:, :, 3) -Aexp(:, :, k) + this.augB(:, :, 3) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k))' * P_sum31 * (this.augA(:, :, 3) -Aexp(:, :, k) + this.augB(:, :, 3) * L(:, :, k) - Bexp(:, :, k) * L(:, :, k)) ...
                        + (this.augA(:, :, 3) + this.augB(:, :, 3) * L(:, :, k))' * P_sum32 * (this.augA(:, :, 3) + this.augB(:, :, 3) * L(:, :, k));
                    PUnder(:, :, 3) = (Aexp(:, :, k) + Bexp(:, :, k) * L(:, :, k) - K(:, :, k) * S_tilde(:, :, k) * this.augC - this.augB(:, :, 3) * L(:, :, k))' * P_sum31 * (Aexp(:, :, k) + Bexp(:, :, k) * L(:, :, k) - K(:, :, k) * S_tilde(:, :, k) * this.augC - this.augB(:, :, 3) * L(:, :, k)) ...
                        + (this.augB(:, :, 3) * L(:, :, k))' * P_sum32 * (this.augB(:, :, 3) * L(:, :, k)) ...
                        + L(:, :, k)' * this.J(:, :, 3)' * this.R * this.J(:, :, 3) * L(:, :, k);
                    omega(3) = omega_sum3 + trace((P_sum31 + P_sum32) * this.augW) ...
                        + trace(P_sum31 *  K(:, :, k)* S_tilde(:, :, k) * this.augV * (K(:, :, k)* S_tilde(:, :, k))');
                end
                
                costs = trace(PUpper(:, :, 1) * (XUpper(:, :, 1, 1) + XUnder(:, :, 1, 1)) + PUnder(:, :, 1) * XUpper(:, :, 1, 1)) ...
                    + trace(PUpper(:, :, 2) * (XUpper(:, :, 2, 1) + XUnder(:, :, 2, 1)) + PUnder(:, :, 2) * XUpper(:, :, 2, 1)) ...
                    + trace(PUpper(:, :, 3) * (XUpper(:, :, 3, 1) + XUnder(:, :, 3, 1)) + PUnder(:, :, 3) * XUpper(:, :, 3, 1)) ...
                    + caModeProbs(1, 1) * omega(1) + caModeProbs(2, 1) * omega(2) + caModeProbs(3, 1) * omega(3);               
                oldCosts = costs;  
            end
        end
        
        %% computeGainsMeasurementAndMode
        function [K, L, Aexp, Bexp] = computeGainsMeasurementAndMode(this)            
            T = this.scDelayTransMat;
            % both measurements are available, so mode is 1 1 -> 3
            % horizon is three            
            
            S_tilde = zeros(2 * this.dimY,2 * this.dimY, this.horizonLength);
            S_tilde(:, :, 1) = this.S(:, :, 4); % initially, mode is 3
            % k+1
            availProbs = zeros(1, 4);
            % we need delay probs at time k (packet arrived at time k, so
            % delay is 0)
            measDelayProbs = [1 0 0 0 0 0];
            % for mode 0 0 (y_k and y_k+1 not available)
            % delay of y_k ~= 1 and delay of y_k+1 ~= 0            
            for j=[2:6]                
                for m=[1 3:6]            
                    availProbs(1) = availProbs(1) + T(m,j)* measDelayProbs(m);            
                end
            end
            % for mode 1 0 (y_k+1 available and y_k not available)
            % delay of y_k ~= 1 and delay of y_k+1 == 0
            for m=[1 3:6]                       
                availProbs(2) = availProbs(2) + T(m, 1)* measDelayProbs(m);
            end
            % for mode 0 1 (y_k+1 not available and y_k available)
            % delay of y_k == 1 and delay of y_k+1 ~= 0
            for j=[2:6]                           
                availProbs(3) = availProbs(3) + T(2,j)* measDelayProbs(2);
            end
            % for mode 1 1 (y_k available and y_k+1 available)
            % delay of y_k == 1 and delay of y_k+1 == 0
            availProbs(4) =  T(2,1)* measDelayProbs(2);            
            S_tilde(:, :, 2) = sum(reshape(availProbs, 1, 1, []) .* this.S, 3);
             
            % k+2
            availProbs = zeros(1, 4);
            % we need delay probs at time k+1 
            measDelayProbs = measDelayProbs * T;
            % for mode 0 0 (y_k and y_k+1 not available)
            % delay of y_k ~= 1 and delay of y_k+1 ~= 0            
            for j=[2:6]                
                for m=[1 3:6]            
                    availProbs(1) = availProbs(1) + T(m,j)* measDelayProbs(m);            
                end
            end
            % for mode 1 0 (y_k+1 available and y_k not available)
            % delay of y_k ~= 1 and delay of y_k+1 == 0
            for m=[1 3:6]                       
                availProbs(2) = availProbs(2) + T(m, 1)* measDelayProbs(m);
            end
            % for mode 0 1 (y_k+1 not available and y_k available)
            % delay of y_k == 1 and delay of y_k+1 ~= 0
            for j=[2:6]                           
                availProbs(3) = availProbs(3) + T(2,j)* measDelayProbs(2);
            end
            % for mode 1 1 (y_k available and y_k+1 available)
            % delay of y_k == 1 and delay of y_k+1 == 0
            availProbs(4) =  T(2,1)* measDelayProbs(2);
            
            S_tilde(:, :, 3) = sum(reshape(availProbs, 1, 1, []) .* this.S, 3);            
               
            % ca Mode: we have observed the previous mode to be the first
            numCaModes = 3; % sequence length + 1
            caModeProbs = zeros(numCaModes, this.horizonLength + 1);
            caModeProbs(:, 1) = this.modeTransitionMatrix' * [1; 0; 0];
            caModeProbs(:, 2) = this.modeTransitionMatrix' * caModeProbs(:, 1);
            caModeProbs(:, 3) = this.modeTransitionMatrix' * caModeProbs(:, 2);
            caModeProbs(:, 4) = this.modeTransitionMatrix' * caModeProbs(:, 3);
            
            [K, L, Aexp, Bexp] = this.computeGains(caModeProbs, S_tilde);
        end
        
        %% computeGainsNoModes
        function [K, L, Aexp, Bexp] = computeGainsNoModes(this)
            T = this.scDelayTransMat;            
            % the current measurement is available, so mode is 1 0 -> 1 as per de2bi
            % horizon is three            
            
            S_tilde = zeros(2 * this.dimY,2 * this.dimY, this.horizonLength);
            S_tilde(:, :, 1) = this.S(:, :, 2); % initially, mode is 1 (1 0)
            % k+1
            availProbs = zeros(1, 4);
            measDelayProbs = [1 0 0 0 0 0];
            % for mode 0 0 (y_k and y_k+1 not available)
            % delay of y_k ~= 1 and delay of y_k+1 ~= 0            
            for j=[2:6]                
                for m=[1 3:6]            
                    availProbs(1) = availProbs(1) + T(m,j)* measDelayProbs(m);            
                end
            end
            % for mode 1 0 (y_k+1 available and y_k not available)
            % delay of y_k ~= 1 and delay of y_k+1 == 0
            for m=[1 3:6]                       
                availProbs(2) = availProbs(2) + T(m, 1)* measDelayProbs(m);
            end
            % for mode 0 1 (y_k+1 not available and y_k available)
            % delay of y_k == 1 and delay of y_k+1 ~= 0
            for j=[2:6]                           
                availProbs(3) = availProbs(3) + T(2,j)* measDelayProbs(2);
            end
            % for mode 1 1 (y_k available and y_k+1 available)
            % delay of y_k == 1 and delay of y_k+1 == 0
            availProbs(4) =  T(2,1)* measDelayProbs(2);            
            S_tilde(:, :, 2) = sum(reshape(availProbs, 1, 1, []) .* this.S, 3);
                        
            % k+2
            availProbs = zeros(1, 4);
            % we need delay probs at time k+1 
            measDelayProbs = measDelayProbs * T;
            % for mode 0 0 (y_k and y_k+1 not available)
            % delay of y_k ~= 1 and delay of y_k+1 ~= 0            
            for j=[2:6]                
                for m=[1 3:6]            
                    availProbs(1) = availProbs(1) + T(m,j)* measDelayProbs(m);            
                end
            end
            % for mode 1 0 (y_k+1 available and y_k not available)
            % delay of y_k ~= 1 and delay of y_k+1 == 0
            for m=[1 3:6]                       
                availProbs(2) = availProbs(2) + T(m, 1)* measDelayProbs(m);
            end
            % for mode 0 1 (y_k+1 not available and y_k available)
            % delay of y_k == 1 and delay of y_k+1 ~= 0
            for j=[2:6]                           
                availProbs(3) = availProbs(3) + T(2,j)* measDelayProbs(2);
            end
            % for mode 1 1 (y_k available and y_k+1 available)
            % delay of y_k == 1 and delay of y_k+1 == 0
            availProbs(4) =  T(2,1)* measDelayProbs(2);            
            S_tilde(:, :, 3) = sum(reshape(availProbs, 1, 1, []) .* this.S, 3); 
            
            % ca Mode: initially we are in last mode, as we have none
            % observed
            numCaModes = 3;
            caModeProbs = zeros(numCaModes, this.horizonLength + 1);
            caModeProbs(:, 1) = [0 0 1]';
            caModeProbs(:, 2) = this.modeTransitionMatrix' * caModeProbs(:, 1);
            caModeProbs(:, 3) = this.modeTransitionMatrix' * caModeProbs(:, 2);
            caModeProbs(:, 4) = this.modeTransitionMatrix' * caModeProbs(:, 3);
            
            [K, L, Aexp, Bexp] = this.computeGains(caModeProbs, S_tilde);
        end
        
        %% computeGainsNoMeasurementsNoModes
        function [K, L, Aexp, Bexp] = computeGainsNoMeasurementsNoModes(this)
            T = this.scDelayTransMat;            
            % no measurements available, so mode is 0 0 -> 0 as per de2bi
            % horizon is three    
            
            S_tilde = zeros(2 * this.dimY,2 * this.dimY, this.horizonLength);
            S_tilde(:, :, 1) = this.S(:, :, 1); % initially, mode is 0 (0 0)
            % k+1
            availProbs = zeros(1, 4);
            % we need delay probs at time k (packet did not arrive, so
            % delay is not known, we only know it's larger than 0)
            measDelayProbs = 1/5 * [0 1 1 1 1 1]; % at time k
            % for mode 0 0 (y_k and y_k+1 not available)
            % delay of y_k ~= 1 and delay of y_k+1 ~= 0            
            for j=[2:6]                
                for m=[1 3:6]            
                    availProbs(1) = availProbs(1) + T(m,j)* measDelayProbs(m);            
                end
            end
            % for mode 1 0 (y_k+1 available and y_k not available)
            % delay of y_k ~= 1 and delay of y_k+1 == 0
            for m=[1 3:6]                       
                availProbs(2) = availProbs(2) + T(m, 1)* measDelayProbs(m);
            end
            % for mode 0 1 (y_k+1 not available and y_k available)
            % delay of y_k == 1 and delay of y_k+1 ~= 0
            for j=[2:6]                           
                availProbs(3) = availProbs(3) + T(2,j)* measDelayProbs(2);
            end
            % for mode 1 1 (y_k available and y_k+1 available)
            % delay of y_k == 1 and delay of y_k+1 == 0
            availProbs(4) =  T(2,1)* measDelayProbs(2);            
            S_tilde(:, :, 2) = sum(reshape(availProbs, 1, 1, []) .* this.S, 3);
                        
            % k+2
            availProbs = zeros(1, 4);
            % we need delay probs at time k+1 
            measDelayProbs = measDelayProbs * T;
            % for mode 0 0 (y_k and y_k+1 not available)
            % delay of y_k ~= 1 and delay of y_k+1 ~= 0            
            for j=[2:6]                
                for m=[1 3:6]            
                    availProbs(1) = availProbs(1) + T(m,j)* measDelayProbs(m);            
                end
            end
            % for mode 1 0 (y_k+1 available and y_k not available)
            % delay of y_k ~= 1 and delay of y_k+1 == 0
            for m=[1 3:6]                       
                availProbs(2) = availProbs(2) + T(m, 1)* measDelayProbs(m);
            end
            % for mode 0 1 (y_k+1 not available and y_k available)
            % delay of y_k == 1 and delay of y_k+1 ~= 0
            for j=[2:6]                           
                availProbs(3) = availProbs(3) + T(2,j)* measDelayProbs(2);
            end
            % for mode 1 1 (y_k available and y_k+1 available)
            % delay of y_k == 1 and delay of y_k+1 == 0
            availProbs(4) =  T(2,1)* measDelayProbs(2);            
            S_tilde(:, :, 3) = sum(reshape(availProbs, 1, 1, []) .* this.S, 3); 
            
            % ca Mode: initially we are in last mode, as we have none
            % observed
            numCaModes = 3;
            caModeProbs = zeros(numCaModes, this.horizonLength + 1);
            caModeProbs(:, 1) = [0 0 1]';
            caModeProbs(:, 2) = this.modeTransitionMatrix' * caModeProbs(:, 1);
            caModeProbs(:, 3) = this.modeTransitionMatrix' * caModeProbs(:, 2);
            caModeProbs(:, 4) = this.modeTransitionMatrix' * caModeProbs(:, 3);
            
            [K, L, Aexp, Bexp] = this.computeGains(caModeProbs, S_tilde);
        end
    end
    
    methods(Test)
        %% testRecedingHorizonUdpLikeControllerInvalidSystemMatrix
        function testRecedingHorizonUdpLikeControllerInvalidSystemMatrix(this)
             expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
             
             invalidSysMatrix = eye(this.dimX, this.dimX + 1); % not square
             this.verifyError(@() RecedingHorizonUdpLikeController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
             
             invalidSysMatrix = eye(this.dimX, this.dimX); % square but not finite
             invalidSysMatrix(1, end) = inf;
             this.verifyError(@() RecedingHorizonUdpLikeController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidInputMatrix
        function testRecedingHorizonUdpLikeControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidInputMatrix = eye(this.dimX +1, this.dimU); % invalid dims
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
             
            invalidInputMatrix = eye(this.dimX, this.dimU); % correct dims, but not finite
            invalidInputMatrix(1, end) = nan;
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidCostMatrices
        function testRecedingHorizonUdpLikeControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = eye(this.dimX); % Q is not symmetric
            invalidQ(1, end) = 1;
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidMeasurementMatrix
        function testRecedingHorizonUdpLikeControllerInvalidMeasurementMatrix(this)
            expectedErrId = 'Validator:ValidateMeasurementMatrix:InvalidMeasMatrix';
            
            invalidMeasMatrix = eye(this.dimY, this.dimX + 1); % invalid dims
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
             
            invalidMeasMatrix = eye(this.dimY, this.dimX); % correct dims, but not finite
            invalidMeasMatrix(1, end) = nan;
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidCaModeTransitionMatrix
        function testRecedingHorizonUdpLikeControllerInvalidCaModeTransitionMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';
            
            invalidModeTransitionMatrix = blkdiag(1, this.modeTransitionMatrix);% invalid dimensions
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidModeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidModeTransitionMatrix = [0 0.1 0.8 0.2];% not a matrix
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidModeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
                     
            invalidModeTransitionMatrix = this.modeTransitionMatrix;
            invalidModeTransitionMatrix(1,1) = -invalidModeTransitionMatrix(1,1); % negative entry
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidModeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidModeTransitionMatrix = this.modeTransitionMatrix;
            invalidModeTransitionMatrix(1,1) = 1.1; % does not sum up to 1
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                invalidModeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidScDelayProbs
        function testRecedingHorizonUdpLikeControllerInvalidScDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, invalidDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidMaxMeasDelay
        function testRecedingHorizonUdpLikeControllerInvalidMaxMeasDelay(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:InvalidMaxMeasDelay';
            
            invalidMaxMeasDelay = [1 2]; % not a scalar
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidMaxMeasDelay = -1; % negative scalar
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidMaxMeasDelay = 1.5; % not an integer
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidMaxMeasDelay = inf; % not finite
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, invalidMaxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidSysNoiseCovariance
        function testRecedingHorizonUdpLikeControllerInvalidSysNoiseCovariance(this)
            expectedErrId = 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCovDim';
            
            invalidW = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, invalidW, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
                       
            invalidW = ones(this.dimU); % W is not pd
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, invalidW, this.V, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId); 
        end
        
        
        %% testRecedingHorizonUdpLikeControllerInvalidMeasNoiseCovariance
        function testRecedingHorizonUdpLikeControllerInvalidMeasNoiseCovariance(this)
            expectedErrId = 'Validator:ValidateMeasNoiseCovarianceMatrix:InvalidCovDim';
            
            invalidV = eye(this.dimY + 1, this.dimY); % not square
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, invalidV, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
                       
            invalidV = ones(this.dimY); % V is not pd
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, invalidV, this.horizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
         %% testRecedingHorizonUdpLikeControllerInvalidHorizonLegth
        function testRecedingHorizonUdpLikeControllerInvalidHorizonLegth(this)
            expectedErrId = 'Validator:ValidateHorizonLength:InvalidHorizonLength';
            
            invalidHorizonLength = [1 2]; % not a scalar
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidHorizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidHorizonLength = -1; % negative scalar
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidHorizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidHorizonLength = 0; 
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidHorizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidHorizonLength = 1.5; % not an integer
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidHorizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
            
            invalidHorizonLength = inf; % not finite
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, invalidHorizonLength, this.x0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidX0
        function testRecedingHorizonUdpLikeControllerInvalidX0(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:InvalidX0';
            
            invalidX0 = ones(this.dimX +1); % not a vector
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, invalidX0, this.x0Cov), ...
                expectedErrId);
            
            invalidX0 = ones(this.dimX +1, 1); % wrong dimensions
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, invalidX0, this.x0Cov), ...
                expectedErrId);
            
            invalidX0 = ones(this.dimX, 1); % not finite
            invalidX0(1) = nan;
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, invalidX0, this.x0Cov), ...
                expectedErrId);
        end
        
        %% testRecedingHorizonUdpLikeControllerInvalidX0Cov
        function testRecedingHorizonUdpLikeControllerInvalidX0Cov(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:InvalidX0Cov';
            
            invalidX0Cov = ones(this.dimX +1, 1); % not a matrix
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, invalidX0Cov), ...
                expectedErrId);
            
            invalidX0Cov = ones(this.dimX +1, this.dimX); % wrong dimensions
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, invalidX0Cov), ...
                expectedErrId);
            
            invalidX0Cov = ones(this.dimX); % not pd
            this.verifyError(@() RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, invalidX0Cov), ...
                expectedErrId);
        end
%%
%%
        %% RecedingHorizonUdpLikeController
        function testRecedingHorizonUdpLikeController(this)
            controller = RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov);
                        
            this.verifyEqual(controller.getControllerPlantState(), this.x0);
            this.verifyEqual(controller.maxMeasurementDelay, this.maxMeasDelay);
            this.verifyEqual(controller.horizonLength, this.horizonLength);
            this.verifyEqual(controller.sequenceLength, this.sequenceLength);
            this.verifyEqual(controller.lastNumIterations, 0);
            this.verifyFalse(controller.requiresExternalStateEstimate); % does not require a filter or state feedback
            
            %by default, we use the mex implementation to obtain the
            %controller gains
            this.verifyTrue(controller.useMexImplementation);
            
            % check if augmented dynamics is properly initialized
            this.verifyEqual(controller.augA, this.augA);
            this.verifyEqual(controller.augB, this.augB);
            this.verifyEqual(controller.augW, this.augW);
            this.verifyEqual(controller.augQ, this.augQ);
            this.verifyEqual(controller.S, this.S);
            this.verifyTrue(issparse(controller.augmentedMeasMatrix));
            this.verifyEqual(full(controller.augmentedMeasMatrix), this.augC);
            
            newMaxMeasDelay = 2;
            controller = RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, newMaxMeasDelay, this.W, this.V, this.horizonLength, this.x0, this.x0Cov);
            % sc delay probs stored by controller have M+2 =4 elements
            
            % check if history of mode transition matrices is created
            % correctly (maxMeasHistory is 2)
            expectedHistory = cat(3, this.modeTransitionMatrix, this.modeTransitionMatrix);
            this.verifyEqual(controller.transitionMatrixCaHistory, expectedHistory);
            
            % check the measurement availability
            expectedStates = [0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0  1 1; 1 1 1];
            this.verifyEqual(controller.measAvailabilityStates, expectedStates);
            
            % check the corresponding indices (index is delay +1)
            % 0 0 0 ~=1 ~=2 ~=3
            % 1 0 0 ==1 ~=2 ~=3
            % 0 1 0 ~=1 ==2 ~=3
            % 1 1 0 ==1 ==2 ~=3
            % 0 0 1 ~=1 ~=2 ==3
            % 1 0 1 ==1 ~=2 ==3
            % 0 1 1 ~=1 ==2 ==3
            % 1 1 1 ==1 ==2 ==3
            this.verifyTrue(iscell(controller.measAvailabilityIdx));
            this.verifySize(controller.measAvailabilityIdx, [3, 8]);
            this.verifyEqual(controller.measAvailabilityIdx{1, 1}, [2:4]);
            this.verifyEqual(controller.measAvailabilityIdx{2, 1}, [1 3 4]);
            this.verifyEqual(controller.measAvailabilityIdx{3, 1}, [1 2 4]);
            this.verifyEqual(controller.measAvailabilityIdx{1, 2}, 1);
            this.verifyEqual(controller.measAvailabilityIdx{2, 2}, [1 3 4]);
            this.verifyEqual(controller.measAvailabilityIdx{3, 2}, [1 2 4]);
            this.verifyEqual(controller.measAvailabilityIdx{1, 3}, [2:4]);
            this.verifyEqual(controller.measAvailabilityIdx{2, 3}, 2);
            this.verifyEqual(controller.measAvailabilityIdx{3, 3}, [1 2 4]);
            this.verifyEqual(controller.measAvailabilityIdx{1, 4}, 1);
            this.verifyEqual(controller.measAvailabilityIdx{2, 4}, 2);
            this.verifyEqual(controller.measAvailabilityIdx{3, 4}, [1 2 4]);
            this.verifyEqual(controller.measAvailabilityIdx{1, 5}, [2:4]);
            this.verifyEqual(controller.measAvailabilityIdx{2, 5}, [1 3 4]);
            this.verifyEqual(controller.measAvailabilityIdx{3, 5}, 3);
            this.verifyEqual(controller.measAvailabilityIdx{1, 6}, 1);
            this.verifyEqual(controller.measAvailabilityIdx{2, 6}, [1 3 4]);
            this.verifyEqual(controller.measAvailabilityIdx{3, 6}, 3);
            this.verifyEqual(controller.measAvailabilityIdx{1, 7}, [2:4]);
            this.verifyEqual(controller.measAvailabilityIdx{2, 7}, 2);
            this.verifyEqual(controller.measAvailabilityIdx{3, 7}, 3);
            this.verifyEqual(controller.measAvailabilityIdx{1, 8}, 1);
            this.verifyEqual(controller.measAvailabilityIdx{2, 8}, 2);
            this.verifyEqual(controller.measAvailabilityIdx{3, 8}, 3);
            
            % now check masurement availability matrices S
            % 0 0 0 ~=1 ~=2 ~=3
            % 1 0 0 ==1 ~=2 ~=3
            % 0 1 0 ~=1 ==2 ~=3
            % 1 1 0 ==1 ==2 ~=3
            % 0 0 1 ~=1 ~=2 ==3
            % 1 0 1 ==1 ~=2 ==3
            % 0 1 1 ~=1 ==2 ==3
            % 1 1 1 ==1 ==2 ==3
            % here 1 indicates that measurement is available for processing
            expS = zeros(3 * this.dimY, 3 * this.dimY, 8);
            expS(:, :, 2) = blkdiag(eye(this.dimY), zeros(this.dimY), zeros(this.dimY));
            expS(:, :, 3) = blkdiag(zeros(this.dimY), eye(this.dimY), zeros(this.dimY));
            expS(:, :, 4) = blkdiag(eye(this.dimY), eye(this.dimY), zeros(this.dimY));
            expS(:, :, 5) = blkdiag(zeros(this.dimY), zeros(this.dimY), eye(this.dimY));
            expS(:, :, 6) = blkdiag(eye(this.dimY), zeros(this.dimY), eye(this.dimY));
            expS(:, :, 7) = blkdiag(zeros(this.dimY), eye(this.dimY), eye(this.dimY));
            expS(:, :, 8) = blkdiag(eye(this.dimY), eye(this.dimY), eye(this.dimY));
            
            this.verifyTrue(issparse(controller.augmentedMeasMatrix));
            this.verifyEqual(controller.S, expS);
        end
%%
%%      
        
        %% testSetInitialControllerGainsInvalidK
        function testSetInitialControllerGainsInvalidK(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:SetInitialControllerGains:InvalidK';
                        
            invalidK = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.setInitialControllerGains(invalidK, this.initialL), ...
                expectedErrId);
            
            invalidK = this.initialL; % wrong dimension
            this.verifyError(@() this.controllerUnderTest.setInitialControllerGains(invalidK, this.initialL), ...
                expectedErrId);
            
            invalidK = this.initialK(:, :, 1:end-1); % wrong number of slices
            this.verifyError(@() this.controllerUnderTest.setInitialControllerGains(invalidK, this.initialL), ...
                expectedErrId);
        end
        
         %% testSetInitialControllerGainsInvalidL
        function testSetInitialControllerGainsInvalidL(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:SetInitialControllerGains:InvalidL';
                        
            invalidL = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.setInitialControllerGains(this.initialK, invalidL), ...
                expectedErrId);
            
            invalidL = this.initialK; % wrong dimension
            this.verifyError(@() this.controllerUnderTest.setInitialControllerGains(this.initialK, invalidL), ...
                expectedErrId);
            
            invalidL = this.initialL(:, :, 1); % wrong number of slices
            this.verifyError(@() this.controllerUnderTest.setInitialControllerGains(this.initialK, invalidL), ...
                expectedErrId);
        end
        
        %% testSetInitialControllerGains
        function testSetInitialControllerGains(this)
            % this call should not crash
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);
            
            this.verifyEqual(this.controllerUnderTest.L, this.initialL);
            this.verifyEqual(this.controllerUnderTest.K, this.initialK);
        end
%%
%%
        %% testSetControllerPlantStateInvalidState
        function testSetControllerPlantStateInvalidState(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:SetControllerPlantState:InvalidState';
            
            invalidState = this; % not a Distribution
            this.verifyError(@() this.controllerUnderTest.setControllerPlantState(invalidState), expectedErrId);
            
            invalidState = Gaussian(0, 1); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.setControllerPlantState(invalidState), expectedErrId);
        end
        
        %% testSetControllerPlantState
        function testSetControllerPlantState(this)
            stateMean = 0.2 * ones(this.dimX, 1);                        
            initialState = Gaussian(stateMean, 42 * eye(this.dimX));
            
            this.controllerUnderTest.setControllerPlantState(initialState);
            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), stateMean);
        end
%%
%%
        %% testChangeCaDelayProbs
        function testChangeCaDelayProbs(this)
            oldTransitionMatrix = this.modeTransitionMatrix;
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);

            this.assertTrue(isa(this.controllerUnderTest, 'CaDelayProbsChangeable'));
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);                       
            this.assertEqual(this.controllerUnderTest.transitionMatrixCa, oldTransitionMatrix);
            this.assertEqual(this.controllerUnderTest.transitionMatrixCaHistory, oldTransitionMatrix);
           
            % now change the ca delay probs
            newCaDelayProbs = ones(1, 10) / 10;        

            this.controllerUnderTest.changeCaDelayProbs(newCaDelayProbs);            
            
            % we observed the previous mode and received two measurements
            previousMode = 1;
            modeDelay = 1;
            receivedMeasurements = [ones(this.dimY, 1) 1.5 * ones(this.dimY, 1)];
            measurementDelays = [1 0];

            % expected controller gains, based on the new transition matrix
            this.modeTransitionMatrix = Utility.calculateDelayTransitionMatrix(...
                    Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.sequenceLength + 1));                         
            [K, L, Aexp, Bexp] = this.computeGainsMeasurementAndMode();
         
            expectedInputSequence = L(:, :, 1) * this.augX0;
            expectedMeasurement = this.augC * this.augX0;
            innovation = [receivedMeasurements(:, 2); receivedMeasurements(:, 1)] - expectedMeasurement;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence(receivedMeasurements, measurementDelays, previousMode, modeDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check the side effect
            this.verifyEqual(this.controllerUnderTest.transitionMatrixCa, this.modeTransitionMatrix);
            this.verifyEqual(this.controllerUnderTest.transitionMatrixCaHistory, this.modeTransitionMatrix);
        end    
%%
%%
        %% testChangeScDelayProbs
        function testChangeScDelayProbs(this)
            oldTransitionMatrix = this.truncatedScDelayTransMat;            
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);

            this.assertTrue(isa(this.controllerUnderTest, 'ScDelayProbsChangeable'));
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);                       
            this.assertEqual(this.controllerUnderTest.scDelayTransitionMatrix, oldTransitionMatrix, 'AbsTol', 1e-15);
            this.assertEqual(this.controllerUnderTest.transitionMatrixScHistory, oldTransitionMatrix, 'AbsTol', 1e-15);
           
            % now change the sc delay probs
            newScDelayProbs = ones(1, 10) / 10;        
                        
            this.controllerUnderTest.changeScDelayProbs(newScDelayProbs);            
            
            % we observed the previous mode and received two measurements
            previousMode = 1;
            modeDelay = 1;
            receivedMeasurements = [ones(this.dimY, 1) 1.5 * ones(this.dimY, 1)];
            measurementDelays = [1 0];

            % expected controller gains, based on the new transition matrix
            newScDelayProbs = Utility.truncateDiscreteProbabilityDistribution(newScDelayProbs, size(this.scDelayTransMat, 1));
            newTruncatedScDelayProbs = Utility.truncateDiscreteProbabilityDistribution(newScDelayProbs, numel(this.truncatedScDelayProbs));
            this.truncatedScDelayTransMat = repmat(reshape(newTruncatedScDelayProbs, 1, []), numel(newTruncatedScDelayProbs), 1);            
            this.scDelayTransMat = repmat(reshape(newScDelayProbs, 1, []), numel(newScDelayProbs), 1);
            
            [K, L, Aexp, Bexp] = this.computeGainsMeasurementAndMode();
         
            expectedInputSequence = L(:, :, 1) * this.augX0;
            expectedMeasurement = this.augC * this.augX0;
            innovation = [receivedMeasurements(:, 2); receivedMeasurements(:, 1)] - expectedMeasurement;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence(receivedMeasurements, measurementDelays, previousMode, modeDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check the side effect
            this.verifyEqual(this.controllerUnderTest.scDelayTransitionMatrix, this.truncatedScDelayTransMat);
            this.verifyEqual(this.controllerUnderTest.transitionMatrixScHistory, this.truncatedScDelayTransMat);
        end    
%%
%%
        %% testChangeModelParametersInvalidAMatrix
        function testChangeModelParametersInvalidAMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidDimensions';
            
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B, this.W), ...
                expectedErrId);
            
            invalidA = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B, this.W), ...
                expectedErrId);
            
            invalidA = eye(this.dimX); % correct dims, but inf
            invalidA(end, end) = inf;
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B, this.W), ...
                expectedErrId);
        end
        
        %% testChangeModelParametersInvalidBMatrix
        function testChangeModelParametersInvalidBMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrixDims';
            
            invalidB = eye(this.dimX, this.dimU + 1); % invalid dims
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB, this.W), ...
                expectedErrId);            
             
            invalidB = eye(this.dimX, this.dimU); % correct dims, but inf
            invalidB(end, end) = inf;
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB, this.W), ...
                expectedErrId);
        end
        
        %% testChangeModelParametersInvalidWMatrix
        function testChangeModelParametersInvalidWMatrix(this)
            expectedErrId = 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCovDim';
            
            invalidW = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, this.B, invalidW), ...
                expectedErrId);
            
            invalidW = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, this.B, invalidW), ...
                expectedErrId);
            
            invalidW = zeros(this.dimX); % matrix is square, but not positive definite
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, this.B, invalidW), ...
                expectedErrId);
        end
        
        %% testChangeModelParametersNewA
        function testChangeModelParametersNewA(this)
            dimEta = 2;
            
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
            
            % now change A, affects augA (A_tilde in the paper) 
            newA = this.A + eye(this.dimX);                        
            
            this.controllerUnderTest.changeModelParameters(newA, this.B, this.W);
            % now process one measurement and compute a control sequence
            % using new system matrix A
            
            [F, ~, H, ~] = Utility.createActuatorMatrices(this.sequenceLength, this.dimU);            
            D = eye(this.dimX);
            E = zeros(this.dimX);
            
            % this is the new augA (A_tilde)
            this.A = newA;
            this.augA = repmat(blkdiag([newA zeros(this.dimX); D E], F), 1,1, this.sequenceLength + 1);
            this.augA(1:this.dimX, end - dimEta + 1:end, :) = mtimesx(this.B, H);

            % one measurement received, without delay
            receivedMeasurement = ones(this.dimY, 1);
            measurementDelay = 0;
            
            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsNoModes();
            
            expectedInputSequence = L(:, :, 1) * this.augX0;
            innovation = [receivedMeasurement; zeros(this.dimY, 1)] - this.S(:, :, 2) * this.augC * this.augX0;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
                        
            actualInputSequence = this.controllerUnderTest.computeControlSequence(receivedMeasurement, measurementDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
        end
        
        %% testChangeModelParametersNewB
        function testChangeModelParametersNewB(this)
            dimEta = 2;
            dimAugState = 2 * this.dimX + dimEta;
            numModes = this.sequenceLength + 1;
            
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
            
            % now change B, affects both augA and augB
            newB = -this.B;
            
            this.controllerUnderTest.changeModelParameters(this.A, newB, this.W);
            % now process one measurement and compute a control sequence
            % using new input matrix B
            
            [F, G, H, ~] = Utility.createActuatorMatrices(this.sequenceLength, this.dimU);           
            D = eye(this.dimX);
            E = zeros(this.dimX);
            
            % this is the new augA (A_tilde)
            this.B = newB;
            this.augA = repmat(blkdiag([this.A zeros(this.dimX); D E], F), 1,1, numModes);
            this.augA(1:this.dimX, end - dimEta + 1:end, :) = mtimesx(newB, H);
            % this is the new augB (B_tilde)
            this.augB = repmat([zeros(dimAugState - size(G, 1), size(G, 2)); G], 1, 1, numModes);
            this.augB(1:this.dimX, :, :) = mtimesx(newB, this.J);
            
            receivedMeasurement = ones(this.dimY, 1);
            measurementDelay = 0;
            
            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsNoModes();
            
            expectedInputSequence = L(:, :, 1) * this.augX0;
            innovation = [receivedMeasurement; zeros(this.dimY, 1)] - this.S(:, :, 2) * this.augC * this.augX0;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
                        
            actualInputSequence = this.controllerUnderTest.computeControlSequence(receivedMeasurement, measurementDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
        end
        
        %% testChangeModelParametersNewW
        function testChangeModelParametersNewW(this)
            dimEta = 2;
            
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
               
            % reduce the noise
            newW = 0.75 * this.W;
            this.augW = blkdiag(newW, zeros(this.dimX), zeros(dimEta));            
            
            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsNoMeasurementsNoModes();
     
            expectedInputSequence = L(:, :, 1) * this.augX0;            
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence;
            
            this.controllerUnderTest.changeModelParameters(this.A, this.B, newW);            
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence();
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol); 
        end
        
%%
%%
        %% testGetLastComputationMeasurementDataNoMeasurements
        function testGetLastComputationMeasurementDataNoMeasurements(this)
            % compute a sequence without using measurements or mode
            % observations
            this.controllerUnderTest.computeControlSequence();
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
            
            % compute a sequence without using measurements but with mode
            % observations
            modeObservation = 1;
            modeDelay = 1;
            this.controllerUnderTest.computeControlSequence([], [], modeObservation, modeDelay);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 0);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
         %% testGetLastComputationMeasurementDataNoMeasDiscarded
        function testGetLastComputationMeasurementDataNoMeasDiscarded(this)
            measurements = ones(this.dimY, 2); % two meas, column-wise arranged
            delays = [0 1];
            
            this.controllerUnderTest.computeControlSequence(measurements, delays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 0);
            
            % now the same but with a single mode observation
            % should not result in a change
            modeObservation = 1;
            modeDelay = 1;
            this.controllerUnderTest.computeControlSequence(measurements, delays, modeObservation, modeDelay);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 0);
        end
        
        %% testGetLastComputationMeasurementData
        function testGetLastComputationMeasurementData(this)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            
            this.applyFixture(...
                SuppressedWarningsFixture('RecedingHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsDelayTooLarge'));
            
            measurements = ones(this.dimY, 4); % four meas, column-wise arranged
            delays = [0 2 1 3]; % two meas are too old
            
            this.controllerUnderTest.computeControlSequence(measurements, delays);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 2);
            
            % now the same but with a single mode observation
            % should not result in a change
            modeObservation = 1;
            modeDelay = 1;
            this.controllerUnderTest.computeControlSequence(measurements, delays, modeObservation, modeDelay);
            
            [actualNumUsedMeas, actualNumDiscardedMeas] ...
                = this.controllerUnderTest.getLastComputationMeasurementData();
            
            this.verifyEqual(actualNumUsedMeas, 2);
            this.verifyEqual(actualNumDiscardedMeas, 2);
        end
%%
%%
        %% testComputeControlSequenceInvalidMeasurements
        function testComputeControlSequenceInvalidMeasurements(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeas';
            
            % two measurements, column-wise arranged, but of wrong dimension
            invalidMeasurements = ones(this.dimY + 1, 2);
            delays = [0 1];            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(invalidMeasurements, delays), ...
                expectedErrId);
            
            % two measurements, column-wise arranged, but of wrong dimension
            invalidMeasurements = ones(this.dimY + 1, 2);
            delays = [0 1];
            % existence of mode observation shoud not change anything
            modeObservation = 2;
            modeDelay = 1;
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(invalidMeasurements, delays, modeObservation, modeDelay), ...
                expectedErrId);
        end

        %% testComputeControlSequenceInvalidMeasDelays
        function testComputeControlSequenceInvalidMeasDelays(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:CheckMeasurementsAndDelays:InvalidMeasDelay';
            
            modeObservation = 2;
            modeDelay = 1;
            measurements = ones(this.dimY, 3); % three measurements, column-wise arranged
            invalidDelays = [0 -1 1]; % negative delay            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 0.5 1]; % fractional delay            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays, modeObservation, modeDelay), ...
                expectedErrId);
            
            invalidDelays = [0 1]; % not enough entries           
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays, modeObservation, modeDelay), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceInvalidMeasDelaysNotUnique
        function testComputeControlSequenceInvalidMeasDelaysNotUnique(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsNotUnique';
            
            measurements = ones(this.dimY, 3); % three measurements, column-wise arranged
            invalidDelays = [0 1 0]; % only one measurement per delay allowed          
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays), ...
                expectedErrId);            
            
            % existence of mode observation shoud not change anything
            modeObservation = 2;
            modeDelay = 1;
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(measurements, invalidDelays, modeObservation, modeDelay), ...
                expectedErrId);  
        end
        
        %% testComputeControlSequenceInvalidModes
        function testComputeControlSequenceInvalidModes(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:CheckModeMeasurementsAndDelays:InvalidModeObservations';
                        
            % two mode observations, column-wise arranged, but of wrong dimension
            invalidModes = ones(2, 2);
            delays = [0 1];        
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([],[], invalidModes, delays), ...
                expectedErrId);
            
            % two mode observations, but out of scope
            invalidModes = [-1 1];
            delays = [0 1];
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], invalidModes, delays), ...
                expectedErrId);
            
            % two mode observations, but out of scope
            invalidModes = [1 4];
            delays = [0 1];
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], invalidModes, delays), ...
                expectedErrId);
            
            % two mode observations, but out of scope
            invalidModes = [0 1];
            delays = [0 1];
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], invalidModes, delays), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceInvalidModeDelays
        function testComputeControlSequenceInvalidModeDelays(this)
            expectedErrId = 'RecedingHorizonUdpLikeController:CheckModeMeasurementsAndDelays:InvalidModeDelays';
                        
            modes = [1 2];
            invalidDelays = 1; % too few        
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([],[], modes, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [-1 0]; % negative
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], modes, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 0.5]; % fractional
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], modes, invalidDelays), ...
                expectedErrId);
            
            invalidDelays = [0 inf]; % not finite
            this.verifyError(@() this.controllerUnderTest.computeControlSequence([], [], modes, invalidDelays), ...
                expectedErrId);
        end
        
        %% testComputeControlSequenceIssueWarningDiscardMeasurement
        function testComputeControlSequenceIssueWarningDiscardMeasurement(this)
            import matlab.unittest.constraints.IsScalar;
            
            expectedWarnId = 'RecedingHorizonUdpLikeController:GetApplicableMeasurements:IgnoringMeasurementsDelayTooLarge';
            
            zeroState = zeros(this.dimX, 1);
            this.controllerUnderTest.setControllerPlantState(Gaussian(zeroState, this.x0Cov));
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
            % the initial state is the origin, and we have a linear control
            % law, so perform a sanity check
            expectedInputSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            measurements = ones(this.dimY, 4); % four measurements, column-wise arranged
            delays = [0 2 1 3]; % two measurements are too old     
            actualInputSequence = this.verifyWarning(@() this.controllerUnderTest.computeControlSequence(measurements, delays), ...
                expectedWarnId);            
     
            this.verifyEqual(actualInputSequence, expectedInputSequence);
            % the controller state should have changed
            this.verifyNotEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(this.controllerUnderTest.lastNumIterations, 0);
            this.verifyLessThanOrEqual(this.controllerUnderTest.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(this.controllerUnderTest.lastNumIterations, IsScalar);
        end
        
        %% testComputeControlSequenceZeroStateNoMeasurementsNoModes
        function testComputeControlSequenceZeroStateNoMeasurementsNoModes(this)
            import matlab.unittest.constraints.IsScalar;

            zeroState = zeros(this.dimX, 1);
            this.controllerUnderTest.setControllerPlantState(Gaussian(zeroState, this.x0Cov));
                        
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
            this.assertTrue(this.controllerUnderTest.useMexImplementation);
            % the initial state is the origin, and we have a linear control
            % law, so perform a sanity check
            expectedInputSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence();
            
            this.verifyEqual(actualInputSequence, expectedInputSequence);
            % the controller state should not have changed, as we have no
            % measurements that could be used to correct it
            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), zeroState);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(this.controllerUnderTest.lastNumIterations, 0);
            this.verifyLessThanOrEqual(this.controllerUnderTest.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(this.controllerUnderTest.lastNumIterations, IsScalar);            
        end
        
        
        %% testComputeControlSequenceZeroStateNoMeasurementsNoModesNoMex
        function testComputeControlSequenceZeroStateNoMeasurementsNoModesNoMex(this)
            import matlab.unittest.constraints.IsScalar;

            % use the matlab implementation to obtain the controller gains
            useMexImplementation = false;
            controller = RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, ...
                this.W, this.V, this.horizonLength, this.x0, this.x0Cov, useMexImplementation);
                        
            zeroState = zeros(this.dimX, 1);
            controller.setControllerPlantState(Gaussian(zeroState, this.x0Cov));
                        
            this.assertEqual(controller.getControllerPlantState(), zeroState);
            this.assertFalse(controller.useMexImplementation);
            % the initial state is the origin, and we have a linear control
            % law, so perform a sanity check
            expectedInputSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            actualInputSequence = controller.computeControlSequence();
            
            this.verifyEqual(actualInputSequence, expectedInputSequence);
            % the controller state should not have changed, as we have no
            % measurements that could be used to correct it
            this.verifyEqual(controller.getControllerPlantState(), zeroState);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(controller.lastNumIterations, 0);
            this.verifyLessThanOrEqual(controller.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(controller.lastNumIterations, IsScalar);            
        end
        
        %% testComputeControlSequenceNoMeasurementsNoModes
        function testComputeControlSequenceNoMeasurementsNoModes(this)
            import matlab.unittest.constraints.IsScalar;
            
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
            this.assertTrue(this.controllerUnderTest.useMexImplementation);

            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsNoMeasurementsNoModes();
     
            expectedInputSequence = L(:, :, 1) * this.augX0;            
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence;
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence();
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(this.controllerUnderTest.lastNumIterations, 0);
            this.verifyLessThanOrEqual(this.controllerUnderTest.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(this.controllerUnderTest.lastNumIterations, IsScalar);
        end
        
        %% testComputeControlSequenceNoMeasurementsNoModesNoMex
        function testComputeControlSequenceNoMeasurementsNoModesNoMex(this)
            import matlab.unittest.constraints.IsScalar;
            
            % use the matlab implementation to obtain the controller gains
            useMexImplementation = false;
            controller = RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, ...
                this.W, this.V, this.horizonLength, this.x0, this.x0Cov, useMexImplementation);
            controller.maxNumIterations = RecedingHorizonUdpLikeControllerTest.numIterations;
            
            controller.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(controller.getControllerPlantState(), this.x0);
            % ensure that matlab implementation is used
            this.assertFalse(controller.useMexImplementation);
            this.assertEqual(controller.K, this.initialK);
            this.assertEqual(controller.L, this.initialL);

            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsNoMeasurementsNoModes();
     
            expectedInputSequence = L(:, :, 1) * this.augX0;            
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence;
            
            actualInputSequence = controller.computeControlSequence();
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(controller.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(controller.lastNumIterations, 0);
            this.verifyLessThanOrEqual(controller.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(controller.lastNumIterations, IsScalar);
        end
        
        %% testComputeControlSequenceNoModes
        function testComputeControlSequenceNoModes(this)
            import matlab.unittest.constraints.IsScalar;
            
            receivedMeasurement = ones(this.dimY, 1);
            measurementDelay = 0;
            
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
            this.assertTrue(this.controllerUnderTest.useMexImplementation);

            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsNoModes();

            expectedInputSequence = L(:, :, 1) * this.augX0;
            innovation = [receivedMeasurement; zeros(this.dimY, 1)] - this.S(:, :, 2) * this.augC * this.augX0;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
            
            actualInputSequence = this.controllerUnderTest.computeControlSequence(receivedMeasurement, measurementDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(this.controllerUnderTest.lastNumIterations, 0);
            this.verifyLessThanOrEqual(this.controllerUnderTest.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(this.controllerUnderTest.lastNumIterations, IsScalar);
        end
        
        %% testComputeControlSequenceNoModesNoMex
        function testComputeControlSequenceNoModesNoMex(this)
            import matlab.unittest.constraints.IsScalar;
            
            % use the matlab implementation to obtain the controller gains
            useMexImplementation = false;
            controller = RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, ...
                this.W, this.V, this.horizonLength, this.x0, this.x0Cov, useMexImplementation);
            
            receivedMeasurement = ones(this.dimY, 1);
            measurementDelay = 0;
            
            controller.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(controller.getControllerPlantState(), this.x0);
            this.assertFalse(controller.useMexImplementation);

            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsNoModes();

            expectedInputSequence = L(:, :, 1) * this.augX0;
            innovation = [receivedMeasurement; zeros(this.dimY, 1)] - this.S(:, :, 2) * this.augC * this.augX0;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
            
            actualInputSequence = controller.computeControlSequence(receivedMeasurement, measurementDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(controller.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(controller.lastNumIterations, 0);
            this.verifyLessThanOrEqual(controller.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(controller.lastNumIterations, IsScalar);
        end
        
        %% testComputeControlSequenceMeasurementsMode
        function testComputeControlSequenceMeasurementsMode(this)
            import matlab.unittest.constraints.IsScalar;
            
            receivedMeasurements = [ones(this.dimY, 1) 1.5 * ones(this.dimY, 1)];
            measurementDelays = [1 0];
            % we observed the previous mode
            previousMode = 1;
            modeDelay = 1;
            
            this.controllerUnderTest.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(this.controllerUnderTest.getControllerPlantState(), this.x0);
            this.assertTrue(this.controllerUnderTest.useMexImplementation);

            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsMeasurementAndMode();

            expectedInputSequence = L(:, :, 1) * this.augX0;
            innovation = [receivedMeasurements(:, 2); receivedMeasurements(:, 1)] - this.S(:, :, 4) * this.augC * this.augX0;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
                        
            actualInputSequence = this.controllerUnderTest.computeControlSequence(receivedMeasurements, measurementDelays, previousMode, modeDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(this.controllerUnderTest.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(this.controllerUnderTest.lastNumIterations, 0);
            this.verifyLessThanOrEqual(this.controllerUnderTest.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(this.controllerUnderTest.lastNumIterations, IsScalar);
        end
        
        %% testComputeControlSequenceMeasurementsModeNoMex
        function testComputeControlSequenceMeasurementsModeNoMex(this)
            import matlab.unittest.constraints.IsScalar;
            
            % use the matlab implementation to obtain the controller gains
            useMexImplementation = false;
            controller = RecedingHorizonUdpLikeController(this.A, this.B, this.C, this.Q, this.R, ...
                this.modeTransitionMatrix, this.scDelayProbs, this.sequenceLength, this.maxMeasDelay, ...
                this.W, this.V, this.horizonLength, this.x0, this.x0Cov, useMexImplementation);
                        
            receivedMeasurements = [ones(this.dimY, 1) 1.5 * ones(this.dimY, 1)];
            measurementDelays = [1 0];
            % we observed the previous mode
            previousMode = 1;
            modeDelay = 1;
            
            controller.setInitialControllerGains(this.initialK, this.initialL);
            % initial control state is non-zero
            this.assertEqual(controller.getControllerPlantState(), this.x0);
            % and we use Matlab to compute the gains
            this.assertFalse(controller.useMexImplementation);

            % expected controller gains                  
            [K, L, Aexp, Bexp] = this.computeGainsMeasurementAndMode();

            expectedInputSequence = L(:, :, 1) * this.augX0;
            innovation = [receivedMeasurements(:, 2); receivedMeasurements(:, 1)] - this.S(:, :, 4) * this.augC * this.augX0;
            expectedNewControllerState = Aexp(:, :, 1) * this.augX0 + Bexp(:, :, 1) * expectedInputSequence + K(:, :, 1) * innovation;
            
            actualInputSequence = controller.computeControlSequence(receivedMeasurements, measurementDelays, previousMode, modeDelay);
            this.verifyEqual(actualInputSequence, expectedInputSequence, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);

            this.verifyEqual(controller.getControllerPlantState(), expectedNewControllerState(1:this.dimX), ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
            
            % check that at least one iteration was performed to obtain the
            % gains
            this.verifyGreaterThan(controller.lastNumIterations, 0);
            this.verifyLessThanOrEqual(controller.lastNumIterations, ...
                RecedingHorizonUdpLikeControllerTest.numIterations);
            this.verifyThat(controller.lastNumIterations, IsScalar);
        end
%%
%%
        %% testDoStageCostsComputation
        function testDoStageCostsComputation(this)
            state = ones(this.dimX, 1);
            input = 1.5 * ones(this.dimU, 1);
            timestep = 1;
            
            expectedStageCosts = state' * this.Q * state + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
        end
%%
%%
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            trajectoryLength = 10;
            inputs = ones(this.dimU, trajectoryLength);
            expectedErrId = 'RecedingHorizonUdpLikeController:DoCostsComputation';
            
            invalidStates = ones(this.dimX, trajectoryLength + 2); %  trajectory too long
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStates, inputs), ...
                expectedErrId);
            
            invalidStates = ones(this.dimX, trajectoryLength); %  trajectory too short
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStates, inputs), ...
                expectedErrId);
        end
        
         %% testDoCostsComputation
        function testDoCostsComputation(this)
            trajectoryLength = 100;  
            inputs = ones(this.dimU, trajectoryLength);
            states = 2.25 * ones(this.dimX, trajectoryLength + 1);
            
            % cost function is simple a quadratic
            expectedCosts = states(:, end)' * this.Q * states(:, end);
            for j=1:trajectoryLength
                expectedCosts = expectedCosts + states(:, j)' * this.Q * states(:, j) ...
                    + inputs(:, j)' * this.R * inputs(:, j);
            end
                        
            actualCosts = this.controllerUnderTest.computeCosts(states, inputs);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', RecedingHorizonUdpLikeControllerTest.absTol);
        end
    end
end

