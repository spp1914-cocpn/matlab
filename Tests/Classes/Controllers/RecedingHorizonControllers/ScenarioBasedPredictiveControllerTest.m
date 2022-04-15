classdef ScenarioBasedPredictiveControllerTest < matlab.unittest.TestCase
    % Test cases for ScenarioBasedPredictiveController.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2022  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        absTol = 1e-7;
    end
    
    properties
        A;
        B;
        Q;
        R;
        
        Z;
        refTrajectory;
        Qref;
        stateTrajectoryRef;
        inputTrajectoryRef;
        
        sequenceLength;
        dimX;
        dimU;    
        
        initialPlantState;
        initialBufferedSequences;
        
        stateTrajectory;
        inputTrajectory;
        horizonLength;       
        
        caDelayProbs;
        modeTransitionMatrix;
        scenarioProbabilities;
        
        controllerUnderTest;
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
                
            this.sequenceLength = 3;
            
            this.caDelayProbs = [1/3 5/12 1/6 1/12 0 0];          
            this.modeTransitionMatrix = Utility.calculateDelayTransitionMatrix( ...
                Utility.truncateDiscreteProbabilityDistribution(this.caDelayProbs, this.sequenceLength + 1));
            
            this.dimX = 2;
            this.dimU = 2;
            this.horizonLength = 5;
            this.initialPlantState = [0.1; 0];
            
            V = diag([0.01, 1]); % in the book: z = V*x
            this.A = diag([0.9512, 0.9048]);
           
            this.B = [4.877 4.877; -1.1895 3.569];
            this.Q = V * diag([50 0.02]) * V; % R_3 in the book
            this.R = diag([1/3, 3]); % R_2 in the book   
            
            % assume that the buffer is not empty
            this.initialBufferedSequences(:, :, 3) = [13 14 15 
                                                      16 17 18];            
            this.initialBufferedSequences(:, :, 2) = [7 8 9 
                                                      10 11 12];
            this.initialBufferedSequences(:, :, 1) = [1 2 5 
                                                      3 4 6];                                                  
            
            this.computeScenarioProbabilities([0 0 0 1]); 
                                           
            % only the first state variable is to be driven to the origin
            this.Z = [1 0];
            this.Qref = this.Q(1,1);
            this.refTrajectory = zeros(1, 2 * this.horizonLength);
        end
    end
    
    methods (Access = private)
        %% setupControllerUnderTest
        function setupControllerUnderTest(this, doTrackRef)
            lazyInitOptimProb = true;
            if doTrackRef
                this.controllerUnderTest = ScenarioBasedPredictiveController(this.A, this.B, this.Qref, this.R, this.sequenceLength, ...
                    this.modeTransitionMatrix, lazyInitOptimProb, this.Z, this.refTrajectory);
            else
                this.controllerUnderTest = ScenarioBasedPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength, ...
                    this.modeTransitionMatrix, lazyInitOptimProb);
            end
            
            this.controllerUnderTest.setBufferedControlSequences(this.initialBufferedSequences);            
            this.controllerUnderTest.changeHorizonLength(this.horizonLength);
        end
        
        %% computeScenarioProbabilities
        function computeScenarioProbabilities(this, inputProbs)            
            alphas = cell(1, this.sequenceLength);            
           
            alphas{1} = inputProbs(:);           
            
            for k=1:this.sequenceLength - 1
                alphas{k+1} = (this.modeTransitionMatrix')^k*alphas{1};
                alphas{k+1} = alphas{k+1}(k+1:end) / sum(alphas{k+1}(k+1:end));
            end
            
            % compute the scenario probs
            % hardcoded
            numScenarios = 24;
            this.scenarioProbabilities = ones(1, numScenarios);
                        
            this.scenarioProbabilities(1) = alphas{1}(1) * alphas{2}(1) * alphas{3}(1); % 1 1 1
            this.scenarioProbabilities(2) = alphas{1}(2) * alphas{2}(1) * alphas{3}(1); % 2 1 1 
            this.scenarioProbabilities(3) = alphas{1}(3) * alphas{2}(1) * alphas{3}(1); % 3 1 1
            this.scenarioProbabilities(4) = alphas{1}(4) * alphas{2}(1) * alphas{3}(1); % 4 1 1
            this.scenarioProbabilities(5) = alphas{1}(1) * alphas{2}(2) * alphas{3}(1); % 1 2 1
            this.scenarioProbabilities(6) = alphas{1}(2) * alphas{2}(2) * alphas{3}(1); % 2 2 1
            this.scenarioProbabilities(7) = alphas{1}(3) * alphas{2}(2) * alphas{3}(1); % 3 2 1
            this.scenarioProbabilities(8) = alphas{1}(4) * alphas{2}(2) * alphas{3}(1); % 4 2 1
            this.scenarioProbabilities(9) = alphas{1}(1) * alphas{2}(3) * alphas{3}(1); % 1 3 1
            this.scenarioProbabilities(10) = alphas{1}(2) * alphas{2}(3) * alphas{3}(1); % 2 3 1
            this.scenarioProbabilities(11) = alphas{1}(3) * alphas{2}(3) * alphas{3}(1); % 3 3 1
            this.scenarioProbabilities(12) = alphas{1}(4) * alphas{2}(3) * alphas{3}(1); % 4 3 1
            this.scenarioProbabilities(13) = alphas{1}(1) * alphas{2}(1) * alphas{3}(2); % 1 1 2
            this.scenarioProbabilities(14) = alphas{1}(2) * alphas{2}(1) * alphas{3}(2); % 2 1 2
            this.scenarioProbabilities(15) = alphas{1}(3) * alphas{2}(1) * alphas{3}(2); % 3 1 2
            this.scenarioProbabilities(16) = alphas{1}(4) * alphas{2}(1) * alphas{3}(2); % 4 1 2
            this.scenarioProbabilities(17) = alphas{1}(1) * alphas{2}(2) * alphas{3}(2); % 1 2 2
            this.scenarioProbabilities(18) = alphas{1}(2) * alphas{2}(2) * alphas{3}(2); % 2 2 2
            this.scenarioProbabilities(19) = alphas{1}(3) * alphas{2}(2) * alphas{3}(2); % 3 2 2
            this.scenarioProbabilities(20) = alphas{1}(4) * alphas{2}(2) * alphas{3}(2); % 4 2 2
            this.scenarioProbabilities(21) = alphas{1}(1) * alphas{2}(3) * alphas{3}(2); % 1 3 2
            this.scenarioProbabilities(22) = alphas{1}(2) * alphas{2}(3) * alphas{3}(2); % 2 3 2
            this.scenarioProbabilities(23) = alphas{1}(3) * alphas{2}(3) * alphas{3}(2); % 3 3 2
            this.scenarioProbabilities(24) = alphas{1}(4) * alphas{2}(3) * alphas{3}(2); % 4 3 2
        end
        
        %% setupScenarioTree
        function scenarioProb = setupScenarioTree(this, doTrackRef)
            if doTrackRef
                [~,P,~] = lqry(ss(this.A, this.B, this.Z, [], -1), this.Qref, this.R);            
                P = P(1,1); % terminal weighting
            else                
                P = idare(this.A, this.B, this.Q, this.R);
            end
            
            defaultInput = zeros(this.dimX, 1);
            % 24 scenarios           
            controls = optimvar('controls', [this.dimU, this.horizonLength]);
            scenarioStates = optimvar('scenarioStates', [this.dimX, this.horizonLength + 1, 24]);
            scenarioInputs = optimvar('scenarioInputs', [this.dimU, this.horizonLength, 24]);
                        
            inputConstraints = optimconstr(this.dimU, this.horizonLength * 24);
            stateConstraints = optimconstr(this.dimX, (this.horizonLength + 1) * 24);
            costs = optimexpr(24 * (this.horizonLength + 1));
            
            % counter for constraints
            inputIdx = 1;
            costIdx = 1;
            stateIdx = 1;
            for s=1:24
                if doTrackRef
                    % terminal costs
                    z_term = this.Z * scenarioStates(:, this.horizonLength + 1, s) - this.refTrajectory(:, this.horizonLength + 1);
                    costs(costIdx) = this.scenarioProbabilities(s) * z_term' * P * z_term;
                else
                    % terminal costs
                    costs(costIdx) = scenarioStates(:, end, s)' * P * scenarioStates(:, end, s) * this.scenarioProbabilities(s);
                end
                % initial state constraint
                stateConstraints(:, stateIdx) = scenarioStates(:, 1, s) == this.initialPlantState;
                stateIdx = stateIdx + 1;
                costIdx = costIdx + 1;
                % first stage
                switch s
                    case {1, 5, 9, 13, 17, 21}                        
                        inputConstraints(:, inputIdx) = scenarioInputs(:, 1, s) == controls(:, 1);
                    case {2, 6, 10, 14, 18, 22}                        
                        inputConstraints(:, inputIdx) = scenarioInputs(:, 1, s) == this.initialBufferedSequences(:, 2, 1);
                    case {3, 7, 11, 15, 19, 23}                        
                        inputConstraints(:, inputIdx) = scenarioInputs(:, 1, s) == this.initialBufferedSequences(:, 3, 2);
                    otherwise                        
                        inputConstraints(:, inputIdx) = scenarioInputs(:, 1, s) == defaultInput;
                end
                if doTrackRef
                    z_1 = this.Z * scenarioStates(:, 1, s) - this.refTrajectory(:, 1);
                    costs(costIdx) = this.scenarioProbabilities(s) * (z_1' * this.Qref * z_1 + scenarioInputs(:, 1, s)' * this.R * scenarioInputs(:, 1, s));
                else
                    costs(costIdx) = this.scenarioProbabilities(s) * (scenarioStates(:, 1, s)' * this.Q * scenarioStates(:, 1, s) ...
                         + scenarioInputs(:, 1, s)' * this.R * scenarioInputs(:, 1, s));                
                end
                stateConstraints(:, stateIdx) = scenarioStates(:, 2, s) == this.A * scenarioStates(:, 1, s) + this.B *  scenarioInputs(:, 1, s);
                
                inputIdx = inputIdx + 1;
                stateIdx = stateIdx + 1;
                costIdx = costIdx + 1;
                
                % second stage
                switch s
                    case {1, 2, 3, 4, 13, 14, 15, 16}                        
                        inputConstraints(:, inputIdx) = scenarioInputs(:, 2, s) == controls(:, 2);
                    case {5, 6, 7, 8, 17, 18, 19, 20}                        
                        inputConstraints(:, inputIdx) = scenarioInputs(:, 2, s) == this.initialBufferedSequences(:, 3, 1);
                    otherwise                        
                        inputConstraints(:, inputIdx) = scenarioInputs(:, 2, s) == defaultInput;
                end
                if doTrackRef
                    z_2 = this.Z * scenarioStates(:, 2, s) - this.refTrajectory(:, 2);
                    costs(costIdx) = this.scenarioProbabilities(s) * (z_2' * this.Qref * z_2 + scenarioInputs(:, 2, s)' * this.R * scenarioInputs(:, 2, s));
                else
                    costs(costIdx) = this.scenarioProbabilities(s) * (scenarioStates(:, 2, s)' * this.Q * scenarioStates(:, 2, s) ...
                         + scenarioInputs(:, 2, s)' * this.R * scenarioInputs(:, 2, s));
                end
                stateConstraints(:, stateIdx) = scenarioStates(:, 3, s) == this.A * scenarioStates(:, 2, s) + this.B *  scenarioInputs(:, 2, s);
                
                inputIdx = inputIdx + 1;
                stateIdx = stateIdx + 1;
                costIdx = costIdx + 1;
                
                % third stage
                if s <= 12
                    inputConstraints(:, inputIdx) = scenarioInputs(:, 3, s) == controls(:, 3);
                else
                    inputConstraints(:, inputIdx) = scenarioInputs(:, 3, s) == defaultInput;
                end
                if doTrackRef
                    z_3 = this.Z * scenarioStates(:, 3, s) - this.refTrajectory(:, 3);
                	costs(costIdx) = this.scenarioProbabilities(s) * (z_3' * this.Qref * z_3 + scenarioInputs(:, 3, s)' * this.R * scenarioInputs(:, 3, s));
                else
                    costs(costIdx) = this.scenarioProbabilities(s) * (scenarioStates(:, 3, s)' * this.Q * scenarioStates(:, 3, s) ...
                         + scenarioInputs(:, 3, s)' * this.R * scenarioInputs(:, 3, s));
                end
                stateConstraints(:, stateIdx) = scenarioStates(:, 4, s) == this.A * scenarioStates(:, 3, s) + this.B *  scenarioInputs(:, 3, s);
                
                inputIdx = inputIdx + 1;
                stateIdx = stateIdx + 1;
                costIdx = costIdx + 1;
                 
                % remaining stages
                for h=4:this.horizonLength
                    inputConstraints(:, inputIdx) = scenarioInputs(:, h, s) == controls(:, h);
                    if doTrackRef
                        z_h = this.Z * scenarioStates(:, h, s) - this.refTrajectory(:, h);
                        costs(costIdx) = this.scenarioProbabilities(s) * (z_h' * this.Qref * z_h + scenarioInputs(:, h, s)' * this.R * scenarioInputs(:, h, s));
                    else
                        costs(costIdx) = this.scenarioProbabilities(s) * (scenarioStates(:, h, s)' * this.Q * scenarioStates(:, h, s) ...
                            + scenarioInputs(:, h, s)' * this.R * scenarioInputs(:, h, s));
                    end
                    stateConstraints(:, stateIdx) = scenarioStates(:, h + 1, s) == this.A * scenarioStates(:, h, s) + this.B *  scenarioInputs(:, h, s);
                    
                    inputIdx = inputIdx + 1;
                    stateIdx = stateIdx + 1;
                    costIdx = costIdx + 1;
                end
            end                   
           
            % define the problem, and solve via quadprog
            scenarioProb = optimproblem;
            scenarioProb.Objective = sum(costs);
            scenarioProb.Constraints.stateConstraints = stateConstraints;
            scenarioProb.Constraints.inputConstraints = inputConstraints;            
        end
        
        %% computeTrajectories
        function [stateTrajectory, inputTrajectory] = computeTrajectories(this)      
            
            % construct the scenario tree, and solve via quadprog
            scenarioProb = this.setupScenarioTree(false);
            
            qp = prob2struct(scenarioProb);           
            [x,~,exitflag,~] = quadprog(qp.H, qp.f, [], [], qp.Aeq, qp.beq);
            assert(exitflag == 1); % problem feasible, solution found
            
            % extract desired variables
            controlsIdx = varindex(scenarioProb, 'controls');
            scenarioStatesIdx = varindex(scenarioProb, 'scenarioStates');
            % we need the states resulting from first scenario
            endIdx = this.dimX * (this.horizonLength +1);
            
            inputTrajectory = reshape(x(controlsIdx), this.dimU, this.horizonLength);
            stateTrajectory = reshape(x(scenarioStatesIdx(1:endIdx)) , this.dimX, this.horizonLength + 1);             
        end
        
        %% computeTrajectoriesRef
        function [stateTrajectoryRef, inputTrajectoryRef] = computeTrajectoriesRef(this)
            
            % construct the scenario tree, and solve via quadprog
            scenarioProb = this.setupScenarioTree(true);
                       
            qp = prob2struct(scenarioProb);            
            [x,~,exitflag,~] = quadprog(qp.H, qp.f, [], [], qp.Aeq, qp.beq);
            assert(exitflag == 1); % problem feasible, solution found                      
            
            % extract desired variables
            controlsIdx = varindex(scenarioProb, 'controls');
            scenarioStatesIdx = varindex(scenarioProb, 'scenarioStates');
            % we need the states resulting from first scenario
            endIdx = this.dimX * (this.horizonLength + 1);
            
            inputTrajectoryRef = reshape(x(controlsIdx), this.dimU, this.horizonLength);
            stateTrajectoryRef = reshape(x(scenarioStatesIdx(1:endIdx)) , this.dimX, this.horizonLength + 1);         
        end
        
        %% computeCosts
        function costs = computeCosts(this)
            % compute the costs in a straightforward manner            
            costs =  this.stateTrajectory(:, end)' * this.Q * this.stateTrajectory(:, end);
            for j = 1:this.horizonLength
                costs = costs + this.stateTrajectory(:, j)' * this.Q * this.stateTrajectory(:, j) ...
                    + this.inputTrajectory(:, j)' * this.R * this.inputTrajectory(:, j);
            end
        end
        
        %% computeCostsRef
        function costsRef = computeCostsRef(this)
            % compute the costs in a straightforward manner
            performance = this.Z * this.stateTrajectoryRef(:, 1:end) - this.refTrajectory(:, 1:this.horizonLength + 1);
            costsRef =  performance(:, end)' * this.Qref * performance(:, end);
            for j = 1:this.horizonLength
                costsRef = costsRef + performance(:, j)' * this.Qref * performance(:, j) ...
                    + this.inputTrajectoryRef(:, j)' * this.R * this.inputTrajectoryRef(:, j);
            end
        end
    end
    
    methods (Test)
        %% testScenarioBasedPredictiveControllerInvalidSysMatrix
        function testScenarioBasedPredictiveControllerInvalidSysMatrix(this)            
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
                        
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() ScenarioBasedPredictiveController(invalidA, this.B, this.Q, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            invalidA = inf(this.dimX, this.dimX); % square but inf
             this.verifyError(@() ScenarioBasedPredictiveController(invalidA, this.B, this.Q, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
        end
        
        %% testScenarioBasedPredictiveControllerInvalidInputMatrix
        function testScenarioBasedPredictiveControllerInvalidInputMatrix(this)            
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidB = eye(this.dimX + 1, this.dimU); % invalid dims
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, invalidB, this.Q, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            invalidB = eye(this.dimX, this.dimU); % correct dims, but nan
            invalidB(1, 1) = nan;
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, invalidB, this.Q, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
        end
        
        %% testScenarioBasedPredictiveControllerInvalidCostMatrices
        function testScenarioBasedPredictiveControllerInvalidCostMatrices(this)            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
             this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
             this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = -eye(this.dimX); % Q is not psd
             this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, invalidQ, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, invalidR, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, invalidR, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, invalidR, ...
                this.sequenceLength, this.modeTransitionMatrix, true), ...
                expectedErrId);
        end
        
        %% testScenarioBasedPredictiveControllerInvalidModeTransitionMatrix
        function testScenarioBasedPredictiveControllerInvalidModeTransitionMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';
            
            invalidModeTransitionMatrix = this.modeTransitionMatrix(2:end, 2:end);% invalid dimensions
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, invalidModeTransitionMatrix, true), expectedErrId);
            
            invalidModeTransitionMatrix = [0 0.1 0.8 0.2];% not a matrix
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, invalidModeTransitionMatrix, true), expectedErrId);
                     
            invalidModeTransitionMatrix = this.modeTransitionMatrix;
            invalidModeTransitionMatrix(1,1) = 1.1; % does not sum up to 1
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, invalidModeTransitionMatrix, true), expectedErrId);
            
            invalidModeTransitionMatrix = this.modeTransitionMatrix;
            invalidModeTransitionMatrix(1,1) = -invalidModeTransitionMatrix(1,1); % negative entry
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, invalidModeTransitionMatrix, true), expectedErrId);
        end
        
         %% testScenarioBasedPredictiveControllerInvalidNumArgs
        function testScenarioBasedPredictiveControllerInvalidNumArgs(this)            
            expectedErrId = 'ScenarioBasedPredictiveController:InvalidNumberOfArguments';
            
            % use 8 arguments, we expect either 7 or 9
            this.verifyError(@() ScenarioBasedPredictiveController(this.A, this.B, this.Q, this.R, ...
                this.sequenceLength, this.modeTransitionMatrix, true, []), ...
                expectedErrId);
        end
%%
%%        
        %% testScenarioBasedPredictiveController
        function testScenarioBasedPredictiveController(this)            
            % should not crash
            controller = ScenarioBasedPredictiveController(this.A, this.B, this.Q, this.R, this.sequenceLength, this.modeTransitionMatrix, false);
              
            this.verifyTrue(controller.requiresExternalStateEstimate); % needs a filter or state feedback
            this.verifyEqual(controller.horizonLength, this.sequenceLength); % horizon is by default equal to sequence length
            
            this.verifyEqual(controller.scenarioProbabilities, this.scenarioProbabilities);
            this.verifyNotEmpty(controller.getScenarioTree());
        end
%%
%%
         %% testChangeModelParametersInvalidSystemMatrix
        function testChangeModelParametersInvalidSystemMatrix(this)
            this.setupControllerUnderTest(false);
            
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
            
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidDimensions';
            
            invalidA = this.A(1,1); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
                        
            invalidA = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
            
            invalidA = this.A(:, 1); % not square
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
            
            invalidA = this.A; % not finite
            invalidA(end) = inf;
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
        end
        
         %% testChangeModelParametersInvalidInputMatrix
        function testChangeModelParametersInvalidInputMatrix(this)
            this.setupControllerUnderTest(false);
            
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
            
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrixDims';
            
            invalidB = this.B(1,1); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
                        
            invalidB = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
                        
            invalidB = this.B; % not finite
            invalidB(end) = inf;
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
        end
        
         %% testChangeModelParameters
        function testChangeModelParameters(this)
            this.setupControllerUnderTest(false);
            % use a Gaussian mixture state
            state = GaussianMixture(repmat(this.initialPlantState, 1, 4), ...
                repmat(eye(this.dimX), 1, 1, 4), [0 0 0 1]); 
            timestep = 1;
            mode = 1; % not needed but has to be passed
                       
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
                                    
            newA = this.A - 0.1 * eye(this.dimX);
            newB = this.B + 2 * eye(this.dimU);
                        
            this.controllerUnderTest.changeModelParameters(newA, newB);
            % check if computation of new sequence is affected
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            
            this.A = newA;
            this.B = newB;
            [~, newInputs] = this.computeTrajectories();            
            
            expectedInputs = newInputs(:, 1:this.sequenceLength);                                 
            this.verifyEqual(actualSequence, expectedInputs(:), 'AbsTol', 1e-6);
        end
        %% testChangeModelParametersNewWIgnored
        function testChangeModelParametersNewWIgnored(this)
            this.setupControllerUnderTest(false);
            
            this.computeScenarioProbabilities([0 0 0 1] * this.modeTransitionMatrix);
            [~, this.inputTrajectory] = this.computeTrajectories();    
            % we pass Gaussian state, so controller internally propagates input probs
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
                       
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
            
            % pass a new noise covariance, has no effect
            newW = eye(this.dimX);
                        
            this.controllerUnderTest.changeModelParameters(this.A, this.B, newW);
            % should not change the input sequence
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);            
                        
            expectedInputs = this.inputTrajectory(:, 1:this.sequenceLength);                                 
            this.verifyEqual(actualSequence, expectedInputs(:), 'AbsTol', 1e-6);
        end
%%
%%
        %% testChangeCaDelayProbs
        function testChangeCaDelayProbs(this)            
            this.setupControllerUnderTest(false);
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            this.assertTrue(isa(this.controllerUnderTest, 'CaDelayProbsChangeable'));
                        
            newDelayProbs = [0 0 1 0 0]; % fixed delay of 2 time steps
            this.caDelayProbs = newDelayProbs;
            this.modeTransitionMatrix = Utility.calculateDelayTransitionMatrix([0 0 1 0]);
            this.computeScenarioProbabilities([0 0 1 0]);
            
            this.controllerUnderTest.changeCaDelayProbs(newDelayProbs);                   
                        
            % check if computation of new sequence is affected
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            
            % check side effect
            this.verifyEqual(this.controllerUnderTest.scenarioProbabilities, this.scenarioProbabilities, 'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
            
            [~, newInputs] = this.computeTrajectories();            
            % only the last element of the computed control sequence is
            % relevant, due to fixed delay of two time steps, the other two
            % will never be applied                        
            this.verifyEqual(actualSequence(end-this.dimX+1:end), newInputs(:, this.sequenceLength), 'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
%%
%%
        %% testDoControlSequenceComputationInvalidTimestep
        function testDoControlSequenceComputationInvalidTimestep(this)
            this.setupControllerUnderTest(false);
            
            expectedErrId = 'ScenarioBasedPredictiveController:DoControlSequenceComputation:InvalidTimestep';
            
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            mode = 1; % not needed but has to be passed
            
            invalidTimestep = -1; % not positive
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(zeroState, mode, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = 1.5; % not an integer            
            this.verifyError(@() this.controllerUnderTest.computeControlSequence(zeroState, mode, invalidTimestep), ...
                expectedErrId);
        end 
        
        %% testDoControlSequenceComputationZeroState
        function testDoControlSequenceComputationZeroState(this)
            this.setupControllerUnderTest(false);
            % clear the buffer
            this.controllerUnderTest.reset();
            
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            % because we emptied the actuator buffer
            % no matter what the scenario probs are
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            expectedSequence = zeros(this.dimU * this.sequenceLength, 1);
            
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState, mode, timestep);
            
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);      
        end
        
        %% testDoControlSequenceComputationZeroStateRef
        function testDoControlSequenceComputationZeroStateRef(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            % clear the buffer
            this.controllerUnderTest.reset();
            timestep = 1;
            
            % assert that initial difference to ref is zero
            state = [this.refTrajectory(:, timestep); 42]; % second component of state arbitrary
            this.assertTrue(this.Z * state == this.refTrajectory(:, timestep));
                        
            mode = 1; % not needed but has to be passed
            expectedSequence = zeros(this.dimU * this.sequenceLength, 1);

            actualSequence = this.controllerUnderTest.computeControlSequence(Gaussian(state, eye(this.dimX)), mode, timestep);
            
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationGaussianMixtureState
        function testDoControlSequenceComputationGaussianMixtureState(this)
            this.setupControllerUnderTest(false);
            
            inputProbs =  [0.3 0.2 0.4 0.1];
            this.computeScenarioProbabilities(inputProbs);
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            expectedInputs = this.inputTrajectory(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);
            
            % use a Gaussian mixture state
            state = GaussianMixture(repmat(this.initialPlantState, 1, 4), ...
                repmat(eye(this.dimX), 1, 1, 4), inputProbs); 
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationGaussianState
        function testDoControlSequenceComputationGaussianState(this)
            this.setupControllerUnderTest(false);
            
            inputProbs = [0 0 0 1] * this.modeTransitionMatrix;
            this.computeScenarioProbabilities(inputProbs);
            [~, this.inputTrajectory] = this.computeTrajectories();         
            
            expectedInputs = this.inputTrajectory(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);
            
            % we pass Gaussian state, so controller internally propagates input probs
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationLongerHorizonGaussianMixtureState
        function testDoControlSequenceComputationLongerHorizonGaussianMixtureState(this)
            this.setupControllerUnderTest(false); 
            
            inputProbs =  [0.3 0.2 0.4 0.1];
            this.computeScenarioProbabilities(inputProbs);           
            
            newHorizonLength = 8;
            this.controllerUnderTest.changeHorizonLength(newHorizonLength);
            this.assertEqual(this.controllerUnderTest.horizonLength, newHorizonLength);
            
            this.horizonLength = newHorizonLength;
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            expectedInputs = this.inputTrajectory(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);            
             
            % use a Gaussian mixture state
            state = GaussianMixture(repmat(this.initialPlantState, 1, 4), ...
                repmat(eye(this.dimX), 1, 1, 4), inputProbs); 
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationLongerHorizonGaussianState
        function testDoControlSequenceComputationLongerHorizonGaussianState(this)
            this.setupControllerUnderTest(false); 
            
          	inputProbs = [0 0 0 1] * this.modeTransitionMatrix;
            this.computeScenarioProbabilities(inputProbs);           
            
            newHorizonLength = 8;
            this.controllerUnderTest.changeHorizonLength(newHorizonLength);
            this.assertEqual(this.controllerUnderTest.horizonLength, newHorizonLength);
            
            this.horizonLength = newHorizonLength;
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            expectedInputs = this.inputTrajectory(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);            
             
            % we pass Gaussian state, so controller internally propagates input probs
            state = Gaussian(this.initialPlantState, eye(this.dimX)); 
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        
        %% testDoControlSequenceComputationRefGaussianMixtureState
        function testDoControlSequenceComputationRefGaussianMixtureState(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            
            inputProbs =  [0.3 0.2 0.4 0.1];
            this.computeScenarioProbabilities(inputProbs);
            
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
                        
            expectedInputs = this.inputTrajectoryRef(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);
            
            % use a Gaussian mixture state
            state = GaussianMixture(repmat(this.initialPlantState, 1, 4), ...
                repmat(eye(this.dimX), 1, 1, 4), inputProbs);
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationRefGaussianState
        function testDoControlSequenceComputationRefGaussianState(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            
            inputProbs = [0 0 0 1] * this.modeTransitionMatrix;
            this.computeScenarioProbabilities(inputProbs);
            
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
                        
            expectedInputs = this.inputTrajectoryRef(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);
            
            % we pass Gaussian state, so controller internally propagates input probs
            state = Gaussian(this.initialPlantState, eye(this.dimX)); 
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationRefLongerHorizonGaussianMixtureState
        function testDoControlSequenceComputationRefLongerHorizonGaussianMixtureState(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            
            inputProbs =  [0.3 0.2 0.4 0.1];
            this.computeScenarioProbabilities(inputProbs);
            
            newHorizonLength = 8;
            this.controllerUnderTest.changeHorizonLength(newHorizonLength);
            this.assertEqual(this.controllerUnderTest.horizonLength, newHorizonLength);
            
            this.horizonLength = newHorizonLength;            
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
                        
            expectedInputs = this.inputTrajectoryRef(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);
            
            % use a Gaussian mixture state
            state = GaussianMixture(repmat(this.initialPlantState, 1, 4), ...
                repmat(eye(this.dimX), 1, 1, 4), inputProbs);
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoControlSequenceComputationRefLongerHorizonGaussianState
        function testDoControlSequenceComputationRefLongerHorizonGaussianState(this)
            % include reference tracking
            this.setupControllerUnderTest(true);
            
            inputProbs = [0 0 0 1] * this.modeTransitionMatrix;
            this.computeScenarioProbabilities(inputProbs);
            
            newHorizonLength = 8;
            this.controllerUnderTest.changeHorizonLength(newHorizonLength);
            this.assertEqual(this.controllerUnderTest.horizonLength, newHorizonLength);
            
            this.horizonLength = newHorizonLength;            
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
                        
            expectedInputs = this.inputTrajectoryRef(:, 1:this.sequenceLength);
            expectedSequence = expectedInputs(:);
            
            % we pass Gaussian state, so controller internally propagates input probs
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            timestep = 1;
            mode = 1; % not needed but has to be passed
            
            actualSequence = this.controllerUnderTest.computeControlSequence(state, mode, timestep);
            this.verifyEqual(actualSequence, expectedSequence, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
%%
%%
        %% testDoStageCostsComputation
        function testDoStageCostsComputation(this)
            this.setupControllerUnderTest(false);
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            state = this.stateTrajectory(:, end-1);
            input = this.inputTrajectory(:, end);
            timestep = size(this.inputTrajectory, 2);
            
            expectedStageCosts = state' * this.Q * state + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoStageCostsComputationRef
        function testDoStageCostsComputationRef(this)
            % include ref tracking
            this.setupControllerUnderTest(true);
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
            
            state = this.stateTrajectoryRef(:, end-1);
            input = this.inputTrajectoryRef(:, end);
            timestep = size(this.inputTrajectoryRef, 2);
            
            expectedStageCosts = (this.Z * state - this.refTrajectory(:, timestep))' * this.Qref * (this.Z * state - this.refTrajectory(:, timestep)) ...
                + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqual(actualStageCosts, expectedStageCosts, ...
                'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
%%
%%
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            % include ref tracking
            this.setupControllerUnderTest(true);
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
            
            expectedErrId = 'ScenarioBasedPredictiveController:DoCostsComputation:InvalidStateTrajectory';
            
            % state trajectory too short
            appliedInputs = this.inputTrajectoryRef;
            invalidStateTrajectory = this.stateTrajectoryRef(:, 1:end-2);
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
            
            % state trajectory too long
            appliedInputs = this.inputTrajectoryRef(:, 1:end-2);
            invalidStateTrajectory = this.stateTrajectoryRef;
            this.verifyError(@() this.controllerUnderTest.computeCosts(invalidStateTrajectory, appliedInputs), ...
                expectedErrId);
        end
        
        %% testDoCostsComputation
        function testDoCostsComputation(this)
            this.setupControllerUnderTest(false);
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            expectedCosts = this.computeCosts();
            actualCosts = this.controllerUnderTest.computeCosts(this.stateTrajectory, this.inputTrajectory);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
        
        %% testDoCostsComputationRef
        function testDoCostsComputationRef(this)
            this.setupControllerUnderTest(true);
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
            
            expectedCosts = this.computeCostsRef();
            actualCosts = this.controllerUnderTest.computeCosts(this.stateTrajectoryRef, this.inputTrajectoryRef);
            
            this.verifyEqual(actualCosts, expectedCosts, 'AbsTol', ScenarioBasedPredictiveControllerTest.absTol);
        end
%%
%%
        %% testDoGetDeviationFromRefForState
        function testDoGetDeviationFromRefForState(this)
            this.setupControllerUnderTest(false);
            [this.stateTrajectory, this.inputTrajectory] = this.computeTrajectories();
            
            state = this.initialPlantState;
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(state, 1);
            
            % we want to reach the origin, so there should be a deviation
            this.verifyEqual(actualDeviation, state);
        end
        
         %% testDoGetDeviationFromRefForStateRef
        function testDoGetDeviationFromRefForStateRef(this)
            this.setupControllerUnderTest(true);
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
                       
            state = this.initialPlantState;
            actualDeviation = this.controllerUnderTest.getDeviationFromRefForState(state, 1);
            % we want to reach the origin with the first variable, so there should be a deviation
            this.verifyEqual(actualDeviation, state(1));
        end
        
        %% testDoGetDeviationFromRefForStateInvalidTimestep
        function testDoGetDeviationFromRefForStateInvalidTimestep(this)
            this.setupControllerUnderTest(true);
            [this.stateTrajectoryRef, this.inputTrajectoryRef] = this.computeTrajectoriesRef();
            
            expectedErrId = 'ScenarioBasedPredictiveController:GetDeviationFromRefForState:InvalidTimestep';
            
            state = this.initialPlantState;
            invalidTimestep = -1; % not positive
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(state, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = 1.5; % not an integer
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(state, invalidTimestep), ...
                expectedErrId);
            
            invalidTimestep = size(this.refTrajectory, 2) + 1; % exceed the length of the reference trajectory
            this.verifyError(@() this.controllerUnderTest.getDeviationFromRefForState(state, invalidTimestep), ...
                expectedErrId);
        end
    end
end

