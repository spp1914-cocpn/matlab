classdef JumpLinearSystemModelTest < matlab.unittest.TestCase
    % Test cases for JumpLinearSystemModel.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Constant)
        seed = 1;
    end
    
    properties
        transitionMatrix;
        A1;
        A2;
        B1;
        B2;
        
        % noise covs
        W1; 
        W2;
        numModes;
        modeModels;
        
        negativeNumModes;
        fractionaNumModes;
        
        modelUnderTest;
        state;
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.A1 = [1 1; 0 1];
            this.A2 = [1.2 1.2; 0 1];
            this.B1 = [0; 1];
            this.B2 = [0.1; 0.8];

            this.W1 = eye(2);
            this.W2 = eye(2);
            this.numModes = 2;
            
            this.transitionMatrix = [0.9 0.1; 0.1 0.9];
            
            mode1Model = LinearPlant(this.A1, this.B1, this.W1);
            mode2Model = LinearPlant(this.A2, this.B2, this.W2);
            
            this.modeModels = {mode1Model; mode2Model};
            
            this.modelUnderTest = JumpLinearSystemModel(this.numModes, this.modeModels);
            
            this.negativeNumModes = -2;
            this.fractionaNumModes = 1.5;
            
            this.state = [1; 1];
        end
    end
    
    methods (Test)
        %% testJumpLinearSystemModel
        function testJumpLinearSystemModelInvalidNumModes(this)
            expectedErrId = 'JumpLinearSystemModel:InvalidNumModes';
            
            this.verifyError(@() JumpLinearSystemModel(this.negativeNumModes, this.modeModels), expectedErrId);
            this.verifyError(@() JumpLinearSystemModel(this.fractionaNumModes, this.modeModels), expectedErrId);
        end
        
        %% testJumpLinearSystemModel
        function testJumpLinearSystemModel(this)
            model = JumpLinearSystemModel(this.numModes, this.modeModels);
            
            this.verifyClass(model.modeSystemModels, ?cell);
            this.verifyNumElements(model.modeSystemModels, this.numModes);
            this.verifyEqual(model.modeSystemModels{1}, this.modeModels{1});
            this.verifyEqual(model.modeSystemModels{2}, this.modeModels{2});
        end
        
        %% testSetModeSystemModelsInvalidModels
        function testSetModeSystemModelsInvalidModels(this)
            expectedErrId = 'JumpLinearSystemModel:InvalidModeSystemModels';
            
            % first, no cell array
            invalidModeModels = eye(2);
            this.verifyError(@() this.modelUnderTest.setModeSystemModels(invalidModeModels), expectedErrId);
            % now , cell array with only 1 element
            invalidModeModels = this.modeModels(1);
            this.verifyError(@() this.modelUnderTest.setModeSystemModels(invalidModeModels), expectedErrId);
            % finally, a cell array with invalid content
            invalidModeModels = {this.modeModels{1}, 42};
            this.verifyError(@() this.modelUnderTest.setModeSystemModels(invalidModeModels), expectedErrId);
        end
        
         %% testSetModeSystemModels
        function testSetModeSystemModels(this)
            
            % first, cell array with <numModes> elements
            expectedNewModeModels = {this.modeModels{2}; this.modeModels{1}};
            this.modelUnderTest.setModeSystemModels(expectedNewModeModels);
            
            actualNewModeModels = this.modelUnderTest.modeSystemModels;
            
            this.verifyClass(actualNewModeModels, ?cell);
            this.verifyNumElements(actualNewModeModels, this.numModes);
            this.verifyEqual(actualNewModeModels{1}, expectedNewModeModels{1});
            this.verifyEqual(actualNewModeModels{2}, expectedNewModeModels{2});
            
            % first, cell array with <numModes +1> elements
            expectedNewModeModels = {this.modeModels{2}, this.modeModels{1}, this.modeModels{2}};
            this.modelUnderTest.setModeSystemModels(expectedNewModeModels);
            actualNewModeModels = this.modelUnderTest.modeSystemModels;
            
            this.verifyClass(actualNewModeModels, ?cell);
            this.verifyNumElements(actualNewModeModels, this.numModes);
            this.verifyEqual(actualNewModeModels{1}, expectedNewModeModels{1});
            this.verifyEqual(actualNewModeModels{2}, expectedNewModeModels{2});
        end
        
        %% testSetSystemInputInvalidInput
        function testSetSystemInputInvalidInput(this)
            inputs = repmat(42, 1, this.numModes);
            invalidInputs = mat2cell(inputs,1, this.numModes); % inputs given as cell array is invalid
            
            expectedErrId = 'JumpLinearSystemModel:InvalidSystemInput';
            this.verifyError(@() this.modelUnderTest.setSystemInput(invalidInputs), expectedErrId);
            
        end
        
        %% testSetSystemInput
        function testSetSystemInput(this)
            % empty input is valid
            emptyInput = [];
            
            this.modelUnderTest.setSystemInput(emptyInput);
            this.verifyEmpty(this.modelUnderTest.getSystemInput());
            
            % single input for all modes is valid
            singleInput = 42;
            this.modelUnderTest.setSystemInput(singleInput);
            actualInputs = this.modelUnderTest.getSystemInput();
            
            % actual inputs should be column-wise arranged
            this.verifySize(actualInputs, [1 this.numModes]);
            arrayfun(@(i) this.verifyEqual(actualInputs(i), singleInput), 1:this.numModes);
            
            % now different inputs for all modes
            inputs = 42 + [1:this.numModes];
            this.modelUnderTest.setSystemInput(inputs);
            actualInputs = this.modelUnderTest.getSystemInput();
            
            % actual inputs should be column-wise arranged
            this.verifySize(actualInputs, [1 this.numModes]);
            arrayfun(@(i) this.verifyEqual(actualInputs(i), inputs(i)), 1:this.numModes);
        end
        
        %% testSetSystemMatrixForModeInvalidMode
        function testSetSystemMatrixForModeInvalidMode(this)
            newSysMatrix = this.A1 + this.A2;
            
            expectedErrId = 'JumpLinearSystemModel:InvalidMode';
            
            % first, mode is negative
            this.verifyError(@() this.modelUnderTest.setSystemMatrixForMode(newSysMatrix, this.negativeNumModes), expectedErrId);
            % now mode is fractional
            this.verifyError(@() this.modelUnderTest.setSystemMatrixForMode(newSysMatrix, this.fractionaNumModes), expectedErrId);
            % now mode is out of bounds
            this.verifyError(@() this.modelUnderTest.setSystemMatrixForMode(newSysMatrix, this.numModes + 1), expectedErrId);
        end
        
        %% testSetSystemMatrixForMode
        function testSetSystemMatrixForMode(this)
            expectedNewSysMatrix = this.A1 + this.A2;
            mode = 1;
            
            this.modelUnderTest.setSystemMatrixForMode(expectedNewSysMatrix, mode);
            
            actualNewSysMatrix = this.modelUnderTest.modeSystemModels{mode}.sysMatrix;
            this.verifyEqual(actualNewSysMatrix, expectedNewSysMatrix);
        end
        
        %% testSetActiveModeInvalidMode
        function testSetActiveModeInvalidMode(this)
            expectedErrId = 'JumpLinearSystemModel:InvalidMode';
             
            % first, mode is negative
            this.verifyError(@() this.modelUnderTest.setActiveMode(this.negativeNumModes), expectedErrId);
            % now mode is fractional
            this.verifyError(@() this.modelUnderTest.setActiveMode(this.fractionaNumModes), expectedErrId);
            % now mode is out of bounds
            this.verifyError(@() this.modelUnderTest.setActiveMode(this.numModes + 1), expectedErrId);
        end
        
        %% testSimulateInvalidState
        function testSimulateInvalidState(this)
            expectedErrId = 'SystemModel:InvalidSystemState';
            
            invalidState = this.state'; % row vector instead of col vector
            this.verifyError(@() this.modelUnderTest.simulate(invalidState), expectedErrId);
        end
        
        %% testSimulate
        function testSimulate(this)
            newActiveMode = 2;
            this.modelUnderTest.setActiveMode(newActiveMode);
            
            % seeding is required to obtain reproducible results
            rng(JumpLinearSystemModelTest.seed);
            expectedPredictedState = this.modeModels{newActiveMode}.simulate(this.state);
            rng(JumpLinearSystemModelTest.seed);
            actualPredictedState = this.modelUnderTest.simulate(this.state);
            
            this.verifyEqual(actualPredictedState, expectedPredictedState);
        end
        
        %% testIsMeanSquareStableInvalidTransitionMatrix
        function testIsMeanSquareStableInvalidTransitionMatrix(this)
            expectedErrId = 'Validator:ValidateTransitionMatrix:InvalidTransitionMatrixDim';
            
            % invalid transition matrix: not square
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.7 0.1];
            this.verifyError(@() this.modelUnderTest.isMeanSquareStable(invalidTransitionMatrix), expectedErrId); 
            
            % invalid transition matrix: square, but row sums not 1
            invalidTransitionMatrix = [0.7 0.3; 0.2 0.5];
            this.verifyError(@() this.modelUnderTest.isMeanSquareStable(invalidTransitionMatrix), expectedErrId);
            
            % invalid transition matrix: square, but wrong dimension
            invalidTransitionMatrix = [0.7 0.2 0.1; 0.2 0.3 0.5];
            this.verifyError(@() this.modelUnderTest.isMeanSquareStable(invalidTransitionMatrix), expectedErrId);
            
            % invalid transition matrix: square, but elements not all in [0,1]
            invalidTransitionMatrix = [0.7 0.3; 1.1 -0.1];
            this.verifyError(@() this.modelUnderTest.isMeanSquareStable(invalidTransitionMatrix), expectedErrId);
        end
        
        %% testIsMeanSquareStable
        function testIsMeanSquareStable(this)
            % if the transition matrix defined above is utilized to create
            % a MJLS, the system becomes mean square instable (cf. section IV) in
            % 
            %   Maxim Dolgov, Christof Chlebek, and Uwe D. Hanebeck,
            %   Dynamic Compensation of Markov Jump Linear Systems without Mode Observation,
            %   Proceedings of the 2016 European Control Conference (ECC 2016),
            %   Aalborg, Denmark, June 2016.
            
            this.verifyFalse(this.modelUnderTest.isMeanSquareStable(this.transitionMatrix));
            
            % now change the transition matrix
            newTransitionMatrix = [0.7 0.3; 0.2 0.8];
            this.verifyFalse(this.modelUnderTest.isMeanSquareStable(newTransitionMatrix));
        end
    end
end

