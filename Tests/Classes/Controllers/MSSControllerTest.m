classdef MSSControllerTest < matlab.unittest.TestCase
    % Test cases for MSSController.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
        A;
        B;
        
        sequenceLength;
        dimX;
        dimU;
        dimEta;
        controllerDelta;
        initialPlantState;
        expectedNumVertices;
                
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)        
            this.applyFixture(matlab.unittest.fixtures.PathFixture('matlab/external/sdpt3/', 'IncludingSubfolders', true));
            
            this.controllerDelta = 0.2;
            % use an unstable scalar system as example
            this.dimX = 1;
            this.dimU = 1;
                        
            this.A = 1.2;
            this.B = 0.1;
            this.sequenceLength = 3; 
            this.dimEta = this.dimU * this.sequenceLength * (this.sequenceLength - 1) / 2;
            this.expectedNumVertices = 36; % this is the expected number of vertices of the transition matrix polytope
            
            this.initialPlantState = 42; % because that's the answer
            this.controllerUnderTest = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta);
        end
    end
    
    methods (Test)
        %% testCtorInvalidSysMatrix
        function testCtorInvalidSysMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
            
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() MSSController(invalidA, this.B, this.sequenceLength, this.controllerDelta), ...
                expectedErrId);
            
            invalidA = inf(this.dimX, this.dimX); % square but inf
            this.verifyError(@() MSSController(invalidA, this.B, this.sequenceLength, this.controllerDelta), ...
                expectedErrId);
        end
        
        %% testCtorInvalidInputMatrix
        function testCtorInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidB = eye(this.dimX + 1, this.dimU); % invalid dims
            this.verifyError(@() MSSController(this.A, invalidB, this.sequenceLength, this.controllerDelta), ...
                expectedErrId);
            
            invalidB = eye(this.dimX, this.dimU); % correct dims, but nan
            invalidB(1, 1) = nan;
            this.verifyError(@() MSSController(this.A, invalidB, this.sequenceLength, this.controllerDelta), ...
                expectedErrId);
        end
        
        %% testCtorInvalidPlant
        function testCtorInvalidPlant(this)
            expectedErrId = 'MSSController:InvalidPlant';
            
            invalidB = 0; % plant runs open loop, thus not controllable
            this.verifyError(@() MSSController(this.A, invalidB, this.sequenceLength, this.controllerDelta), ...
                expectedErrId);
        end
        
        %% testCtorSolverNotFound
        function testCtorSolverNotFound(this)           
            expectedWarnId = 'MSSController:SolverNotFound';

            % prepare
            rmpath('matlab/external/sdpt3/');
            this.assertFalse(exist('sdpt3', 'file') == 2);
              
            this.verifyWarning(@() MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta), ...
                expectedWarnId);
            
            %restore
            addpath('matlab/external/sdpt3/');
        end
        
        %% testCtor
        function testCtor(this)
            % validate the side effecs: flag useLmiLab is set properly
            % vertices of transition matrix polytope are computed corrected
            this.verifyFalse(this.controllerUnderTest.useLmiLab); % false by default
            
            % check the transition matrix polytope
            expectedSize = [this.sequenceLength + 1, this.sequenceLength + 1, this.expectedNumVertices];
            this.verifySize(this.controllerUnderTest.vertices, expectedSize);            
            
            % check if some particular vertices are present
            P1 = zeros(this.sequenceLength + 1);
            P2 = zeros(this.sequenceLength + 1);
            P3 = zeros(this.sequenceLength + 1);
            P4 = diag([1 1 1], 1);
            P1(:, 1) = 1; % first column = 1
            P2(:, 2) = 1; % second column = 1
            P3(1, 2) = 1;
            P3(2:end, 3) = 1;
            P4(end) = 1;
            
            this.verifyEqual(sum(ismember(P1, this.controllerUnderTest.vertices), 'all'), (this.sequenceLength + 1)^2);
            this.verifyEqual(sum(ismember(P2, this.controllerUnderTest.vertices), 'all'), (this.sequenceLength + 1)^2);
            this.verifyEqual(sum(ismember(P3, this.controllerUnderTest.vertices), 'all'), (this.sequenceLength + 1)^2);
            this.verifyEqual(sum(ismember(P4, this.controllerUnderTest.vertices), 'all'), (this.sequenceLength + 1)^2);
            
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, this.expectedNumVertices-1);            
            T = sdpvar(this.sequenceLength + 1, this.sequenceLength + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:this.expectedNumVertices -1
                sumMat = sumMat + alpha(i) * this.controllerUnderTest.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * this.controllerUnderTest.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytop
            T_feas = Utility.calculateDelayTransitionMatrix([0.1 0.4 0.4 0.1]);
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: T_feas2 is in the polytop
            T_feas2 = Utility.calculateDelayTransitionMatrix([[0.4 0.4 0.1 0.1]', ...
                [1e-20 0.7 1-0.7-2e-20 1e-20]', [0.1 0.4 0.4 0.1]']); % time-varying delay probs
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
            
            % now check for the infeasible T
            T_infeas = T_feas;
            T_infeas(end, end-1) = T_feas(end, end);
            T_infeas(end, end) = T_feas(end, end-1); % last two elements in last row shifted, so last two rows are no longer equal
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % finally check if flag is set according to given param
            controller = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, true);
            this.verifyTrue(controller.useLmiLab);
        end
%%
%%
        %% testSetEtaStateInvalidEta
        function testSetEtaStateInvalidEta(this)
            expectedErrId = 'MSSController:SetEtaState:InvalidEta';
            
            invalidEta = this; % not a vector
            this.verifyError(@() this.controllerUnderTest.setEtaState(invalidEta), expectedErrId);
            
            invalidEta = ones(this.dimEta + 2, 1); % invalid dims
            this.verifyError(@() this.controllerUnderTest.setEtaState(invalidEta), expectedErrId);
            
            invalidEta = ones(this.dimEta, 1); % not finite
            invalidEta(2) = inf;
            this.verifyError(@() this.controllerUnderTest.setEtaState(invalidEta), expectedErrId);
        end
%%
%%
        %% testIsMeanSquareStable
        function testIsMeanSquareStable(this)
            import matlab.unittest.constraints.IsLessThan
            import matlab.unittest.constraints.IsGreaterThanOrEqualTo
            
            [maxRho, isMss] = this.controllerUnderTest.isMeanSquareStable();
            
            this.verifyTrue(isMss);
            this.verifyThat(maxRho, IsGreaterThanOrEqualTo(0) & IsLessThan(1));
        end
        
        %% testIsMeanSquareStableLmiLab
        function testIsMeanSquareStableLmiLab(this)
            import matlab.unittest.constraints.IsLessThan
            import matlab.unittest.constraints.IsGreaterThanOrEqualTo
            
            useLmiLab = true;
            controller = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, useLmiLab);
            
            [maxRho, isMss] = controller.isMeanSquareStable();
            
            this.verifyTrue(isMss);
            this.verifyThat(maxRho, IsGreaterThanOrEqualTo(0) & IsLessThan(1));
        end
%%
%%
        %% testDoControlSequenceComputationZeroState
        function testDoControlSequenceComputationZeroState(this)
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));

            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            this.verifyEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
            
            % sanity check: change eta portion of state, so that augmented
            % state is no longer zero
            newEta = ones(this.dimEta, 1);
            this.controllerUnderTest.setEtaState(newEta);
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            this.verifyNotEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
        end
        
        %% testDoControlSequenceComputationZeroStateLmiLab
        function testDoControlSequenceComputationZeroStateLmiLab(this)
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            useLmiLab = true;
            controller = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, useLmiLab);
            
            actualSequence = controller.computeControlSequence(zeroState);
            
            this.verifyEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
            
            % sanity check: change eta portion of state, so that augmented
            % state is no longer zero
            newEta = ones(this.dimEta, 1);
            this.controllerUnderTest.setEtaState(newEta);
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            this.verifyNotEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
        end
%%
%%
        %% testDoControlSequenceComputation
        function testDoControlSequenceComputation(this)
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            actualSequence = this.controllerUnderTest.computeControlSequence(state);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];            
            
            this.verifySize(actualSequence, expectedSize);            
        end
        
        %% testDoControlSequenceComputationLmiLab
        function testDoControlSequenceComputationLmiLab(this)
            useLmiLab = true;
            controller = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, useLmiLab);
            
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            actualSequence = controller.computeControlSequence(state);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];            
            
            this.verifySize(actualSequence, expectedSize);            
        end
    end
end

