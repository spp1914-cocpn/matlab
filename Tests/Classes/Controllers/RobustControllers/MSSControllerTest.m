classdef MSSControllerTest < matlab.unittest.TestCase
    % Test cases for MSSController.
    
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
    
    properties (Constant, Access = private)
        numConvexCombinations = 5000;
        absTol = 1e-12;
    end
    
    properties (Access = private)
        A;
        B;
        
        sequenceLength;
        dimX;
        dimU;
        dimEta;
        controllerDelta;
        initialPlantState; 
                
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
                        
            this.initialPlantState = 42; % because that's the answer
            this.controllerUnderTest = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, false, false);
        end
    end
    
    methods (Test)
        %% testCtorInvalidSysMatrix
        function testCtorInvalidSysMatrix(this)
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
            
            invalidA = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() MSSController(invalidA, this.B, this.sequenceLength, this.controllerDelta, true, false), ...
                expectedErrId);
            
            invalidA = inf(this.dimX, this.dimX); % square but inf
            this.verifyError(@() MSSController(invalidA, this.B, this.sequenceLength, this.controllerDelta, true, false), ...
                expectedErrId);
        end
        
        %% testCtorInvalidInputMatrix
        function testCtorInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidB = eye(this.dimX + 1, this.dimU); % invalid dims
            this.verifyError(@() MSSController(this.A, invalidB, this.sequenceLength, this.controllerDelta, true, false), ...
                expectedErrId);
            
            invalidB = eye(this.dimX, this.dimU); % correct dims, but nan
            invalidB(1, 1) = nan;
            this.verifyError(@() MSSController(this.A, invalidB, this.sequenceLength, this.controllerDelta, true, false), ...
                expectedErrId);
        end
        
        %% testCtorInvalidPlant
        function testCtorInvalidPlant(this)
            expectedErrId = 'MSSController:InvalidPlant';
            
            invalidB = 0; % plant runs open loop, thus not controllable
            this.verifyError(@() MSSController(this.A, invalidB, this.sequenceLength, this.controllerDelta, true, false), ...
                expectedErrId);
        end
        
        %% testCtorSolverNotFound
        function testCtorSolverNotFound(this)           
            expectedWarnId = 'MSSController:SolverNotFound';

            % prepare
            rmpath('matlab/external/sdpt3/');
            this.assertFalse(exist('sdpt3', 'file') == 2);
              
            this.verifyWarning(@() MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, true, false, false), ...
                expectedWarnId);
            
            %restore
            addpath('matlab/external/sdpt3/');
        end
        
        %% testCtor
        function testCtor(this)
            % validate the side effecs: flag useLmiLab is set properly
            % vertices of transition matrix polytope are computed corrected
            this.verifyFalse(this.controllerUnderTest.useLmiLab); % false by default
            this.verifyFalse(this.controllerUnderTest.assumeCorrelatedDelays); % false by default
            % check if gains struct is properly setup
            this.verifyTrue(isstruct(this.controllerUnderTest.gains));
            this.verifySize(this.controllerUnderTest.gains, [1 1]);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'A'));
            this.verifyEqual(this.controllerUnderTest.gains(1).A, this.A);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'B'));
            this.verifyEqual(this.controllerUnderTest.gains(1).B, this.B);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'L'));
            this.verifyEqual(this.controllerUnderTest.gains(1).L, this.controllerUnderTest.getControllerGain());
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'delta'));
            this.verifyEqual(this.controllerUnderTest.gains(1).delta, this.controllerDelta);
            
            % finally check if flag is set according to given param
            controller = MSSController(this.A, this.B, 2, this.controllerDelta, false, true, true);
            this.verifyTrue(controller.useLmiLab);
            this.verifyTrue(controller.assumeCorrelatedDelays);
            
            % gain is not completely zero and has correct dimensions
            this.verifySize(controller.getControllerGain(), [2, 2])
            this.verifyNotEqual(controller.getControllerGain(), zeros(2,2))
            
            lazyInitGain = true;
            controller = MSSController(this.A, this.B, 2, this.controllerDelta, lazyInitGain, false);
            % side effect: controller gain not yet computed
            this.verifyEmpty(controller.getControllerGain());
            this.verifyEmpty(controller.gains);
        end
%%
%%
        %% testPolytopeSeqOne
        function testPolytopeSeqOne(this)
            seqLenth = 1;
            expectedNumVertices = 2; % this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, false);
            
             % check the transition matrix polytope if sequence length is 2
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1]');
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present            
            P1 = [1-this.controllerDelta, this.controllerDelta;
                  1-this.controllerDelta, this.controllerDelta];
                                    
            this.verifyEqual(sum(ismember(P1, controller.vertices), 'all'), (seqLenth + 1)^2);
            
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = Utility.calculateDelayTransitionMatrix([0.9 0.1]);
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);            
                 
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
         
            % now check for the infeasible T
            T_infeas = T_feas;
            T_infeas(end, end-1) = T_feas(end, end);
            T_infeas(end, end) = T_feas(end, end-1); % last two elements in last row shifted, so last two rows are no longer equal
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(1, 1) = T_feas(1, 2);
            T_infeas(1, 2) = T_feas(1, 1); % shift the elements in first row, first column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible        
            
            % another check for the infeasible T
            T_infeas = Utility.calculateDelayTransitionMatrix([0.6 0.4]); % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % check the column and row constraints     
                this.verifyEqual(T_res(1, :), T_res(2, :));
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(2, 2), this.controllerDelta);   
            end             
        end
        
        %% testPolytopeSeqOneCorrelatedDelays
        function testPolytopeSeqOneCorrelatedDelays(this)
            seqLenth = 1;
            expectedNumVertices = 4; % this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, true);
            
            % check the transition matrix polytope if sequence length is 1
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1]', 'AbsTol', 1e-15);
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present            
            P1 = [1-this.controllerDelta, this.controllerDelta;
                  1-this.controllerDelta, this.controllerDelta];           
            this.verifyEqual(sum(ismembertol(P1, controller.vertices, 1e-15), 'all'), (seqLenth + 1)^2);
            
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = [0.9 0.1; 0.9 0.1];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);            
                 
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
         
            % now check for the infeasible T            
            T_infeas = Utility.calculateDelayTransitionMatrix([0.6 0.4]); % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);           
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(2, 2), this.controllerDelta);   
            end
        end
%%
%%
        %% testPolytopeSeqTwo
        function testPolytopeSeqTwo(this)
            seqLenth = 2;
            expectedNumVertices = 4; % this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, false);
            
            % check the transition matrix polytope if sequence length is 2
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1 1]');
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present
            P1 = [0, 1, 0;                  
                  0, 1-this.controllerDelta, this.controllerDelta;
                  0, 1-this.controllerDelta, this.controllerDelta];
            P2 = [1-this.controllerDelta, this.controllerDelta, 0;
                  1-this.controllerDelta, 0, this.controllerDelta;
                  1-this.controllerDelta, 0, this.controllerDelta];
                                    
            this.verifyEqual(sum(ismember(P1, controller.vertices), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismember(P2, controller.vertices), 'all'), (seqLenth + 1)^2);
                        
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = Utility.calculateDelayTransitionMatrix([0.1 0.8 0.1]);
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: T_feas2 is in the polytope
            T_feas2 = Utility.calculateDelayTransitionMatrix([[0.8 0.1 0.1]', ...
                [0.1 0.8 0.1]']); % time-varying delay probs
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 0] = 1
            T_feas3 = [1 0 0;
                       1 0 0;
                       1 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas3);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 1] = 1
            T_feas4 = [0 1 0;
                       0 1 0;
                       0 1 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas4);            
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
         
            % now check for the infeasible T
            T_infeas = T_feas;
            T_infeas(end, end-1) = T_feas(end, end);
            T_infeas(end, end) = T_feas(end, end-1); % last two elements in last row shifted, so last two rows are no longer equal
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(1, 1) = T_feas(1, 2);
            T_infeas(1, 2) = T_feas(1, 1); % shift the elements in first row, first column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % corner case: Pr[delay = 2] = 1
            T_infeas = [0 1 0;
                        0 0 1;
                        0 0 1];
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % not a hessenberg matrix, infeasible as well
            T_infeas = abs(normalize(gallery('moler', 3), 2, 'norm', 1));
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = Utility.calculateDelayTransitionMatrix([0.2 0.4 0.4]); % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
             
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % check the column and row constraints     
                this.verifyEqual(T_res(:, 1), [T_res(1,1);T_res(1,1);T_res(1,1)]);
                this.verifyEqual(T_res(2:end, 2), [T_res(2,2); T_res(2,2)]);
                this.verifyEqual(T_res(3, :), T_res(2, :));
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(3, 3), this.controllerDelta);
                % check if zero entries are as expected
                this.verifyEqual(T_res(1,3), 0);     
            end            
        end

        %% testPolytopeSeqTwoCorrelatedDelays
        function testPolytopeSeqTwoCorrelatedDelays(this)
            seqLenth = 2;
            expectedNumVertices = 2 * seqLenth * factorial(seqLenth + 1); % this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, true);
            
            % check the transition matrix polytope if sequence length is 1
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1 1]', 'AbsTol', 1e-15);
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present
            P1 = [0, 1, 0;                  
                  0, 0, 1;
                  0, 1-this.controllerDelta, this.controllerDelta];
            P2 = [0, 1, 0;                  
                  0, 0, 1;
                  1, 0, 0];
            P3 = [1, 0, 0;                  
                  0, 0, 1;
                  0.8, 0, 0.2];
            this.verifyEqual(sum(ismembertol(P1, controller.vertices, 1e-15), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismembertol(P2, controller.vertices, 1e-15), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismembertol(P3, controller.vertices, 1e-15), 'all'), (seqLenth + 1)^2);
            
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = [0.43 0.57 0; 0.15 0.08 0.77; 0.19 0.8 0.01];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);            
            
             % this problem is supposed to be feasible: T_feas2 is in the polytope
             % white delays are a special case of correlated delays
            T_feas2 = [0.1 0.9 0; 0.1 0.7 0.2; 0.1 0.7 0.2];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 0] = 1
            T_feas3 = [1 0 0;
                       1 0 0;
                       1 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas3);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 1] = 1
            T_feas4 = [0 1 0;
                       0 1 0;
                       0 1 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas4);            
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
            
            % now check for the infeasible T            
            T_infeas = [0.07 0.93 0; 0.2 0.5 0.3; 0.16 0.19 0.65]; % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % validate the Hessenberg structure
                this.verifyEqual(T_res(1, end), 0);                
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(end, end), this.controllerDelta);
                % check the row sums
                this.verifyEqual(sum(T_res, 2), [1 1 1]', 'AbsTol', 1e-15);
            end
        end
%%
%%
        %% testPolytopeSeqThree
        function testPolytopeSeqThree(this)
            seqLenth = 3;
            expectedNumVertices = 6; % this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, false);
            
            % check the transition matrix polytope if sequence length is  3
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1 1 1]');
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present
            P1 = [0, 1, 0, 0;
                  0, 0, 1, 0;
                  0, 0, 1-this.controllerDelta, this.controllerDelta;
                  0, 0, 1-this.controllerDelta, this.controllerDelta];
            P2 = [0, 1, 0, 0;
                  0, 1-this.controllerDelta, this.controllerDelta, 0;
                  0, 1-this.controllerDelta, this.controllerDelta, 0;
                  0, 1-this.controllerDelta, this.controllerDelta, 0];
                        
            this.verifyEqual(sum(ismember(P1, controller.vertices), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismember(P2, controller.vertices), 'all'), (seqLenth + 1)^2);
            
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = Utility.calculateDelayTransitionMatrix([0.1 0.4 0.4 0.1]);
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: T_feas2 is in the polytope
            T_feas2 = Utility.calculateDelayTransitionMatrix([[0.4 0.4 0.1 0.1]', ...
                [1e-20 0.7 1-0.7-2e-20 1e-20]', [0.1 0.4 0.4 0.1]']); % time-varying delay probs
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 0] = 1
            T_feas3 = [1 0 0 0;
                       1 0 0 0;
                       1 0 0 0;
                       1 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas3);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 1] = 1
            T_feas4 = [0 1 0 0;
                       0 1 0 0;
                       0 1 0 0;
                       0 1 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas4);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 2] = 1
            T_feas5 = [0 1 0 0;
                       0 0 1 0;
                       0 0 1 0;
                       0 0 1 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas5);            
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
         
            % now check for the infeasible T
            T_infeas = T_feas;
            T_infeas(end, end-1) = T_feas(end, end);
            T_infeas(end, end) = T_feas(end, end-1); % last two elements in last row shifted, so last two rows are no longer equal
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(1, 1) = T_feas(1, 2);
            T_infeas(1, 2) = T_feas(1, 1); % shift the elements in first row, first column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(2, 2) = T_feas(2, 3);
            T_infeas(2, 3) = T_feas(2, 2); % shift the elements in second row, second column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % corner case: Pr[delay = 3] = 1
            T_infeas = [0 1 0 0;
                        0 0 1 0;
                        0 0 0 1;
                        0 0 0 1];
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % not a hessenberg matrix, infeasible as well
            T_infeas = abs(normalize(gallery('moler', 4), 2, 'norm', 1));
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = Utility.calculateDelayTransitionMatrix([0.1 0.4 0.1 0.4]); % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
             
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % check the column constraints     
                this.verifyEqual(repmat(T_res(1, 1), seqLenth, 1),T_res(2:end, 1));
                this.verifyEqual(repmat(T_res(2, 2), seqLenth-1, 1),T_res(3:end, 2));
                this.verifyEqual(T_res(3,3),T_res(4, 3));
                this.verifyEqual(T_res(3,4),T_res(4, 4));
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(4, 4), this.controllerDelta);
                % check if zero entries are as expected
                this.verifyEqual(T_res(1,3:4), [0 0]);
                this.verifyEqual(T_res(2,4), 0);
            end
            
        end
        
        
        %% testPolytopeSeqThreeCorrelatedDelays
        function testPolytopeSeqThreeCorrelatedDelays(this)
            seqLenth = 3;
            expectedNumVertices = 2 * seqLenth * factorial(seqLenth + 1); % 144, this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, true);
            
            % check the transition matrix polytope if sequence length is 1
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1 1 1]', 'AbsTol', 1e-14);
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present
            P1 = [0, 1, 0, 0;                  
                  0, 0, 1, 0;
                  0, 0, 0, 1
                  0, 0, 1-this.controllerDelta, this.controllerDelta];
            P2 = [0, 1, 0, 0;                  
                  0, 0, 1, 0;
                  0, 0, 1, 0
                  0, 1-this.controllerDelta, 0, this.controllerDelta];
            P3 = [0, 1, 0, 0;                  
                  0, 0, 1, 0;
                  0, 0, 0, 1
                  1-this.controllerDelta, 0, 0, this.controllerDelta];
            this.verifyEqual(sum(ismembertol(P1, controller.vertices, 1e-15), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismembertol(P2, controller.vertices, 1e-15), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismembertol(P3, controller.vertices, 1e-15), 'all'), (seqLenth + 1)^2);
            
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = [0.4294 0.5706 0 0; 
                      0.0657 0.0374 0.8969 0; 
                      0.1225 0.1913 0.1993 0.4869
                      0.1216 0.4885 0.2497 0.1402];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);            
            
            % this problem is supposed to be feasible: T_feas2 is in the polytope
            % white delays are a special case of correlated delays
            T_feas2 = [0 1 0 0; 
                       0 0.1 0.9 0;
                       0 0.1 0.7 0.2
                       0 0.1 0.7 0.2];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
            
           % corner case: Pr[delay = 0] = 1
            T_feas3 = [1 0 0 0;
                       1 0 0 0;
                       1 0 0 0;
                       1 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas3);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 1] = 1
            T_feas4 = [0 1 0 0;
                       0 1 0 0;
                       0 1 0 0;
                       0 1 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas4);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 2] = 1
            T_feas5 = [0 1 0 0;
                       0 0 1 0;
                       0 0 1 0;
                       0 0 1 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas5);            
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
            
            % now check for the infeasible T
            % last entry is larger then delta
            T_infeas = [0.4294 0.5706 0 0; 
                        0.0657 0.0374 0.8969 0; 
                        0.1225 0.1913 0.1993 0.4869
                        0.1216 0.1885 0.1497 0.5402]; 
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % validate the Hessenberg structure
                this.verifyEqual(T_res(1, end), 0);
                this.verifyEqual(T_res(1, end-1:end), [0 0]); 
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(end, end), this.controllerDelta);
                % check the row sums
                this.verifyEqual(sum(T_res, 2), [1 1 1 1]', 'AbsTol', 1e-14);
            end
        end
%%
%%
        %% testPolytopeSeqFour
        function testPolytopeSeqFour(this)
            seqLenth = 4;
            expectedNumVertices = 8; % this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, false);
            
            % check the transition matrix polytope if sequence length is  4
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1 1 1 1]');
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present
            P1 = [0, 1, 0, 0, 0;
                  0, 0, 1, 0, 0;
                  0, 0, 0, 1, 0
                  0, 0, 0, 1-this.controllerDelta, this.controllerDelta;
                  0, 0, 0, 1-this.controllerDelta, this.controllerDelta];
            P2 = [1-this.controllerDelta, this.controllerDelta, 0, 0, 0;
                  1-this.controllerDelta, 0, this.controllerDelta, 0, 0;
                  1-this.controllerDelta, 0, 0, this.controllerDelta, 0;
                  1-this.controllerDelta, 0,0,0, this.controllerDelta;
                  1-this.controllerDelta, 0,0,0, this.controllerDelta];
                        
            this.verifyEqual(sum(ismember(P1, controller.vertices), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismember(P2, controller.vertices), 'all'), (seqLenth + 1)^2);
            
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = Utility.calculateDelayTransitionMatrix([0.1 0.4 0.2 0.2 0.1]);
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: T_feas2 is in the polytope
            T_feas2 = Utility.calculateDelayTransitionMatrix([[0.1 0.4 0.2 0.2 0.1]', ...
                [0.2 0.2 0.4 0.1 0.1]', [0.4 0.1 0.2 0.2 0.1]', [0.2 0.2 0.2 0.3 0.1]']); % time-varying delay probs
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 0] = 1
            T_feas3 = [1 0 0 0 0;
                       1 0 0 0 0;
                       1 0 0 0 0;
                       1 0 0 0 0;
                       1 0 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas3);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 1] = 1
            T_feas4 = [0 1 0 0 0;
                       0 1 0 0 0;
                       0 1 0 0 0;
                       0 1 0 0 0;
                       0 1 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas4);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 2] = 1
            T_feas5 = [0 1 0 0 0;
                       0 0 1 0 0;
                       0 0 1 0 0;
                       0 0 1 0 0;
                       0 0 1 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas5);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 3] = 1
            T_feas6 = [0 1 0 0 0;
                       0 0 1 0 0;
                       0 0 0 1 0;
                       0 0 0 1 0;
                       0 0 0 1 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas6);            
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
         
            % now check for the infeasible T
            T_infeas = T_feas;
            T_infeas(end, end-1) = T_feas(end, end);
            T_infeas(end, end) = T_feas(end, end-1); % last two elements in last row shifted, so last two rows are no longer equal
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(1, 1) = T_feas(1, 2);
            T_infeas(1, 2) = T_feas(1, 1); % shift the elements in first row, first column has different values
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(2, 2) = T_feas(2, 3);
            T_infeas(2, 3) = T_feas(2, 2); % shift the elements in second row, second column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(3, 3) = T_feas(3, 4);
            T_infeas(3, 4) = T_feas(3, 3); % shift the elements in third row, third column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % corner case: Pr[delay = 4] = 1
            T_infeas = [0 1 0 0 0;
                        0 0 1 0 0;
                        0 0 0 1 0;
                        0 0 0 0 1;
                        0 0 0 0 1];
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % not a hessenberg matrix, infeasible as well
            T_infeas = abs(normalize(gallery('moler', 5), 2, 'norm', 1));
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = Utility.calculateDelayTransitionMatrix([0.1 0.2 0.1 0.2 0.4]); % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
             
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % check the column constraints     
                this.verifyEqual(repmat(T_res(1, 1), seqLenth, 1),T_res(2:end, 1));
                this.verifyEqual(repmat(T_res(2, 2), seqLenth-1, 1),T_res(3:end, 2));
                this.verifyEqual(repmat(T_res(3, 3), seqLenth-2, 1),T_res(4:end, 3));
                this.verifyEqual(T_res(4, :),T_res(5, :));
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(5, 5), this.controllerDelta);
                % check if zero entries are as expected
                this.verifyEqual(T_res(1,3:end), [0 0 0]);
                this.verifyEqual(T_res(2, 4:end), [0 0]);
                this.verifyEqual(T_res(3, end), 0);
            end            
        end
        
        %% testPolytopeSeqFive
        function testPolytopeSeqFive(this)
            seqLenth = 5;
            expectedNumVertices = 10; % this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, false);
            
            % check the transition matrix polytope if sequence length is  5
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1 1 1 1 1]');
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if some particular vertices are present
            P1 = [0, 1, 0, 0, 0, 0;
                  0, 0, 1, 0, 0, 0;
                  0, 0, 0, 1, 0, 0;
                  0, 0, 0, 0, 1, 0;
                  0, 0, 0, 0, 1-this.controllerDelta, this.controllerDelta;
                  0, 0, 0, 0, 1-this.controllerDelta, this.controllerDelta];
            P2 = [0, 1, 0, 0, 0, 0;
                  0, 1-this.controllerDelta, this.controllerDelta, 0, 0, 0;
                  0, 1-this.controllerDelta, 0, this.controllerDelta, 0, 0;
                  0, 1-this.controllerDelta, 0, 0, this.controllerDelta, 0;
                  0, 1-this.controllerDelta, 0, 0, 0, this.controllerDelta;
                  0, 1-this.controllerDelta, 0, 0, 0, this.controllerDelta];
            P3 = [0, 1, 0, 0, 0, 0;
                  0, 0, 1, 0, 0, 0;
                  0, 0, 0, 1, 0, 0;
                  0, 0, 0, 1-this.controllerDelta, this.controllerDelta, 0;
                  0, 0, 0, 1-this.controllerDelta, 0, this.controllerDelta;
                  0, 0, 0, 1-this.controllerDelta, 0, this.controllerDelta];
              
              
            this.verifyEqual(sum(ismember(P1, controller.vertices), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismember(P2, controller.vertices), 'all'), (seqLenth + 1)^2);
            this.verifyEqual(sum(ismember(P3, controller.vertices), 'all'), (seqLenth + 1)^2);
               
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = Utility.calculateDelayTransitionMatrix([0.1 0.4 0.2 0.2 0.06 0.04]);
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: T_feas2 is in the polytope
            T_feas2 = Utility.calculateDelayTransitionMatrix([[0.1 0.4 0.2 0.2 0.06 0.04]', ...
                [0.1 0.4 0.2 0.2 0.04 0.06]', [0.4 0.1 0.2 0.2 0.06 0.04]', ...
                [0.2 0.2 0.2 0.3 0.05 0.05]', [0.1 0.1 0.1 0.3 0.2 0.2]']); % time-varying delay probs
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 0] = 1
            T_feas3 = [1 0 0 0 0 0;
                       1 0 0 0 0 0;
                       1 0 0 0 0 0;
                       1 0 0 0 0 0;
                       1 0 0 0 0 0;
                       1 0 0 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas3);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 1] = 1
            T_feas4 = [0 1 0 0 0 0;
                       0 1 0 0 0 0;
                       0 1 0 0 0 0;
                       0 1 0 0 0 0;
                       0 1 0 0 0 0;
                       0 1 0 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas4);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 2] = 1
            T_feas5 = [0 1 0 0 0 0;
                       0 0 1 0 0 0;
                       0 0 1 0 0 0;
                       0 0 1 0 0 0;
                       0 0 1 0 0 0;
                       0 0 1 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas5);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 3] = 1
            T_feas6 = [0 1 0 0 0 0;
                       0 0 1 0 0 0;
                       0 0 0 1 0 0;
                       0 0 0 1 0 0;
                       0 0 0 1 0 0;
                       0 0 0 1 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas6);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 4] = 1
            T_feas7 = [0 1 0 0 0 0;
                       0 0 1 0 0 0;
                       0 0 0 1 0 0;
                       0 0 0 0 1 0;
                       0 0 0 0 1 0;
                       0 0 0 0 1 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas7);            
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
         
            % now check for the infeasible T
            T_infeas = T_feas;
            T_infeas(end, end-1) = T_feas(end, end);
            T_infeas(end, end) = T_feas(end, end-1); % last two elements in last row shifted, so last two rows are no longer equal
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(1, 1) = T_feas(1, 2);
            T_infeas(1, 2) = T_feas(1, 1); % shift the elements in first row, first column has different values
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(2, 2) = T_feas(2, 3);
            T_infeas(2, 3) = T_feas(2, 2); % shift the elements in second row, second column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(3, 3) = T_feas(3, 4);
            T_infeas(3, 4) = T_feas(3, 3); % shift the elements in third row, third column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % corner case: Pr[delay = 5] = 1
            T_infeas = [0 1 0 0 0 0;
                        0 0 1 0 0 0;
                        0 0 0 1 0 0;
                        0 0 0 0 1 0;
                        0 0 0 0 0 1;
                        0 0 0 0 0 1];
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % not a hessenberg matrix, infeasible as well
            T_infeas = abs(normalize(gallery('moler', 6), 2, 'norm', 1));
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = Utility.calculateDelayTransitionMatrix([0.05 0.2 0.1 0.2 0.05 0.4]); % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
             
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % check the column constraints     
                this.verifyEqual(repmat(T_res(1, 1), seqLenth, 1),T_res(2:end, 1));
                this.verifyEqual(repmat(T_res(2, 2), seqLenth-1, 1),T_res(3:end, 2));
                this.verifyEqual(repmat(T_res(3, 3), seqLenth-2, 1),T_res(4:end, 3));
                this.verifyEqual(repmat(T_res(4, 4), seqLenth-3, 1),T_res(5:end, 4));
                this.verifyEqual(T_res(5, :),T_res(6, :));
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(6, 6), this.controllerDelta);
                % check if zero entries are as expected
                this.verifyEqual(T_res(1,3:end), [0 0 0 0]);
                this.verifyEqual(T_res(2, 4:end), [0 0 0]);
                this.verifyEqual(T_res(3, 5:end), [0 0]);
                this.verifyEqual(T_res(4, end), 0);
            end            
        end
        
        %% testPolytopeSeqSix
        function testPolytopeSeqSix(this)
            seqLenth = 6;
            expectedNumVertices = 12; % 2*6 this is the expected number of vertices of the transition matrix polytope
            controller = MSSController(this.A, this.B, seqLenth, this.controllerDelta, true, false);
            
            % check the transition matrix polytope if sequence length is  6
            expectedSize = [seqLenth + 1 seqLenth + 1 expectedNumVertices];
            this.verifySize(controller.vertices, expectedSize);            
            % all vertices are different and valid transition matrices
            for i=1:expectedNumVertices
                this.verifyEqual(sum(controller.vertices(:, :, i), 2), [1 1 1 1 1 1 1]');
                for j=i+1:expectedNumVertices                    
                    this.verifyNotEqual(controller.vertices(:, :, j), controller.vertices(:, :, i));                    
                end
            end
            
            % check if a particular vertex is present
            P1 = [0, 1, 0, 0, 0, 0, 0;
                  0, 0, 1, 0, 0, 0, 0;
                  0, 0, 0, 1, 0, 0, 0;
                  0, 0, 0, 0, 1, 0, 0;
                  0, 0, 0, 0, 0, 1, 0;
                  0, 0, 0, 0, 0, 1-this.controllerDelta, this.controllerDelta;
                  0, 0, 0, 0, 0, 1-this.controllerDelta, this.controllerDelta];              
              
            this.verifyEqual(sum(ismember(P1, controller.vertices), 'all'), (seqLenth + 1)^2); 
%             
            % check for a hessenberg transition matrix          
            % use yalmip to solve the problem
            alpha = sdpvar(1, expectedNumVertices-1);            
            T = sdpvar(seqLenth + 1, seqLenth + 1, 'full'); % parameter
            sumMat = 0;
            for i=1:expectedNumVertices -1
                sumMat = sumMat + alpha(i) * controller.vertices(:, :, i); 
            end
            sumMat = sumMat + (1-sum(alpha)) * controller.vertices(:, :, end); % last element
            constraints = [0<=alpha, sum(alpha) <= 1, sumMat == T];       
            options = sdpsettings('solver', 'sdpt3');
            options.verbose = 0; % mute the solver
            solver = optimizer(constraints, [], options, T, alpha);
            
            % this problem is supposed to be feasible: T_feas is in the polytope
            T_feas = Utility.calculateDelayTransitionMatrix([0.1 0.4 0.2 0.15 0.05 0.06 0.04]);
            % spanned by the vertices
            [~, errorcode] = solver(T_feas);
            this.verifyTrue(errorcode == 0);
            
            % this problem is supposed to be feasible: T_feas2 is in the polytope
            T_feas2 = Utility.calculateDelayTransitionMatrix([[0.1 0.4 0.2 0.15 0.05 0.06 0.04]', ...
                [0.1 0.4 0.2 0.1 0.1 0.04 0.06]', [0.05 0.35 0.1 0.2 0.2 0.06 0.04]', ...
                [0.2 0.2 0.2 0.1 0.2 0.05 0.05]', [0.1 0.1 0.1 0.15 0.15 0.2 0.2]', ...
                [0.15 0.15 0.1 0.1 0.1 0.2 0.2]']); % time-varying delay probs
            % spanned by the vertices
            [~, errorcode] = solver(T_feas2);
            this.verifyTrue(errorcode == 0);
                        
            % corner case: Pr[delay = 0] = 1
            T_feas3 = [1 0 0 0 0 0 0;
                       1 0 0 0 0 0 0;
                       1 0 0 0 0 0 0;
                       1 0 0 0 0 0 0;
                       1 0 0 0 0 0 0;
                       1 0 0 0 0 0 0;
                       1 0 0 0 0 0 0;];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas3);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 1] = 1
            T_feas4 = [0 1 0 0 0 0 0;
                       0 1 0 0 0 0 0;
                       0 1 0 0 0 0 0;
                       0 1 0 0 0 0 0;
                       0 1 0 0 0 0 0;
                       0 1 0 0 0 0 0;
                       0 1 0 0 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas4);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 2] = 1
            T_feas5 = [0 1 0 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 1 0 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas5);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 3] = 1
            T_feas6 = [0 1 0 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 0 1 0 0 0;
                       0 0 0 1 0 0 0;
                       0 0 0 1 0 0 0;
                       0 0 0 1 0 0 0;
                       0 0 0 1 0 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas6);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 4] = 1
            T_feas7 = [0 1 0 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 0 1 0 0 0;
                       0 0 0 0 1 0 0;
                       0 0 0 0 1 0 0;
                       0 0 0 0 1 0 0;
                       0 0 0 0 1 0 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas7);            
            this.verifyTrue(errorcode == 0);
            
            % corner case: Pr[delay = 5] = 1
            T_feas8 = [0 1 0 0 0 0 0;
                       0 0 1 0 0 0 0;
                       0 0 0 1 0 0 0;
                       0 0 0 0 1 0 0;
                       0 0 0 0 0 1 0;
                       0 0 0 0 0 1 0;
                       0 0 0 0 0 1 0];
            % spanned by the vertices
            [~, errorcode] = solver(T_feas8);            
            this.verifyTrue(errorcode == 0);
                 
            % this problem is supposed to be feasible: each vertex is in
            % the polytope
            for j=1:expectedNumVertices                
                [~, errorcode] = solver(controller.vertices(:, :, j));
                this.verifyTrue(errorcode == 0);
            end
         
            % now check for the infeasible T
            T_infeas = T_feas;
            T_infeas(end, end-1) = T_feas(end, end);
            T_infeas(end, end) = T_feas(end, end-1); % last two elements in last row shifted, so last two rows are no longer equal
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(1, 1) = T_feas(1, 2);
            T_infeas(1, 2) = T_feas(1, 1); % shift the elements in first row, first column has different values
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(2, 2) = T_feas(2, 3);
            T_infeas(2, 3) = T_feas(2, 2); % shift the elements in second row, second column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = T_feas;
            T_infeas(3, 3) = T_feas(3, 4);
            T_infeas(3, 4) = T_feas(3, 3); % shift the elements in third row, third column has different values
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % corner case: Pr[delay = 6] = 1
            T_infeas = [0 1 0 0 0 0 0;
                        0 0 1 0 0 0 0;
                        0 0 0 1 0 0 0;
                        0 0 0 0 1 0 0;
                        0 0 0 0 0 1 0;
                        0 0 0 0 0 0 1;
                        0 0 0 0 0 0 1];
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % not a hessenberg matrix, infeasible as well
            T_infeas = abs(normalize(gallery('moler', seqLenth + 1), 2, 'norm', 1));
            [~, errorcode] = solver(T_infeas);            
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
            
            % another check for the infeasible T
            T_infeas = Utility.calculateDelayTransitionMatrix([0.05 0.2 0.1 0.18 0.02 0.05 0.4]); % last entry is larger then delta
            [~, errorcode] = solver(T_infeas);
            this.verifyTrue(errorcode == 1); % supposed to be infeasible
             
            % another sanity check: randomly generate convex combination of
            % vertices and check structure of resulting transition matrix
            weights = randfixedsum(expectedNumVertices, MSSControllerTest.numConvexCombinations, 1, 0, 1);
            
            for j=1:MSSControllerTest.numConvexCombinations                
                T_res = sum(controller.vertices .* reshape(weights(:, j), 1, 1, expectedNumVertices), 3);
                % check the column constraints     
                this.verifyEqual(repmat(T_res(1, 1), seqLenth, 1),T_res(2:end, 1));
                this.verifyEqual(repmat(T_res(2, 2), seqLenth-1, 1),T_res(3:end, 2));
                this.verifyEqual(repmat(T_res(3, 3), seqLenth-2, 1),T_res(4:end, 3));
                this.verifyEqual(repmat(T_res(4, 4), seqLenth-3, 1),T_res(5:end, 4));
                this.verifyEqual(repmat(T_res(5, 5), seqLenth-4, 1),T_res(6:end, 5));
                this.verifyEqual(T_res(6, :),T_res(7, :), 'AbsTol', MSSControllerTest.absTol); % last two columns
                % check if last entry is less than or equal to specified delta
                this.verifyLessThanOrEqual(T_res(7, 7), this.controllerDelta);
                % check if zero entries are as expected
                this.verifyEqual(T_res(1,3:end), [0 0 0 0 0]);
                this.verifyEqual(T_res(2, 4:end), [0 0 0 0]);
                this.verifyEqual(T_res(3, 5:end), [0 0 0]);
                this.verifyEqual(T_res(4, 6:end), [0 0]);
                this.verifyEqual(T_res(5, end), 0);
             end            
        end
%%
%%
        %% testChangeCaDelayProbs
        function testChangeCaDelayProbs(this)
            this.assertTrue(isa(this.controllerUnderTest, 'CaDelayProbsChangeable'));
            
            L_old = this.controllerUnderTest.getControllerGain();
            
            % new delay probs yield a new controller delta and thus new vertices
            newDelayProbs = [0 0.25 0.5 0.2 0.05];         
            expectedNewControllerDelta = newDelayProbs(end) + newDelayProbs(end-1);
            
            this.controllerUnderTest.changeCaDelayProbs(newDelayProbs);
            
            actualNewControllerDelta = this.controllerUnderTest.controllerDelta;            
            this.verifyEqual(actualNewControllerDelta, expectedNewControllerDelta);
            
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            % this call triggers the recomputation of the gain
            this.controllerUnderTest.computeControlSequence(zeroState);
           
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            L_new = this.controllerUnderTest.getControllerGain();
            
            this.verifyEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
            
            % gain should have changed
            this.verifyNotEqual(L_new, L_old);
            
            % sanity check: change eta portion of state, so that augmented
            % state is no longer zero
            newEta = ones(this.dimEta, 1);
            this.controllerUnderTest.setEtaState(newEta);
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            this.verifyNotEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
            
            % check if gains struct is properly updated
            this.verifyTrue(isstruct(this.controllerUnderTest.gains));
            this.verifySize(this.controllerUnderTest.gains, [1 2]);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'A'));
            this.verifyEqual(this.controllerUnderTest.gains(1).A, this.A);
            this.verifyEqual(this.controllerUnderTest.gains(2).A, this.A); % A did not change
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'B'));
            this.verifyEqual(this.controllerUnderTest.gains(1).B, this.B);
            this.verifyEqual(this.controllerUnderTest.gains(2).B, this.B); % B did not change
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'L'));
            this.verifyEqual(this.controllerUnderTest.gains(1).L, L_old);
            this.verifyEqual(this.controllerUnderTest.gains(2).L, L_new);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'delta'));
            this.verifyEqual(this.controllerUnderTest.gains(1).delta, this.controllerDelta);
            this.verifyEqual(this.controllerUnderTest.gains(2).delta, expectedNewControllerDelta);
            
        end
        
        %% testChangeCaDelayProbsSameDelta
        function testChangeCaDelayProbsSameDelta(this)
           this.assertTrue(isa(this.controllerUnderTest, 'CaDelayProbsChangeable'));
            
            L_old = this.controllerUnderTest.getControllerGain();
            
            % new delay probs but no new controller delta and no new vertices
            newDelayProbs = [0 0.2 0.6 0.2];         
            expectedNewControllerDelta = this.controllerDelta;
            
            this.controllerUnderTest.changeCaDelayProbs(newDelayProbs);
            
            actualNewControllerDelta = this.controllerUnderTest.controllerDelta;            
            this.verifyEqual(actualNewControllerDelta, expectedNewControllerDelta);                        
            
            % this call triggers the recomputation of the gain
            this.controllerUnderTest.computeControlSequence(Gaussian(zeros(this.dimX, 1), eye(this.dimX)));            
            L_new = this.controllerUnderTest.getControllerGain();
            
            % no change expected
            this.verifyEqual(L_new, L_old, 'AbsTol', 1e-8);
                        
            % accordingly, gains struct should be the same as before
            this.verifyTrue(isstruct(this.controllerUnderTest.gains));
            this.verifySize(this.controllerUnderTest.gains, [1 1]);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'A'));
            this.verifyEqual(this.controllerUnderTest.gains(1).A, this.A);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'B'));
            this.verifyEqual(this.controllerUnderTest.gains(1).B, this.B);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'L'));
            this.verifyEqual(this.controllerUnderTest.gains(1).L, L_old);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'delta'));
            this.verifyEqual(this.controllerUnderTest.gains(1).delta, this.controllerDelta);            
        end
%%
%%
        %% testChangeModelParametersInvalidSystemMatrix
        function testChangeModelParametersInvalidSystemMatrix(this)
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
            
            expectedErrId = 'Validator:ValidateSystemMatrix:InvalidDimensions';
            
            invalidA = eye(4); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
                        
            invalidA = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
            
            invalidA = eye(2,3); % not square
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
            
            invalidA = inf; % not finite
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(invalidA, this.B), expectedErrId);
        end
        
        %% testChangeModelParametersInvalidInputMatrix
        function testChangeModelParametersInvalidInputMatrix(this)
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
            
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrixDims';
            
            invalidB = ones(2,2); % wrong dimension
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
                        
            invalidB = this; % not a matrix
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
                        
            invalidB = inf; % not finite      
            this.verifyError(@() this.controllerUnderTest.changeModelParameters(this.A, invalidB), expectedErrId);
        end
        
        %% testChangeModelParameters
        function testChangeModelParameters(this)
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
            
            L_old = this.controllerUnderTest.getControllerGain();
            
            % change A and B
            newA = -this.A;
            newB = -this.B;
                        
            this.controllerUnderTest.changeModelParameters(newA, newB);
            
            % perform a sanity check: given state is origin, so computed control sequence should be also the zero vector
            % due to the underlying linear control law
            zeroState = Gaussian(zeros(this.dimX, 1), eye(this.dimX));
            % this call triggers the recomputation of the gain
            this.controllerUnderTest.computeControlSequence(zeroState);
           
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            L_new = this.controllerUnderTest.getControllerGain();
            
            this.verifyEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
            
            % gain should have changed
            this.verifyNotEqual(L_new, L_old);
            
            % sanity check: change eta portion of state, so that augmented
            % state is no longer zero
            newEta = ones(this.dimEta, 1);
            this.controllerUnderTest.setEtaState(newEta);
            actualSequence = this.controllerUnderTest.computeControlSequence(zeroState);
            
            this.verifyNotEqual(actualSequence, zeros(this.dimU * this.sequenceLength, 1));
                        
            % check if gains struct is properly updated
            this.verifyTrue(isstruct(this.controllerUnderTest.gains));
            this.verifySize(this.controllerUnderTest.gains, [1 2]);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'A'));
            this.verifyEqual(this.controllerUnderTest.gains(1).A, this.A);
            this.verifyEqual(this.controllerUnderTest.gains(2).A, newA); 
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'B'));
            this.verifyEqual(this.controllerUnderTest.gains(1).B, this.B);
            this.verifyEqual(this.controllerUnderTest.gains(2).B, newB);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'L'));
            this.verifyEqual(this.controllerUnderTest.gains(1).L, L_old);
            this.verifyEqual(this.controllerUnderTest.gains(2).L, L_new); % changed
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'delta'));
            this.verifyEqual(this.controllerUnderTest.gains(1).delta, this.controllerDelta);
            this.verifyEqual(this.controllerUnderTest.gains(2).delta, this.controllerDelta); % did not change   
        end
        
        %% testChangeModelParametersNewWIgnored
        function testChangeModelParametersNewWIgnored(this)
            this.assertTrue(isa(this.controllerUnderTest, 'ModelParamsChangeable'));
            
            L_old = this.controllerUnderTest.getControllerGain();           
            
            % pass a new noise covariance, has no effect
            newW = eye(this.dimX);
                        
            this.controllerUnderTest.changeModelParameters(this.A, this.B, newW);
            % this call triggers the recomputation of the gain
            % but existing gain should be reused
            this.controllerUnderTest.computeControlSequence(Gaussian(zeros(this.dimX, 1), eye(this.dimX)));            
            L_new = this.controllerUnderTest.getControllerGain();
                        
            this.verifyEqual(L_new, L_old);
            
            % accordingly, gains struct should be the same as before
            this.verifyTrue(isstruct(this.controllerUnderTest.gains));
            this.verifySize(this.controllerUnderTest.gains, [1 1]);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'A'));
            this.verifyEqual(this.controllerUnderTest.gains(1).A, this.A);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'B'));
            this.verifyEqual(this.controllerUnderTest.gains(1).B, this.B);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'L'));
            this.verifyEqual(this.controllerUnderTest.gains(1).L, L_old);
            this.verifyTrue(isfield(this.controllerUnderTest.gains, 'delta'));
            this.verifyEqual(this.controllerUnderTest.gains(1).delta, this.controllerDelta);
            
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
            controller = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, false, useLmiLab);
            
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
            controller = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, false, false, useLmiLab);
            
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
            controller = MSSController(this.A, this.B, this.sequenceLength, this.controllerDelta, false, false, useLmiLab);
            
            state = Gaussian(this.initialPlantState, eye(this.dimX));
            actualSequence = controller.computeControlSequence(state);
            % we expect the whole sequence as a stacked column vector
            expectedSize = [this.dimU * this.sequenceLength 1];            
            
            this.verifySize(actualSequence, expectedSize);            
        end
    end
end

