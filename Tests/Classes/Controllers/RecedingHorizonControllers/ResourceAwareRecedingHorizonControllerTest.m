classdef ResourceAwareRecedingHorizonControllerTest < matlab.unittest.TestCase
    % Test cases for ResourceAwareRecedingHorizonController.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (Constant, Access = private)
        absTol = 1e-8;
    end
    
    properties (Access = private)
        A;
        B;
        dimU;
        dimX;
        dimY;
        C;
        Q;
        R;
        alpha;
        
        W;
        V;
        
        x0;
        x0Cov;
        
        maxMeasDelay;
        horizonLength;
        sequenceLength;
        caDelayProbs;
        
        F;
        G;
        H;
        J;
        
        sendingCostFunc;
        neighborFunc;
        initScheduleFunc;
        startScheduleFunc;
        
        refTrajectory;
        
        controllerUnderTest;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            rng(42);
            % use (noise-free) stirred tank example (Example 6.15, p. 500-501) from
            %
            % Huibert Kwakernaak, and Raphael Sivan, 
            % Linear Optimal Control Systems,
            % Wiley-Interscience, New York, 1972.
            %            
            this.dimX = 2;
            this.dimU = 2;
            this.dimY = 2;
            
            this.A = diag([0.9512, 0.9048]);           
            this.B = [4.877 4.877; -1.1895 3.569];
            this.C = 2 * ones(this.dimY, this.dimX);
            
            this.Q = diag([0.01, 1]) * diag([50 0.02]) * diag([0.01, 1]); % R_3 in the book
            this.R = diag([1/3, 3]); % R_2 in the book
            this.alpha = 0.1;
            
            this.W = eye(this.dimX) * 0.01;
            this.V = eye(this.dimY) * 0.001;
            
            this.x0 = ones(this.dimX,1);
            this.x0Cov = eye(this.dimX) * 0.2;
            
            this.maxMeasDelay = 3;
            this.horizonLength = 15;
            
            this.sequenceLength = 3; % so we have four modes
            this.caDelayProbs = ones(1, 5) / 5;
            
            this.initAdditionalProperties();
            
            this.controllerUnderTest = ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc);
        end

    end
    
    methods (Access = private)
        
        %% initAdditionalProperties
        function initAdditionalProperties(this)
            [this.F, this.G, this.H, this.J] = Utility.createActuatorMatrices(this.sequenceLength, this.dimU);
            
            this.neighborFunc = @(schedule) RecedingHorizonUdpLikeControllerTest.standardGetNeighbor(schedule);
            this.sendingCostFunc = @(eventSchedule, eventHistory) sum(eventSchedule) * 0;
            this.initScheduleFunc = @(horizon) ones(1, horizon);
            this.startScheduleFunc = @(lastSchedule) [lastSchedule(2:end), randi(2) - 1];
            
            this.refTrajectory = ones(this.dimX, 2 * this.horizonLength);
        end        
        
        %% verifyEqualWithAbsTol
        function verifyEqualWithAbsTol(this, actual, expected)
            this.verifyEqual(actual, expected, 'AbsTol', ResourceAwareRecedingHorizonControllerTest.absTol);
        end

    end

    methods(Test)
        
        %% testResourceAwareRecedingHorizonControllerInvalidNumArgs
        function testResourceAwareRecedingHorizonControllerInvalidNumArgs(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidNumberOfArguments';
            
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidSystemMatrix
        function testResourceAwareRecedingHorizonControllerInvalidSystemMatrix(this)
             expectedErrId = 'Validator:ValidateSystemMatrix:InvalidMatrix';
             
             invalidSysMatrix = eye(this.dimX, this.dimX + 1); % not square
             this.verifyError(@() ResourceAwareRecedingHorizonController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
             
             invalidSysMatrix = eye(this.dimX, this.dimX); % square but not finite
             invalidSysMatrix(1, end) = inf;
             this.verifyError(@() ResourceAwareRecedingHorizonController(invalidSysMatrix, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidInputMatrix
        function testResourceAwareRecedingHorizonControllerInvalidInputMatrix(this)
            expectedErrId = 'Validator:ValidateInputMatrix:InvalidInputMatrix';
            
            invalidInputMatrix = eye(this.dimX +1, this.dimU); % invalid dims
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
             
            invalidInputMatrix = eye(this.dimX, this.dimU); % correct dims, but not finite
            invalidInputMatrix(1, end) = nan;
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, invalidInputMatrix, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidMeasurementMatrix
        function testResourceAwareRecedingHorizonControllerInvalidMeasurementMatrix(this)
            expectedErrId = 'Validator:ValidateMeasurementMatrix:InvalidMeasMatrix';
            
            invalidMeasMatrix = eye(this.dimY, this.dimX + 1); % invalid dims
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
             
            invalidMeasMatrix = eye(this.dimY, this.dimX); % correct dims, but not finite
            invalidMeasMatrix(1, end) = nan;
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, invalidMeasMatrix, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidCostMatrices
        function testResourceAwareRecedingHorizonControllerInvalidCostMatrices(this)
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrix';
  
            invalidQ = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidQ = eye(this.dimX + 1); % matrix is square, but of wrong dimension
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidQ = eye(this.dimX); % correct dims, but inf
            invalidQ(end, end) = inf;
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidQMatrixPSD';
            invalidQ = eye(this.dimX); % Q is not symmetric
            invalidQ(1, end) = 1;
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidQ = -eye(this.dimX); % Q is not psd
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, invalidQ, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            % now test for the R matrix
            expectedErrId = 'Validator:ValidateCostMatrices:InvalidRMatrix';
            
            invalidR = eye(this.dimU + 1, this.dimU); % not square
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidR = eye(this.dimU); % correct dims, but inf
            invalidR(1,1) = inf;
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidR = ones(this.dimU); % R is not pd
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, invalidR, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidCaDelayProbs
        function testResourceAwareRecedingHorizonControllerInvalidCaDelayProbs(this)
            expectedErrId = 'Validator:ValidateDiscreteProbabilityDistribution:InvalidProbs';
            
            invalidDelayProbs = [-0.1 0.1 0.8 0.2]; % negative entry
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, invalidDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidDelayProbs = [inf 0.1 0.8 0.2];% inf entry
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, invalidDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
                     
            invalidDelayProbs = [0.06 0.05 0.8 0.1];% does not sum up to 1
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, invalidDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidHorizonLegth
        function testResourceAwareRecedingHorizonControllerInvalidHorizonLegth(this)
            expectedErrId = 'Validator:ValidateHorizonLength:InvalidHorizonLength';
            
            invalidHorizonLength = [1 2]; % not a scalar
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, invalidHorizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidHorizonLength = -1; % negative scalar
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, invalidHorizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidHorizonLength = 0; 
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, invalidHorizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidHorizonLength = 1.5; % not an integer
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, invalidHorizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidHorizonLength = inf; % not finite
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, invalidHorizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidMaxMeasDelay
        function testResourceAwareRecedingHorizonControllerInvalidMaxMeasDelay(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidMaxMeasDelay';
            
            invalidMaxMeasDelay = [1 2]; % not a scalar
            this.verifyError(@()  ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, invalidMaxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidMaxMeasDelay = -1; % negative scalar
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, invalidMaxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidMaxMeasDelay = 1.5; % not an integer
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, invalidMaxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidMaxMeasDelay = inf; % not finite
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, invalidMaxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidSysNoiseCovariance
        function testResourceAwareRecedingHorizonControllerInvalidSysNoiseCovariance(this)
            expectedErrId = 'Validator:ValidateSysNoiseCovarianceMatrix:InvalidCovDim';
            
            invalidW = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                invalidW, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
                       
            invalidW = ones(this.dimU); % W is not pd
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                invalidW, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidMeasNoiseCov
        function testResourceAwareRecedingHorizonControllerInvalidMeasNoiseCov(this)
            expectedErrId = 'Validator:ValidateMeasNoiseCovarianceMatrix:InvalidCovDim';
            
           invalidV = eye(this.dimX + 1, this.dimX); % not square
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, invalidV, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
                       
            invalidV = ones(this.dimY); % W is not pd
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, invalidV, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidX0
        function testResourceAwareRecedingHorizonControllerInvalidX0(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidX0';
            
            invalidX0 = ones(this.dimX +1); % not a vector
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, invalidX0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidX0 = ones(this.dimX +1, 1); % wrong dimensions
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, invalidX0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidX0 = ones(this.dimX, 1); % not finite
            invalidX0(1) = nan;
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, invalidX0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidX0Cov
        function testResourceAwareRecedingHorizonControllerInvalidX0Cov(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidX0Cov';
            
            invalidX0Cov = ones(this.dimX +1, 1); % not a matrix
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, invalidX0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidX0Cov = ones(this.dimX +1, this.dimX); % wrong dimensions
             this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, invalidX0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidX0Cov = ones(this.dimX); % not pd
             this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, invalidX0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidSendingCostFun
        function testResourceAwareRecedingHorizonControllerInvalidSendingCostFun(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidSendingCostFunction';
            
            invalidFun = this; % not a function
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                invalidFun, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidFun = @(a,b,c) a; % too many arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                invalidFun, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidFun = @(a) a; % too few arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                invalidFun, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidSendingCostFun
        function testResourceAwareRecedingHorizonControllerInvalidNeighborFun(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidNeighborFunction';
            
            invalidFun = this; % not a function
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, invalidFun, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidFun = @(a,b,c) a; % too many arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, invalidFun, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
            
            invalidFun = @() this; % too few arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, invalidFun, ...
                this.initScheduleFunc, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidInitScheduleFun
        function testResourceAwareRecedingHorizonControllerInvalidInitScheduleFun(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidInitScheduleFunction';
            
            invalidFun = this; % not a function
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                invalidFun, this.startScheduleFunc), expectedErrId);
            
            invalidFun = @(a,b,c) a; % too many arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                invalidFun, this.startScheduleFunc), expectedErrId);
            
            invalidFun = @() this; % too few arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                invalidFun, this.startScheduleFunc), expectedErrId);
        end
        
        %% testResourceAwareRecedingHorizonControllerInvalidStartScheduleFun
        function testResourceAwareRecedingHorizonControllerInvalidStartScheduleFun(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:InvalidStartScheduleFunction';
            
            invalidFun = this; % not a function
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, invalidFun), expectedErrId);
            
            invalidFun = @(a,b,c) a; % too many arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, invalidFun), expectedErrId);
            
            invalidFun = @() this; % too few arguments
            this.verifyError(@() ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, this.horizonLength, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, invalidFun), expectedErrId);
        end
%%
%%
        %% testResourceAwareRecedingHorizonController
        function testResourceAwareRecedingHorizonController(this)
            % check the side effects (creation of F and G) and of matrices K            
            dimEta = size(this.F, 1);
            dimInputSeq = size(this.G, 2);
            
            horizon = 1;
            expectedF = eye(dimEta);
            expectedG = zeros(size(this.G));
            controller = ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, horizon, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc);
            
            this.verifyEqual(controller.rollF, expectedF);
            this.verifyEqual(controller.rollG, expectedG);
            
            % now check K
            actualKs = controller.K;
            this.verifyEqual(size(actualKs, 3), 1);
            expectedK = this.B' * this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 1), expectedK);
            
            horizon = 3;
            expectedF = [expectedF; this.F; this.F^2];       
            expectedG = [zeros(dimEta, dimInputSeq*3); this.G, zeros(dimEta, dimInputSeq*2);this.F*this.G, this.G, zeros(dimEta, dimInputSeq)];
            
            controller = ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, horizon, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc);
            
            this.verifyEqual(controller.rollF, expectedF);
            this.verifyEqual(controller.rollG, expectedG);
            
            % now check K           
            actualKs = controller.K;
            this.verifyEqual(size(actualKs, 3), 3);
            expectedK = this.B' * this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 3), expectedK);
            expectedK = this.B' * this.A' * this.A* this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 2), expectedK);            
            expectedK = this.B' * this.A' * this.A * this.A' * this.A * this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 1), expectedK);
      
            
            horizon = 5;
            expectedF = [expectedF; this.F^3; this.F^4];       
            mZero = zeros(size(this.G));
            expectedG = [mZero, mZero, mZero, mZero, mZero; 
                this.G, mZero, mZero, mZero, mZero; 
                this.F * this.G, this.G, mZero, mZero, mZero; 
                this.F * this.F * this.G, this.F * this.G, this.G, mZero, mZero;
                this.F * this.F * this.F*this.G, this.F * this.F * this.G, this.F * this.G, this.G, mZero
                ];
            controller = ResourceAwareRecedingHorizonController(this.A, this.B, this.C, this.Q, this.R, ...
                this.alpha, this.caDelayProbs, this.sequenceLength, horizon, this.maxMeasDelay, ...
                this.W, this.V, this.x0, this.x0Cov, ...
                this.sendingCostFunc, this.neighborFunc, ...
                this.initScheduleFunc, this.startScheduleFunc);
            
            this.verifyEqual(controller.rollF, expectedF);
            this.verifyEqual(controller.rollG, expectedG);
            
            % now check K           
            actualKs = controller.K;
            this.verifyEqual(size(actualKs, 3), 5);
            expectedK = this.B' * this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 5), expectedK);
            expectedK = this.B' * this.A' * this.A* this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 4), expectedK);            
            expectedK = this.B' * this.A' * this.A * this.A' * this.A * this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 3), expectedK);
            expectedK = this.B' * this.A' * this.A * this.A' * this.A * this.A' * this.A * this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 2), expectedK);
            expectedK = this.B' * this.A' * this.A * this.A' * this.A * this.A' * this.A * this.A' * this.A * this.B;
            this.verifyEqualWithAbsTol(actualKs(:, :, 1), expectedK);
                             
            % horizon is 15
            expectedF = [eye(dimEta); this.F; this.F^2; this.F^3; this.F^4; this.F^5; this.F^6; this.F^7; this.F^8;
                        this.F^9; this.F^10; this.F^11; this.F^12; this.F^13; this.F^14];               
            expectedG = [repmat(mZero, 1, 15);
                this.G, repmat(mZero, 1, 14); 
                this.F * this.G, this.G ,repmat(mZero, 1, 13); 
                this.F ^2 * this.G, this.F * this.G, this.G, repmat(mZero, 1, 12);
                this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 11);
                this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 10);
                this.F^5 * this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 9);
                this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 8);
                this.F^7 * this.G, this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 7);
                this.F^8 * this.G, this.F^7 * this.G, this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 6);
                this.F^9 * this.G,this.F^8 * this.G, this.F^7 * this.G, this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 5);
                this.F^10 * this.G,this.F^9 * this.G,this.F^8 * this.G, this.F^7 * this.G, this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 4);
                this.F^11 * this.G,this.F^10 * this.G,this.F^9 * this.G,this.F^8 * this.G, this.F^7 * this.G, this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 3);
                this.F^12 * this.G,this.F^11 * this.G,this.F^10 * this.G,this.F^9 * this.G,this.F^8 * this.G, this.F^7 * this.G, this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, repmat(mZero, 1, 2);
                this.F^13 * this.G, this.F^12 * this.G,this.F^11 * this.G,this.F^10 * this.G,this.F^9 * this.G,this.F^8 * this.G, this.F^7 * this.G, this.F^6 * this.G,this.F^5*this.G,this.F ^4 *this.G, this.F ^3 *this.G, this.F ^ 2 *this.G, this.F * this.G, this.G, mZero;
                ];
            this.verifyEqual(this.controllerUnderTest.rollF, expectedF);
            this.verifyEqual(this.controllerUnderTest.rollG, expectedG);
            
            actualKs = this.controllerUnderTest.K;
            
            for i = 0:this.horizonLength - 1
                expectedK = this.B' * (this.A^(this.horizonLength - 1 - i))' * this.A^(this.horizonLength - 1 - i) * this.B;
                this.verifyEqualWithAbsTol(actualKs(:, :, i + 1), expectedK);
            end
            
            % validate that history is initialized properly
            this.verifySize(this.controllerUnderTest.stateHistory, [1 this.maxMeasDelay + 1]);
            this.verifyEqual(this.controllerUnderTest.stateHistory{1}, Gaussian(this.x0, this.x0Cov));
            for j=2:this.maxMeasDelay + 1
                this.verifyEmpty(this.controllerUnderTest.stateHistory{j});
            end
            
            this.verifySize(this.controllerUnderTest.measurementHistory, [1 this.maxMeasDelay]);
            for j=1:this.maxMeasDelay
                this.verifyEmpty(this.controllerUnderTest.measurementHistory{j});
            end
            
            % possible inputs are all zero
            this.verifyEqual(this.controllerUnderTest.inputHistory, zeros(this.dimU, this.sequenceLength, this.maxMeasDelay));
            expectedInputProbs = [0 0 0 1; 
                                  0 0 0 1;
                                  0 0 0 1];
            this.verifyEqual(this.controllerUnderTest.inputProbsHistory, expectedInputProbs);
                              
            this.verifyEqual(this.controllerUnderTest.eventHistory, zeros(1, this.sequenceLength - 1));
                        
            expectedInitSchedule = this.initScheduleFunc(this.horizonLength);           
            this.verifyEqual(this.controllerUnderTest.lastSchedule, expectedInitSchedule);      
            
            truncatedProbs = [this.caDelayProbs(1:3)'; this.caDelayProbs(4) + this.caDelayProbs(5)];
            expectedDelayProbMat = [truncatedProbs truncatedProbs truncatedProbs];
            this.verifyEqualWithAbsTol(this.controllerUnderTest.delayProbs, expectedDelayProbMat);
            
            % check the solver options
            saOptions = this.controllerUnderTest.saOptions;
            this.verifyEqual(saOptions.InitialTemperature, this.horizonLength);
            this.verifyEqual(saOptions.MaxStallIterations, 10);
            this.verifyEqual(saOptions.ReannealInterval, 800);
            this.verifyEqual(saOptions.Display, 'off');
        end
%%
%%
        %% testSetControllerPlantStateInvalidState
        function testSetControllerPlantStateInvalidState(this)
            expectedErrId = 'ResourceAwareRecedingHorizonController:SetControllerPlantState:InvalidState';
            
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
        %% testDoStageCostsComputation
        function testDoStageCostsComputation(this)
            state = ones(this.dimX, 1);
            input = 1.5 * ones(this.dimU, 1);
            timestep = 1;
            
            expectedStageCosts = state' * this.Q * state + input' * this.R * input;
            actualStageCosts = this.controllerUnderTest.computeStageCosts(state, input, timestep);
            
            this.verifyEqualWithAbsTol(actualStageCosts, expectedStageCosts);
        end
%%
%%
        %% testDoCostsComputationInvalidStateTrajectory
        function testDoCostsComputationInvalidStateTrajectory(this)
            trajectoryLength = 10;
            inputs = ones(this.dimU, trajectoryLength);
            expectedErrId = 'ResourceAwareRecedingHorizonController:DoCostsComputation:InvalidStateTrajectory';
            
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
            
            this.verifyEqualWithAbsTol(actualCosts, expectedCosts);
        end
    end
    
    methods (Access = private, Static)
        %% standardGetNeighbor
        function neighbor = standardGetNeighbor(schedule)
            % Standard getNeighbor implementation.
            % This function changes a random bit in the given schedule.
            %
            % Parameters:
            %   >> schedule (Binary vector)
            %      The current schedule.
            %
            % Returns:
            %   >> neighbor (Binary vector)
            %      The given schedule with a randomly changed bit.
            n = length(schedule);
            randomIndex = randi(n);
            neighbor = schedule;

            if neighbor(randomIndex) == 1
                neighbor(randomIndex) = 0;
            else
                neighbor(randomIndex) = 1;
            end
        end
    end
end

