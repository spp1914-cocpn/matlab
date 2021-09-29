classdef ProjectionTest < matlab.unittest.TestCase
    properties (Constant, Access = private)
        name = "Projection Test";
        epsilon = 0.01;
        scale = 10;
    end
    
    properties (Access = private)
        delayDistribution;
        transitionMatrix;
        dimension;
        project;
    end
    
    methods (TestMethodSetup)
        %% initProperties
        function initProperties(this)
            this.delayDistribution = [0.0100 0.5000 0.4900];
            this.transitionMatrix = Utility.calculateDelayTransitionMatrix(...
                this.delayDistribution);
            this.dimension = size(this.delayDistribution, 2);
            this.project = @RKLUpdater.project;
        end
    end
    
    methods (Test)
        function testWithValidTransitionMatrix(this)
            projection = this.project(this.transitionMatrix);
            this.verifyEqual(projection, this.transitionMatrix);
        end
        
        function testWithValidButNotHessenbergTransitionMatrix(this)
            transitionMatrix = (1 / this.dimension) * ones(this.dimension);
            projection = this.project(transitionMatrix);
            Validator.validateTransitionMatrix(projection);
        end
        
        function testWithPartlyValidTransitionMatrix(this)
            transitionMatrix = this.transitionMatrix;
            transitionMatrix(1, :) = ones(1, this.dimension);
            projection = this.project(transitionMatrix);
            Validator.validateTransitionMatrix(projection);
        end
        
        function testWithSumGreaterOne(this)
            transitionMatrix = this.transitionMatrix ...
                + this.epsilon * ones(this.dimension);
            projection = this.project(transitionMatrix);
            Validator.validateTransitionMatrix(projection);
        end
        
        function testWithSumMuchGreaterOne(this)
            transitionMatrix = this.scale * this.transitionMatrix;
            projection = this.project(transitionMatrix);
            Validator.validateTransitionMatrix(projection);
        end
        
        function testWithSumSmallerOne(this)
            transitionMatrix = this.epsilon * this.transitionMatrix;
            projection = this.project(transitionMatrix);
            Validator.validateTransitionMatrix(projection);
        end
        
        function testWithValuesSmallerZero(this)
            transitionMatrix = this.transitionMatrix ...
                - eye(this.dimension);
            projection = this.project(transitionMatrix);
            Validator.validateTransitionMatrix(projection);
        end
    end
end