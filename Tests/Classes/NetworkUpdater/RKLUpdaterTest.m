classdef RKLUpdaterTest < BaseUpdaterTest
    properties (Access = private, Constant)
        name = 'Tests for RKL-based Network Updater';
    end
    
    properties (Access = private)
        dimension = 1;
        %
        %
        currentProbabilityMatrixEstimate = [];
        
        %
        %
        currentModeWeights = [];
        
        % scale factor, needs to fullfill the following contraints (to be
        % implemented):
        %   infinite sum of every step size needs to be infinite
        %   infinite sum of every squared step size needs to be smaller
        %   than infinity
        %   => stepSize needs to decrease with every recursion step
        %       this decreasing needs to be slow for the algorithm to
        %       converge. 
        stepSize;
        
        %
        %
        currentMeasurementMean = [];
        
        %
        %
        currentMeasurementCovariance = [];
        
        %
        %
        stateEstimate = [];
        
        %
        %
        measurementMatrix = [];
        
        %
        %
        systemMatrix = [];
        
        %
        %
        stateCovariance = [];
        
        %
        %
        measurementCovariance = [];
        
        distribution = [];
        updater = [];
    end
    
    methods (Access = private)
        function updater = createStandardRKLUpdater(this)
            weights = (1/ 5) * ones(1, 5);
            config = BuildDoubleIntegratorConfig();
            
            updater = RKLUpdater(weights, this.stepSize, ...
                config.C, config.A, config.B, ...
                config.W, config.V);          
        end
    end
    
    methods (TestMethodSetup)
        function initProperties(this)
            this.dimension = 2;
            this.stepSize = 0.02;
            initProperties@BaseUpdaterTest(this);
        end
    end
    
    methods (Test)        
        
    end
end