classdef InvertedPendulumTest < matlab.unittest.TestCase
    %INVERTEDPENDULUMTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        massCart; % mass of cart
        massPendulum; % mass of pendulum
        friction; % friction of the cart
        inertia; % moment of inertia of the pendulum
    
        length; % length of pendulum
        samplingInterval;
        
        pendulumUnderTest;
        
        upwardEquilibrium;
        downwardEquilibrium;
    end
    
     methods (TestMethodSetup)
        function initProperties(this)
            this.massCart = 0.5;
            this.massPendulum = 0.5;
            this.friction = 0.1;
            this.inertia = 0.006;
            this.length = 0.3;
            
            this.samplingInterval = 0.01;
            
            this.pendulumUnderTest = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.friction, this.inertia, this.length, this.samplingInterval);
            
            this.upwardEquilibrium = [0 0 pi 0]';
            this.downwardEquilibrium = [2 0 0 0]'; % different position though
        end
     end
    
    methods (Test)
        %% testNonlinearDynamics
        function testNonlinearDynamics(this)
            
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyEqual(actualState, this.upwardEquilibrium, 'AbsTol', 1e-8);
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyEqual(actualState, this.downwardEquilibrium, 'AbsTol', 1e-8);
        end
    end
end

