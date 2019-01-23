classdef InvertedPendulumTest < matlab.unittest.TestCase
    % Test cases for InvertedPendulum.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2017-2018  Florian Rosenthal <florian.rosenthal@kit.edu>
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
            this.inertia = 0.015;
            this.length = 0.3;
            
            this.samplingInterval = 0.01;
            
            this.pendulumUnderTest = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval);
            
            this.upwardEquilibrium = [0 0 pi 0]';
            this.downwardEquilibrium = [2 0 0 0]'; % different position though
        end
     end
    
    methods (Test)
        %% testInvertedPendulum
        function testInvertedPendulum(this)
            pendulum = InvertedPendulum(this.massCart, this.massPendulum, ...
                this.length, this.friction, this.samplingInterval);
            % check if moment of inertia is computed as expected (J = ml²/3)
            this.verifyEqual(pendulum.inertia, this.inertia);
        end
        %% testNonlinearDynamics
        function testNonlinearDynamics(this)
            
            actualState = this.pendulumUnderTest.simulate(this.upwardEquilibrium);
            this.verifyEqual(actualState, this.upwardEquilibrium, 'AbsTol', 1e-8);
            
            actualState = this.pendulumUnderTest.simulate(this.downwardEquilibrium);
            this.verifyEqual(actualState, this.downwardEquilibrium, 'AbsTol', 1e-8);
        end
    end
end

