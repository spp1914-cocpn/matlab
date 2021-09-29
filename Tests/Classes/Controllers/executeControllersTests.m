function results = executeControllersTests(runInParallel)
    % Function to run all test cases in matlab/Tests/Classes/Controllers.
    %
    % Parameters:
    %   >> runInParallel (Flag, (i.e., a logical scalar), Optional)
    %      A flag to indicate whether the test cases shall be ran in
    %      parallel. If left out, the default value <true> is used.
    
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
    
    arguments
        runInParallel(1,1) logical = true;
    end
    
    import matlab.unittest.TestSuite;
    import matlab.unittest.TestRunner;
    
    tests = [
        TestSuite.fromClass(?SequenceBasedControllerTest) ...
        TestSuite.fromClass(?SequenceBasedTrackingControllerTest) ...
        TestSuite.fromClass(?FiniteHorizonControllerTest) ...
        TestSuite.fromClass(?FiniteHorizonTrackingControllerTest) ...
        TestSuite.fromClass(?InfiniteHorizonControllerTest) ...
        TestSuite.fromClass(?InfiniteHorizonUdpLikeControllerTest) ...
        TestSuite.fromClass(?RecedingHorizonUdpLikeControllerTest) ...
        TestSuite.fromClass(?IMMBasedRecedingHorizonControllerTest) ...
        TestSuite.fromClass(?ResourceAwareRecedingHorizonControllerTest) ...
        TestSuite.fromClass(?LinearlyConstrainedPredictiveControllerTest) ...
        TestSuite.fromClass(?ExpectedInputPredictiveControllerTest) ...
        TestSuite.fromClass(?NominalPredictiveControllerTest) ...
        TestSuite.fromClass(?PolePlacementPredictiveControllerTest) ...
        TestSuite.fromClass(?MSSControllerTest)];

    if runInParallel
        results = TestRunner.withTextOutput.runInParallel(tests);
    else
        results = TestRunner.withTextOutput.run(tests);
    end
end

