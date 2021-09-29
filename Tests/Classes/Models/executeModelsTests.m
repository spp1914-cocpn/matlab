function result = executeModelsTests()
    % Function to run all test cases in matlab/Tests/Classes/Models.
    
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

    import matlab.unittest.TestSuite;
    import matlab.unittest.TestRunner;
    
    tests = [
        TestSuite.fromClass(?JumpLinearSystemModelTest) ...
        TestSuite.fromClass(?LinearPlantTest) ...
        TestSuite.fromClass(?NonlinearPlantTest) ...
        TestSuite.fromClass(?InvertedPendulumTest) ...
        TestSuite.fromClass(?DoubleInvertedPendulumTest)
        ];
    
    result = TestRunner.withTextOutput.run(tests);
end

