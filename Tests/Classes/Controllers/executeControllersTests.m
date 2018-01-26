function result = executeControllersTests()
    % Function to run all test cases in matlab/Tests/Classes/Controllers.
    
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
    
    tests = [
        TestSuite.fromClass(?SequenceBasedControllerTest) ...
        TestSuite.fromClass(?FiniteHorizonControllerTest) ...
        TestSuite.fromClass(?FiniteHorizonTrackingControllerTest) ...
        TestSuite.fromClass(?InfiniteHorizonControllerTest) ...
        TestSuite.fromClass(?NominalPredictiveControllerTest)];

    result = tests.run();
end

