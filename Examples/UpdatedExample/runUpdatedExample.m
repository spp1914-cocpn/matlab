clear; % clears matlab variable cache so no old values are used
rng(42); % Intialization of random generator

%% Initialize simulation meta data
simTime = 1000; % number of timesteps to be simulated
numRuns = 1;

flipIntervals = 200;
distributionUpdateIntervals = 150;%simTime-1;

updaterName = 'Histogram';
    
folderName = sprintf(strcat('Examples/UpdatedExample/Results/%d_Runs/%s',...
        'Updater/'), numRuns, updaterName);
    
if ~(7 == exist(folderName, 'dir'))
    mkdir (folderName);
end

compare = zeros(numRuns, 15);

transitionMatrices = [];
transitionMatrices.GT = [];
transitionMatrices.Fixed = [];
transitionMatrices.Variable = [];

%% Intialize printer for output
printer = SimulationInfoPrinter('Run', 1, numRuns);
printer.printSimulationStart();

%% Simulation
for i = 1:numRuns
    printer.printProgress(1, i);

    if strcmp(updaterName, 'Histogram')
        [groundTruthResults, fixedResults, variableResults] = ...
            UpdatedDistNCSExample.simulate(simTime, flipIntervals, ...
            distributionUpdateIntervals);
    elseif strcmp(updaterName, 'ConvexOptimization')
         [groundTruthResults, fixedResults, variableResults] = ...
            UpdatedNCSExample.simulate(simTime, flipIntervals, updaterName);
    else
        [groundTruthResults, fixedResults, variableResults] = ...
            UpdatedNCSExample.simulate(simTime, flipIntervals, updaterName);
    end

    compareGT = calculateErrors(groundTruthResults.estimatedTrajectory, ...
        groundTruthResults.realTrajectory);
    compareFixed = calculateErrors(fixedResults.estimatedTrajectory, ...
        groundTruthResults.realTrajectory);
    compareVariable = calculateErrors(variableResults.estimatedTrajectory, ...
        groundTruthResults.realTrajectory);

    compare(i,:) = [compareGT (groundTruthResults.time / simTime) ...
        compareFixed (fixedResults.time / simTime) ...
        compareVariable (variableResults.time / simTime)]; 

    % TODO: Iterate over struct
    lengthTransitionMatrix = length(groundTruthResults.lastTransitionMatrix);
    transitionMatrices.GT = [transitionMatrices.GT; ...
        groundTruthResults.lastTransitionMatrix(lengthTransitionMatrix, :)];
    transitionMatrices.Fixed = [transitionMatrices.Fixed; ...
        fixedResults.lastTransitionMatrix(lengthTransitionMatrix, :)];
    transitionMatrices.Variable = [transitionMatrices.Variable; ...
        variableResults.lastTransitionMatrix(lengthTransitionMatrix, :)];
end
printer.printSimulationEnd();

gt = transitionMatrices.GT
var = transitionMatrices.Variable
filename = sprintf(strcat('%s%d_step_flipsAfter%d_updateIntervals%d_%s_3intervals.csv'), ...
        folderName, simTime, flipIntervals, distributionUpdateIntervals, datetime('today'));

labels = {'GTRMSEX1','GTRMSEX2','GTMAEX1','GTMAEX2', 'GTTimePerStep', ...
    'FixedRMSEX1','FixedRMSEX2','FixedMAEX1','FixedMAEX2', 'FixedTimePerStep', ...
    'VariableRMSEX1','VariableRMSEX2','VariableMAEX1','VariableMAEX2', 'VariableTimePerStep'};
compareTable = array2table(compare, 'VariableNames', labels);
writetable(compareTable, filename);

%writeTransitionMatrices('groundtruth', filename, transitionMatrices.GT);
%writeTransitionMatrices('fixed', filename, transitionMatrices.Fixed);
%writeTransitionMatrices('variable', filename, transitionMatrices.Variable);

%% Error calculation
function err = rimmse(vec1, vec2)
    % Calculates the root mean squared error of two vectors with the same
    % length.
    %
    % Parameters:
    %   >> vec1: first input vector
    %   >> vec2: second input vector
    % Returns:
    %   << err: root mean squared error of vec1 and vec2
    err = realsqrt(immse(vec1, vec2));
end

function err = maerr(vec1,vec2)
    % Calculates the mean absolute error of two vectors with the same
    % length.
    %
    % Parameters:
    %   >> vec1: first input vector
    %   >> vec2: second input vector
    % Returns:
    %   << err: mean absolute error of vec1 and vec2
    err = sum(abs(vec1 - vec2)) / length(vec1);
end

function errors = calculateErrors(vec1,vec2)
    % Creates an array of different errors calculated from two vectors.
    %
    % Parameters:
    %   >> vec1: first input vector
    %   >> vec2: second input vector
    % Returns:
    %   << errors: array of calculated root mean squared errors and mean
    %   absolute errors
    errors = [rimmse(vec1(1,:), vec2(1,:)), ...      % RMSE of state difference
                rimmse(vec1(2,:), vec2(2,:)), ...
                maerr(vec1(1,:), vec2(1,:)),  ...    % Mean Absolute Error of state difference
                maerr(vec1(2,:), vec2(2,:))]; ...
end

%%
function writeTransitionMatrices(estimatortype, filename, matrix)
    % Writes a matrix to a file.
    %
    % Parameters:
    %   >> estimatorType: estimator type which created the matrix
    %   >> filename: name of comma seperated values which stores the
    %   transition matrix
    %   >> matrix: transition probability matrix (TPM) to be stored
    %   >> estimator: name of the actual algorithm which produced the TPM
    name = strsplit(filename, '.');
    matrixfilename = sprintf('%s_%s.csv', name{1}, estimatortype);
    matrix = array2table(matrix);
    writetable(matrix, matrixfilename);
end
