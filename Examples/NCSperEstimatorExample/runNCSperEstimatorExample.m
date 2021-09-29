clear;
rng(42); % Intialization of random generator

%% Initialize simulation meta data
simTime = 3000; % number of timesteps to be simulated
numRuns = 1;

flipIntervals = 1200;
distributionUpdateIntervals = 150;%simTime-1;

updaterNames = ["Histogram", "ConvexOptimization", "Random", "Super", "RKL"];

%parpool(2);
for updaterIdx = 1:length(updaterNames)
    updaterName = updaterNames(updaterIdx)
    folderName = sprintf(strcat('Examples/NCSperEstimatorExample/Results/%d_Runs/%s',...
            'Updater/'), numRuns, updaterName);

    if ~(7 == exist(folderName, 'dir'))
        mkdir (folderName);
    end

    compare = zeros(numRuns, 15);
    timeFix = zeros(numRuns, simTime+1);
    timeVar = zeros(numRuns, simTime+1);
    
    GT = zeros(numRuns, 5);
    Fixed = zeros(numRuns, 5);
    Variable = zeros(numRuns, 5);

    %% Intialize printer for output
    printer = SimulationInfoPrinter('Run', 1, numRuns);
    printer.turnOff;
    printer.printSimulationStart();

    %% Simulation
    for i = 1:numRuns
        %printer.printProgress(1, i);

        setUpRoutine = SetUpScenarios(simTime, updaterName);
        setUpUpdate = ChangeUpScenarios(simTime, 1, updaterName);

        if strcmp(updaterName, 'Histogram')
            [groundTruthResults, fixedResults, variableResults] = ...
                UpdatedDistNCSExample.simulate(simTime, flipIntervals, ...
                distributionUpdateIntervals);
        else
            variableResults = ...
                NCSTemplate.simulate(simTime, flipIntervals, 'variable', ...
                setUpRoutine, setUpUpdate);
        end

        timeVar(i,:) = variableResults.time;

        groundTruthResults = ...
                NCSTemplate.simulate(simTime, flipIntervals, 'groundtruth', ...
                setUpRoutine, setUpUpdate);

        fixedResults = ...
            NCSTemplate.simulate(simTime, flipIntervals, 'fixed', ...
            setUpRoutine, setUpUpdate);

        compareGT = calculateErrors(groundTruthResults.estimatedTrajectory, ...
            groundTruthResults.realTrajectory);
        compareFixed = calculateErrors(fixedResults.estimatedTrajectory, ...
            groundTruthResults.realTrajectory);
        compareVariable = calculateErrors(variableResults.estimatedTrajectory, ...
            groundTruthResults.realTrajectory);

        compare(i,:) = [compareGT 0 ...
            compareFixed 0 ...
            compareVariable 0]; 

        timeFix(i,:) = groundTruthResults.time;
        
        % TODO: Iterate over struct
        %lengthTransitionMatrix = length(groundTruthResults.lastTransitionMatrix);
        GT(i,:) = groundTruthResults.lastTransitionMatrix(5, :);
        Fixed(i,:) = fixedResults.lastTransitionMatrix(5, :);
        Variable(i,:) = variableResults.lastTransitionMatrix(5, :);
    end
    printer.printSimulationEnd();

    gt = GT
    var = Variable
    fix = Fixed
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename = sprintf(strcat('%s%d_step_flipsAfter%d_updateIntervals%d_%s_2intervals.csv'), ...
            folderName, simTime, flipIntervals, distributionUpdateIntervals, datetime('today'));

    labels = {'GTRMSEX1','GTRMSEX2','GTMAEX1','GTMAEX2', 'GTTimePerStep', ...
        'FixedRMSEX1','FixedRMSEX2','FixedMAEX1','FixedMAEX2', 'FixedTimePerStep', ...
        'VariableRMSEX1','VariableRMSEX2','VariableMAEX1','VariableMAEX2', 'VariableTimePerStep'};
    compareTable = array2table(compare, 'VariableNames', labels);
    writetable(compareTable, filename);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename = sprintf(strcat('%s%d_step_flipsAfter%d_updateIntervals%d_%s_2intervals_runtimevar.csv'), ...
            folderName, simTime, flipIntervals, distributionUpdateIntervals, datetime('today'));
    compareTable = array2table(timeVar);
    writetable(compareTable, filename);
    
    filename = sprintf(strcat('%s%d_step_flipsAfter%d_updateIntervals%d_%s_2intervals_runtimefix.csv'), ...
            folderName, simTime, flipIntervals, distributionUpdateIntervals, datetime('today'));
    compareTable = array2table(timeFix);
    writetable(compareTable, filename);
    
    writeTransitionMatrices('groundtruth', filename, GT, updaterName);
    %writeTransitionMatrices('fixed', filename, Fixed);
    writeTransitionMatrices('variable', filename, Variable, updaterName);

end

save('june10');

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
function writeTransitionMatrices(estimatortype, filename, matrix, estimator)
    % Writes a matrix to a file.
    %
    % Parameters:
    %   >> estimatorType: estimator type which created the matrix
    %   >> filename: name of comma seperated values which stores the
    %   transition matrix
    %   >> matrix: transition probability matrix (TPM) to be stored
    %   >> estimator: name of the actual algorithm which produced the TPM
    name = strsplit(filename, '.');
    matrixfilename = sprintf('%s_%s_%s.csv', name{1}, estimatortype, estimator);
    matrix = array2table(matrix);
    writetable(matrix, matrixfilename);
end
