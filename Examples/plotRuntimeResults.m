function plotRuntimeResults()
    numRuns = 100;
    simTime = 30000;
    flipIntervals = 100000;
    updateIntervals = 150;
    date = '26-Mar-2019';
    filename = sprintf('Examples/NCSPerEstimatorExample/Results/%d_Runs/%s', ...
        numRuns, sprintf('%sUpdater/%d_step_flipsAfter%d_updateIntervals%d_%s_%dintervals.csv', ...
        'RKL', simTime, flipIntervals, updateIntervals, date, 3) );
    filename = sprintf('Examples/NCSPerEstimatorExample/Results/%d_Runs/%s', ...
        numRuns, sprintf('%sUpdater/30000_step_flipsAfter3000000_updateIntervals150_29-Apr-2019_2intervals_runtimevar.csv', ...
        'Super') );
    T = readtable(filename);
    n = 10;
    allIdx = 1:(simTime+1);
    allIdx = allIdx(1 : n : end);

    mean = false;
    T = table2array(T);
    T = T(1, :);
    
    T = T(1 : n : end);  % => 1 4 7 10
    %T = T(n : n : end);  % => 3 6 9
    size(T)
    sum(T) / (simTime+1)
    p = scatter(allIdx, T, 'filled', ...
        'MarkerEdgeColor',[0 .5 .5],...
        'MarkerFaceColor',[0 .7 .7],...
        'DisplayName', 'Ohne Updater' , 'LineWidth', 0.5);
    hold on
    if mean
        plot([1,simTime+1], [sum(T) / (simTime+1), sum(T) / (simTime+1)], 'Color',[.7 .7 0], ...
            'DisplayName', 'Mean runtime', 'LineWidth', 2);
        hold on
    end
    axis([0 (simTime+1) 0 0.05])
    grid on

    xlabel('Zeitschritt in der Simulation');
    
    %if error_id <= 1
    %    ylabel('RMSE');
    %elseif error_id == 2 || error_id == 3
    %    ylabel('MAE');
    %else
        ylabel('Laufzeit pro Schritt');
    %end
    
    lgnd = legend('show', 'Location', 'NorthWest');
end