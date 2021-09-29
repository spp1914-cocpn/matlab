function plotResults()
    numRuns = 100;
    simTime = 10000;
    flipIntervals = 100000;
    updateIntervals = 150;
    date = '26-Mar-2019';
    filename = sprintf('Examples/UpdatedExample/Results/%d_Runs/%s', ...
        numRuns, sprintf('%sUpdater/%d_step_flipsAfter%d_updateIntervals%d_%s_%dintervals.csv', ...
        'RKL', simTime, flipIntervals, updateIntervals, date, 3) );
    filename = sprintf('Examples/NCSPerEstimatorExample/Results/%d_Runs/%s', ...
        numRuns, sprintf('%sUpdater/30000_step_flipsAfter1200_updateIntervals150_01-May-2019_2intervals.csv', ...
        'Random') );
    %filename = '30000_step_flipsAfter1200_updateIntervals150_26-Apr-2019_2intervals_runtimevar.csv';
    T = readtable(filename);
    allIdx = 1:numRuns;

    mean = true;
    T = table2array(T);
    %T = T(1, :);
    
    error_id = 2;
    T(:,error_id+5);
    p = plot(allIdx, T(:,error_id), 'Color', 'k', ...
        'DisplayName', 'GT-Schätzer', 'LineWidth', 1.5);
    hold on
    if mean
        plot([1,numRuns], [sum(T(:,error_id)) / numRuns, sum(T(:,error_id)) / numRuns], 'Color', 'k', ...
            'DisplayName', 'GT Durchschnitt', 'LineWidth', 2);
        hold on
    end
    plot(allIdx, T(:,error_id+10), 'Color', [255, 128, 0] / 255, ...
        'LineStyle', '--', ... 
        'DisplayName', 'Variable-Schätzer','LineWidth', 1.5);
    hold on
    if mean
        plot([1,numRuns], [sum(T(:,error_id+10)) / numRuns, sum(T(:,error_id+10)) / numRuns], ...
            'Color', [255, 128, 0] / 255, ...
            'LineStyle', '--', ... 
            'DisplayName', 'Variable Durchschnitt', 'LineWidth', 2);
    hold on
    end
    plot(allIdx, T(:,error_id+5), 'Color', [0, 128, 255] / 255, ...
        'LineStyle', '-.', ...
        'DisplayName', 'Fixed-Schätzer','LineWidth', 1.5);
    if mean
        hold on
        plot([1,numRuns], [sum(T(:,error_id+5)) / numRuns, sum(T(:,error_id+5)) / numRuns], ...
            'Color', [0, 128, 255] / 255, ...
            'LineStyle', '-.', ... 
            'DisplayName', 'Fixed Durchschnitt', 'LineWidth', 2);
    end
    grid on

    xlabel('Durchgang');
    
    if error_id <= 2
        ylabel('RMSE');
    elseif error_id == 2 || error_id == 3
        ylabel('MAE');
    else
        ylabel('Average Runtime per Step');
    end
    
    lgnd = legend('show', 'Location', 'West');
    axis([1 numRuns 6 11])
end