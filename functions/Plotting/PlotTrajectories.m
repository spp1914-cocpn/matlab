function PlotTrajectories(filters, trueTrajectories, estimatedTrajectories, ...
    referenceTrajectory, numRun, yLabel, lineWidth, fontSize)
    figure();
    hold on;
    
    ax = gca;
    ax.FontSize = fontSize;
    ax.Box = 'On';
    set(gca, 'DefaultLineLineWidth', lineWidth);
        
    xlabel('Time step');
    ylabel(yLabel);
    
    numFilters = filters.getNumFilters();
    trueStates = trueTrajectories(:, numRun);
    estimatedStates = estimatedTrajectories(:, :, numRun);
    allIdx = 1:numel(referenceTrajectory);
    plot(allIdx, trueStates, 'Color', 'k', ...
        'DisplayName', 'True Trajectory');
    plot(allIdx, referenceTrajectory, 'Color', 'r', ...
        'DisplayName', 'Reference Trajectory');
    
    for i = 1:numFilters
        f = filters.get(i);
        color = f.getColor();
        name  = f.getName();
        
        try
            % Plot estimate
            plot(allIdx, estimatedStates(i, :), color{:},  ...
                 'DisplayName', name);
        catch
            plot(allIdx, estimatedStates(i, :), 'Color', 'b', ...
                'DisplayName', name);
        end
    end
    
    legend('show', 'Location', 'East');
end

