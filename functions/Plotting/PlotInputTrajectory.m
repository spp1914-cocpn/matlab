function PlotInputTrajectory(appliedInputs, numRun, lineWidth, fontSize)
    figure();
    hold on;
    
    set(gca, 'FontSize', fontSize);
    set(gca, 'Box', 'On');
    
    xlabel('Time step');
    ylabel(sprintf('Applied Input'));
    
    inputs = appliedInputs(1, :, numRun);
    
    plot(0:numel(inputs)-1, inputs, 'Color', 'b', 'LineWidth', lineWidth, ...
                         'DisplayName', 'Applied Inputs');
end

