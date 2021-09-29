function PlotFilterResults(filters, resultData, yLabel, lineWidth, fontSize)
    numFilters = filters.getNumFilters();
    figure();
    hold on;
    
    set(gca, 'FontSize', fontSize);
    set(gca, 'Box', 'On');
    
    xlabel('Time step');
    ylabel(yLabel);
    
    hPlots = zeros(numFilters, 1);
    
    for i = 1:numFilters
        % Select filter
        f = filters.get(i);
        color = f.getColor();
        name  = f.getName();
        data = resultData(i, :, :);
        try
            hPlots(i) = plot(1:numel(data), data, color{:}, 'LineWidth', lineWidth, ...
                         'DisplayName', name);
        catch
             hPlots(i) = plot(1:numel(data), data, 'Color', 'b', 'LineWidth', lineWidth, ...
                         'DisplayName', name);
        end
    end
        
    legend('show', 'Location', 'East');

end

