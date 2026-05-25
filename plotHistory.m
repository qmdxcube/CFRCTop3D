function plotHistory(compliance, volumefrac, contentfracHis)
    iter = 1:numel(compliance);
    figure('Color', 'w', 'Name', 'Iteration History'); 

    % --- Compliance ---
    yyaxis left;
    h1 = plot(iter, compliance, 'r-', 'LineWidth', 2);
    ylabel('Compliance');
    ax = gca; 
    ax.YAxis(1).Color = 'r'; 

    % --- Volume fractions ---
    yyaxis right;
    h2 = plot(iter, volumefrac, 'b-', 'LineWidth', 2); 
    hold on;
    h3 = plot(iter, contentfracHis, 'b--', 'LineWidth', 2); 
    ylabel('Volume Fractions');
    ylim([0 1]);
    ax.YAxis(2).Color = 'b'; 

    grid on; box on;
    xlabel('Optimization Step');
    set(gca, 'FontSize', 10);
    
    legend([h1, h2, h3], {'Compliance', 'Volume Fraction (Material)', 'Volume Fraction (Fiber)'}, ...
        'Location', 'northeast');
end