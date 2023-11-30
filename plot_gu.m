function plot_gu(warden, data, radius)
    figure(990);
    hold on;

    % Plot the region
    viscircles([0, 0], radius, 'color', 'k', 'linestyle', '-', 'LineWidth', 1);

    % Plot GUs
    for i = 1:size(data, 1)
        plot(data(i, 1), data(i, 2), 'color', 'k', 'Marker', 'o', 'Markerfacecolor', 'k', 'MarkerSize', 3, 'LineStyle', 'none', 'DisplayName', '   GU');
    end

    % Plot Willie's position
    plot(warden(:, 1), warden(:, 2), 'color', 'k', 'Marker', '^', 'Markerfacecolor', 'k', 'MarkerSize', 8, 'LineStyle', 'none', 'DisplayName', '   Willie');

    axis equal; % Maintain consistent axis ratios to draw perfect circles
    box on;
    % Reduce white space in exported images
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1), pos(2) - 3, 6, 6]);
    xlim([-1000 1000]);
    ylim([-1000 1000]);

    xlabel('x (m)');
    ylabel('y (m)');

    % Create a legend
    s = plot(NaN, NaN, 'DisplayName', ' ', 'color', 'none', 'Markerfacecolor', 'none');
    h = plot(NaN, NaN, 'DisplayName', '   Willie', 'color', 'none', 'Markerfacecolor', 'none', 'LineStyle', 'none');
    legend([s, h], 'Orientation', 'horizontal', 'Location', 'northoutside', 'FontSize', 12);

    hold off;
end
