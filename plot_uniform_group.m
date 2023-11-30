function plot_uniform_group(warden, data, ctr, wdx, r_w, radius, uniformRadius)
    % Check if uniformRadius is not provided
    if nargin == 6
        uniformRadius = 0;
    end

    % Display the clustered results
    figure(992);
    hold on;

    % Plot the region
    viscircles([0, 0], radius, 'color', 'k', 'linestyle', '-', 'LineWidth', 1);

    m = size(data, 1);
    k = size(ctr, 1);
    w = size(wdx, 1);

    % Plot cluster centers and groups
    for i = 1:k
        viscircles(ctr(i, 1:2), ctr(i, 3), 'color', 'k', 'linestyle', '--', 'LineWidth', 1);
        text(ctr(i, 1), ctr(i, 2), num2str(i + w), 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    % Plot warden's position and guard zone
    plot(warden(:, 1), warden(:, 2), 'color', 'k', 'Marker', '^', 'Markerfacecolor', 'k', 'MarkerSize', 8);
    viscircles(warden, r_w, 'color', 'k', 'linestyle', '-', 'LineWidth', 2);

    % Plot GUs
    for i = 1:m
        plot(data(i, 1), data(i, 2), 'color', 'k', 'Marker', 'o', 'Markerfacecolor', 'k', 'MarkerSize', 3);
    end

    % Plot GUs around Willie
    for i = 1:w
        plot(wdx(i, 1), wdx(i, 2), 'color', 'k', 'Marker', 'o', 'Markerfacecolor', 'k', 'MarkerSize', 3);
        text(wdx(i, 1), wdx(i, 2), num2str(i), 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end

    axis equal; % Maintain consistent axis ratios to draw perfect circles
    box on;
    % Reduce white space in exported images
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) + 1, pos(2) - 2, 6, 6]);
    xlim([-1000 1000]);
    ylim([-1000 1000]);

    xlabel('x (m)');
    ylabel('y (m)');

    % Create a legend
    h = zeros(3, 1);
    h(1) = plot(NaN, NaN, 'o', 'DisplayName', ' MG', 'color', 'none', 'Markerfacecolor', 'none');
    h(2) = plot(NaN, NaN, 'DisplayName', '', 'color', 'none', 'Markerfacecolor', 'none');
    h(3) = plot(NaN, NaN, 'o ', 'DisplayName', ' Guard Zone', 'color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
    hLegend = legend(h, 'Orientation', 'horizontal', 'Location', 'northoutside', 'FontSize', 16);

    % Add a dashed circular annotation at the legend position
    annotation('ellipse', [hLegend.Position(1) + 0.03, hLegend.Position(2) + 0.01, 0.03, 0.03], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

    hold off;
end
