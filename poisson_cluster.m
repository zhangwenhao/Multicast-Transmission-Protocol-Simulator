function [random_points, num_clusters] = poisson_cluster(N, num_clusters, clusteringProbability, radius)

    % Generate N random uniformly distributed points
    random_points = generate_points_in_circle(N, radius); % rand(N, 2) * 2000 - 1000; % Generate coordinates in the range [-1000, 1000]
    
    if num_clusters > 0
        % Select cluster centers randomly from the generated points
        cluster_centers = random_points(randperm(N, num_clusters), :);
        
        % Get inner GUs (Ground Units)
        innerGUs = setdiff(random_points, cluster_centers, 'rows');
        
        innerGUsNum = size(innerGUs, 1);
        
        % Scale the remaining points to the nearest cluster's coordinates
        for i = 1:innerGUsNum
            % If the point does not belong to the selected clusters
            if ~ismember(innerGUs(i, :), cluster_centers, 'rows')
                % Calculate the distance to the nearest cluster
                distances = sqrt((cluster_centers(:, 1) - innerGUs(i, 1)).^2 + (cluster_centers(:, 2) - innerGUs(i, 2)).^2);
                % Find the index of the nearest cluster
                [~, idx] = min(distances);
                % Scale the point's coordinates to the nearest cluster's coordinates
                try
                    innerGUs(i, :) = (1 - clusteringProbability) * innerGUs(i, :) + clusteringProbability * cluster_centers(idx, :);
                catch exception
                    disp(exception.message);
                end
            end
        end
        
        % Combine cluster centers and inner GUs
        random_points = [cluster_centers; innerGUs];
    end
end


function points = generate_points_in_circle(N, radius)
    % Generate N random points within a circle of given radius
    angles = 2 * pi * rand(N, 1);
    radii = radius * sqrt(rand(N, 1));
    points = [radii .* cos(angles), radii .* sin(angles)];
end