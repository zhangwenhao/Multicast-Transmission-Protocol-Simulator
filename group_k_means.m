function [ idx, sortedCtr, sortedWdx] = group_k_means( data, warden, r_w, k)
    % Perform k-means clustering with guards around Willie
    if nargin == 3
        k = 0;
    end

    loop = true;
    while loop
        k= k + 1;
        loop = false;
        
        % Perform k-means clustering
        [idx, ctr, wdx] = k_means(data, warden, r_w, k, 1000);
        
        % Check if any cluster center is within the guard zone
        c_s = size(ctr, 1);
        for i=1:c_s
            distances_g_w = norm(ctr(i, 1:2) - warden);
            if distances_g_w < r_w+ctr(i,3)
                loop=true;
                break;
            end
        end 
    end

    % Sort GUs within the guard zone by distance to Willie
    sortedWdx=[];
    if size(wdx,1)>0
        wdxDistances = sqrt((wdx(:, 1) - warden(1)).^2 + (wdx(:, 2) - warden(2)).^2);
        [~, sortedWdxIndices] = sort(wdxDistances);
        sortedWdx = wdx(sortedWdxIndices, :);
    end

    % Sort GUs out of the guard zone by distance to Willie
    sortedCtr=[];
    if size(ctr,1)>0
        ctrDistances = sqrt((ctr(:, 1) - warden(1)).^2 + (ctr(:, 2) - warden(2)).^2);
        [~, sortedCtrIndices] = sort(ctrDistances);
        sortedCtr = ctr(sortedCtrIndices, :);
    end    
end



function [ idx, ctr, wdx] = k_means( data, warden, r_w, k, iterations )

idx=[];
ctr=[];
wdx=[];
% First step: Remove data points around Willie
distances = sqrt(sum((data - warden).^2, 2));

% Find the data points around Willie, and add a column of radius 0 (to facilitate subsequent calculations)
wdx = data(distances <= r_w, :);
wdx = [wdx,zeros(size(wdx,1),1)];
data = data(distances > r_w, :);


[m, n] = size(data);

    if m>0
        ctr = [];
        ctr(1, :) = data(randi(m), :);
        
        if nargin == 4
            iterations = 0;
        end
        
        % Choose initial cluster centers using k-means++ algorithm
        for i = 2 : k
            % Calculate the shortest distance from each data point to the selected cluster center
            min_distances = min(pdist2(data, ctr(1:i-1, :)), [], 2);
            if min_distances==0
                break;
            end
            probabilities = min_distances.^2 / sum(min_distances.^2);
            next_center_idx = randsample(size(data, 1), 1, true, probabilities);
            ctr(i, :) = data(next_center_idx, :);
        end 
        
        % Add an extra column for radius
        ctr = [ctr, ones(size(ctr, 1), 1)];
        
        %n+1：One more dimension to store categories
        idx = zeros(m, n+1);
        idx(:,1:2) = data;
        
        %  Iterate to update cluster centers
        for iteration = 1:iterations
            % Calculate the distance of each data point to the cluster center
            distances = pdist2(data, ctr(:, 1:2));
            
            % Assign data points to the nearest cluster center
            [~, idx(:,3)] = min(distances, [], 2);
            
            % Update cluster centers
            for i = 1:k
               cluster_points= data(idx(:,3) == i, :);
                if ~isempty(cluster_points) 
                    [ctr(i, 1:2), ctr(i, 3)] = min_circular(cluster_points);
                end
            end
            
            % Check if iteration should stop
            if iteration > 1 && isequal(ctr, prev_ctr)
                break;
            end
            
            prev_ctr = ctr;
        end
    end
end