function [sortedMBSLocations, finalRadius, sortedWdx] = group_uniform_radius(GULocations, warden, r_w, initialClusterRadius)
    finalRadius = 0;
    
    % Find GUs around Willie and add a column for radius (for later calculation)
    distances = sqrt(sum((GULocations - warden).^2, 2));
    wdx = GULocations(distances <= r_w, :);
    wdx = [wdx, zeros(size(wdx, 1), 1)];
    GULocations = GULocations(distances > r_w, :);

    % Initialize cluster centers
    MBSLocations = [];
   
    if size(GULocations, 1) > 0
        % Iterate until Willie's conditions are satisfied
        while true
            % Cluster GUs
            MBSLocations = spiralMBSPlacement(GULocations, initialClusterRadius);

            % Check if Willie is in any group
            WInAnyGroup = checkOverLapped(MBSLocations, r_w, warden);

            % If groups and guard zone do not overlap
            if ~WInAnyGroup
                break;
            end

            % Decrease clustering radius
            initialClusterRadius = initialClusterRadius - 10;  % Adjust step size based on requirements
        end
        
        % Set final clustering radius
        finalRadius = initialClusterRadius;
    end
    
    % Sort GUs within the guard zone by distance to Willie
    sortedWdx = [];
    if size(wdx, 1) > 0
        wdxDistances = sqrt((wdx(:, 1) - warden(1)).^2 + (wdx(:, 2) - warden(2)).^2);
        [~, sortedWdxIndices] = sort(wdxDistances);
        sortedWdx = wdx(sortedWdxIndices, :);
    end
    
    % Sort GUs out of the guard zone by distance to Willie
    sortedMBSLocations = [];
    if size(MBSLocations, 1) > 0
        ctrDistances = sqrt((MBSLocations(:, 1) - warden(1)).^2 + (MBSLocations(:, 2) - warden(2)).^2);
        [~, sortedCtrIndices] = sort(ctrDistances);
        sortedMBSLocations = MBSLocations(sortedCtrIndices, :);
    end 
end

function [MBSLocations] = spiralMBSPlacement(GULocations, r)
    % Input: GULocations - Matrix of GU coordinates, each row is a GU's coordinates
    %        r - MBS coverage radius
    % Output: MBSLocations - Matrix of MBS coordinates and radius

    % Initialize set of uncovered GUs
    uncoveredGUs = GULocations;

    % Initialize MBS and coverage information
    MBSLocations = [];

    while ~isempty(uncoveredGUs)
        % Check if uncovered GUs quantity is greater than or equal to 3
        if size(unique(uncoveredGUs, 'rows'), 1) < 3
            % If uncovered GUs quantity is less than 3, connect these points directly
            MBSLocations = [MBSLocations; [unique(uncoveredGUs, 'rows'), zeros(size(unique(uncoveredGUs, 'rows'), 1), 1)]];
            uncoveredGUs = [];  % All GUs are covered
        else
            % Calculate convex hull of uncovered GUs
            convexHull = convhull(uncoveredGUs(:, 1), uncoveredGUs(:, 2));
            
            % Get GUs on the convex hull
            boundaryGUs = uncoveredGUs(convexHull, :);
            
            % Get GUs inside the convex hull
            innerGUs = setdiff(uncoveredGUs, boundaryGUs, 'rows');

            % Randomly choose a boundary GU as the starting point
            k0 = boundaryGUs(randi([1, size(boundaryGUs, 1)]),:);

            [u, Pprio] = LocalCover(k0, k0, setdiff(boundaryGUs, k0, 'rows'), r);

            [~, Pprio] = LocalCover(u, Pprio, innerGUs, r);

            % Add the current MBS location to the result set
            [center, radius] = min_circular(Pprio);
            temp_u = [center, radius];
            MBSLocations = [MBSLocations; temp_u];

            % Update the set of uncovered GUs
            uncoveredGUs = setdiff(uncoveredGUs, Pprio, 'rows');
        end
    end
end

function [u, Pprio] = LocalCover(u, Pprio, Psec, r)
    Presec = [];
    while ~isempty(Psec) && ~isequal(Presec, Psec)
        Presec = Psec;
        % Exclude GUs more than 2r away from any GU in Pprio
        distances = min(pdist2(Psec, Pprio), [], 2);
        Psec(distances > 2*r, :) = [];
        
        % Include/exclude GUs within distance r to u
        distances = min(pdist2(Psec, u), [], 2);
        Pprio = [Pprio; Psec(distances <= r, :)];
        Psec(distances <= r, :) = [];
      
        % Find GU with the shortest distance to u
        distances = min(pdist2(Psec, u), [], 2);
        [~, idx] = min(distances);
        k1 = Psec(idx, :);

        % Add or remove GU from Pprio or Psec
        [Pprio, Psec, u] = refineSets(u, k1, Pprio, Psec, r);
    end
end

function [Pprio, Psec, u] = refineSets(u, k1, Pprio, Psec, r)
    % Add or remove GU from Pprio or Psec
    % based on the minicircle
    [covered, center] = isCovered(Pprio, k1, r);
    if covered
        u = center;
        Pprio = [Pprio; k1];
        Psec = setdiff(Psec, k1, 'rows');
    end
end

function [covered, center] = isCovered(Pprio, k1, r)
    % Check if GU k can be covered by refining MBS location
    Pprio = [Pprio; k1];
    [center, distance] = min_circular(Pprio);
    covered = distance <= r;
end

function overLapped = checkOverLapped(clusterCenters, r_w, Willie)
    % Check if Willie is in any group
    overLapped = false;
    for i = 1:size(clusterCenters, 1)
        if norm(Willie - clusterCenters(i, 1:2)) <= r_w + clusterCenters(i, 3)
            overLapped = true;
            break;
        end
    end
end
