function [center, radius] = min_circular(points)
    n = size(points, 1);

    % If the point set is empty, return invalid center and radius
    if isempty(points)
        center = NaN;
        radius = NaN;
        return;
    elseif numel(points) == 1
        % If there is only one point, return that point as the center with radius 0
        center = points(1, :);
        radius = 0;
        return;
    elseif size(points, 1) == 2
        % If there are two points, use their midpoint as the center and half of the distance as the radius
        center = mean(points, 1);
        radius = pdist2(points(1, :), points(2, :)) / 2;
        return;
    end

    % Initialize the center and radius
    center = points(1, :);
    radius = 0;

    % Iterate through points to update the minimum enclosing circle
    for i = 2:n
        if norm(center - points(i, :)) - radius > 0 % If the point is not inside the circle
            center = points(i, :);
            radius = 0;

            for j = 1:i-1
                if norm(center - points(j, :)) - radius > 0
                    center = (points(i, :) + points(j, :)) / 2.0;
                    radius = norm(center - points(i, :));

                    for k = 1:j-1
                        if norm(center - points(k, :)) - radius > 0
                            center = get_c(points(i, :), points(j, :), points(k, :));
                            radius = norm(center - points(i, :));
                        end
                    end
                end
            end
        end
    end
end

function c = get_c(a, b, c)
    p = (a + b) / 2;
    q = (a + c) / 2;
    v = rotate(b - a, pi/2.0);
    w = rotate(c - a, pi/2.0);

    if abs(cross(v, w)) < eps
        if abs(norm(a - b) + norm(b - c) - norm(a - c)) < eps
            c = (a + c) / 2;
        elseif abs(norm(b - a) + norm(a - c) - norm(b - c)) < eps
            c = (b + c) / 2;
        elseif abs(norm(a - c) + norm(c - b) - norm(a - b)) < eps
            c = (a + b) / 2;
        end
    else
        c = jiao(p, v, q, w);
    end
end

function result = rotate(a, t)
    result = [a(1) * cos(t) - a(2) * sin(t), a(1) * sin(t) + a(2) * cos(t)];
end

function result = cross(a, b)
    result = a(1) * b(2) - a(2) * b(1);
end

function result = jiao(p, v, q, w)
    u = p - q;
    t = cross(w, u) / cross(v, w);
    result = p + v * t;
end
