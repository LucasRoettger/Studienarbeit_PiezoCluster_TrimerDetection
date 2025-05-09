function intersect = checkTriangleIntersect3D(T1, T2)
% Input: T1 and T2 are 3x3 Matrices with [x y z] coordinates in each line
% Prerequisites: Sympolic Math Toolbox
    
    intersect = false;
    
    %% Get distances form T1 to pi2 and T2 to pi1
    Triangle = struct();
    Triangle(1).Vertices = T1;      Triangle(2).Vertices = T2;
    Triangle(1).dist = zeros(3,1);  Triangle(2).dist = zeros(3,1);
    Triangle(1).p_V = zeros(3,1);   Triangle(2).p_V = zeros(3,1);
    
    for i = 1:2
        k = 3-i;
        [Triangle(k).pi, Triangle(k).N, Triangle(k).d] = getPlane(Triangle(k).Vertices);
        
        for ii = 1:3
            Triangle(i).dist(ii) = dot(Triangle(k).N, Triangle(i).Vertices(ii,:)) + Triangle(k).d;
        end
        %%
        Triangle(i).S = sign(Triangle(i).dist);
        
        if all(Triangle(i).S > 0) || all(Triangle(i).S < 0) % All above or all below
            % reject
            % disp("Reject")
            return
        % elseif any(Triangle(i).S == 0) && ~all(Triangle(i).S == 0)
        %     %intersecting
        %     intersect = true;
        %     return
        elseif all(Triangle(i).dist == 0)
            % coplanar
            % disp("Co-planar")
            Triangle(i).coplanar = true;
        else
            % continue
            % disp("Not co-planar")
            Triangle(i).coplanar = false;
        end
    end
    
    if Triangle(1).coplanar && Triangle(2).coplanar
        % Coplanar continuation
        % Project both triangles onto an axis-aligned plane where they have the area surface
        % Find plane to project to
        N0 = Triangle(1).N/norm(Triangle(1).N); % Normalize N of one of the planes
        
        for i = 1:2
            Triangle(i).Vertices(:,find(N0 == max(N0), 1)) = []; % Delete coordinates in the direction in witch the normal vector is longest.
        end
    
        intersect = checkTriangleIntersect2D(Triangle(1).Vertices, Triangle(2).Vertices);
    
    elseif ~Triangle(1).coplanar && ~Triangle(2).coplanar
        % Non-coplanar continuation
        
        % Calculate directional vector D of intersection line L
        D = cross(Triangle(1).N, Triangle(2).N);
        
        for i = 1:2
            % Project each vertex onto L
            for ii = 1:3
                Triangle(i).p_V(ii) = Triangle(i).Vertices(ii,find(D == max(D),1));
            end
    
            % Calculate the t for the intersections of each triangle with L
            Triangle(i).t = zeros(2,1); 
            corners = find(Triangle(i).S == mode(Triangle(i).S));
            switch num2str(corners')
                case '1  3'
                    base = 2;
                case '1  2'
                    base = 3;
                case '2  3'
                    base = 1;
                otherwise
                    disp("CAVE! Edge case: One point is part of the intersection line!")
                    base = find(Triangle(i).S ~= 0, 1);
                    corners = [find(Triangle(i).S == 0) find(Triangle(i).S == sign(Triangle(i).dist(base))*-1)];
            end
    
            for ii = 1:2
                p_corner = Triangle(i).p_V(corners(ii));
                p_base = Triangle(i).p_V(base);
                d_corner = Triangle(i).dist(corners(ii));
                d_base = Triangle(i).dist(base);
                Triangle(i).t(ii) = p_corner+(p_base-p_corner)*d_corner/(d_corner-d_base);
            end
            Triangle(i).t = sort(Triangle(i).t);
        end
        
        if (Triangle(1).t(1) < Triangle(2).t(2)) && (Triangle(2).t(1) < Triangle(1).t(2))
            intersect = true;
        else
            intersect = false;
        end
    
        else
        error("Error: One triangle is coplanar, one is not!")
    end 
    %===========%
    
    
    function isIntersecting = checkTriangleIntersect2D(triangle1, triangle2)
        triangles{1} = triangle1; triangles{2} = triangle2;
        % triangle1 and triangle2 are 3x2 matrices where each row is a vertex [x, y]
        isIntersecting = false;
        for j = 1:2
            % Check for intersect using the Separating Axis Theorem\    
            for jj = 1:3
                % Get the edge of the triangle
                edge = triangles{j}(jj, :) - triangles{j}(mod(jj, 3) + 1, :);
                axis = [-edge(2), edge(1)]; % Perpendicular axis
        
                % Project both triangles onto the axis
                [min1, max1] = projectTriangle(triangles{1}, axis);
                [min2, max2] = projectTriangle(triangles{2}, axis);
        
                % Check for intersect
                if max1 < min2 || max2 < min1
                    return; % No intersect
                end
            end
        
            isIntersecting = true; % intersect detected
        end
    
        function [minProj, maxProj] = projectTriangle(triangle, axis)
            % Project the triangle onto the given axis
            projections = triangle * axis';
            minProj = min(projections);
            maxProj = max(projections);
        end
    end
%===============
end
