function [plane, N, d] = getPlane(T)
% Input: T is a 3x3 array with [x y z] coordinates in each row
% Prerequisites: Symbolic Math Toolbox

    % Verify Input
    if width(T) ~= 3 || height(T) ~= 3
        error("T must be a 3x3 matrix!");
    end
    
    % Calculate Plane
    V1 = T(1,:)'; V2 = T(2,:)'; V3 = T(3,:)';
    N = cross(V2-V1, V3-V1);
    d = dot(-N, V1);
    
    syms x [3 1]
    plane = dot(N,x) + d == 0;

end
