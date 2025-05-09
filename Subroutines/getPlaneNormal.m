function normalVector = getPlaneNormal(p1, p2, p3)
    % Input: p1, p2, p3 - 3D points as 1x3 vectors
    % Output: normalVector - Normal vector of the plane defined by the points

    % Create vectors from the points
    v1 = p2 - p1;
    v2 = p3 - p1;

    % Calculate the normal vector using the cross product
    normalVector = cross(v1, v2);

    % Normalize the normal vector
    normalVector = normalVector / norm(normalVector);
end