function [N, P, M] = getFitNormal(cluster, fitobj, point)

    % Input:    cluster - nx3 vector of 3D coordinates,
    %           fit - sfit-object
    %           point - 1x3 vector of 3D coordinates
    % Output:   N: Normal of the fitted surface closest to point,
    %           P: Point on the fitted surface that is closest to point
    
    numP = 256; % Increase for higher precision
    tic
    disp("Fetching fit normal vector")
    % Get closest point on fit (P)
    
    [gridX, gridY] = meshgrid(linspace(min(cluster(:,1)), max(cluster(:,1)), numP), linspace(min(cluster(:,2)), max(cluster(:,2)), numP));
    gridZ = fitobj(gridX, gridY);
    
    M = zeros(numP, numP, 3);
    M(:,:,1) = gridX;
    M(:,:,2) = gridY;
    M(:,:,3) = gridZ;
    
    fitPoints = [gridX(:), gridY(:), gridZ(:)];

    % surface(fitPoints(:,1), fitPoints(:,2), fitPoints(:,3), 10);
    
    [D, IDX] = pdist2(fitPoints, point, "euclidean", "Smallest", 3);

    closestPoints = [fitPoints(IDX(1), :); fitPoints(IDX(2), :); fitPoints(IDX(3), :)];

    P = closestPoints(1,:);

    % Get normal vector (V)
    N = getPlaneNormal(closestPoints(1,:), closestPoints(2,:), closestPoints(3,:)); %Three-point-method
    
    toc
end


