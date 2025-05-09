function clusterAndFit = getClusterFit(cluster)
    % Input: cluster: nx3 array of 3D coordinates
    % Outputs: clusterAndFit: struct with the fields
    %           - .I: Circumcenter of the cluster
    %           - .r: radius of the cluster' circumcircle
    %

    tic;
    disp("Fitting cluster");

    clusterAndFit = struct();    
    cluster(:,3) = cluster(:,3)-min(cluster(:,3));
    
    % Get circumcircle
    I = [min(cluster(:,1))+(max(cluster(:,1))-min(cluster(:,1)))/2, min(cluster(:,2))+(max(cluster(:,2))-min(cluster(:,2)))/2, 0];   
    r = max(getDistancesToPoint(cluster(:,1:2), I(1:2)));
    
    % Move cluster center to xy [0 0]
    cluster(:,1) = cluster(:,1)-I(1);
    cluster(:,2) = cluster(:,2)-I(2);
    
    % Populate struct
    clusterAndFit.cluster = cluster;
    clusterAndFit.I = I;
    clusterAndFit.r = r;
    
    % Fit 3D spline
    K0 = min(cluster(:,3)); % top
    K1 = max(cluster(:,3)); % cluster depth
    % K2 = r/3; % cluster width
    % K3 = 30; % slope
    K4 = 0; % Center X
    K5 = 0; % Center Y
 
    % Compute model
    % pitModel = @(K2, K3, x, y) K0 + K1 - K1 ./ (1 + exp(-((((x - K4).^2) - K2^2) ./ K3^2) - ((((y - K5).^2) - K2^2) ./ K3^2)));    
    % clusterAndFit.Fit = fit([cluster(:,1), cluster(:,2)], cluster(:,3), pitModel, 'StartPoint', [r/3 30], 'Algorithm', 'Levenberg-Marquardt');

    K2 = r/3; % cluster width
    K3 = 30; % slope
    clusterAndFit.Fit = @(x, y) K0 + K1 - K1 ./ (1 + exp(-((((x - K4).^2) - K2^2) ./ K3^2) - ((((y - K5).^2) - K2^2) ./ K3^2)));
    
    toc
end

function distances = getDistancesToPoint(matrix, point)
    distances = zeros(height(matrix), 1);
    for ii = 1:height(matrix)
        distances(ii) = abs(norm(matrix(ii,:)-point));
    end
end
