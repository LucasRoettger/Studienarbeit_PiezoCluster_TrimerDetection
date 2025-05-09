%% Creates struct with dummy clusters to run the trimer detection script with

% Variables
nClusters = 3; % Enter a number of clusters, that you want to generate
szCluster = 100; % Enter a size value for the clusters
nTraces = 30; % Enter a number of traces per cluster you want to generate
jitter = 0.3; % Amount of jitter you want to add to the data


K0 = 0; % top
K1 = 120; % cluster depth
K2 = 10; % cluster width
K3 = 50; % slope
K4 = szCluster/2; % Center X
K5 = szCluster/2; % Center Y

% Functions
fitModel = @(x, y) K0 + K1 - K1 ./ (1 + exp(-((((x - K4).^2) - K2^2) ./ K3^2) - ((((y - K5).^2) - K2^2) ./ K3^2)));

% Code block ----------------------
fig = figure;
tly = tiledlayout(fig, "flow");

IndivClusters_OPEN = struct;
for i = 1:nClusters
 [gridX, gridY] = meshgrid(linspace(1, szCluster, round(sqrt(nTraces))), linspace(1, szCluster, round(sqrt(nTraces))));
 gridZ = fitModel(gridX, gridY);

 % Jitter data
 for ii = 1: numel(gridX)
     gridX(ii) = gridX(ii) + (rand * jitter * szCluster);
     gridY(ii) = gridY(ii) + (rand * jitter * szCluster);
     gridZ(ii) = gridZ(ii) + (rand * jitter * K1);
 end

 field = strcat("dummyCluster_", num2str(i));
 IndivClusters_OPEN.(field) = [gridX(:), gridY(:), gridZ(:)];

 nexttile(tly)
 scatter3(gridX(:), gridY(:), gridZ(:))
 axis equal
end

save("ExampleData\modelClusters.mat", "IndivClusters_OPEN");
