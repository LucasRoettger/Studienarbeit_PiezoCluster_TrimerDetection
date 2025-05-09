%% TrimerDetectionInClusters
% Input:    Struct with clusters (Named: IndivClusters_OPEN)
% Prerequisites: Symbolic Math Toolbox, Curve Fitting Toolbox

clc, clear

% Variables
threshold = 40; % Threshold for maximum distance between vertices of potential trimers

NNexclusion = 0; % Set 0 if you want to exclude all trimers with another point nearer to the circumcenter than it's vertices. Set 1 if you want to include them
weightSOD = 3; % Sum of distances
weightMIA = 1; % Max inner angle
weightSVA = 1; % Surface-vector-angle

vP = [330, 45]; % view parameters
%% Step 1: Find all potential clusters based on tace distances

tStart = tic;
disp("Step 1")
% Load and Process data
[file, path] = uigetfile(".mat", "Select cluster data to load", "ExampleData/");
load(fullfile(path,file));

% Evaluate input
if ~exist("IndivClusters_OPEN", "var")
    error("Loaded file '" + file + "' does not contain a struct named IndivClusters_OPEN\nPlease rename your data or change reference in lines 33 and 41. If you do, delete lines 24 - 27 and please proceed with care");
end

% Generate random model data
rng("default"); 
% traceMeans.RandomModelData = rand(100,3)*200;

fields = fieldnames(IndivClusters_OPEN);
figWB = waitbar(0, "Finding all potential trimers");
clusterStruct = struct();

for i = 1:length(fields)
    waitbar(i/length(fields),figWB, "Finding all potential trimers");

    disp([num2str(i) "/" length(fields)]);
    stru = getClusterStruct(IndivClusters_OPEN.(fields{i}));

    % Find density factor (median of neighbor traces within threshold)
    
    D = pdist2(stru.allTraces(:,1:3), stru.allTraces(:,1:3), "euclidean", "Smallest", 100); % 100 is set to cut down on processing time
    
    densityFactor = max(sum(D<threshold,1));
    
    if densityFactor < 3
        disp("No potential trimers detected in cluster " + fields{i});
        continue
    end

    % Run pdist2 again
    [D, I] = pdist2(stru.allTraces(:,1:3), stru.allTraces(:,1:3), "euclidean", "Smallest", densityFactor);

    % Clean up all values > threshold

    for traceNR = 1: width(D)
        cutoff = find(D(:,traceNR)>threshold);
        D(cutoff, traceNR) = NaN;
        I(cutoff, traceNR) = NaN;
    end
    
    % Find potential Trimers

    potTrimers = cell(ceil(length(I)/3), 2);

    k = 0;
    for traceNR = 1:width(I) % Loop through every trace
        for lineNR = 2:length(I(:,traceNR)) % Loop through every neighbor (LineNR = 1 is the trace itself)
            neighborIDX = I(lineNR,traceNR);
            
            if isnan(neighborIDX)
                continue
            end

            C = intersect(I(:,traceNR), I(:,neighborIDX));
            
            if length(C) < 3 % Skip
                continue
            elseif length(C) == 3 % Add to potential Clusters
                if ~any(cellfun(@(x) isequal(x, C'), potTrimers)) % Check for duplicates
                    trimerIDX = char(strjoin([num2str(traceNR) "-" num2str(neighborIDX)], ''));
                    k                = k+1;
                    potTrimers{k, 1} = trimerIDX;
                    potTrimers{k, 2} = C';
                end
            elseif length(C) > 3 % Find all sub clusters an add to potential clusters
                trimerIDX = string([num2str(traceNR) "-" num2str(neighborIDX)]);
                potThirdTraceList = C(C~= traceNR & C~= neighborIDX);

                for subTrimerNR = 1:length(potThirdTraceList)                    
                    subTrimerIDX = char(strjoin([trimerIDX "-" num2str(potThirdTraceList(subTrimerNR))], ''));
                    subC = sort([traceNR, neighborIDX, potThirdTraceList(subTrimerNR)]);

                    if ~any(cellfun(@(x) isequal(x, subC), potTrimers)) % Check for duplicates
                        k                = k+1;
                        potTrimers{k, 1} = subTrimerIDX;
                        potTrimers{k, 2} = subC;
                    end
                end
            end
        end
    end 
    stru.potTrimers = cell2table(potTrimers, "VariableNames", ["Index" "Traces"]);
    
    if iscell(stru.potTrimers.Traces) && any(cell2mat(cellfun(@(x) isempty(x), stru.potTrimers.Traces, "UniformOutput",false)))
        disp("No potential trimers detected in cluster " + fields{i});
        continue
    else
        % Filter out all traces without trimers
        allIndices = 1:height(stru.allTraces);
    
        tracesWithoutTrimers = allIndices(~ismember(allIndices(:), stru.potTrimers.Traces));
        stru.tracesWithoutTrimer = stru.allTraces(tracesWithoutTrimers, :);
    
        % Filter out all traces within trimers
        tracesInTrimers = allIndices(ismember(allIndices(:), stru.potTrimers.Traces));
        stru.tracesInTrimers = stru.allTraces(tracesInTrimers, :);    
    end

    subFields = fieldnames(stru);
    for traceNR = 1:length(subFields)
        clusterStruct(i).Name = fields{i};
        clusterStruct(i).(subFields{traceNR}) = stru.(subFields{traceNR});
    end
end
close(figWB);
toc(tStart)

%% Step 2: Evaluation and selection of (potential) trimers

% Plotting options
plOptionNames = ["plotPotTrimers" "plotVerices" "plotFit" "plotVertexVectors" "plotTrimerCenters" "plotNormalVectors" "plotClusterCenter" "plotCenterVectors" "plotlikelihoods"];
plOptions = table('Size', [1 length(plOptionNames)], 'VariableTypes', repelem("logical",length(plOptionNames)), 'VariableNames', plOptionNames);
plOptions.plotPotTrimers = true;
plOptions.plotVertices = false;
plOptions.plotFit = true;
plOptions.plotVertexVectors = false;
plOptions.plotTrimerCenters = false;
plOptions.plotNormalVectors = false;
plOptions.plotSurfaceNormalVectors = false;
plOptions.plotlikelihoods = true;
clear plOptionNames

% Create waitbar figure
figWB = uifigure("Name", "Please wait trimer detection in progress.", "Position",[500 500 500 400]);
wbFields = uiprogressdlg(figWB, "Title","Processing clusters");
wbCluster = uiprogressdlg(figWB, "Title","Processing trimers", "ShowPercentage","on");

for i = 1:length(clusterStruct)
    Step2 = tic;
    wbFields.Message = "Processing cluster (" + num2str(i) + "/" + num2str(height(fields)) + ")";
    wbFields.Value = i/height(fields);

    disp("Processing cluster: " + fields{i});

    nPotTrimers = height(clusterStruct(i).potTrimers);
    clusterStruct(i).vecAB = {nPotTrimers, 1};
    clusterStruct(i).vecAC = {nPotTrimers, 1};
    clusterStruct(i).vecBC = {nPotTrimers, 1};
    clusterStruct(i).COG = {nPotTrimers, 1};
    clusterStruct(i).distSum = {nPotTrimers, 1};
    clusterStruct(i).maxAng = {nPotTrimers, 1};
    clusterStruct(i).normVector = {nPotTrimers, 1};
    clusterStruct(i).surfaceNormalVector = {nPotTrimers, 1};
    clusterStruct(i).surfaceVectorAngle = {nPotTrimers, 1};

    clustersAndFits(i) = getClusterFit(clusterStruct(i).allTraces(:,1:3));
    % Move cluster to origin
    clusterStruct(i).allTraces(:, 1) = clusterStruct(i).allTraces(:, 1) - clustersAndFits(i).I(1);
    clusterStruct(i).allTraces(:, 2) = clusterStruct(i).allTraces(:, 2) - clustersAndFits(i).I(2);
    clusterStruct(i).allTraces(:, 3) = clusterStruct(i).allTraces(:, 3) - min(clusterStruct(i).allTraces(:, 3));

    if plOptions.plotPotTrimers
        % Generate a figure for each cluster
        triFig = figure("Name", num2str(i), "Position", [400 500 1200 400]);
        if plOptions.plotFit
            triFigL = tiledlayout(1,3, "TileSpacing", "compact");
        else
            triFigL = tiledlayout(1,2, "TileSpacing", "compact");   
        end
        cTab = rand(nPotTrimers, 3);        
    end
    
    %% a) Generate qualitative values
        wbCluster.Title = "Processing potential trimers. (Step 1/3)" ;
    for ii = 1:nPotTrimers
        wbCluster.Message = "Generating qualitative values ("+ num2str(ii) + "/" + num2str(nPotTrimers)+")";
        wbCluster.Value = ii/nPotTrimers;

        % Get the points
        pointA = clusterStruct(i).allTraces(clusterStruct(i).potTrimers.Traces(ii,1),1:3);
        pointB = clusterStruct(i).allTraces(clusterStruct(i).potTrimers.Traces(ii,2),1:3);
        pointC = clusterStruct(i).allTraces(clusterStruct(i).potTrimers.Traces(ii,3),1:3);
        trimer = [pointA;pointB;pointC];

        clusterStruct(i).vecAB{ii} = pointB-pointA;
        clusterStruct(i).vecAC{ii} = pointC-pointA;
        clusterStruct(i).vecBC{ii} = pointC-pointB;
    
        % 1st: Total distance between all Points
        clusterStruct(i).distSum{ii} = norm(clusterStruct(i).vecAB{ii}) + norm(clusterStruct(i).vecAC{ii}) + norm(clusterStruct(i).vecBC{ii});
    
        % 2nd: Biggest angle of the trimer
        alfa = getVectorAngle(clusterStruct(i).vecAB{ii}, clusterStruct(i).vecAC{ii});
        beta = getVectorAngle(clusterStruct(i).vecBC{ii}, clusterStruct(i).vecAB{ii}*-1);
        gamma = getVectorAngle(clusterStruct(i).vecAC{ii}, clusterStruct(i).vecBC{ii});
    
        clusterStruct(i).maxAng{ii} = max([alfa beta gamma]);
    
        % 3rd: Vector to cluster center and angles to the normal vector
        clusterStruct(i).normVector{ii} = getPlaneNormal(pointA, pointB, pointC);
        
        
        trimerCOG = mean(trimer, 1); % Trimer's center of gravity
        
        [N,P,M] = getFitNormal(clusterStruct(i).allTraces(:,1:3), clustersAndFits(i).Fit, trimerCOG);

        clusterStruct(i).COG{ii} = trimerCOG;
        clusterStruct(i).surfaceNormalVector{ii} = N;        
        SNVangle = getVectorAngle(clusterStruct(i).normVector{ii}, clusterStruct(i).surfaceNormalVector{ii});
        
        if SNVangle>90 % Check if normal Vector points outward
            clusterStruct(i).normVector{ii} = clusterStruct(i).normVector{ii}*-1;
            clusterStruct(i).surfaceVectorAngle{ii} = getVectorAngle(clusterStruct(i).normVector{ii}, clusterStruct(i).surfaceNormalVector{ii});                
        else
            clusterStruct(i).surfaceVectorAngle{ii} = SNVangle;
        end
        
        % 4th: Distance from circumcenter to nearest point
        [~, trimerCC, radius] = triangle_circumcircle(pointA', pointB', pointC'); % Trimer's circumcenter
        CCdistances = zeros(height(clusterStruct(i).allTraces),1);

        for iii = 1:numel(CCdistances)
            CCdistances(iii) = norm(clusterStruct(i).allTraces(iii,1:3)-trimerCC');
        end
        
        clusterStruct(i).nearestNeighbour(ii,:) = [radius, min(CCdistances)];

        % Plot potential trimers     
             
        if plOptions.plotPotTrimers
            ax1 = nexttile(triFigL,1);
            hold on
            fill3(trimer(:,1), trimer(:,2), trimer(:,3), cTab(ii));
            xlabel("X (nm)"); ylabel("Y (nm)"); zlabel("Z (nm)");
            title("Potential Trimers", "Interpreter","none");
            colormap("parula");
            axis equal
            grid(ax1, "on");
            view(vP);
            if plOptions.plotVerices
                scatter3(trimer(:,1), trimer(:,2), trimer(:,3), 10, "cyan", "filled");            
                text(trimer(1,1), trimer(1,2), trimer(1,3), "A", "HorizontalAlignment","right");
                text(trimer(2,1), trimer(2,2), trimer(2,3), "B", "HorizontalAlignment","right");
                text(trimer(3,1), trimer(3,2), trimer(3,3), "C", "HorizontalAlignment","right");
            end
            
            if plOptions.plotVertexVectors
                quiver3(pointA(1), pointA(2), pointA(3), clusterStruct(i).vecAB{ii}(1), clusterStruct(i).vecAB{ii}(2), clusterStruct(i).vecAB{ii}(3), "off", "filled", "Color","k", "LineWidth",3);
                quiver3(pointA(1), pointA(2), pointA(3), clusterStruct(i).vecAC{ii}(1), clusterStruct(i).vecAC{ii}(2), clusterStruct(i).vecAC{ii}(3), "off", "filled", "Color","k", "LineWidth",3);
                quiver3(pointB(1), pointB(2), pointB(3), clusterStruct(i).vecBC{ii}(1), clusterStruct(i).vecBC{ii}(2), clusterStruct(i).vecBC{ii}(3), "off", "filled", "Color","k", "LineWidth",3);
        
            end
            if plOptions.plotTrimerCenters
                scatter3(trimerCOG(1), trimerCOG(2), trimerCOG(3), 40, 'red', 'filled');
            end
            if plOptions.plotNormalVectors
                magn = 5;
                quiver3(trimerCOG(1), trimerCOG(2), trimerCOG(3), clusterStruct(i).normVector{ii}(1)*magn, clusterStruct(i).normVector{ii}(2)*magn, clusterStruct(i).normVector{ii}(3)*magn, "off", "filled", "Color","g", "LineWidth",3);
            end
            if plOptions.plotSurfaceNormalVectors
                quiver3(trimerCOG(1), trimerCOG(2), trimerCOG(3), clusterStruct(i).surfaceNormalVector{ii}(1), clusterStruct(i).surfaceNormalVector{ii}(2), clusterStruct(i).surfaceNormalVector{ii}(3), "off", "filled", "Color","m", "LineWidth",1);
            end
            
            hold off
            Lims = [xlim; ylim; zlim];
          
        end        
    end

    if plOptions.plotFit            
        ax2 = nexttile(triFigL,2);
        scatter3(clusterStruct(i).allTraces(:,1), clusterStruct(i).allTraces(:,2), clusterStruct(i).allTraces(:,3), "filled", "red")
        surface(M(1,:,1), M(:,1,2), M(:,:,3), 'LineStyle', ':', 'EdgeAlpha', 0.8)

        xlabel("X (nm)"); ylabel("Y (nm)"); zlabel("Z (nm)");
        title("Cluster Fit", "Interpreter","none");
        colormap("parula");
        view(vP); axis equal
    end
    %% b) Calculate likelihood values to each trimer
    wbCluster.Title = "Processing potential trimers. (Step 2/3)" ;

    likeTable = cell(nPotTrimers, 5);

    for ii = 1:nPotTrimers
        wbCluster.Message = "Calculating likelihood estimations ("+ num2str(ii) + "/" + num2str(nPotTrimers)+")";
        wbCluster.Value = ii/nPotTrimers;

        disp("Processing potential trimer # " + num2str(ii));
        % 1: Total distance between all Points
         % Left border: 24  % Right border: 120 % Optimum: 72       
         likeTable{ii,2} = getGaussianLike(clusterStruct(i).distSum{ii}, 24, 120);
    
        % 2: Maximum internal angle angle
         % Optimum: 60°   % Worst: 178°
         pFactorIntAng = 1/(178-60);
         likeTable{ii,3} = 1-(clusterStruct(i).maxAng{ii}-60)*pFactorIntAng;
         
        % 3: Angle between fit surface and trimer
         % Optimum: 0°  % Worst: 90°
         pFactorCentAngle = 1/90;
         likeTable{ii,4} = 1-clusterStruct(i).surfaceVectorAngle{ii}*pFactorCentAngle;
        
        % 4: Nearest neighbor analysis
        if clusterStruct(i).nearestNeighbour(ii,1) > clusterStruct(i).nearestNeighbour(ii,2)
            likeTable{ii,5} = 1;
        else
            likeTable{ii,5} = NNexclusion; 
        end

        % 5: Total likelihood
        likeTable{ii,1} = (likeTable{ii,2}*weightSOD+likeTable{ii,3}*weightMIA+likeTable{ii,3}*weightSVA) / (weightSOD+weightSVA+weightMIA); % Additive
    end

    clusterStruct(i).distSum(2,:) = likeTable(:,2)';
    clusterStruct(i).maxAng(2,:) = likeTable(:,3)';
    clusterStruct(i).surfaceVectorAngle(2,:) = likeTable(:,4)';
    clusterStruct(i).potTrimers.likelihood = likeTable(:,1);
    
    % Plot likelihoods for evaluation
    if plOptions.plotlikelihoods
        if ~exist("likeFig", "var")
            likeFig = figure("Name", "likelihoods", "Position", [400 200 1200 600]);
            tL = tiledlayout(likeFig, 2, 4);
            
            axDistSum = nexttile(tL, 1);    
            axMaxAng = nexttile(tL, 2);
            axNormAngle = nexttile(tL, 3);

            histDistSum = nexttile(tL, 5);
            histMaxAngle = nexttile(tL, 6);
            histNormAngle = nexttile(tL, 7);

            allLikes = nexttile(tL, 4, [2,1]);
        end
        
        hold(axDistSum, "on")
        scatter(axDistSum, cell2mat(clusterStruct(i).distSum(1,:)), cell2mat(clusterStruct(i).distSum(2,:)), 10, [0 0.4470 0.7410] ,'filled', 'MarkerFaceAlpha', 0.8);
        title(axDistSum,"Sum of Distances");
        xlabel(axDistSum, "Distance (nm)"); ylabel(axDistSum, "likelihood"); axDistSum.YLim = [0 1]; axDistSum.XLim = [0 150];
        hold(axDistSum, "off")

        hold(axMaxAng, "on")
        scatter(axMaxAng, cell2mat(clusterStruct(i).maxAng(1,:)), cell2mat(clusterStruct(i).maxAng(2,:)), 10, [0.9290 0.6940 0.1250] , 'filled', 'MarkerFaceAlpha', 0.8);
        title(axMaxAng, "Maximum internal angle");
        xlabel(axMaxAng, "Angle (deg)"); ylabel(axMaxAng, "likelihood"); axMaxAng.YLim = [0 1]; axMaxAng.XLim = [60 180];
        hold(axMaxAng, "off")

        hold(axNormAngle, "on")
        scatter(axNormAngle, cell2mat(clusterStruct(i).surfaceVectorAngle(1,:)), cell2mat(clusterStruct(i).surfaceVectorAngle(2,:)), 10, [0.8500 0.3250 0.0980],'filled', 'MarkerFaceAlpha', 0.8);
        title(axNormAngle, "Angle to center vector")
        xlabel(axNormAngle, "Angle (deg)"); ylabel(axNormAngle, "likelihood"); axNormAngle.YLim = [0 1]; axNormAngle.XLim = [0 90];
        hold(axNormAngle, "off")

        hold(histDistSum, "on")
        histogram(histDistSum, cell2mat(clusterStruct(i).distSum(1,:)), 20, "FaceAlpha", 0.2, "FaceColor", [0 0.4470 0.7410], "EdgeColor", "none");
        xlabel(axDistSum, "Distance (nm)"); histDistSum.XLim = [0 150]; ylabel(histDistSum, "Occurrences");
        hold(histDistSum, "off")
        hold(histMaxAngle, "on")
        histogram(histMaxAngle, cell2mat(clusterStruct(i).maxAng(1,:)), 20, "FaceAlpha", 0.2, "FaceColor", [0.9290 0.6940 0.1250], "EdgeColor", "none");
        xlabel(axMaxAng, "Angle (deg)"); histMaxAngle.XLim = [60 180]; ylabel(histMaxAngle, "Occurrences");
        hold(histMaxAngle, "off")
        hold(histNormAngle, "on")
        histogram(histNormAngle, cell2mat(clusterStruct(i).surfaceVectorAngle(1,:)), 20, "FaceAlpha", 0.2, "FaceColor", [0.8500 0.3250 0.0980], "EdgeColor", "none");
        xlabel(axNormAngle, "Angle (deg)"); histNormAngle.XLim = [0 90]; ylabel(histNormAngle, "Occurrences");
        hold(histNormAngle, "off")
    end

    %% c) Evaluate clusters by likelihood and select most likely
    wbCluster.Title = "Evaluatiung likelihood (Step 3/3)" ;

    disp("Evaluating likelihoods!")
    clusterStruct(i).definiteTrimers = table();

    % Make the working table
    trimerTable = array2table(zeros(0));
    for ii = 1:height(clusterStruct(i).tracesInTrimers)
        trimerTable.(num2str(clusterStruct(i).tracesInTrimers(ii,4))) = sortrows(getMatchingLines(clusterStruct(i).potTrimers.Traces, clusterStruct(i).tracesInTrimers(ii,4), clusterStruct(i).potTrimers.likelihood)',2, "descend")';
    end

    % Iterate through table to find the highest likelihood trimer and move it to the definite trimers
    k = 0;
    while length(trimerTable.Properties.VariableNames) > 3
        % Endless-Loop-Check
        k = k+1;
        if k>100
            error("While loop stuck!")
        end
        
        % Find cluster with highest likelihood
        maxLikeTrimers = cell(length(trimerTable.Properties.VariableNames), 2);
        
        for ii = 1:length(trimerTable.Properties.VariableNames)
            if isempty(trimerTable.(cell2mat(trimerTable.Properties.VariableNames(ii))))
                continue
            end
            breakout = trimerTable.(cell2mat(trimerTable.Properties.VariableNames(ii)))';
            maxLikeTrimers(ii, :) = breakout(1,:); 
        end
        pArray = cell2mat(maxLikeTrimers(:,2));
        maxLike = max(pArray);

        if sum(maxLike == pArray) == 1 || sum(maxLike == pArray) == 3
            foundTrimer = maxLikeTrimers{find(pArray == maxLike, 1),1};
        elseif sum(maxLike == pArray) == 2 || sum(maxLike == pArray) > 3
            possibleTraces = maxLikeTrimers{find(pArray == maxLike, 1),1};
            
            % Find the second most likely trimer for each trace
            subtrimerTable = array2table(zeros(0));
            for iii = 1:length(possibleTraces)
                subtrimerTable.(num2str(possibleTraces(iii))) = trimerTable.(num2str(possibleTraces(iii)))(:,2);
            end

            error("Weiter schreiben") % If this case ever occurs
        else
            error("404: " + sum(maxLike == pArray) + " trimers found")

        end
        
        % Add trimer to definite trimers
        foundTrimerIDX = find(sum(clusterStruct(i).potTrimers.Traces == foundTrimer,2)==3, 1);
        clusterStruct(i).definiteTrimers.Index(k) = clusterStruct(i).potTrimers.Index(foundTrimerIDX);
        clusterStruct(i).definiteTrimers.Traces(k, :) = clusterStruct(i).potTrimers.Traces(foundTrimerIDX, :);
        clusterStruct(i).definiteTrimers.likelihood(k, :) = clusterStruct(i).potTrimers.likelihood{foundTrimerIDX};

        % Delete variables from table
    	trimerTable = removevars(trimerTable, string(foundTrimer));
        trimerTableVarnames = trimerTable.Properties.VariableNames;

        % Remove all trimers sharing traces with the found trimer
        for ii = 1:length(trimerTable.Properties.VariableNames)
            varName = trimerTable.(trimerTableVarnames{ii});
            delList = zeros(1,width(varName));
            for iii = 1:width(varName)
                if any(ismember(foundTrimer, varName{1,iii}))
                    delList(iii) = 1;
                end                
            end
            varName(:,logical(delList)) = [];
            if isempty(varName)
                trimerTable = removevars(trimerTable, trimerTableVarnames{ii});
            else
                trimerTable.(trimerTableVarnames{ii}) = varName;
            end
        end
    end

    % Check if trimers intersect
    disp("Checking for intersections!");
    tempDefTrimers = clusterStruct(i).definiteTrimers;
    delList = zeros(height(tempDefTrimers), 1);
    idxList = tempDefTrimers.Index;
    for ii = 1:height(clusterStruct(i).definiteTrimers)        
        T1 = [  clusterStruct(i).allTraces(tempDefTrimers.Traces(ismember(tempDefTrimers.Index,idxList{1}),1),1:3);...
                clusterStruct(i).allTraces(tempDefTrimers.Traces(ismember(tempDefTrimers.Index,idxList{1}),2),1:3);...
                clusterStruct(i).allTraces(tempDefTrimers.Traces(ismember(tempDefTrimers.Index,idxList{1}),3),1:3)];
        idxList(1, :) = [];
        for iii = 1:numel(idxList)
            % disp("Comparing " + tempDefTrimers.Index{ii} + " and " + idxList{iii});
            T2 = [  clusterStruct(i).allTraces(tempDefTrimers.Traces(ismember(tempDefTrimers.Index,idxList{iii}),1),1:3);...
                    clusterStruct(i).allTraces(tempDefTrimers.Traces(ismember(tempDefTrimers.Index,idxList{iii}),2),1:3);...
                    clusterStruct(i).allTraces(tempDefTrimers.Traces(ismember(tempDefTrimers.Index,idxList{iii}),3),1:3)];
            if checkTriangleIntersect3D(T1, T2)
                % disp("Triangles " + tempDefTrimers.Index{ii} + " and " + idxList(iii) + " are intersecting!")

                % Compare likelihoods of intersecting triangles and delete least likely
                pT1 = tempDefTrimers.likelihood(ii);
                pT2 = tempDefTrimers.likelihood(find(ismember(tempDefTrimers.Index,idxList{iii})));
                
                if pT1 > pT2
                    delList(ismember(tempDefTrimers.Index,idxList{iii})) = 1;
                elseif pT1 < pT2
                    delList(ii) = 1;
                elseif pT1 == pT2
                    delList(ii) = 1; % delete either
                else
                    error("likelihoods of intersecting clusters could not be compared!")
                end
            end
        end
    end
    tempDefTrimers(logical(delList),:) = [];
    clusterStruct(i).definiteTrimers = tempDefTrimers;
    clearvars("tempDefTrimers", "idxList", "delList");

    if plOptions.plotPotTrimers

        % Generate a figure for each cluster
        if plOptions.plotFit
            k = 3;
        else
            k = 2;
        end
         ax3 = nexttile(triFigL, k);      
        hold on
        for ii = 1:height(clusterStruct(i).definiteTrimers)  

            % Plot definite trimers
            pointA = clusterStruct(i).allTraces(clusterStruct(i).definiteTrimers.Traces(ii,1),1:3);
            pointB = clusterStruct(i).allTraces(clusterStruct(i).definiteTrimers.Traces(ii,2),1:3);
            pointC = clusterStruct(i).allTraces(clusterStruct(i).definiteTrimers.Traces(ii,3),1:3);
            trimer = [pointA;pointB;pointC];
            fill3(trimer(:,1), trimer(:,2), trimer(:,3), clusterStruct(i).definiteTrimers.likelihood(ii));        
        end
        colormap("parula");
        title("Definite Trimers", "Interpreter","none");
        xlabel("X (nm)"); ylabel("Y (nm)"); zlabel("Z (nm)");
        xlim = Lims(1,:);
        ylim = Lims(2,:);
        zlim = Lims(3,:);
        colorbar(ax3, "eastoutside");
        view(vP);
        grid(ax3, "on");
        axis equal
    end

    h = height(clusterStruct(i).definiteTrimers);

    % Plot likelihood
    if plOptions.plotlikelihoods
        X = (1:h)/h;
        hold(allLikes, "on")
        plot(allLikes, X, clusterStruct(i).definiteTrimers.likelihood, "LineWidth", 1)
        colormap("jet");
        allLikes.YLim = [0,1];
        ylabel(allLikes, "likelihoods", 'Interpreter','none'); %allLikes.YAxis.Label = text("likelihood", 'Interpreter', 'None');
        allLikes.XAxis.Visible = "off";
        title(allLikes,"Sorted likelihoods per cluster")
        hold(allLikes, "off")
    end

    %% Generate metadata

    % Mean inter-blade distance
    mibD = cell(h,1);
    ibD = zeros(h, 3);
    
    for ii = 1:h
        IDX = find(ismember(clusterStruct(i).potTrimers.Index, clusterStruct(i).definiteTrimers.Index{ii}));
        mibD{ii} = mean([norm(clusterStruct(i).vecAB{IDX}), norm(clusterStruct(i).vecAC{IDX}), norm(clusterStruct(i).vecBC{IDX})]);
        ibD(ii,:) = [norm(clusterStruct(i).vecAB{IDX}), norm(clusterStruct(i).vecAC{IDX}), norm(clusterStruct(i).vecBC{IDX})];
    end
    ibD = ibD(:);

    % Max internal angle
    maxA = cell(h,1);
    for ii = 1:h
        maxA{ii} = clusterStruct(i).maxAng{1,ismember(clusterStruct(i).potTrimers.Index, clusterStruct(i).definiteTrimers.Index{ii})};
    end
    
    % Center of gravity
    COG = cell(h,1);
    for ii = 1:h
        COG{ii} = clusterStruct(i).COG{1,ismember(clusterStruct(i).potTrimers.Index, clusterStruct(i).definiteTrimers.Index{ii})};
    end
    % Add to clusterStruct
    clusterStruct(i).definiteTrimers.MeanInterbladeDistances = mibD;
    clusterStruct(i).definiteTrimers.MaxInternalAngle = maxA;
    clusterStruct(i).definiteTrimers.CenterOfGravity = COG;
    clusterStruct(i).allInterbladeDistances = ibD;

    toc(Step2)
end
close(figWB);

disp("Done! Overall Time:")
toc(tStart)

