function potClusters = getPotentialTrimers(clusterStruct, threshold)
%GETPOTENTIALTRIMERS Generates an array of all potential clusters solely
%based on distances

    %threshold = 40;
    
    allTraces = clusterStruct.tracesWithoutTrimer;
    [D, I] = pdist2(allTraces(:,1:3), "euclidean");

    allTraces = allTraces(:, :) <= threshold;
    
    potClusters = cell(length(allTraces), 1);
    
    for i = 1:width(allTraces)
        potClusters{i} = find(allTraces(:,i));
    
    end

end
