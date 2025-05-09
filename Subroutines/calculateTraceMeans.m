function [traces_means] = calculateTraceMeans(traces)
    fns = fieldnames(traces);

    for k=1:size(fns,1)
    % read cluster data from structure
    ClusXYZRAW = traces.(fns{k});
    ClusXYZ = ClusXYZRAW(:,1:4);

    % calculate trace center
    [uv_tid, ~, id_tid] = unique(ClusXYZ(:,4)); 
    ClusterTraces_AVG = [accumarray(id_tid,ClusXYZ(:,1),[],@mean) accumarray(id_tid,ClusXYZ(:,2),[],@mean) accumarray(id_tid,ClusXYZ(:,3),[],@mean)];

    % run DBSCAN to find signals from the same fluorophore – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ClusterIDs = dbscan(ClusterTraces_AVG(:,1:2),25,2);
    
    % add column with cluster ID to main coordinates table
    ClusXYZ = [ClusterTraces_AVG ClusterIDs];
    
    % split ClusXYZ into 2 submatrices - one with single signals and one with multiple signals per channel 
    singleloc = ClusXYZ(ClusXYZ(:,4) == -1,:);
    multiloc = ClusXYZ(ClusXYZ(:,4) > -1,:);
    
    % calculate averages of multiloc clusters (i.e. signals from one channel)
    [uv, ~, id] = unique(multiloc(:,4));
    MyMeanX = [accumarray(id,multiloc(:,1),[],@mean)];
    MyMeanY = [accumarray(id,multiloc(:,2),[],@mean)];
    MyMeanZ = [accumarray(id,multiloc(:,3),[],@mean)];
    ClID = [accumarray(id,multiloc(:,4),[],@mean)];
    multiloc = [MyMeanX, MyMeanY, MyMeanZ, ClID];
    
    % merge submatrices into single matrix
    ClusterTraces_AVG = cat(1,singleloc(:,1:4),multiloc(:,1:4));
    traces_means.(fns{k}) = ClusterTraces_AVG(:,1:3);

    end
end

