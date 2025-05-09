function outputStruct = getClusterStruct(cluster)
%GETCLUSTERSTRUCT Constructs the struct-variable for the input cluster that
%is used in the further analysis
%INPUT : 2D-array with mean 3D data for all traces in a potential cluster in
%colums 1 - 3 and additional data in columns 4 onward

outputStruct = struct();
outputStruct.RAWinput = cluster;

% Make a table with the cluster coordinates and indivudual indentifyers for
% all traces


outputStruct.allTraces = cluster(:,1:3);
outputStruct.allTraces(:,4) =  1:height(cluster);

end

