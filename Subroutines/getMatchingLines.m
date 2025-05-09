function lineCell = getMatchingLines(array, n, prob)
    %GETMATCHINGLINES Takes an array and retuns all lines that contain n as
    %a cell array
    lineCell = {};
    for i = 1:length(array)
        if ismember(n, array(i,:))
            k = width(lineCell)+1;
            lineCell{1,k} = array(i,:);
            lineCell{2,k} = cell2mat(prob(i));
        end
    end
end