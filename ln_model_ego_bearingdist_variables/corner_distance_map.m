function [cornerdistgrid,cdistVec] = corner_distance_map(cornerdist,binsize)

% binedges of corner distance
xAxis = 0:binsize:ceil(max(cornerdist + binsize)); 
[N,~,cornerdistVec] = histcounts(cornerdist, xAxis);
cornerdistgrid = zeros(numel(cornerdist),size(N,2));

for ii = 1:numel(cornerdist)

    % figure out the cornerdist index
    idx = cornerdistVec(ii);
    cornerdistgrid(ii,idx) = 1;

end

cdistVec = xAxis(2:end);
return