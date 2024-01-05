function [posgrid,xnbins,ynbins] = pos_map(pos, binsize)

if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end

xAxis = 0:binsize:ceil(max(pos(:,1)+binsize)); %binedges of x position
yAxis = 0:binsize:ceil(max(pos(:,2)+binsize)); %binedges of y position
[N,~,~,binY,binX] = histcounts2(pos(:,2), pos(:,1), yAxis, xAxis);

position_vec = 1:(size(N,1)*size(N,2));
position_mat = reshape(position_vec, [size(N,1),size(N,2)]);

% store grid
xnbins = size(N,1);
ynbins = size(N,2);
posgrid = zeros(length(pos), xnbins*ynbins);

% loop over positions
for idx = 1:length(pos)
    posVec = position_mat(binY(idx), binX(idx));
    posgrid(idx, posVec) = 1;
    
end

end