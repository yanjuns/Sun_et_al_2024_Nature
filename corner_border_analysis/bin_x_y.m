function [xbin,ybin] = bin_x_y(x,y,binsize)
% this function aims to bin the vector of x and meanwhile average the value
% of y in each bined x
% yanjuns@stanford.edu 11/22/22
if ~exist('binsize','var')||isempty(binsize)
    binsize = 1;
end

edges = [0:binsize:max(x)+1];
[xbin,~,bin] = histcounts(x,edges);
ybin = NaN(1,length(xbin));
for ii = 1:length(xbin)
    ybin(ii) = nanmean(y(bin ==ii));
end

end

