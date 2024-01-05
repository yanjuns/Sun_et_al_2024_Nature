function [N,Ns,xAxis,yAxis] = anat_density(C, xAxis, yAxis, sigma, binsize)
%this function aims to calculate the anatomical density of neurons from
%miniscope imaging. 
%C is centroid of neuron
%N is the raw map, Ns is the smoothed map
if ~exist('binsize','var') || isempty(binsize)
    binsize = 25;
end
if ~exist('sigma','var') || isempty(sigma)
    sigma = 1;
end
if ~exist('yAxis','var') || isempty(yAxis)
    yAxis = round(min(C(:,2)-binsize)):binsize:ceil(max(C(:,2)+binsize)); %binedges of y position
end
if ~exist('xAxis','var') || isempty(xAxis)
    xAxis = round(min(C(:,1)-binsize)):binsize:ceil(max(C(:,1)+binsize)); %binedges of x position
end

[N,~,~] = histcounts2(C(:,1), C(:,2), xAxis, yAxis);
Ns = imgaussfilt(N,sigma);

end

