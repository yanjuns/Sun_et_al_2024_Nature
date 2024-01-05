function [firingrateAll,countAll,countTime,smoccu] = calculate_firing_ratemap(neuron,behav,thresh,binsize,smoothocc,timeInterval,dataType)
% This function aims to calculate the unsmoothed raw firing rate map for a
% miniscope recording session, it works similar to the previous function
% named 'calculatingCellSpatialForSingleData', but with improved algorithm
% and efficiency. 
% Input arg:
% neuron: a source2D variable extracted from CNMF-E
% behav: s struct contains animals behavior information
% thresh: threshold for determine a spike
% binsize: actuall size for each bin, cm
% dataType, neuron.S or neuron.trace
% Output:
% countAll: spike counts in each bin
% countTime: animal's occupancy time in each bin
% firingrateAll: countAll./countTime
% Yanjun Sun, Stanford University, 9/9/2019

if ~exist('dataType','var') || isempty(dataType)
    dataType = 'S';
end
if ~exist('timeInterval','var') || isempty(timeInterval)
    timeInterval = median(diff(behav.time))/1000;
end
if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end
if ~exist('smoothocc','var') || isempty(smoothocc)
    smoothocc = false;
end
%% calculate occupancy time
xAxis = 0:binsize:ceil(max(behav.position(:,1)+binsize)); %binedges of x position
yAxis = 0:binsize:ceil(max(behav.position(:,2)+binsize)); %binedges of y position
[timeframe,~,~] = histcounts2(behav.position(:,1), behav.position(:,2), xAxis, yAxis);
countTime = (timeframe*timeInterval)';

%% calculate spike counts
neuron.pos = interp1(behav.time,behav.position,neuron.time); 
numcells = size(neuron.S,1);
countAll = cell(1,numcells);
for ii = 1:numcells;
    if strcmpi(dataType,'S')
        idx = find(neuron.S(ii,:)> thresh(ii));
    elseif strcmpi(dataType,'trace')
        idx = find(neuron.trace(ii,:)> thresh(ii));
    end
    [count,~,~] = histcounts2(neuron.pos(idx,1), neuron.pos(idx,2), xAxis, yAxis);
    countAll{ii} = count';
end
%% firing rate maps
if smoothocc
    smoccu = smooth_occupancy(countTime);
else
    smoccu = countTime;
end

firingrateAll = cell(1,numcells);
for jj = 1:numcells
    firingrateAll{jj} = countAll{jj}./smoccu;
end

end