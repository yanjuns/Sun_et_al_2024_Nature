function [neuron1_mch, behav1_mch, neuron2_mch, behav2_mch] = downsample_match_position(neuron1,behav1,neuron2,behav2,binsize)
% This function aims to matche two sessions based on behav position
% based on Hardcastle et al (2017) and Kiah Hardcastle's code
% Yanjun Sun, Stanford University, 9/11/2019
if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end
if isempty(neuron1.pos)
    neuron1.pos = interp1(behav1.time, behav1.position, neuron1.time);
end
if isempty(neuron2.pos)
    neuron2.pos = interp1(behav2.time, behav2.position, neuron2.time);
end
%% get parameters and initial the matrices
x1 = behav1.position(:,1); y1 = behav1.position(:,2);
x2 = behav2.position(:,1); y2 = behav2.position(:,2);
xmax = max([max(x1);max(x2)]);
ymax = max([max(y1);max(y2)]);
xAxis = 0:binsize:ceil(xmax+binsize);
yAxis = 0:binsize:ceil(ymax+binsize);
%initialize matrices
map1 = zeros(numel(yAxis), numel(xAxis));
map2 = zeros(numel(yAxis), numel(xAxis));
dsidx1 = [];
dsidx2 = [];
%% Determine number of observations in each bin in session1
for i = 1:numel(xAxis)-1
    for j = 1:numel(yAxis)-1
        ind1 = find(x1 >= xAxis(i) & x1 < xAxis(i+1) & y1 >= yAxis(j) & y1 < yAxis(j+1));
        map1(j,i)= numel(ind1);
        %same thing for map 2
        ind2 = find(x2 >= xAxis(i) & x2 < xAxis(i+1) & y2 >= yAxis(j) & y2 < yAxis(j+1));
        map2(i,j) = numel(ind2);
        % Number of observations to keep in each map
        numToKeep = min(numel(ind1),numel(ind2));        
        dsidx1 = [dsidx1; datasample(ind1, numToKeep, 'Replace', false)];
        dsidx2 = [dsidx2; datasample(ind2, numToKeep, 'Replace', false)];
    end
end
dsidx1 = sort(dsidx1);
dsidx2 = sort(dsidx2);
%% result output
neuron1_mch = struct; behav1_mch = struct;
neuron2_mch = struct; behav2_mch = struct;

behav1_mch.position = behav1.position(dsidx1,:);
behav1_mch.time = behav1.time(dsidx1,:);
behav1_mch.speed = behav1.speed(:,dsidx1);

behav2_mch.position = behav2.position(dsidx2,:);
behav2_mch.time = behav2.time(dsidx2,:);
behav2_mch.speed = behav2.speed(:,dsidx2);

behavtime1 = behav1.time;
idx1 = knnsearch(behavtime1, neuron1.time);
neurontime1 = behavtime1(idx1);
idxn1 = ismember(neurontime1,behav1_mch.time);
neuron1_mch.pos = neuron1.pos(idxn1,:);
neuron1_mch.time = neuron1.time(idxn1,:);
neuron1_mch.S = neuron1.S(:,idxn1);
neuron1_mch.trace = neuron1.trace(:,idxn1);

behavtime2 = behav2.time;
idx2 = knnsearch(behavtime2, neuron2.time);
neurontime2 = behavtime2(idx2);
idxn2 = ismember(neurontime2,behav2_mch.time);
neuron2_mch.pos = neuron2.pos(idxn2,:);
neuron2_mch.time = neuron2.time(idxn2,:);
neuron2_mch.S = neuron2.S(:,idxn2);
neuron2_mch.trace = neuron2.trace(:,idxn2);

end