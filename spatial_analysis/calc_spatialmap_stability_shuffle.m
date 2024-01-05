function Rxy_shuffle = calc_spatialmap_stability_shuffle(neuron,behav,thresh,binsize,nboot,deltaTall)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~exist('nboot','var') || isempty(nboot)
    nboot = 1000;
end
if ~exist('deltaTall','var') || isempty(deltaTall)
    deltaTall = randi([round(0.05*length(neuron.time)),length(neuron.time)-round(0.05*length(neuron.time))],nboot,1); %msec
end
if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end

Rxy_shuffle = zeros(length(thresh),nboot);
for nE = 1:nboot
    deltaT = deltaTall(nE);
    neuronboot = neuron;
    timeboot = (1:length(neuronboot.time)) + deltaT;
    idx = timeboot > length(neuronboot.time);
    timeboot(idx) = timeboot(idx) - length(neuronboot.time);
    [~,index] = sort(timeboot);
    neuronboot.S = neuronboot.S(:,index);
    neuronboot.trace = neuronboot.trace(:,index);
    Rxy = calc_spatialmap_stability(neuronboot,behav,thresh,binsize);
    Rxy_shuffle(:,nE) = Rxy;
end