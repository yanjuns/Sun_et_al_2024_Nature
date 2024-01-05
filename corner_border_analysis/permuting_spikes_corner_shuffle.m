function cscore_shuffle = permuting_spikes_corner_shuffle(neuron,behav,thresh,env_coor,mask,threshlevel,smoothocc,deltaTall,nboot,binsize)
% This function helps to identify corner cells by generating a shuffled
% corner score distribution. The final output will be a matrix of numCells
% x number of shuffles with shuffled corner score at each entry. 
% Input arguments:
% neuron, behav, thresh: general data
% smoothocc: logical, whether smooth occupancy map or not. 
% nboot: times of shuffling
% binsize: 2 by default
% deltaTall: circular shift amount
% Yanjun Sun 10/12/22; yanjuns@stanford.edu
if ~exist('mask','var') || isempty(mask)
    mask = false;
end
if ~exist('threshlevel','var')||isempty(threshlevel)
    threshlevel = 0.3;
end
if ~exist('smoothocc','var') || isempty(smoothocc)
    smoothocc = true;
end
if ~exist('nboot','var') || isempty(nboot)
    nboot = 1000;
end
if ~exist('deltaTall','var') || isempty(deltaTall)
    deltaTall = randi([round(0.05*length(neuron.time)),length(neuron.time)-round(0.05*length(neuron.time))],nboot,1); %msec
end
if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end

cscore_shuffle = zeros(size(neuron.S,1),nboot);
% time0 = neuron.time;
% S0 = neuron.S;
% trace0 = neuron.trace;
for nE = 1:nboot
    deltaT = deltaTall(nE);
%     neuron.time = time0; neuron.S = S0; neuron.trace = trace0;
    neuronboot = neuron;
    timeboot = (1:length(neuronboot.time)) + deltaT;
    idx = timeboot > length(neuronboot.time);
    timeboot(idx) = timeboot(idx) - length(neuronboot.time);
    [~,index] = sort(timeboot);
    neuronboot.S = neuronboot.S(:,index);
    neuronboot.trace = neuronboot.trace(:,index);
    [firingrateAll,~,~,~] = calculate_firing_ratemap(neuronboot,behav,thresh,binsize,smoothocc);
    ratemap = cellfun(@(x) filter2DMatrices(x,1), firingrateAll, 'uni', 0);
    cmetrics = compute_corner_score(ratemap,env_coor,mask,false,threshlevel);
    cscore_shuffle(:,nE) = cmetrics.cornerscore;
end

end
