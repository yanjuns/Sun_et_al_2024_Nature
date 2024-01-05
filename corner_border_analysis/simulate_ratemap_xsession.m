function [neuronSim,ratemapSim] = simulate_ratemap_xsession(datapath)
%This function aims to simulate the activity of a neuron that fires
%constantly at the overall mean firing rate of the session. The simulated
%activity is 1 and other entries are 0. 
%Yanjun Sun, 12/18/2022

if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end
cd(datapath)
%simulate the activity of a neuron
load('neuronIndividualsf.mat','neuronIndividualsf')
load('behavIndividualsf.mat','behavIndividualsf')
load('firingrateAll.mat', 'meanFRindividuals')
meanFR = mean(cellfun(@mean, meanFRindividuals));
neuronSim = cell(1,length(neuronIndividualsf));
firingrateAll = cell(1,length(neuronIndividualsf));
countAll = cell(1,length(neuronIndividualsf));
countTime = cell(1,length(neuronIndividualsf));
smoccu = cell(1,length(neuronIndividualsf)); %smoothed countTime
ratemapSim = cell(1,length(neuronIndividualsf));
for ii = 1:length(neuronIndividualsf)
    neuron = neuronIndividualsf{ii};
    neuronSim_ea = struct;
    neuronSim_ea.S = zeros(1,size(neuron.S,2));
    frq = round(15/meanFR);
    neuronSim_ea.S(1:frq:end) = 1;
    neuronSim_ea.time = neuron.time;
    neuronSim_ea.pos = neuron.pos;
    neuronSim{ii} = neuronSim_ea;
end
%compute the ratemap of the simulated neuron
thresh = 0;
for ii = 1:length(neuronIndividualsf)
    [firingrateAll{ii},countAll{ii},countTime{ii},smoccu{ii}] = ...
        calculate_firing_ratemap(neuronSim{ii},behavIndividualsf{ii},thresh,2,true);
end
for ii = 1:length(firingrateAll)
    fr = firingrateAll{ii};
    % to remove potential inf in the unsmoothed ratemap
    idx = cellfun(@isinf, fr,'UniformOutput',false);
    for jj = 1:length(fr)
        fr{jj}(idx{jj}) = 0;
    end
    % smoothed ratemap
    ratemapSim{ii} = cellfun(@(x) filter2DMatrices(x,1), fr, 'uni', 0);
end

end

